import { readFile, writeFile, mkdir } from 'node:fs/promises';
import { join, dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { performance } from 'node:perf_hooks';
import { spawnSync } from 'node:child_process';
import { calculateSquare } from './calculo/canonico.mjs';
import { perfil, nystrom } from './calculo/nystrom.mjs';
import { curva } from './calculo/actividades.mjs';
import { todosAjustes, extrapolar } from './calculo/ajustes.mjs';
import { comprobar } from './comprobar.mjs';

const base = dirname(fileURLToPath(import.meta.url));
const destino = resolve(process.argv[2] || join(base, 'corrida_' + new Date().toISOString().replace(/[:-]/g,'')));
await mkdir(destino, {recursive:false});
const leer = nombre => readFile(join(base,nombre),'utf8');
const guardar = (nombre,valor) => writeFile(join(destino,nombre),valor,{flag:'wx'});
const json = data => JSON.stringify(data,null,2)+'\n';
const csv = data => [Object.keys(data[0]).join(','),...data.map(r=>Object.values(r).join(','))].join('\n')+'\n';
const controles = async nombre => {
  const [cab,...datos]=(await leer(nombre)).trim().split(/\r?\n/);
  const nombres=cab.split(',');
  return datos.map(l=>Object.fromEntries(l.split(',').map((x,i)=>[nombres[i],Number(x)])));
};
const inicio=performance.now(), tiempos={};
function marca(nombre,desde) {tiempos[nombre]=(performance.now()-desde)/1000;}
const ref=JSON.parse(await leer('referencias/canonico_archivado.json'));
const completos=ref.map(({n})=>calculateSquare(n));
for(const [i,r] of completos.entries()) {
  if(r.logQ.toFixed(7)!==ref[i].q || r.gap.toFixed(5)!==ref[i].g) throw Error('Predictor o brecha en n='+r.n);
  if(!(r.pMin>0 && r.pMax<.25 && r.hMin>0 && r.gapLower>0 && r.stationarity<1e-10 && r.maxMassError<1e-10)) throw Error('Interioridad o residuo en n='+r.n);
}
await guardar('canonico.csv',csv(completos));
const can=completos.map(r=>({n:r.n,logQpred:r.logQ,gap:r.gap,pmin:r.pMin,pmax:r.pMax,hmin:r.hMin}));
const fits=todosAjustes(can);
marca('canonico_y_ajustes',inicio);
const ti=performance.now();
const rejillas=await controles('ny_resultados.csv'), valores=[];
let ne=0, p;
for(const esperado of rejillas) {
  if(ne!==esperado.perfil) {ne=esperado.perfil;p=perfil(ne);}
  const r=nystrom(esperado.m,p);
  for(const [campo,v] of Object.entries(r)) if(Math.abs(v-esperado[campo])>2e-10*Math.max(1,Math.abs(v))) throw Error('Cuadratura '+r.m+','+ne+' '+campo);
  valores.push(r);
}
const punto=m=>valores.find(r=>r.m===m && r.perfil===m);
const ext=(i,j,campo)=>extrapolar(punto(i),punto(j),campo);
const nys={I2:ext(192,384,'I2'),det2:Math.exp(ext(192,384,'logdet2')),par:ext(192,384,'C'),impar:ext(193,385,'C'),mallas:[[192,384],[193,385]]};
await guardar('ny_resultados.csv',csv(valores));
marca('nystrom',ti);
const ta=performance.now();
const esperadoA=await controles('actividades_793.csv');
const actividades=[...new Set(esperadoA.map(x=>x.n))].flatMap(n=>curva(n));
if(actividades.length!==793) throw Error('Número de actividades inesperado');
for(const [i,r] of actividades.entries()) for(const [campo,v] of Object.entries(r)) if(Math.abs(v-esperadoA[i][campo])>2e-10*Math.max(1,Math.abs(v))) throw Error('Actividad '+i+' '+campo);
await guardar('actividades_793.csv',csv(actividades));
marca('actividad_finita',ta);
const tb=performance.now(), executable=process.env.PYTHON || (process.platform==='win32'?'python':'python3');
const py=spawnSync(executable,[join(base,'barrera_4_5.py')],{encoding:'utf8',env:{...process.env,PYTHONIOENCODING:'utf-8',OPENBLAS_NUM_THREADS:'1',OMP_NUM_THREADS:'1'}});
if(py.error || py.status!==0) throw Error('Python/NumPy: '+(py.error?.message||py.stderr));
await guardar('barrera_4_5_salida.txt',py.stdout);
const barrera={};
for(const ln of py.stdout.split(/\r?\n/)){
  const m=ln.match(/^n=(4|5): Λ = 10\^\s*([-\d.]+) → máx π̂ = ([\d.]+),.*en el borde/);
  if(m) barrera[m[1]]={Lambda:10**Number(m[2]),maxPi:Number(m[3])};
}
marca('barrera',tb);
const exact=JSON.parse(await leer('referencias/conteos.json'));
const canonico5=JSON.parse(await leer('referencias/canonico_n5.json')).maxPi;
const anunciados=JSON.parse(await leer('referencias/anunciados.json'));
const datos={can,exact,fits,nys,actividades,alpha:[1.944000752,1.944001082],alpha0:1.944001,canonico5,barrera,anunciados};
await guardar('verificacion_datos.json',json(datos));
const nombres=['00_resumen.tex','01a_introduccion.tex','01b_introduccion.tex','09_numerico.tex','10_discusion.tex'];
const contenido=Object.fromEntries(await Promise.all(nombres.map(async n=>[n,await leer('../articulo/secciones/'+n)])));
const informe=comprobar(contenido,datos);
await guardar('verificacion.json',json(informe));
marca('total',inicio);
await guardar('tiempos_segundos.json',json(tiempos));
if(!informe.ok) throw Error(informe.errores.join('\n'));
console.log(informe.comprobaciones.length+' verificaciones correctas. Archivos en '+destino);
