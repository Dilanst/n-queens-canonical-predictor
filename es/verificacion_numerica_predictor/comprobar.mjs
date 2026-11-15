// Verificación reproducible de las afirmaciones numéricas del artículo B. Node.js >= 18.
export function comprobar(contenidos, datos) {
  const errores = [], comprobaciones = [];
  const registrar = (id, bien, detalle = '') => {
    comprobaciones.push(id);
    if (!bien) errores.push(id + (detalle ? ' — ' + detalle : ''));
  };
  const resultado = () => ({ok: errores.length === 0, comprobaciones, errores});
  const decimal = (v, cifras) => Number(v).toFixed(cifras);
  const cerca = (x, y, tol) => Number.isFinite(x) && Number.isFinite(y) && Math.abs(x-y) <= tol;
  const finito = v => typeof v === 'number' && Number.isFinite(v);
  const lf = n => { let s = 0; for (let j=2; j<=n; j++) s += Math.log(j); return s; };
  const nombres = ['00_resumen.tex','01a_introduccion.tex','01b_introduccion.tex','09_numerico.tex','10_discusion.tex'];
  const [res, intro, relacionados, num, disc] = nombres;
  const textos = {}, cubiertos = {};
  try {
    for (const nombre of nombres) {
      registrar(nombre + ' presente', typeof contenidos[nombre] === 'string');
      // Conservar posiciones en el texto normalizado. Ignorar citas y etiquetas al buscar cifras.
      textos[nombre] = (contenidos[nombre] ?? '').replace(/(?<!\\)%[^\r\n]*/g,'')
        .replace(/\\(?:cite|eqref|ref)\{[^}]*\}/g,'REF')
        .replace(/\{,\}/g,'.').replace(/(?<=\d),(?=\d)/g,'.').replace(/\$/g,'')
        .replace(/\\[,;!]/g,' ').replace(/~/g,' ').replace(/\s+/g,' ').trim();
      cubiertos[nombre] = new Set();
    }
    const R = String.raw;
    // Cada @ representa una cifra empírica cuyo valor se calcula independientemente.
    const patron = "(?<![\\d.])(-?\\d+(?:\\.\\d+)?)(?!\\d|\\.\\d)";
    function cotejar(nombre, id, plantilla, esperados) {
      const re = new RegExp(plantilla.replaceAll('@', patron), 'gd');
      const encontrados = [...textos[nombre].matchAll(re)];
      const etiqueta = nombre + ' · ' + id;
      if (encontrados.length !== 1) {
        registrar(etiqueta, false, 'se esperaba una declaración, encontradas ' + encontrados.length);
        return;
      }
      const m = encontrados[0];
      const vistos = m.slice(1);
      for (const [i] of m.indices.slice(1)) cubiertos[nombre].add(i);
      registrar(etiqueta, vistos.length === esperados.length && vistos.every((v,i) => v === String(esperados[i])),
        'artículo [' + vistos.join(', ') + '], referencia [' + esperados.join(', ') + ']');
    }
    const formas = (nombre, id, re) => registrar(nombre + ' · ' + id, re.test(textos[nombre]));
    const {can, exact, fits, nys, actividades, alpha, alpha0, anunciados} = datos;
    const ns = [20,27,40,60,80,100,130,160,200,240];
    const ne = Array.from({length:12}, (_,i)=>i+16);
    const na = [6,7,8,10,12,16,20,27,40,60,80,120,160];
    const porN = new Map(can.map(f=>[f.n,f]));
    registrar('Canónico · diez tamaños distintos en la serie tabulada', ns.every(n=>can.filter(f=>f.n===n).length===1));
    registrar('Canónico · condición de interioridad y brecha', ns.every(n=>{
      const f = porN.get(n);
      return f && [f.logQpred,f.pmin,f.pmax,f.hmin,f.gap].every(finito) &&
        f.pmin>0 && f.pmin<=f.pmax && f.pmax<0.25 && f.hmin>0 && f.gap>0;
    }));
    registrar('Conteos · enteros positivos en el intervalo exacto', ne.every(n=>typeof exact[n]==='string' && /^[1-9]\d*$/.test(exact[n])));
    registrar('Referencia · intervalo y punto de referencia', Array.isArray(alpha) && alpha.length===2 &&
      [alpha0,...alpha].every(finito) && alpha[0]<alpha0 && alpha0<alpha[1]);
    if (errores.length) return resultado();
    const cp = (n,a) => porN.get(n).logQpred - lf(n) - n*(1-a);
    const ce = n => Math.log(Number(exact[n])) - lf(n) - n*(1-alpha0);
    const diferencia = n => Math.log(Number(exact[n])) - porN.get(n).logQpred;
    // Ordenadas al origen con QR modificado, columnas (100/n)^j.
    function intercepto(tamanos, grado, a) {
      const q=[], r=[], proy=[];
      const y = tamanos.map(n=>cp(n,a));
      const punto = (x,z)=>x.reduce((s,e,i)=>s+e*z[i],0);
      for (let j=0; j<=grado; j++) {
        const v=tamanos.map(n=>(100/n)**j);
        for (let k=0;k<j;k++) {
          r[k][j]=punto(q[k],v);
          for (let i=0;i<v.length;i++) v[i]-=r[k][j]*q[k][i];
        }
        r[j]=[]; r[j][j]=Math.hypot(...v);
        q[j]=v.map(e=>e/r[j][j]); proy[j]=punto(q[j],y);
      }
      for (let j=grado;j>=0;j--) {
        for (let k=j+1;k<=grado;k++) proy[j]-=r[j][k]*proy[k];
        proy[j]/=r[j][j];
      }
      return proy[0];
    }
    const fcalc = fits.map(f=>{
      const tamanos=ns.filter(n=>n>=f.lo && n<=f.hi);
      if (![2,3].includes(f.grado) || tamanos.length<=f.grado || !tamanos.includes(f.lo) || !tamanos.includes(f.hi))
        throw Error('Dominio o grado de ajuste desconocido: '+[f.lo,f.hi,f.grado]);
      const [b0,b1,b2]=[alpha0,...alpha].map(a=>intercepto(tamanos,f.grado,a));
      const g={...f,central:b0,min:Math.min(b1,b2),max:Math.max(b1,b2)};
      registrar('Ajuste '+[f.lo,f.hi,f.grado].join('/')+' · coeficientes frente a los datos canónicos',
        ['central','min','max'].every(k=>cerca(g[k],f[k],5e-8)));
      return g;
    });
    registrar('Sensibilidad · ocho ajustes distintos', fits.length===8 && new Set(fits.map(f=>[f.lo,f.hi,f.grado].join('/'))).size===8);
    const central = g => fcalc.find(f=>f.lo===40 && f.hi===240 && f.grado===g).central;
    const banda = grupo => {
      if (!grupo.length) throw Error('Falta un grupo de ajustes de sensibilidad');
      return [decimal(Math.min(...grupo.map(f=>f.min)),4),decimal(Math.max(...grupo.map(f=>f.max)),4)];
    };
    const tardia = banda(fcalc.filter(f=>[40,60].includes(f.lo)));
    const amplia = banda(fcalc.filter(f=>f.lo===20 && f.hi>=27));
    // Estimación a cuatro decimales a partir del centro de la banda publicada.
    const estimacion = decimal((Number(tardia[0])+Number(tardia[1]))/2,4);
    const C = Math.log(2)/2 - 0.5 - 1.5*nys.I2 - Math.log(nys.det2)/2;
    registrar('Nyström · parámetros y valor de la fórmula', [nys.I2,nys.det2,nys.par,nys.impar,C].every(finito) &&
      nys.I2>0 && nys.det2>0 && cerca(nys.par,C,5e-7) && cerca(nys.impar,C,5e-7));
    registrar('Nyström · parejas de mallas', JSON.stringify(nys.mallas)==='[[192,384],[193,385]]');
    const bounds = alpha.map(v=>decimal(v,9));
    const exactos=ne.map(ce);
    const exterior=[decimal(Math.floor(Math.min(...exactos)*1000)/1000,3),decimal(Math.ceil(Math.max(...exactos)*1000)/1000,3)];
    cotejar(res,'estimación agregada, extrapolación y fórmula directa',
      R`Numéricamente, C_\{\\rm pred\}\\approx\s*@\s*: la extrapolación del predictor en tableros finitos da @ y una evaluación directa de la fórmula da @`,
      [estimacion,estimacion,decimal(C,5)]);
    cotejar(intro,'intervalo de Nobel, Agrawal y Boyd',R`acotaron \\alpha en \[@,\s*@\]`,bounds);
    cotejar(intro,'estimación y banda de los ajustes 40/60',
      R`Numéricamente, C_\{\\rm pred\}\\approx\s*@\s*: en todo el intervalo de \\alpha citado, los ajustes del predictor desde n\s*=\s*@ o @ dan entre @ y @, y una evaluación directa de la fórmula da @`,
      [estimacion,40,60,...tardia,decimal(C,5)]);
    cotejar(intro,'banda exterior para los conteos exactos y referencia Phi0',
      R`Con la referencia \\Phi_0 de esa sección, los conteos exactos para @\\le\s*n\\le\s*@ dan \\log Q\(n\)-\\log n!-n\\Phi_0 entre @ y @`,
      [ne[0],ne.at(-1),...exterior]);
    cotejar(num,'intervalo de Nobel, Agrawal y Boyd',R`dan \\alpha\\in\[@,\s*@\]`,bounds);
    cotejar(num,'referencia de ambas tablas',R`como referencia \\alpha_0\s*=\s*@ y \\Phi_0\s*=\s*1-\\alpha_0`,[decimal(alpha0,6)]);
    formas(num,'definición del predictor centrado',/C'_\{\{\\rm pred\},0\}\(n\)=\\log Q_\{\\rm pred\}\(n\)-\\log n!-n\\Phi_0/);
    formas(num,'interioridad en los diez tamaños, en punto flotante',/interioridad en los diez tamaños tabulados, en aritmética de punto flotante/);
    cotejar(num,'ordenadas de grado dos y tres con alpha0',
      R`Con \\alpha_0 y @\\le\s*n\\le\s*@, el ajuste C\+a/n\+b/n\^2 da C\s*=\s*@ y el que añade c/n\^3 da C\s*=\s*@`,
      [40,240,decimal(central(2),5),decimal(central(3),5)]);
    cotejar(num,'sensibilidad con tamaños iniciales 40 o 60',
      R`los ajustes que empiezan en @ o @ dan extremos redondeados de @ y @`,[40,60,...tardia]);
    cotejar(num,'sensibilidad al incluir 20 y 27',
      R`Al incluir @ y @ los extremos redondeados son @ y @`,[20,27,...amplia]);
    cotejar(num,'estimación global de los ajustes',R`Así se obtiene la estimación C_\{\\rm pred\}\\approx\s*@`,[estimacion]);
    cotejar(num,'integral del cuadrado, determinante regularizado y fórmula',
      R`Se obtiene I_2\\approx\s*@ y \\det_2\(I\+J\)\\approx\s*@, y la fórmula del Teorema\s*REF da C_\{\\rm pred\}\\approx\s*@`,
      [decimal(nys.I2,8),decimal(nys.det2,9),decimal(C,6)]);
    cotejar(num,'mallas emparejadas por paridad y extrapolaciones',
      R`En las parejas de mallas de tamaños @ y @ \(pares\), y @ y @ \(impares\), las extrapolaciones dieron, respectivamente, @ y @`,
      [...nys.mallas.flat(),decimal(nys.par,11),decimal(nys.impar,11)]);
    cotejar(num,'alcance de la referencia de los conteos exactos',
      R`Los valores exactos de Q\(n\) se conocen hasta n\s*=\s*@\s*REF`,[Math.max(...Object.keys(exact).map(Number))]);
    formas(num,'definición del conteo exacto centrado',/C'_0\(n\)=\\log Q\(n\)-\\log n!-n\\Phi_0/);
    // Valor anunciado de Q(28): dato informativo, fuera de las tablas y de la figura.
    const q28 = anunciados?.['28']?.valor;
    registrar('Conteos · valor anunciado de Q(28), entero positivo', typeof q28==='string' && /^[1-9]\d*$/.test(q28));
    if (typeof q28==='string') cotejar(num,'valor anunciado de Q(28), sólo informativo',
      R`resulta C'_0\(@\)\\approx\s*@`,[28,decimal(Math.log(Number(q28))-lf(28)-28*(1-alpha0),4)]);
    cotejar(num,'diferencias exacto menos predictor',
      R`En n\s*=\s*@ y n\s*=\s*@ se tiene \\log Q\(n\)-\\log Q_\{\\rm pred\}\(n\)\s*=\s*@ y @`,
      [20,27,decimal(diferencia(20),4),decimal(diferencia(27),4)]);
    cotejar(disc,'residuos a tres decimales',R`R_Q\(@\)\\approx\s*@ y R_Q\(@\)\\approx\s*@`,
      [20,decimal(diferencia(20),3),27,decimal(diferencia(27),3)]);
    function tabla(etiqueta, tamanos, valor, expr) {
      const seccion = textos[num], id = num+' · '+etiqueta;
      const tablas = [...seccion.matchAll(/\\begin\{table\}(?:\[[^\]]*\])?[\s\S]*?\\end\{table\}/g)];
      const candidatas = tablas.filter(m=>m[0].includes('\\label{'+etiqueta+'}'));
      if (candidatas.length!==1) {registrar(id,false,'etiqueta ausente o duplicada');return;}
      const m=candidatas[0], cuerpo=/\\begin\{tabular\}\{[^}]+\}([\s\S]*?)\\end\{tabular\}/d.exec(m[0]);
      if (!cuerpo) {registrar(id,false,'falta el entorno tabular');return;}
      const filas=cuerpo[1].split('\\\\').map(f=>f.replace(/\\(?:toprule|midrule|bottomrule)\b/g,'').trim()).filter(Boolean);
      const conjuntos=filas.map(f=>f.split('&').map(x=>x.trim()));
      const tramos=tamanos.length/2;
      registrar(id+' · estructura de dos grupos', conjuntos.length===4 && conjuntos.every(a=>a.length===tramos+1));
      for (let j=0;j<2;j++) {
        const h=conjuntos[2*j] || [], v=conjuntos[2*j+1] || [];
        registrar(id+' · encabezados del grupo '+(j+1),h[0]==='n' && v[0]===expr);
        for (let k=0;k<tramos;k++) {
          const n=tamanos[j*tramos+k], calculado=decimal(valor(n),4);
          registrar(id+' · n='+n,h[k+1]===String(n) && v[k+1]===calculado,
            'artículo n='+h[k+1]+', '+v[k+1]+', referencia '+calculado);
        }
      }
      const inicio=m.index+cuerpo.indices[1][0];
      for (const d of cuerpo[1].matchAll(/\d+\.\d+/g)) cubiertos[num].add(inicio+d.index);
    }
    tabla('tab:predictor',ns,n=>cp(n,alpha0),"C'_{{\\rm pred},0}(n)");
    tabla('tab:exacto',ne,ce,"C'_0(n)");
    registrar('Actividad finita · 793 observaciones', actividades.length===13*61);
    for (const n of na) {
      const serie=actividades.filter(f=>f.n===n).sort((a,b)=>a.Lambda-b.Lambda);
      registrar('Actividad finita · n='+n+', 61 actividades y rama pequeña',serie.length===61 &&
        serie.every((f,i)=>[f.Lambda,f.pi_max,f.a_rel,f.m_diag,f.residuo].every(finito) &&
          cerca(Math.log10(f.Lambda),-3+i/4,1e-9) && f.pi_max>0 && f.pi_max<0.25 &&
          f.a_rel>0 && f.a_rel<1 && f.m_diag>0 && f.m_diag<1 && f.residuo>=0 && f.residuo<=5e-13));
    }
    cotejar(disc,'umbral de la conjetura sobre actividad finita',
      R`Para n\\ge\s*@, que la rama pequeña exista para todo \\lambda>0 es una conjetura`,[Math.min(...na)]);
    const lista = na.map((_,i)=>i===0?'n\\s*=\\s*@':i===na.length-1?' y @':', @').join('');
    cotejar(disc,'evidencia en las 13 series de 61 actividades',
      R`la encuentra en las @ actividades ensayadas, espaciadas logarítmicamente entre \\Lambda\s*=\s*\\lambda\(n-1\)\s*=\s*10\^\{@\} y 10\^\{@\}, para `+lista, [61,-3,12,...na]);
    cotejar(disc,'n=4, consecuencia exacta de la restricción y de la suma de fila',
      R`Para n\s*=\s*@, K_@ se reduce al punto \\pi\\equiv\s*@/@, porque cada fila tiene cuatro casillas con p_c\\le\s*@/@ y suma uno; por tanto \(H_\\infty\) no se cumple`,
      [4,4,1,4,1,4]);
    cotejar(disc,'n=5, cálculo del canónico a masa uno',
      R`Para n\s*=\s*@, el cálculo indica \\max_c\\pi\^\{\(@\)\}_c\s*=\s*@/@`,[5,5,1,4]);
    registrar('n=5 · referencia canónica a masa uno',cerca(datos.canonico5,1/4,1e-12));
    formas(disc,'n=4,5, calificación empírica del diagnóstico de barrera',
      /En ambos casos, un optimizador de barrera sugiere que la rama deja de existir para \\Lambda grande, pero no está demostrado/);
    for (const [n,L] of [[4,100],[5,10000]]) {
      const f=datos.barrera?.[n];
      registrar('n='+n+' · actividad finita, diagnóstico de barrera',f && f.Lambda===L && cerca(f.maxPi,1/4,1e-8));
    }
    // La lista blanca está ligada a posiciones y a afirmaciones. Detecta incluso
    // una cifra ya conocida cuando se añade fuera de su declaración documentada.
    for (const nombre of nombres) {
      const sinIndices = textos[nombre]
        .replace(/\\begin\{(teorema|lema|corolario|proposicion|conjetura|definicion|proof)\}[\s\S]*?\\end\{\1\}/g,s=>' '.repeat(s.length))
        .replace(/\\label\{[^}]*\}|\b(?:Teorema|Lema|Corolario|Proposición|Sección|Apéndice|Ecuación)\s+\d+(?:\.\d+)+|\([A-Z]+\d+(?:\.\d+)+\)/gi,s=>' '.repeat(s.length));
      for (const m of sinIndices.matchAll(/\d+\.\d+|(?<![\w\d])\.\d+/g)) {
        if (!cubiertos[nombre].has(m.index))
          registrar(nombre+' · decimal empírico sin asignar', false, m[0]+' en «'+textos[nombre].slice(Math.max(0,m.index-28),m.index+m[0].length+28)+'»');
      }
    }
  } catch(e) {
    registrar('Datos o estructura incompletos',false,e instanceof Error ? e.message : String(e));
  }
  return resultado();
}

// Uso: node comprobar.mjs [directorio_articulo] [ruta_json]
// Las rutas que el usuario pasa son relativas al directorio de trabajo.
async function ejecutarCLI() {
  const {readFile} = await import('node:fs/promises');
  const {resolve,join,dirname} = await import('node:path');
  const {fileURLToPath} = await import('node:url');
  const base = dirname(fileURLToPath(import.meta.url));
  const articulo=process.argv[2] ? resolve(process.argv[2]) : resolve(base,'../articulo');
  const fuente=process.argv[3] ? resolve(process.argv[3]) : resolve(base,'verificacion_datos.json');
  const nombres=['00_resumen.tex','01a_introduccion.tex','01b_introduccion.tex','09_numerico.tex','10_discusion.tex'];
  const contenidos=Object.fromEntries(await Promise.all(nombres.map(async n=>
    [n,await readFile(join(articulo,'secciones',n),'utf8')])));
  const datos=JSON.parse((await readFile(fuente,'utf8')).replace(/^\uFEFF/,''));
  const r=comprobar(contenidos,datos);
  if (r.ok) console.log('Correcto: '+r.comprobaciones.length+' comprobaciones.');
  else {
    console.error(r.errores.join('\n'));
    console.error(r.errores.length+' fallos en '+r.comprobaciones.length+' comprobaciones.');
    process.exitCode=1;
  }
}
if (typeof process!=='undefined' && process.argv?.[1]) {
  const {pathToFileURL}=await import('node:url');
  if (import.meta.url===pathToFileURL(process.argv[1]).href)
    ejecutarCLI().catch(e=>{console.error(e.message);process.exitCode=1;});
}
