
// Variables duales identificadas por las simetrías del cuadrado.
// Cada órbita acumula los residuos de todas sus líneas.
function curva(n) {
 const nr=Math.ceil(n/2),dim=nr+n, nc=n*n, mult=new Float64Array(dim);
 const cells=[];
 for(let i=0;i<n;i++) mult[Math.min(i,n-1-i)]+=2;
 for(let d=0;d<n;d++)mult[nr+d]=d===0?2:4;
 for(let i=0;i<n;i++)for(let j=0;j<n;j++){
   const all=[Math.min(i,n-1-i),Math.min(j,n-1-j),nr+Math.abs(i-j),nr+Math.abs(i+j-(n-1))];
   cells.push(all);
 }
 function small(a){if(!(a>=0&&a<27/256))return NaN;
   let v=a;for(let i=0;i<30;i++){const y=1-v;const step=(v*y*y*y-a)/(y*y*(1-4*v));v-=step;if(Math.abs(step)<3e-16)return v;}let lo=0,hi=.25;for(let j=0;j<60;j++){const mid=(lo+hi)/2;if(mid*(1-mid)**3<a)lo=mid;else hi=mid;}return (lo+hi)/2;
 }
 function ev(z,lam,jac){
   const occ=new Float64Array(dim),vac=Float64Array.from(z,(s,i)=>Math.exp(s)/(i<nr?Math.sqrt(lam):1));
   if(vac.some(t=>!Number.isFinite(t)))return null;
   const exp=Float64Array.from(z,Math.exp),M=jac?Array.from({length:dim},()=>new Float64Array(dim)):null;
   if(jac)for(let i=0;i<dim;i++)M[i][i]=mult[i]*vac[i];
   let mx=0,ma=0;
   for(const ids of cells){const [a,b,c,d]=ids,t=exp[a]*exp[b]*exp[c]*exp[d],v=small(t);if(!Number.isFinite(v))return null;
     mx=Math.max(mx,v);ma=Math.max(ma,t/(27/256));
     occ[a]+=v;occ[b]+=v;occ[c]+=v;occ[d]+=v;
     if(jac){const w=v*(1-v)/(1-4*v);for(const i of ids){const row=M[i];for(const j of ids)row[j]+=w;}}
   }
   let residual=0,md=0;const F=new Float64Array(dim);
   for(let i=0;i<dim;i++){F[i]=occ[i]+mult[i]*(vac[i]-1);residual=Math.max(residual,Math.abs(F[i]/mult[i]));if(i>=nr)md=Math.max(md,occ[i]/mult[i]);}
   return {F,M,residual,mx,ma,md};
 }
 function solve(M,F){const k=F.length,b=Float64Array.from(F,v=>-v);
   for(let i=0;i<k;i++){for(let j=0;j<=i;j++){let a=M[i][j];for(let l=0;l<j;l++)a-=M[i][l]*M[j][l];if(i===j&&a<=0)throw Error('no positiva');M[i][j]=i===j?Math.sqrt(a):a/M[j][j];}
     for(let j=0;j<i;j++)b[i]-=M[i][j]*b[j];b[i]/=M[i][i];}
   for(let i=k-1;i>=0;i--){for(let j=i+1;j<k;j++)b[i]-=M[j][i]*b[j];b[i]/=M[i][i];}return b;
 }
 let z=new Float64Array(dim),rows=[];
 for(let t=0;t<=60;t++){let L=10**(-3+t/4),lam=L/(n-1);
   if(t===0)for(let i=0;i<nr;i++)z[i]=Math.log(lam)/2;
   let out;for(let it=0;it<80;it++){out=ev(z,lam,true);if(!out)throw Error('rama en '+n+' '+t);if(out.residual<5e-13)break;
     let dx=solve(out.M,out.F),pass=false;
     for(let l=0;l<35;l++){const zn=Float64Array.from(z,(v,j)=>v+dx[j]*2**(-l));const a=ev(zn,lam,false);if(a&&a.residual<out.residual){z=zn;pass=true;break;}}
     if(!pass)throw Error('descenso '+n+' '+t+' '+out.residual);
   }
   if(!out||!(out.residual<5e-13&&out.mx<.25))throw Error('residuo');
   rows.push({n,Lambda:L,pi_max:out.mx,a_rel:out.ma,m_diag:out.md,residuo:out.residual});
 }
 return rows;
}

export {curva};
