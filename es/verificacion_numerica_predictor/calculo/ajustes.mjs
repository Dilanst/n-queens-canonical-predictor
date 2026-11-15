
export function logfactorial(n) {
 let s=0;for(let i=2;i<=n;i++)s+=Math.log(i);
 return s;
}
export function ajustar(rows,lo,hi,degree,alpha) {
 const a=rows.filter(r=>r.n>=lo&&r.n<=hi).map(({n,logQpred})=>{
  const powers=Array.from({length:degree+1},(_,i)=>(100/n)**i);
  return [...powers,logQpred-logfactorial(n)-n*(1-alpha)];
 });
 const d=degree+1;if(a.length<d)throw Error("Faltan observaciones");
 for(let j=0;j<d;j++){
  const v=a.slice(j).map(r=>r[j]),norm=Math.hypot(...v);
  if(norm===0)throw Error("Diseño singular");
  v[0]+=(v[0]<0?-1:1)*norm;
  const den=v.reduce((s,x)=>s+x*x,0);
  for(let k=j;k<=d;k++){
   let dot=0;for(let t=0;t<v.length;t++)dot+=v[t]*a[t+j][k];dot*=2/den;
   for(let t=0;t<v.length;t++)a[t+j][k]-=dot*v[t];
  }
 }
 const coef=Array(d).fill(0);
 for(let i=d-1;i>=0;i--){let x=a[i][d];for(let j=i+1;j<d;j++)x-=a[i][j]*coef[j];coef[i]=x/a[i][i];}
 return coef[0];
}
export function todosAjustes(rows) {
 return [[40,130],[40,240],[60,240],[20,240]].flatMap(([lo,hi])=>[2,3].map(grado=>({
  lo,hi,grado,min:ajustar(rows,lo,hi,grado,1.944000752),central:ajustar(rows,lo,hi,grado,1.944001),max:ajustar(rows,lo,hi,grado,1.944001082)})));
}
export function extrapolar(a,b,key) {return b[key]+(b[key]-a[key])/((b.m/a.m)**2-1);}
