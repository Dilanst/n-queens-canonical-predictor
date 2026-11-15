
// Matrices y cuadratura en memoria. Puntos medios de las casillas, pesos 1/m.
// La paridad de i-j y de i+j+1 es opuesta, incluso en m impar.
function perfil(n) {
  let r=new Float64Array(n).fill(1),u=new Float64Array(2*n-1).fill(0.5),b=4;
  const rows=new Float64Array(n),ds=new Float64Array(2*n-1),old=new Float64Array(n);
  let res=0,it=0;
  for(;it<4000;it++){
    rows.fill(0);ds.fill(0);
    for(let i=0;i<n;i++)for(let j=0;j<n;j++){
      const d=i-j+n-1, a=i+j, v=b*r[i]*r[j]*u[d]*u[a]/n;
      rows[i]+=v;ds[d]+=v;
    }
    res=Math.max(...rows.map(v=>Math.abs(v-1)),...ds.map((v,i)=>Math.abs(v+u[i]-1)));
    if(res<2e-13)break;
    for(let d=0;d<u.length;d++)u[d]/=u[d]+ds[d];
    old.set(r);let mean=0;
    for(let i=0;i<n;i++){let k=0;for(let j=0;j<n;j++)k+=old[j]*u[i-j+n-1]*u[i+j];r[i]=Math.sqrt(old[i]*n/(b*k));mean+=r[i]/n;}
    for(let i=0;i<n;i++)r[i]/=mean;
    b*=mean*mean;
  }
  if(it===4000)throw Error('Sin convergencia en perfil');
  let i2=0,chk=0; const cols=new Float64Array(n),as=new Float64Array(2*n-1);
  for(let i=0;i<n;i++)for(let j=0;j<n;j++){let v=b*r[i]*r[j]*u[i-j+n-1]*u[i+j];i2+=v*v/(n*n);cols[j]+=v/n;as[i+j]+=v/n;}
  for(let j=0;j<n;j++)chk=Math.max(chk,Math.abs(cols[j]-1));
  for(let a=0;a<as.length;a++)chk=Math.max(chk,Math.abs(as[a]+u[a]-1));
  return {n,r,u,b,it,res:Math.max(res,chk),i2};
}
function interpolar(v,index){const j=Math.min(Math.max(0,Math.floor(index)),v.length-2);return v[j]+(index-j)*(v[j+1]-v[j]);}
function densidad(p,x,y){const n=p.n;return p.b*interpolar(p.r,n*x-.5)*interpolar(p.r,n*y-.5)*interpolar(p.u,n*(x-y)+n-1)*interpolar(p.u,n*(x+y)-1);}
function nystrom(m,p) {
 const l=6*m-2, A=Array.from({length:l},()=>new Float64Array(l));
 for(let i=0;i<l;i++)A[i][i]=1;
 for(let i=0;i<2*m;i++)for(let j=0;j<=i;j++)A[i][j]+=((i<m)===(j<m)?1:-1)/(2*m);
 let s=0, gauge2=0;
 for(let i=0;i<m;i++)for(let j=0;j<m;j++){
   const d=i-j, a=i+j+1;if((Math.abs(d)%2)===a%2)throw Error('Paridad');
   const ids=[i,m+j,2*m+d+m-1,4*m-1+a-1],w=densidad(p,(i+.5)/m,(j+.5)/m)/m;
   for(let f=1;f<4;f++)for(let g=0;g<f;g++)A[ids[f]][ids[g]]+=w;
   s+=w*w;
 }
 let minPivot=Infinity,logdet=0,trace=0;
 for(let i=0;i<l;i++)trace+=A[i][i]-1;
 for(let i=0;i<l;i++){let val=0;for(let j=0;j<2*m;j++){const z=(i>=j?A[i][j]:A[j][i]);val+=(j<m?1:-1)*z;}gauge2+=(val-(i<2*m?(i<m?1:-1):0))**2/(2*m);}
 for(let i=0;i<l;i++){
   const x=A[i];for(let j=0;j<=i;j++){const y=A[j];let s=x[j];for(let t=0;t<j;t++)s-=x[t]*y[t];
     if(i===j){if(!(s>0))throw Error('Pivote '+i);x[j]=Math.sqrt(s);minPivot=Math.min(minPivot,s);logdet+=Math.log(s);}
     else x[j]=s/y[j];}
 }
 let eig=0,vec=Float64Array.from({length:l},(_,i)=>Math.cos(0.087*i)+Math.sin(0.039*i)+0.25);
 for(let z=0;z<36;z++){let norm=Math.hypot(...vec);for(let i=0;i<l;i++)vec[i]/=norm;
   let old=vec.slice();for(let i=0;i<l;i++){const x=A[i];let s=vec[i];for(let j=0;j<i;j++)s-=x[j]*vec[j];vec[i]=s/x[i];}
   for(let i=l-1;i>=0;i--){let s=vec[i];for(let j=i+1;j<l;j++)s-=A[j][i]*vec[j];vec[i]=s/A[i][i];}
   let ip=0;for(let i=0;i<l;i++)ip+=old[i]*vec[i];eig=1/ip;
 }
 const ld=logdet-trace, c=Math.log(2)/2-.5-1.5*p.i2-.5*ld;
 return {m,perfil:p.n,iter:p.it,res:p.res,I2:p.i2,I2_oper:s,logdet2:ld,det2:Math.exp(ld),C:c,minPivot,eig_est:eig,trace,defecto_modo:Math.sqrt(gauge2)};
}

export {perfil,nystrom};
