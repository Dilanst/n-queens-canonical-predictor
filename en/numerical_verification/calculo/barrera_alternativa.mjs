function finiteBarrier(n) {
  const cap=0.25, width=6*n-2, groupIds=new Map(), orbits=[];
  const lineCells=Array.from({length:width},()=>[]);
  for(let r=0;r<n;++r) for(let c=0;c<n;++c) {
    const a=Math.min(r,n-1-r), b=Math.min(c,n-1-c);
    const key=Math.min(a,b)+','+Math.max(a,b);
    if(!groupIds.has(key)) {groupIds.set(key,orbits.length);orbits.push({s:0,lines:new Float64Array(width)});}
    const group=groupIds.get(key); orbits[group].s++;
    for(const l of [r,n+c,2*n+r-c+n-1,4*n-1+r+c]) {
      orbits[group].lines[l]++; lineCells[l].push(group);
    }
  }
  const length=orbits.length;
  function state(gap, logLambda, mu) {
    const p=gap.map(v=>cap-v);
    if(p.some((v,i)=>!(v>0&&gap[i]>0))) return null;
    const h=lineCells.map(cells=>cells.reduce((val,g)=>val-p[g],1));
    if(h.some(v=>!(v>0))) return null;
    const gradient=new Float64Array(length), C=Array.from({length},()=>new Float64Array(length));
    let obj=0;
    for(let g=0;g<length;++g) {
      const q=p[g], u=gap[g], size=orbits[g].s;
      obj+=size*(logLambda*q-q*Math.log(q)+3*(1-q)*Math.log1p(-q)+mu*Math.log(u));
      gradient[g]=size*(logLambda-Math.log(q)-3*Math.log1p(-q)-mu/u);
      C[g][g]=size*((1-4*q)/(q*(1-q))+mu/(u*u));
    }
    for(let l=0;l<width;++l) {
      const hole=h[l], lh=Math.log(hole);
      obj-=hole*lh;
      for(let i=0;i<length;++i) {
        const d=orbits[i].lines[l];
        gradient[i]+=d*lh;
        for(let j=0;j<=i;++j) C[i][j]+=d*orbits[j].lines[l]/hole;
      }
    }
    let norm=0; for(let g=0;g<length;++g) norm=Math.max(norm,Math.abs(gradient[g])/orbits[g].s);
    return {gradient,C,obj,norm,p,h};
  }
  function step(C,gradient) {
    const s=new Float64Array(gradient);
    for(let i=0;i<length;++i) {
      for(let j=0;j<=i;++j) {
        let v=C[i][j]; for(let k=0;k<j;++k) v-=C[i][k]*C[j][k];
        C[i][j]=i===j ? Math.sqrt(v) : v/C[j][j];
      }
      for(let j=0;j<i;++j) s[i]-=C[i][j]*s[j];
      s[i]/=C[i][i];
    }
    for(let i=length-1;i>=0;--i) {
      for(let j=i+1;j<length;++j) s[i]-=C[j][i]*s[j];
      s[i]/=C[i][i];
    }
    return s;
  }
  function optimize(gap,logLambda,mu) {
    let s=state(gap,logLambda,mu);
    let iter=0;
    for(;iter<150 && s.norm>5e-10;++iter) {
      const delta=step(s.C,s.gradient);
      let a=1, slope=delta.reduce((t,v,i)=>t+v*s.gradient[i],0);
      for(let i=0;i<length;++i) {
        if(delta[i]>0) a=Math.min(a,0.995*gap[i]/delta[i]);
        else if(delta[i]<0) a=Math.min(a,0.995*s.p[i]/(-delta[i]));
      }
      for(let l=0;l<width;++l) {
        const dh=lineCells[l].reduce((v,i)=>v+delta[i],0);
        if(dh>0) a=Math.min(a,0.995*s.h[l]/dh);
      }
      let passed=false;
      for(let i=0;i<70;++i) {
        const next=gap.map((x,g)=>x-a*delta[g]), t=state(next,logLambda,mu);
        if(t && (t.norm<s.norm*(1-0.01*a) || (t.obj>s.obj && t.obj>=s.obj+1e-4*a*slope))) {
          gap=next;s=t;passed=true;break;
        }
        a*=0.5;
      }
      if(!passed) break;
    }
    if(s.norm>2e-8) throw Error(`Residuo ${s.norm}, n=${n} logLambda=${logLambda} mu=${mu}`);
    return {gap,norm:s.norm,iter};
  }
  let warm=new Float64Array(length).fill(cap-0.5/n), worst=0, first=null;
  const rows=[];
  for(let j=-4;j<=36;++j) {
    const k=j/4, Lambda=10**k, logLambda=Math.log(Lambda/(n-1));
    let gap=warm,norm=0, maxIter=0;
    for(let m=2;m<=14;++m) {
      const z=optimize(gap,logLambda,10**(-m));gap=z.gap;norm=z.norm;
      worst=Math.max(worst,norm);maxIter=Math.max(maxIter,z.iter);
      if(m===2) warm=gap;
    }
    const minGap=Math.min(...gap), pMax=cap-minGap;
    const row={k,Lambda,pMax,minGap,norm,maxIter}; rows.push(row);
    if(pMax>=cap-1e-5 && first===null) first=row;
  }
  return {n,first,last:rows[rows.length-1],prev:first?rows[rows.indexOf(first)-1]:null,worst,rows};
}
export { finiteBarrier };
