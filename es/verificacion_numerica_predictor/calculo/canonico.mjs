function calculateSquare(n) {
  if (!Number.isSafeInteger(n) || n < 10) throw RangeError("Se requiere un entero n ≥ 10");
  const families = Math.ceil(n / 2), groups = families + n, lines = 6 * n - 2;
  const rows = n, down = 2 * n, up = 4 * n - 1, middle = n - 1;
  const count = n * n, four = 4 * count;
  const local = new Int32Array(four), global = new Int32Array(four);
  const sizes = new Int32Array(groups);
  const diagonalGroup = k => families + n - 1 - Math.abs(k - middle);
  for (let i = 0; i < n; ++i) sizes[Math.min(i,n-1-i)] += 2;
  for (let k = 0; k < 2*n-1; ++k) sizes[diagonalGroup(k)] += 2;
  for (let r = 0; r < n; ++r) {
    for (let c = 0; c < n; ++c) {
      const k = 4*(r*n+c), d = r-c+middle, a = r+c;
      local[k] = Math.min(r,n-1-r); local[k+1] = Math.min(c,n-1-c);
      local[k+2] = diagonalGroup(d); local[k+3] = diagonalGroup(a);
      global[k] = r; global[k+1] = rows+c; global[k+2] = down+d; global[k+3] = up+a;
    }
  }

  function factor(matrix, dimension) {
    for (let i = 0; i < dimension; ++i) {
      const row = i * dimension;
      for (let j = 0; j <= i; ++j) {
        const prev = j * dimension;
        let s = matrix[row+j];
        for (let k = 0; k < j; ++k) s -= matrix[row+k]*matrix[prev+k];
        if (i === j) {
          if (!(s > 0)) throw Error("Matriz sin pivote positivo");
          matrix[row+i] = Math.sqrt(s);
        } else {
          matrix[row+j] = s/matrix[prev+j];
        }
      }
    }
    return matrix;
  }

  function solveLowerProduct(matrix, rhs) {
    const x = new Float64Array(rhs), len = x.length;
    for (let i = 0; i < len; ++i) {
      const row = i*len;
      let v = x[i];
      for (let j = 0; j < i; ++j) v -= matrix[row+j]*x[j];
      x[i] = v/matrix[row+i];
    }
    for (let i = len-1; i >= 0; --i) {
      const row = i*len;
      x[i] /= matrix[row+i];
      for (let j = 0; j < i; ++j) x[j] -= matrix[row+j]*x[i];
    }
    return x;
  }

  function evaluate(z, curvature) {
    const residual = new Float64Array(groups);
    const matrix = curvature ? new Float64Array(groups*groups) : null;
    const p = new Float64Array(count), h = new Float64Array(n);
    for (let g = 0; g < groups; ++g) {
      const hole = g < families ? 0 : Math.exp(z[g]);
      if (g >= families) h[g-families] = hole;
      residual[g] = sizes[g]*(hole-1);
      if (matrix) matrix[g*groups+g] = sizes[g]*hole;
    }
    for (let k = 0; k < four; k += 4) {
      const b = local[k], c = local[k+1], d = local[k+2], a = local[k+3];
      const activity = Math.exp(z[b]+z[c]+z[d]+z[a]);
      if (!(activity < 27/256) || activity === 0) return null;
      let q = activity;
      for (let it = 0; it < 40; ++it) {
        const s = 1-q;
        const step = (activity-q*s*s*s)/(s*s*(1-4*q));
        q += step;
        if (Math.abs(step) <= 3e-16*q) break;
      }
      p[k/4] = q;
      residual[b] += q; residual[c] += q; residual[d] += q; residual[a] += q;
      if (matrix) {
        const t = q*(1-q)/(1-4*q);
        for (let i = k; i < k+4; ++i) {
          const offset = local[i]*groups;
          for (let j = k; j < k+4; ++j) matrix[offset+local[j]] += t;
        }
      }
    }
    let norm = 0;
    for (let i = 0; i < groups; ++i) norm = Math.max(norm,Math.abs(residual[i])/sizes[i]);
    return {residual, norm, matrix, p, h};
  }

  let z = new Float64Array(groups);
  z.fill(-Math.log(n)/2,0,families);
  let solution = evaluate(z,true), iterations = 0;
  while (solution.norm > 3e-14 && iterations < 60) {
    const delta = solveLowerProduct(factor(solution.matrix,groups),solution.residual);
    let accepted = false;
    for (let a = 1; a >= 2**-30; a /= 2) {
      const trial = z.map((v,i) => v - a*delta[i]);
      const proposed = evaluate(trial,false);
      if (proposed && proposed.norm < solution.norm) {
        z = trial; accepted = true; break;
      }
    }
    if (!accepted) throw Error("El paso de Newton no reduce el residuo");
    ++iterations;
    solution = evaluate(z,true);
  }
  if (solution.norm > 3e-14) throw Error("Newton excedió el límite de iteraciones");

  const {p,h} = solution, s = new Float64Array(lines*lines), individual = new Float64Array(lines);
  const t = p.map(q => q*(1-q)/(1-4*q));
  for (let k = 0; k < 2*n-1; ++k) {
    const hole = h[diagonalGroup(k)-families];
    s[(down+k)*(lines+1)] = hole; s[(up+k)*(lines+1)] = hole;
    individual[down+k] = hole; individual[up+k] = hole;
  }
  let entropy = 0, cellLog = 0, pMin = 1, pMax = 0, hMin = 1, stationarity = 0;
  let corrEntropy = 0, corrCell = 0;
  for (let k = 0; k < four; k += 4) {
    const q = p[k/4], weight = t[k/4], empty = Math.log1p(-q);
    pMin = Math.min(pMin,q); pMax = Math.max(pMax,q);
    const residual = Math.log(q)+3*empty -(z[local[k]]+z[local[k+1]]+z[local[k+2]]+z[local[k+3]]);
    stationarity = Math.max(stationarity,Math.abs(residual));
    const dH = -q*Math.log(q)+3*(1-q)*empty-corrEntropy, nextH = entropy+dH;
    corrEntropy = (nextH-entropy)-dH; entropy = nextH;
    const dC = Math.log1p(-4*q)-4*empty-corrCell, nextC = cellLog+dC;
    corrCell = (nextC-cellLog)-dC; cellLog = nextC;
    for (let i = k; i < k+4; ++i) {
      const row = global[i]*lines;
      individual[global[i]] += q;
      for (let j = k; j < k+4; ++j) s[row+global[j]] += weight;
    }
  }
  for (let i = 0; i < n; ++i) {
    const hole = h[i], g = i+families;
    entropy -= sizes[g]*hole*Math.log(hole);
    hMin = Math.min(hMin,hole);
  }
  let maxMassError = 0;
  for (const m of individual) maxMassError = Math.max(maxMassError,Math.abs(m-1));
  const projected = new Float64Array(s);
  for (let i = 0; i < down; ++i) {
    for (let j = 0; j < down; ++j)
      projected[i*lines+j] += (i<n)===(j<n) ? 1/(2*n) : -1/(2*n);
  }
  const chol = factor(projected,lines);
  let logDet = 0;
  for (let i = 0; i < lines; ++i) logDet += 2*Math.log(chol[i*(lines+1)]);
  const logQ = entropy-0.5*(cellLog+logDet)+0.5*Math.log(4*Math.PI*n);

  function project(x) {
    let a = 0;
    for (let i = 0; i < n; ++i) a += x[i]-x[n+i];
    a /= 2*n;
    for (let i = 0; i < n; ++i) { x[i] -= a; x[n+i] += a; }
  }
  function dot(x,y) {
    let sum = 0;
    for (let i = 0; i < x.length; ++i) sum += x[i]*y[i];
    return sum;
  }
  function multiply(x) {
    const y = new Float64Array(lines);
    for (let k = 0; k < 2*n-1; ++k) {
      const hole = h[diagonalGroup(k)-families];
      y[down+k] = hole*x[down+k]; y[up+k] = hole*x[up+k];
    }
    for (let k = 0; k < four; k += 4) {
      const b=global[k], c=global[k+1], d=global[k+2], a=global[k+3];
      const u=t[k/4]*(x[b]+x[c]+x[d]+x[a]);
      y[b]+=u; y[c]+=u; y[d]+=u; y[a]+=u;
    }
    return y;
  }
  let seed = 19260817;
  let x = new Float64Array(lines);
  for (let i = 0; i < lines; ++i) {
    seed = (Math.imul(1664525,seed)+1013904223) >>> 0;
    x[i] = seed/4294967296-0.5;
  }
  project(x);
  let gap = 0, gapResidual = Infinity, inverseIterations = 0;
  while (gapResidual > 1e-11 && inverseIterations < 1500) {
    x = solveLowerProduct(chol,x);
    project(x);
    const norm = Math.sqrt(dot(x,x));
    for (let i = 0; i < lines; ++i) x[i] /= norm;
    const y = multiply(x);
    gap = dot(x,y);
    for (let i = 0; i < lines; ++i) y[i] -= gap*x[i];
    gapResidual = Math.sqrt(dot(y,y));
    ++inverseIterations;
  }
  if (gapResidual > 1e-11) throw Error("La iteración inversa excedió el límite");
  // Si esta matriz es definida positiva, toda la parte no nula de S supera el extremo inferior.
  const gapLower = gap - 1e-8;
  for (let i = 0; i < down; ++i) {
    for (let j = 0; j <= i; ++j) s[i*lines+j] += (i<n)===(j<n) ? 1/(2*n) : -1/(2*n);
  }
  for (let i = 0; i < lines; ++i) s[i*(lines+1)] -= gapLower;
  factor(s,lines);
  const gapUpper = gap + gapResidual;
  return {n, logQ, entropy, cellLog, logDet, iterations, gapLower, gapUpper, pMin, pMax, branchMargin:0.25-pMax, hMin, gap, gapResidual, maxMassError, stationarity, inverseIterations};
}
export { calculateSquare };
