"""Diagnóstico de barrera a actividad finita. Python 3 y NumPy."""
import numpy as np
from math import log

def reinas(n):
    """A (n²×(6n−2)), lista de tipos: 0 fila, 1 columna, 2 diagonal, 3 antidiagonal"""
    lines = []; tipo = []
    for i in range(n): lines.append([i * n + j for j in range(n)]); tipo.append(0)
    for j in range(n): lines.append([i * n + j for i in range(n)]); tipo.append(1)
    for k in range(-(n - 1), n): lines.append([i * n + (i - k) for i in range(n) if 0 <= i - k < n]); tipo.append(2)
    for s in range(0, 2 * n - 1): lines.append([i * n + (s - i) for i in range(n) if 0 <= s - i < n]); tipo.append(3)
    A = np.zeros((n * n, len(lines)))
    for Li, L in enumerate(lines): A[L, Li] = 1
    return A, np.array(tipo), lines

def xlogx(x):
    x = np.asarray(x, float); out = np.zeros_like(x); p = x > 0; out[p] = x[p] * np.log(x[p]); return out

def Phi(lam, p, A, d):
    m = A.T @ p; h = 1 - m
    return log(lam) * p.sum() - xlogx(p).sum() + (d - 1) * xlogx(1 - p).sum() - xlogx(h).sum()

def grad_hess(lam, p, A, d, mu=0.0):
    m = A.T @ p; h = 1 - m
    g = log(lam) - np.log(p) - (d - 1) * np.log(1 - p) + A @ np.log(h)
    Hs = np.diag(-1 / p + (d - 1) / (1 - p)) - A @ np.diag(1 / h) @ A.T
    if mu > 0:
        b = 1.0 / d - p; g = g - mu / b; Hs = Hs - np.diag(mu / b ** 2)
    return g, Hs

def factible(p, A, d, estricto_d=True):
    if np.any(p <= 0): return False
    if estricto_d and np.any(p >= 1.0 / d): return False
    if np.any(A.T @ p >= 1): return False
    return True

def maximiza(lam, A, d, p0, mu=0.0, tol=1e-13, maxit=500):
    """Newton amortiguado sobre la función cóncava Φ_λ + μ Σ log(1/d − π) en Ω°."""
    p = p0.copy()
    def obj(q):
        v = Phi(lam, q, A, d)
        if mu > 0: v += mu * np.sum(np.log(1.0 / d - q))
        return v
    f0 = obj(p)
    for it in range(maxit):
        g, Hs = grad_hess(lam, p, A, d, mu)
        dx = np.linalg.solve(Hs, -g)
        dec = -g @ dx                              # decremento de Newton
        if np.max(np.abs(g)) < 1e-12 or abs(dec) < 1e-28: return p, it, np.max(np.abs(g))
        t = 1.0
        while True:
            q = p + t * dx
            if factible(q, A, d, estricto_d=True):
                fq = obj(q)
                if fq >= f0 + 1e-4 * t * (g @ dx) - 1e-15 * abs(f0): break
            t *= 0.5
            if t < 1e-14: return p, it, np.max(np.abs(g))
        p, f0 = q, fq
    return p, maxit, np.max(np.abs(grad_hess(lam, p, A, d, mu)[0]))

def max_barrera(lam, A, d, p0):
    p = p0.copy()
    for mu in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14]:
        p, it, gr = maximiza(lam, A, d, p, mu=mu)
    return p

for n in (4, 5, 6):
    A, tipo, lines = reinas(n); d = 4
    p = np.full(n * n, 1e-3)
    prev = None
    for k in np.arange(-1, 9.01, 0.25):
        Lam = 10 ** k; lam = Lam / (n - 1)
        p = max_barrera(lam, A, d, np.minimum(p, 0.2499))
        mx = p.max()
        # gradiente de Φ (sin barrera) en las casillas con π máximo: > 0 indica restricción activa π = 1/4
        m = A.T @ p; g = log(lam) - np.log(p) - 3 * np.log(1 - p) + A @ np.log(1 - m)
        diagnostico = 'rama pequeña' if mx < 0.25 - 1e-5 else 'π̂ en el borde 1/4'
        if prev != diagnostico: print(f"n={n}: Λ = 10^{k:5.2f} → máx π̂ = {mx:.7f}, ∂Φ en argmax = {g[np.argmax(p)]:+.2e}: {diagnostico}")
        prev = diagnostico
    print(f"n={n}: último Λ = 10^{k:.2f}: máx π̂ = {mx:.7f} ({diagnostico})")
