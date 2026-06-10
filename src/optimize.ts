/**
 * Small numeric solvers standing in for scipy.optimize.root / minimize.
 * Gear geometry only ever solves 1-3 dimensional, smooth problems.
 */

export interface SolveResult {
  x: number[];
  fun: number;
  success: boolean;
}

const EPS_J = 1e-7;

function solveLinear(A: number[][], b: number[]): number[] | null {
  // Gaussian elimination with partial pivoting; returns null if singular.
  const n = b.length;
  const M = A.map((row, i) => [...row, b[i]]);
  for (let col = 0; col < n; col++) {
    let piv = col;
    for (let r = col + 1; r < n; r++) {
      if (Math.abs(M[r][col]) > Math.abs(M[piv][col])) piv = r;
    }
    if (Math.abs(M[piv][col]) < 1e-14) return null;
    [M[col], M[piv]] = [M[piv], M[col]];
    for (let r = 0; r < n; r++) {
      if (r === col) continue;
      const f = M[r][col] / M[col][col];
      for (let c = col; c <= n; c++) M[r][c] -= f * M[col][c];
    }
  }
  return M.map((row, i) => row[n] / M[i][i]);
}

/**
 * Find a root of a vector function f: R^n -> R^m using Levenberg-Marquardt
 * (minimizing ||f||^2 until it vanishes). Handles rank-deficient Jacobians,
 * including unused variables (mirrors scipy.optimize.root 'hybr' usage here).
 */
export function rootFind(
  f: (x: number[]) => number[],
  x0: number[],
  { tol = 1e-12, maxIter = 200 }: { tol?: number; maxIter?: number } = {},
): SolveResult {
  let x = [...x0];
  const n = x.length;
  let fx = f(x);
  let cost = fx.reduce((acc, v) => acc + v * v, 0);
  let lambda = 1e-3;

  for (let iter = 0; iter < maxIter; iter++) {
    if (Math.sqrt(cost) < tol) break;
    const m = fx.length;
    // numeric Jacobian m x n
    const J: number[][] = [];
    for (let i = 0; i < m; i++) J.push(new Array(n).fill(0));
    for (let j = 0; j < n; j++) {
      const xp = [...x];
      const h = EPS_J * Math.max(1, Math.abs(x[j]));
      xp[j] += h;
      const fp = f(xp);
      for (let i = 0; i < m; i++) J[i][j] = (fp[i] - fx[i]) / h;
    }
    // normal equations: (J^T J + lambda diag) dx = -J^T f
    const JtJ: number[][] = [];
    const Jtf: number[] = new Array(n).fill(0);
    for (let i = 0; i < n; i++) {
      JtJ.push(new Array(n).fill(0));
      for (let j = 0; j < n; j++) {
        let s = 0;
        for (let k = 0; k < m; k++) s += J[k][i] * J[k][j];
        JtJ[i][j] = s;
      }
      let s = 0;
      for (let k = 0; k < m; k++) s += J[k][i] * fx[k];
      Jtf[i] = s;
    }
    const diagScale = Math.max(...JtJ.map((row, i) => row[i]), 1e-12);
    let improved = false;
    for (let attempt = 0; attempt < 25; attempt++) {
      const A = JtJ.map((row, i) =>
        row.map((v, j) => (i === j ? v + lambda * diagScale : v)),
      );
      const dx = solveLinear(A, Jtf.map((v) => -v));
      if (dx) {
        const xNew = x.map((v, i) => v + dx[i]);
        const fNew = f(xNew);
        const costNew = fNew.reduce((acc, v) => acc + v * v, 0);
        if (costNew < cost) {
          x = xNew;
          fx = fNew;
          cost = costNew;
          lambda = Math.max(lambda * 0.3, 1e-12);
          improved = true;
          break;
        }
      }
      lambda *= 10;
      if (lambda > 1e12) break;
    }
    if (!improved) break;
  }
  return { x, fun: Math.sqrt(cost), success: Math.sqrt(cost) < 1e-8 };
}

function bisect(f: (t: number) => number, a: number, b: number): number {
  let fa = f(a);
  for (let i = 0; i < 100; i++) {
    const m = (a + b) / 2;
    const fm = f(m);
    if (fm === 0 || (b - a) / 2 < 1e-15 * Math.max(1, Math.abs(m))) return m;
    if ((fa < 0) === (fm < 0)) {
      a = m;
      fa = fm;
    } else {
      b = m;
    }
  }
  return (a + b) / 2;
}

/**
 * Scalar root find. Tries Levenberg-Marquardt first; if it stalls at a local
 * minimum, scans an expanding interval around the guess for the nearest sign
 * change and bisects (scipy's hybr wanders much further than LM does).
 */
export function rootFind1D(
  f: (t: number) => number,
  guess: number,
  opts?: { tol?: number; maxIter?: number },
): SolveResult {
  const lm = rootFind((x) => [f(x[0])], [guess], opts);
  if (lm.success) return lm;

  // expanding ring scan for the sign change nearest to the guess
  const f0 = f(guess);
  if (f0 === 0) return { x: [guess], fun: 0, success: true };
  let bestRoot: number | null = null;
  for (let radius = 0.5; radius <= 2048 && bestRoot === null; radius *= 2) {
    const lo = radius === 0.5 ? guess - radius : guess - radius;
    const hi = guess + radius;
    const innerLo = radius === 0.5 ? guess : guess - radius / 2;
    const innerHi = radius === 0.5 ? guess : guess + radius / 2;
    const N = 128;
    const candidates: number[] = [];
    // scan only the newly-added shells (plus full interval on first pass)
    const segments: [number, number][] =
      radius === 0.5 ? [[lo, hi]] : [[lo, innerLo], [innerHi, hi]];
    for (const [a, b] of segments) {
      let prevT = a;
      let prevF = f(a);
      for (let i = 1; i <= N; i++) {
        const t = a + ((b - a) * i) / N;
        const ft = f(t);
        if ((prevF < 0) !== (ft < 0)) {
          candidates.push(bisect(f, prevT, t));
        }
        prevT = t;
        prevF = ft;
      }
    }
    if (candidates.length > 0) {
      bestRoot = candidates.reduce((best, c) =>
        Math.abs(c - guess) < Math.abs(best - guess) ? c : best,
      );
    }
  }
  if (bestRoot !== null) {
    return { x: [bestRoot], fun: Math.abs(f(bestRoot)), success: true };
  }
  return lm;
}

/**
 * Derivative-free minimization (Nelder-Mead) of f: R^n -> R.
 * Mirrors the scipy.optimize.minimize fallback paths in py_gearworks.
 */
export function minimizeND(
  f: (x: number[]) => number,
  x0: number[],
  { tol = 1e-14, maxIter = 600, initialStep = 0.05 }: {
    tol?: number; maxIter?: number; initialStep?: number;
  } = {},
): SolveResult {
  const n = x0.length;
  // initial simplex
  let simplex: { x: number[]; f: number }[] = [{ x: [...x0], f: f(x0) }];
  for (let i = 0; i < n; i++) {
    const x = [...x0];
    x[i] += x[i] !== 0 ? initialStep * Math.abs(x[i]) + initialStep : initialStep;
    simplex.push({ x, f: f(x) });
  }
  const alpha = 1, gamma = 2, rho = 0.5, sigma = 0.5;

  for (let iter = 0; iter < maxIter; iter++) {
    simplex.sort((a, b) => a.f - b.f);
    const best = simplex[0];
    const worst = simplex[n];
    if (Math.abs(worst.f - best.f) < tol * (Math.abs(best.f) + tol)) break;

    const centroid = new Array(n).fill(0);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) centroid[j] += simplex[i].x[j] / n;
    }
    const xr = centroid.map((c, j) => c + alpha * (c - worst.x[j]));
    const fr = f(xr);
    if (fr < best.f) {
      const xe = centroid.map((c, j) => c + gamma * (xr[j] - c));
      const fe = f(xe);
      simplex[n] = fe < fr ? { x: xe, f: fe } : { x: xr, f: fr };
    } else if (fr < simplex[n - 1].f) {
      simplex[n] = { x: xr, f: fr };
    } else {
      const xc = centroid.map((c, j) => c + rho * (worst.x[j] - c));
      const fc = f(xc);
      if (fc < worst.f) {
        simplex[n] = { x: xc, f: fc };
      } else {
        simplex = simplex.map((s, i) =>
          i === 0
            ? s
            : (() => {
                const x = best.x.map((b, j) => b + sigma * (s.x[j] - b));
                return { x, f: f(x) };
              })(),
        );
      }
    }
  }
  simplex.sort((a, b) => a.f - b.f);
  return { x: simplex[0].x, fun: simplex[0].f, success: true };
}
