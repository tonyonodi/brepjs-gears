/**
 * Port of py_gearworks.curve — the curve framework used to build tooth profiles.
 *
 * Curves are callable via `.at(s)` where s in [0,1] is an arc-length-equalized
 * parameter backed by a sampled lookup table (21 samples by default), exactly
 * mirroring the Python implementation so that fixture parametrizations match.
 */
import { DELTA, ORIGIN, OUT, PI, RIGHT, UP, Vec3 } from './defs.js';
import {
  angleBetweenVectors, interpolateLin, involuteSphere, octoidFn,
} from './funcGen.js';
import {
  fromEuler, fromRotvec, Mat3, mmul, mvmul, mtranspose,
  vadd, vcross, vdot, vnorm, vnormalize, vscale, vsub, vclone, mclone,
} from './vec.js';
import { minimizeND, rootFind, rootFind1D, SolveResult } from './optimize.js';

/** np.interp equivalent: linear interpolation with clamping outside range. */
function interp(x: number, xp: number[], fp: number[]): number {
  const n = xp.length;
  if (x <= xp[0]) return fp[0];
  if (x >= xp[n - 1]) return fp[n - 1];
  // binary search for the right interval
  let lo = 0;
  let hi = n - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (xp[mid] <= x) lo = mid;
    else hi = mid;
  }
  if (xp[hi] === xp[lo]) return fp[lo];
  return fp[lo] + ((fp[hi] - fp[lo]) * (x - xp[lo])) / (xp[hi] - xp[lo]);
}

/** np.searchsorted(a, v, side='left'): first index i where a[i] >= v. */
function searchsortedLeft(a: number[], v: number): number {
  let lo = 0;
  let hi = a.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (a[mid] < v) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

export abstract class Curve {
  protected _active = true;
  t0: number;
  t1: number;
  protected _length = 0;
  lenApproxN: number;
  lookupT: number[] = [-1e6, 1e6];
  lookupS: number[] = [-1e6, 1e6];
  label = '';

  constructor(t0 = 0, t1 = 1, active = true, lenApproxN = 21) {
    this.t0 = t0;
    this.t1 = t1;
    this._active = active;
    this.lenApproxN = lenApproxN;
  }

  get active(): boolean {
    return this._active;
  }

  set active(v: boolean) {
    this._active = v;
  }

  get length(): number {
    return this._length;
  }

  set length(v: number) {
    this._length = v;
  }

  /** The underlying curve function in its natural parameter t. */
  abstract evalT(t: number): Vec3;

  /** Deep copy. */
  abstract clone(): this;

  /** Initialize lengths — must be called at end of subclass constructors. */
  protected init(): void {
    if (this.active) this.updateLengths();
  }

  at(s: number): Vec3 {
    return this.evalT(this.s2t(s));
  }

  s2t(s: number): number {
    let ret = interp(s, this.lookupS, this.lookupT);
    if (s < 0) {
      ret = interpolateLin(s, this.lookupS[0], this.lookupS[1], this.lookupT[0], this.lookupT[1]);
    } else if (s > 1) {
      const n = this.lookupS.length;
      ret = interpolateLin(
        s, this.lookupS[n - 1], this.lookupS[n - 2], this.lookupT[n - 1], this.lookupT[n - 2],
      );
    }
    return ret;
  }

  t2s(t: number): number {
    let ret = interp(t, this.lookupT, this.lookupS);
    if (t < this.t0) {
      ret = interpolateLin(t, this.lookupT[0], this.lookupT[1], this.lookupS[0], this.lookupS[1]);
    } else if (t > this.t1) {
      const n = this.lookupT.length;
      ret = interpolateLin(
        t, this.lookupT[n - 1], this.lookupT[n - 2], this.lookupS[n - 1], this.lookupS[n - 2],
      );
    }
    return ret;
  }

  updateLengths(): void {
    const n = this.lenApproxN;
    const tRange: number[] = [];
    for (let i = 0; i < n; i++) tRange.push(this.t0 + ((this.t1 - this.t0) * i) / (n - 1));
    const pts = tRange.map((t) => this.evalT(t));
    const cum: number[] = [0];
    for (let i = 1; i < n; i++) cum.push(cum[i - 1] + vnorm(vsub(pts[i], pts[i - 1])));
    this.length = cum[n - 1];
    if (this.length === 0) {
      this.lookupS = tRange.map((_, i) => i / (n - 1));
    } else {
      this.lookupS = cum.map((c) => c / this.length);
    }
    this.lookupT = tRange;
  }

  reverse(): this {
    const tmp = this.t0;
    this.t0 = this.t1;
    this.t1 = tmp;
    this.lookupT = [...this.lookupT].reverse();
    return this;
  }

  derivative(s: number, direction = 0, delta = DELTA): Vec3 {
    if (direction === 0) {
      return vscale(vsub(this.at(s + delta), this.at(s - delta)), 1 / (2 * delta));
    }
    if (direction > 0) {
      return vscale(vsub(this.at(s + delta), this.at(s)), 1 / delta);
    }
    return vscale(vsub(this.at(s), this.at(s - delta)), 1 / delta);
  }

  setStartAndEndOn(s0: number, s1: number): void {
    const newT0 = this.s2t(s0);
    const newT1 = this.s2t(s1);
    this.t0 = newT0;
    this.t1 = newT1;
    this.updateLengths();
  }

  setStartOn(s0: number): void {
    this.t0 = this.s2t(s0);
    this.updateLengths();
  }

  setEndOn(s1: number): void {
    this.t1 = this.s2t(s1);
    this.updateLengths();
  }

  getCurves(): Curve[] {
    return [this];
  }
}

/** A generic function-backed curve (port of plain Curve(func, ...) usage). */
export class FunctionCurve extends Curve {
  fn: (t: number) => Vec3;

  constructor(fn: (t: number) => Vec3, t0 = 0, t1 = 1, active = true, lenApproxN = 21) {
    super(t0, t1, active, lenApproxN);
    this.fn = fn;
    this.init();
  }

  evalT(t: number): Vec3 {
    return this.fn(t);
  }

  clone(): this {
    const c = new FunctionCurve(this.fn, this.t0, this.t1, this.active, this.lenApproxN);
    c.lookupT = [...this.lookupT];
    c.lookupS = [...this.lookupS];
    c.length = this.length;
    c.label = this.label;
    return c as this;
  }
}

export class LineCurve extends Curve {
  p0: Vec3;
  p1: Vec3;

  constructor(p0: Vec3 = ORIGIN, p1: Vec3 = ORIGIN, active = true) {
    super(0, 1, active);
    this.p0 = vclone(p0);
    this.p1 = vclone(p1);
    this.init();
  }

  evalT(t: number): Vec3 {
    return vadd(vscale(this.p0, 1 - t), vscale(this.p1, t));
  }

  override updateLengths(): void {
    this.length = vnorm(vsub(this.p1, this.p0)) * Math.abs(this.t1 - this.t0);
    this.lookupS = [0, 1];
    this.lookupT = [this.t0, this.t1];
  }

  override setStartAndEndOn(s0: number, s1: number): void {
    const a = this.at(s0);
    const b = this.at(s1);
    this.p0 = a;
    this.p1 = b;
    this.updateLengths();
  }

  override setStartOn(s0: number): void {
    this.p0 = this.at(s0);
    this.updateLengths();
  }

  override setEndOn(s1: number): void {
    this.p1 = this.at(s1);
    this.updateLengths();
  }

  transform(fn: (p: Vec3) => Vec3): LineCurve {
    return new LineCurve(fn(this.p0), fn(this.p1), this.active);
  }

  clone(): this {
    const c = new LineCurve(this.p0, this.p1, this.active);
    c.t0 = this.t0;
    c.t1 = this.t1;
    c.label = this.label;
    c.updateLengths();
    return c as this;
  }
}

export class ArcCurve extends Curve {
  radius: number;
  angle: number;
  center: Vec3;
  rotmat: Mat3;

  constructor(
    {
      radius = 1.0,
      angle = PI / 2,
      center = ORIGIN,
      yaw = 0.0,
      pitch = 0.0,
      roll = 0.0,
      rotmat,
      active = true,
    }: {
      radius?: number; angle?: number; center?: Vec3;
      yaw?: number; pitch?: number; roll?: number; rotmat?: Mat3; active?: boolean;
    } = {},
  ) {
    super(0, 1, active);
    this.radius = radius;
    this.angle = angle;
    this.center = vclone(center);
    // scipy from_euler('zyx', [yaw, pitch, roll]) extrinsic
    this.rotmat = rotmat ? mclone(rotmat) : fromEuler('zyx', [yaw, pitch, roll]);
    this.init();
  }

  evalT(t: number): Vec3 {
    const rotArc = fromEuler('z', this.angle * t);
    return vadd(mvmul(mmul(this.rotmat, rotArc), vscale(RIGHT, this.radius)), this.center);
  }

  override updateLengths(): void {
    this.length = this.radius * this.angle;
    this.lookupS = [0, 1];
    this.lookupT = [this.t0, this.t1];
  }

  get p0(): Vec3 {
    return this.at(0);
  }

  get p1(): Vec3 {
    return this.at(1);
  }

  get axis(): Vec3 {
    return mvmul(this.rotmat, OUT);
  }

  static from2PointCenter(
    p0: Vec3, p1: Vec3, center: Vec3, revolutions = 0, active = true,
  ): ArcCurve {
    const r = vnorm(vsub(p0, center));
    const x = vnormalize(vsub(p0, center));
    const z = vnormalize(vcross(vsub(p0, center), vsub(p1, center)));
    const y = vcross(z, x);
    const rotmat: Mat3 = mtranspose([x, y, z]);
    return new ArcCurve({
      radius: r,
      angle: angleBetweenVectors(vsub(p0, center), vsub(p1, center)) + revolutions * PI * 2,
      center,
      rotmat,
      active,
    });
  }

  static fromPointCenterAngle(
    p0: Vec3, center: Vec3, angle: number, axis: Vec3 = OUT, active = true,
  ): ArcCurve {
    const r = vnorm(vsub(p0, center));
    const x = vnormalize(vsub(p0, center));
    const y = vnormalize(vcross(axis, x));
    const z = vcross(x, y);
    const rotmat: Mat3 = mtranspose([x, y, z]);
    return new ArcCurve({ radius: r, angle, center, rotmat, active });
  }

  transform(fn: (p: Vec3) => Vec3): ArcCurve {
    const p0 = fn(this.at(0));
    const center0 = fn(this.center);
    const center2 = fn(vadd(this.center, this.axis));
    const axis2 = vnormalize(vsub(center2, center0));
    return ArcCurve.fromPointCenterAngle(p0, center0, this.angle, axis2, this.active);
  }

  clone(): this {
    const c = new ArcCurve({
      radius: this.radius,
      angle: this.angle,
      center: this.center,
      rotmat: this.rotmat,
      active: this.active,
    });
    c.t0 = this.t0;
    c.t1 = this.t1;
    c.label = this.label;
    c.updateLengths();
    return c as this;
  }
}

export class SphericalInvoluteCurve extends Curve {
  r: number;
  angle: number;
  cSphere: number;
  vOffs: Vec3;
  zOffs: number;

  constructor(
    {
      r = 1, t0 = 0, t1 = 1, angle = 0, cSphere = 1,
      vOffs = ORIGIN, zOffs = 0, active = true,
    }: {
      r?: number; t0?: number; t1?: number; angle?: number; cSphere?: number;
      vOffs?: Vec3; zOffs?: number; active?: boolean;
    } = {},
  ) {
    super(t0, t1, active);
    this.r = r;
    this.angle = angle;
    this.cSphere = cSphere;
    this.vOffs = vclone(vOffs);
    this.zOffs = zOffs;
    this.init();
  }

  evalT(t: number): Vec3 {
    return involuteSphere(t, this.r, this.cSphere, this.angle, this.vOffs, this.zOffs);
  }

  get R(): number {
    return 1 / this.cSphere;
  }

  get centerSphere(): Vec3 {
    return vscale(
      OUT,
      Math.sqrt(this.R * this.R - this.r * this.r) * Math.sign(this.cSphere) + this.zOffs,
    );
  }

  get baseRadius(): number {
    return this.r;
  }

  clone(): this {
    const c = new SphericalInvoluteCurve({
      r: this.r, t0: this.t0, t1: this.t1, angle: this.angle,
      cSphere: this.cSphere, vOffs: this.vOffs, zOffs: this.zOffs, active: this.active,
    });
    c.lookupT = [...this.lookupT];
    c.lookupS = [...this.lookupS];
    c.length = this.length;
    c.label = this.label;
    return c as this;
  }
}

export class OctoidCurve extends Curve {
  r: number;
  angle: number;
  cSphere: number;
  alpha: number;

  constructor(
    {
      r = 1, t0 = 0, t1 = 1, angle = 0, cSphere = 1,
      alpha = (20 * PI) / 180, active = true,
    }: {
      r?: number; t0?: number; t1?: number; angle?: number;
      cSphere?: number; alpha?: number; active?: boolean;
    } = {},
  ) {
    super(t0, t1, active);
    this.r = r;
    this.angle = angle;
    this.cSphere = cSphere;
    this.alpha = alpha;
    this.init();
  }

  evalT(t: number): Vec3 {
    return octoidFn(t, this.r, 1 / this.cSphere, this.alpha, this.angle);
  }

  get R(): number {
    return 1 / this.cSphere;
  }

  get centerSphere(): Vec3 {
    return vscale(
      OUT,
      Math.sqrt(this.R * this.R - this.r * this.r) * Math.sign(this.cSphere),
    );
  }

  get baseRadius(): number {
    return this.r;
  }

  clone(): this {
    const c = new OctoidCurve({
      r: this.r, t0: this.t0, t1: this.t1, angle: this.angle,
      cSphere: this.cSphere, alpha: this.alpha, active: this.active,
    });
    c.lookupT = [...this.lookupT];
    c.lookupS = [...this.lookupS];
    c.length = this.length;
    c.label = this.label;
    return c as this;
  }
}

export class TransformedCurve extends Curve {
  targetCurve: Curve;
  transformFn: (p: Vec3) => Vec3;

  constructor(transformFn: (p: Vec3) => Vec3, curve: Curve, t0 = 0, t1 = 1) {
    super(t0, t1, curve.active);
    this.targetCurve = curve;
    this.transformFn = transformFn;
    this.init();
  }

  evalT(t: number): Vec3 {
    return this.transformFn(this.targetCurve.at(t));
  }

  clone(): this {
    const c = new TransformedCurve(this.transformFn, this.targetCurve.clone(), this.t0, this.t1);
    c.lookupT = [...this.lookupT];
    c.lookupS = [...this.lookupS];
    c.length = this.length;
    c.label = this.label;
    return c as this;
  }
}

export class MirroredCurve extends TransformedCurve {
  constructor(curve: Curve, planeNormal: Vec3 = RIGHT, center: Vec3 = ORIGIN) {
    const normal = vnormalize(planeNormal);
    const mirrorFn = (p: Vec3): Vec3 => {
      const p2 = vsub(p, center);
      const h = vdot(p2, normal);
      return vadd(vsub(p2, vscale(normal, 2 * h)), center);
    };
    super(mirrorFn, curve);
  }
}

export class RotatedCurve extends TransformedCurve {
  constructor(curve: Curve, angle = 0, axis: Vec3 = OUT, center: Vec3 = ORIGIN) {
    const naxis = vnormalize(axis);
    const rotFn = (p: Vec3): Vec3 => {
      const p2 = vsub(p, center);
      return vadd(mvmul(fromRotvec(vscale(naxis, angle)), p2), center);
    };
    super(rotFn, curve);
  }
}

export class CurveChain extends Curve {
  curves: Curve[];
  private _chainActive = true;

  constructor(...curves: Curve[]) {
    super(0, 1, true);
    this.curves = curves;
    this.updateLengths();
  }

  override get active(): boolean {
    return this._chainActive && this.curves.some((c) => c.active);
  }

  override set active(v: boolean) {
    this._chainActive = v;
  }

  override updateLengths(): void {
    for (const c of this.curves ?? []) {
      if (c.active) c.updateLengths();
    }
  }

  get lengthArray(): number[] {
    return this.curves.map((c) => (c.active ? c.length : 0));
  }

  override get length(): number {
    return this.lengthArray.reduce((a, b) => a + b, 0);
  }

  override set length(_v: number) {
    /* chain length is derived */
  }

  get numCurves(): number {
    return this.curves.length;
  }

  override s2t(s: number): number {
    return s;
  }

  override t2s(t: number): number {
    return t;
  }

  get idxActiveMin(): number {
    const idx = this.curves.findIndex((c) => c.active);
    return idx === -1 ? this.curves.length : idx;
  }

  get idxActiveMax(): number {
    for (let i = this.curves.length - 1; i >= 0; i--) {
      if (this.curves[i].active) return i;
    }
    return -1;
  }

  getLengthPortions(): number[] {
    const total = this.length;
    const portions = [0];
    let acc = 0;
    for (const l of this.lengthArray) {
      acc += l;
      portions.push(acc / total);
    }
    return portions;
  }

  getSIndex(s: number): [number, number] {
    const portions = this.getLengthPortions();
    let idx = searchsortedLeft(portions, s);
    if (idx > this.idxActiveMax + 1) idx = this.idxActiveMax + 1;
    if (idx < this.idxActiveMin + 1) idx = this.idxActiveMin + 1;
    let sIdx: number;
    if (portions[idx] - portions[idx - 1] !== 0) {
      sIdx = (s - portions[idx - 1]) / (portions[idx] - portions[idx - 1]);
    } else {
      sIdx = 0.5;
    }
    return [idx - 1, sIdx];
  }

  getTForIndex(idx: number): [number, number] {
    const portions = this.getLengthPortions();
    return [portions[idx], portions[idx + 1]];
  }

  override evalT(t: number): Vec3 {
    // not used; chain evaluates via at()
    return this.at(t);
  }

  override at(s: number): Vec3 {
    const [idx, s2] = this.getSIndex(s);
    return this.curves[idx].at(s2);
  }

  override getCurves(): Curve[] {
    const out: Curve[] = [];
    for (const c of this.curves) {
      if (c instanceof CurveChain) out.push(...c.getCurves());
      else out.push(c);
    }
    return out;
  }

  append(c: Curve): void {
    this.curves.push(c);
    this.updateLengths();
  }

  insert(index: number, c: Curve): void {
    this.curves.splice(index, 0, c);
    this.updateLengths();
  }

  override reverse(): this {
    this.curves.reverse();
    for (const c of this.curves) c.reverse();
    return this;
  }

  override setStartAndEndOn(s0: number, s1: number): void {
    this.setStartOn(s0);
    this.setEndOn(s1);
  }

  override setStartOn(s0: number): void {
    const [idx, s2] = this.getSIndex(s0);
    this.curves = this.curves.slice(idx);
    this.curves[0].setStartOn(s2);
    this.updateLengths();
  }

  override setEndOn(s1: number): void {
    const [idx, s2] = this.getSIndex(s1);
    this.curves = this.curves.slice(0, idx + 1);
    this.curves[this.curves.length - 1].setEndOn(s2);
    this.updateLengths();
  }

  clone(): this {
    const c = new CurveChain(...this.curves.map((cv) => cv.clone()));
    c._chainActive = this._chainActive;
    c.label = this.label;
    return c as this;
  }
}

export enum IntersectMethod {
  EQUALITY = 1,
  MINDISTANCE = 2,
}

export function findCurveIntersect(
  curve1: Curve,
  curve2: Curve,
  guess: [number, number] = [0.5, 0.5],
  method: IntersectMethod = IntersectMethod.EQUALITY,
): SolveResult {
  if (method === IntersectMethod.EQUALITY) {
    return rootFind(
      (x) => {
        const d = vsub(curve1.at(x[0]), curve2.at(x[1]));
        return [d[0], d[1], d[2]];
      },
      [guess[0], guess[1], 0],
    );
  }
  return minimizeND(
    (x) => {
      const d = vscale(vsub(curve1.at(x[0]), curve2.at(x[1])), 1 / DELTA);
      return vdot(d, d);
    },
    [guess[0], guess[1]],
  );
}

export function findCurvePlaneIntersect(
  curve: Curve,
  {
    planeNormal = OUT,
    offset = ORIGIN,
    guess = 0,
    method = IntersectMethod.EQUALITY,
  }: {
    planeNormal?: Vec3; offset?: Vec3; guess?: number; method?: IntersectMethod;
  } = {},
): SolveResult {
  const target = (t: number) => vdot(vsub(curve.at(t), offset), planeNormal);
  if (method === IntersectMethod.EQUALITY) {
    return rootFind1D(target, guess);
  }
  return minimizeND((x) => target(x[0]) ** 2, [guess]);
}

export function findCurveNearestPoint(
  curve: Curve, point: Vec3, guesses: number[] = [0.5],
): number {
  let best: SolveResult | null = null;
  for (const g of guesses) {
    const res = minimizeND(
      (x) => {
        const d = vsub(curve.at(x[0]), point);
        return vdot(d, d);
      },
      [g],
    );
    if (!best || res.fun < best.fun) best = res;
  }
  return best!.x[0];
}

export interface TangentArcResult {
  arc: ArcCurve;
  t1: number;
  t2: number;
  sol: SolveResult;
}

export function calcTangentArc(
  curve1: Curve,
  curve2: Curve,
  radius: number,
  startLocations: [number, number] = [1, 0],
  method: IntersectMethod = IntersectMethod.EQUALITY,
): TangentArcResult {
  const calcCenters = (t1: number, t2: number): [Vec3, Vec3] => {
    const p1 = curve1.at(t1);
    const p2 = curve2.at(t2);
    const tan1 = vnormalize(curve1.derivative(t1));
    const tan2 = vnormalize(curve2.derivative(t2));
    let arcAxis = vcross(tan1, tan2);
    const angle = vnorm(arcAxis);
    arcAxis = vscale(arcAxis, 1 / angle);
    const normal1 = vcross(tan1, arcAxis);
    const normal2 = vcross(tan2, arcAxis);
    return [vsub(p1, vscale(normal1, radius)), vsub(p2, vscale(normal2, radius))];
  };

  let sol: SolveResult;
  if (method === IntersectMethod.EQUALITY) {
    sol = rootFind(
      (x) => {
        const [c1, c2] = calcCenters(startLocations[0] - x[0], startLocations[1] + x[1]);
        const d = vsub(c1, c2);
        return [d[0], d[1], d[2]];
      },
      [0, 0, 0],
    );
  } else {
    sol = minimizeND(
      (x) => {
        const [c1, c2] = calcCenters(startLocations[0] - x[0], startLocations[1] + x[1]);
        const d = vsub(c1, c2);
        return vdot(d, d);
      },
      [0, 0],
    );
  }

  const t1 = startLocations[0] - sol.x[0];
  const t2 = startLocations[1] + sol.x[1];
  const [c1, c2] = calcCenters(t1, t2);
  const center = vscale(vadd(c1, c2), 0.5);
  const arc = ArcCurve.from2PointCenter(curve1.at(t1), curve2.at(t2), center);
  return { arc, t1, t2, sol };
}
