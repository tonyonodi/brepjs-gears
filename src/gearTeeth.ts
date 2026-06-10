/**
 * Port of py_gearworks.gearteeth (octoid bevel teeth only) and the conic
 * tooth-generator base behavior from base_classes.py.
 *
 * Only the cone_angle != 0 code paths are ported — this library targets bevel gears.
 */
import { ORIGIN, OUT, PI, RIGHT, UP, Vec3 } from './defs.js';
import { angleBetweenVectorAndPlane } from './funcGen.js';
import {
  ArcCurve, Curve, CurveChain, findCurveIntersect, findCurvePlaneIntersect,
  IntersectMethod, OctoidCurve, SphericalInvoluteCurve, FunctionCurve,
} from './curve.js';
import { fromEuler, fromRotvec, mmul, mvmul, vadd, vcross, vdot, vnorm, vnormalize, vscale, vsub } from './vec.js';
import { rootFind1D } from './optimize.js';

export const LABEL_INVOLUTE_FLANK = 'involute_flank';

export interface ToothLimitParam {
  hA: number;
  hD: number;
  hO: number;
}

export interface ConicDataParams {
  coneAngle: number;
  baseRadius: number;
}

export function conicGamma(c: ConicDataParams): number {
  return c.coneAngle / 2;
}

export function conicHeight(c: ConicDataParams): number {
  return c.baseRadius / Math.tan(conicGamma(c));
}

/** Spherical radius |r / sin(gamma)|, always positive. */
export function conicSphericalRadius(c: ConicDataParams): number {
  return Math.abs(c.baseRadius / Math.sin(conicGamma(c)));
}

export interface OctoidToothParams {
  pitchIntersectAngle: number;
  pitchRadius: number;
  coneAngle: number;
  pressureAngle: number;
  refLimits: ToothLimitParam;
  pitchAngle: number;
}

/** Port of GearToothConicGenerator.check_lower_curve_limit (conic branch). */
export function checkLowerCurveLimit(
  toothCurve: CurveChain | Curve,
  params: OctoidToothParams,
): CurveChain | Curve {
  if (params.coneAngle === 0) throw new Error('bevel-only port: coneAngle must be != 0');
  const cone: ConicDataParams = {
    coneAngle: params.coneAngle,
    baseRadius: params.pitchRadius,
  };
  const gamma = conicGamma(cone);
  const R = conicSphericalRadius(cone);
  const hDAngle = PI / 2 - gamma + params.refLimits.hD / R;
  const zCenter = R * Math.cos(gamma);
  const planeNormal = mvmul(fromEuler('y', hDAngle), OUT);
  let sol = findCurvePlaneIntersect(toothCurve, {
    planeNormal,
    offset: vscale(OUT, zCenter),
    guess: 0,
  });
  if (!sol.success) {
    sol = findCurvePlaneIntersect(toothCurve, {
      planeNormal,
      offset: vscale(OUT, zCenter),
      guess: 0,
      method: IntersectMethod.MINDISTANCE,
    });
    toothCurve.setStartOn(sol.x[0]);
    const bottomAngle = Math.abs(
      angleBetweenVectorAndPlane(toothCurve.at(0), OUT) - hDAngle,
    );
    const connector = ArcCurve.fromPointCenterAngle(
      toothCurve.at(0),
      vscale(OUT, R * Math.cos(gamma)),
      bottomAngle,
      UP,
    );
    connector.reverse();
    return new CurveChain(connector, toothCurve);
  }
  toothCurve.setStartOn(sol.x[0]);
  return toothCurve;
}

/** Port of generate_flat_rack_curve (conic branch). */
export function generateFlatRackCurve(
  pitchRadius: number,
  pitchIntersectAngle: number,
  refLimits: ToothLimitParam,
  pressureAngle: number,
  coneAngle: number,
): ArcCurve {
  if (coneAngle === 0) throw new Error('bevel-only port: coneAngle must be != 0');
  const cone: ConicDataParams = { coneAngle, baseRadius: pitchRadius };
  const r = cone.baseRadius;
  const R = conicSphericalRadius(cone);

  const rotAlpha = fromEuler('x', PI / 2 - pressureAngle);
  const rotPitch = fromEuler('z', (-pitchIntersectAngle * r) / R);
  const rotAll = mmul(rotPitch, rotAlpha);
  const arc = new ArcCurve({
    radius: R,
    center: ORIGIN,
    angle: (refLimits.hA + refLimits.hD) / R / Math.cos(pressureAngle),
    yaw: -refLimits.hD / R / Math.cos(pressureAngle),
  });
  arc.rotmat = mmul(rotAll, arc.rotmat);
  return arc;
}

/** Port of generate_undercut_curve (conic branch). */
export function generateUndercutCurve(
  pitchRadius: number,
  coneAngle: number,
  undercutRefPoint: Vec3,
): SphericalInvoluteCurve {
  if (coneAngle === 0) throw new Error('bevel-only port: coneAngle must be != 0');
  const gamma = coneAngle / 2;
  const R = Math.abs(pitchRadius / Math.sin(gamma));
  const cSph = 1 / R;
  const rotated = mvmul(fromEuler('y', PI / 2), vsub(undercutRefPoint, vscale(RIGHT, R)));
  const vOffs: Vec3 = [rotated[0], rotated[1], rotated[2] * Math.sign(gamma)];
  const undercut = new SphericalInvoluteCurve({
    r: pitchRadius,
    cSphere: cSph * Math.sign(gamma),
    vOffs,
    t0: 0,
    t1: -1,
  });
  const sol = rootFind1D(
    (t) => {
      const p = undercut.at(t);
      return vdot(p, p) - pitchRadius * pitchRadius;
    },
    1,
  );
  if (sol.x[0] > 0) undercut.setEndOn(sol.x[0]);
  return undercut;
}

/**
 * Port of trim_involute_undercut.
 *
 * The tooth flank and the undercut trochoid can intersect more than once
 * (the trochoid loop crosses the flank twice). py_gearworks' solver ends up
 * on the intersection that keeps the undercut longest — the last crossing of
 * the trochoid along the flank, which is the correct envelope boundary. We
 * find all intersections via a multi-start search and pick that one
 * deterministically.
 */
export function trimInvoluteUndercut(
  toothCurve: CurveChain | Curve,
  undercutCurve: Curve,
): CurveChain {
  const minToothS = toothCurve instanceof CurveChain
    ? toothCurve.getLengthPortions()[1]
    : -0.5;

  const candidates: { st: number; su: number }[] = [];
  const NG = 12;
  for (let i = 0; i < NG; i++) {
    for (let j = 0; j < NG; j++) {
      const g1 = 0.05 + (1.15 * i) / (NG - 1);
      const g2 = 0.05 + (1.15 * j) / (NG - 1);
      const sol = findCurveIntersect(toothCurve, undercutCurve, [g1, g2]);
      const d = vnorm(vsub(toothCurve.at(sol.x[0]), undercutCurve.at(sol.x[1])));
      if (d > 1e-9) continue;
      const st = sol.x[0];
      const su = sol.x[1];
      // reject extrapolated pseudo-intersections outside either curve's domain
      if (st < minToothS || st > 1.05 || su < -0.001 || su > 1.001) continue;
      if (candidates.some((c) => Math.abs(c.st - st) < 1e-3 && Math.abs(c.su - su) < 1e-3)) {
        continue;
      }
      candidates.push({ st, su });
    }
  }

  let sol: { st: number; su: number };
  if (candidates.length > 0) {
    sol = candidates.reduce((a, b) => (b.su > a.su ? b : a));
  } else {
    // fall back to the single-shot solve like the original
    let res = findCurveIntersect(toothCurve, undercutCurve, [1, 1]);
    if (!res.success) {
      res = findCurveIntersect(toothCurve, undercutCurve, [1, 1], IntersectMethod.MINDISTANCE);
    }
    sol = { st: res.x[0], su: res.x[1] };
  }

  toothCurve.setStartOn(sol.st);
  undercutCurve.setEndOn(sol.su);
  return new CurveChain(undercutCurve, toothCurve);
}

/** Port of OctoidTooth / OctoidUndercutTooth generate_octoid_curve (conic branch). */
export function generateOctoidCurve(params: OctoidToothParams): OctoidCurve {
  const { pitchRadius: rp, pressureAngle: alpha, coneAngle } = params;
  if (coneAngle === 0) throw new Error('bevel-only port: coneAngle must be != 0');
  const gamma = coneAngle / 2;
  const R = Math.abs(rp / Math.sin(gamma));
  const cSph = 1 / R;

  const octoidCurve = new OctoidCurve({ r: rp, cSphere: cSph, angle: 0, alpha });
  octoidCurve.angle = -params.pitchIntersectAngle;
  octoidCurve.updateLengths();
  const sol = findCurvePlaneIntersect(octoidCurve, {
    offset: ORIGIN,
    planeNormal: UP,
    guess: 1,
  });
  octoidCurve.setEndOn(sol.x[0]);
  octoidCurve.label = LABEL_INVOLUTE_FLANK;
  return octoidCurve;
}

/** Port of OctoidTooth (no undercut). */
export class OctoidTooth {
  params: OctoidToothParams;

  constructor(params: OctoidToothParams) {
    this.params = params;
  }

  generateToothCurve(): CurveChain | Curve {
    return checkLowerCurveLimit(generateOctoidCurve(this.params), this.params);
  }
}

/** Port of OctoidUndercutTooth. */
export class OctoidUndercutTooth {
  params: OctoidToothParams;

  constructor(params: OctoidToothParams) {
    this.params = params;
  }

  generateToothCurve(): CurveChain | Curve {
    const p = this.params;
    const toothCurve = checkLowerCurveLimit(generateOctoidCurve(p), p);
    const undercutRefPoint = this.getDefaultUndercutRefPoint();

    const baseAngle = Math.asin(Math.sin(p.coneAngle / 2) * Math.cos(p.pressureAngle));
    const cone: ConicDataParams = { coneAngle: p.coneAngle, baseRadius: p.pitchRadius };
    const undercutRefAngle = p.refLimits.hD / conicSphericalRadius(cone);
    const deltaAngle = p.coneAngle / 2 - baseAngle;
    const undercutEnable = undercutRefAngle > deltaAngle;

    if (!undercutEnable) return toothCurve;

    const undercut = generateUndercutCurve(p.pitchRadius, p.coneAngle, undercutRefPoint);
    return trimInvoluteUndercut(toothCurve, undercut);
  }

  getDefaultUndercutRefPoint(): Vec3 {
    const p = this.params;
    const rackCurve = generateFlatRackCurve(
      p.pitchRadius, p.pitchIntersectAngle, p.refLimits, p.pressureAngle, p.coneAngle,
    );
    // pitch angle is on the gear; equivalent rack rotation scales by r/R = sin(gamma)
    const planeNormal = mvmul(
      fromEuler('z', (-p.pitchAngle / 2) * Math.abs(Math.sin(p.coneAngle / 2))),
      UP,
    );
    const sol = findCurvePlaneIntersect(rackCurve, { planeNormal, offset: ORIGIN, guess: 0 });
    if (sol.x[0] > 0) return rackCurve.at(sol.x[0]);
    return rackCurve.at(0);
  }
}

export type ToothGenerator = OctoidTooth | OctoidUndercutTooth;
