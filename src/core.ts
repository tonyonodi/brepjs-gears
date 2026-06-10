/**
 * Port of py_gearworks.core — reference circles, profile trimming/assembly,
 * tip reduction, fillets, and the Gear manager (conic/bevel paths only).
 */
import { DELTA, ORIGIN, OUT, PI, RIGHT, UP, Vec3 } from './defs.js';
import {
  cylindricalToXyz, rotateVectorZ, sphericalToXyz, xyzToCylindrical, xyzToSpherical,
} from './funcGen.js';
import {
  ArcCurve, calcTangentArc, Curve, CurveChain, findCurveIntersect,
  findCurvePlaneIntersect, IntersectMethod, MirroredCurve, RotatedCurve,
} from './curve.js';
import { minimizeND } from './optimize.js';
import { vadd, vcross, vdot, vnorm, vnormalize, vscale, vsub } from './vec.js';
import { GearTransform, Transform } from './transform.js';
import {
  ConicDataParams, conicGamma, conicSphericalRadius, OctoidToothParams,
  OctoidTooth, OctoidUndercutTooth, ToothGenerator, ToothLimitParam,
} from './gearTeeth.js';

export interface FilletParam {
  tipFillet: number;
  rootFillet: number;
  tipReduction: number;
}

/** Polar transform used during profile generation (conic branch only). */
export class GearPolarTransform {
  coneAngle: number;
  baseRadius: number;

  constructor(coneAngle: number, baseRadius: number) {
    if (coneAngle === 0) throw new Error('bevel-only port: coneAngle must be != 0');
    this.coneAngle = coneAngle;
    this.baseRadius = baseRadius;
  }

  get gamma(): number {
    return this.coneAngle / 2;
  }

  /** Spherical radius, always positive. */
  get R(): number {
    return Math.abs(this.baseRadius / Math.sin(this.gamma));
  }

  /** Spherical center (tip side) of the cone, untransformed. */
  get center(): Vec3 {
    return vscale(OUT, this.baseRadius / Math.tan(this.gamma));
  }

  polarTransform(point: Vec3): Vec3 {
    const sign = Math.sign(this.coneAngle);
    const zComp: Vec3 = [1, 1, sign];
    const p: Vec3 = [point[0] * zComp[0], point[1] * zComp[1], point[2] * zComp[2]];
    const c: Vec3 = [
      this.center[0] * zComp[0], this.center[1] * zComp[1], this.center[2] * zComp[2],
    ];
    const sph = xyzToSpherical(p, c);
    return [(PI / 2 - sph[2]) * this.R, sph[1], this.R - sph[0]];
  }

  inversePolarTransform(point: Vec3): Vec3 {
    const sign = Math.sign(this.coneAngle);
    const zComp: Vec3 = [1, 1, sign];
    const c: Vec3 = [
      this.center[0] * zComp[0], this.center[1] * zComp[1], this.center[2] * zComp[2],
    ];
    const sph: Vec3 = [this.R - point[2], point[1], PI / 2 - point[0] / this.R];
    const xyz = sphericalToXyz(sph, c);
    return [xyz[0] * zComp[0], xyz[1] * zComp[1], xyz[2] * zComp[2]];
  }
}

export interface GearRefCircles {
  rACurve: ArcCurve;
  rPCurve: ArcCurve;
  rDCurve: ArcCurve;
  rOCurve: ArcCurve;
}

export function generateReferenceCircles(
  pitchRadius: number,
  limits: ToothLimitParam,
  coneparam: GearPolarTransform,
): GearRefCircles {
  const p0 = vscale(RIGHT, pitchRadius);
  const hO = limits.hO > pitchRadius - DELTA ? pitchRadius - DELTA : limits.hO;

  const pa = coneparam.inversePolarTransform(
    vadd(coneparam.polarTransform(p0), [limits.hA, 0, 0]),
  );
  const pd = coneparam.inversePolarTransform(
    vadd(coneparam.polarTransform(p0), [-limits.hD, 0, 0]),
  );
  const po = coneparam.inversePolarTransform(
    vadd(coneparam.polarTransform(p0), [-hO, 0, 0]),
  );

  const mkCircle = (p: Vec3) =>
    ArcCurve.fromPointCenterAngle(p, vscale(OUT, p[2]), 2 * PI);

  return {
    rPCurve: mkCircle(p0),
    rACurve: mkCircle(pa),
    rDCurve: mkCircle(pd),
    rOCurve: mkCircle(po),
  };
}

/** Port of apply_tip_reduction — returns the (possibly reduced) addendum radius. */
export function applyTipReduction(
  toothCurve: Curve,
  addendumHeight: number,
  dedendumHeight: number,
  tipReduction: number,
  polar: GearPolarTransform,
): number {
  const solData: { x: number[]; rSol: number }[] = [];
  const rah = addendumHeight;
  const rdh = dedendumHeight;
  let rOut = rah;
  for (const guess of [0.1, 0.1 + 0.8 / 3, 0.1 + 1.6 / 3, 0.9]) {
    const sol = findCurvePlaneIntersect(toothCurve, { planeNormal: UP, guess });
    const rSol = polar.polarTransform(toothCurve.at(sol.x[0]))[0];
    if (sol.success && rSol > rdh) {
      solData.push({ x: sol.x, rSol });
    }
  }
  if (solData.length > 0) {
    const best = solData.reduce((a, b) => (b.rSol < a.rSol ? b : a));
    if (best.rSol - tipReduction < rah) {
      rOut = tipReduction > 0 ? best.rSol - tipReduction : best.rSol;
    }
  }
  return rOut;
}

/** Port of apply_fillet (root: direction=1, tip: direction=-1). */
export function applyFillet(
  toothCurve: CurveChain,
  pitchAngle: number,
  targetCircle: ArcCurve,
  filletRadius: number,
  direction: 1 | -1,
): CurveChain {
  const angleCheck = (p: Vec3): boolean => {
    const angle = -Math.atan2(p[1], p[0]);
    return 0 < angle && angle < pitchAngle / 2;
  };

  const refCircleGuess = (-pitchAngle / (2 * PI) / 4) * 1.01;
  const sol1 = findCurveIntersect(toothCurve, targetCircle, [0.5, refCircleGuess]);

  let sharpRoot = false;
  if (sol1.success && angleCheck(targetCircle.at(sol1.x[1]))) {
    let arcRes: ReturnType<typeof calcTangentArc> | null = null;
    const guesses = [0.5, 1, 1.5].map((g) => g * filletRadius);
    if (direction === 1) {
      for (const guess of guesses) {
        const startLocations: [number, number] = [
          sol1.x[1] - guess / targetCircle.length,
          sol1.x[0] + guess / toothCurve.length,
        ];
        arcRes = calcTangentArc(targetCircle, toothCurve, filletRadius, startLocations);
        if (arcRes.sol.success) break;
      }
    } else {
      for (const guess of guesses) {
        const startLocations: [number, number] = [
          sol1.x[0] - guess / toothCurve.length,
          sol1.x[1] + guess / targetCircle.length,
        ];
        arcRes = calcTangentArc(toothCurve, targetCircle, filletRadius, startLocations);
        if (arcRes.sol.success) break;
      }
    }
    const { arc, t1, t2 } = arcRes!;
    if (angleCheck(arc.at(0)) && angleCheck(arc.at(1))) {
      if (direction === 1) {
        toothCurve.setStartOn(t2);
        toothCurve.insert(0, arc);
      } else {
        toothCurve.setEndOn(t1);
        toothCurve.append(arc);
      }
    } else {
      sharpRoot = true;
    }
  } else {
    sharpRoot = true;
  }

  if (sharpRoot) {
    const planeNormal = direction === 1 ? rotateVectorZ(UP, -pitchAngle / 2) : UP;
    const mirrorCurve = new MirroredCurve(toothCurve, planeNormal);
    mirrorCurve.reverse();
    const startLocations: [number, number] = [
      1 - filletRadius / toothCurve.length,
      0 + filletRadius / toothCurve.length,
    ];
    if (direction === 1) {
      const { arc, t2 } = calcTangentArc(mirrorCurve, toothCurve, filletRadius, startLocations);
      arc.setStartOn(0.5);
      toothCurve.setStartOn(t2);
      toothCurve.insert(0, arc);
    } else {
      const { arc, t1 } = calcTangentArc(toothCurve, mirrorCurve, filletRadius, startLocations);
      arc.setEndOn(0.5);
      toothCurve.setEndOn(t1);
      toothCurve.append(arc);
    }
  }

  return toothCurve;
}

export interface GearRefProfile {
  raCurve: ArcCurve;
  rdCurve: ArcCurve;
  roCurve: ArcCurve;
  toothCurve: Curve;
  toothCurveMirror: Curve;
  pitchAngle: number;
  transform: GearTransform;
}

/** profile = rd -> tooth -> ra -> tooth_mirror (one pitch of the gear boundary). */
export function profileChain(p: GearRefProfile): CurveChain {
  return new CurveChain(p.rdCurve, p.toothCurve, p.raCurve, p.toothCurveMirror);
}

export interface GearRefProfileExtended extends GearRefProfile {
  roConnector0: Curve;
  roConnector1: Curve;
  roConnector2: Curve;
  rdConnector: Curve;
  raConnector: Curve;
  roCurveTooth: Curve;
  roCurveDedendum: Curve;
  toothCenterline: Curve;
}

/**
 * profile_closed = rd, tooth, ra, tooth_mirror, ro_connector_2 reversed,
 * ro reversed, ro_connector_0 — a closed loop including the outer ring.
 */
export function profileClosedChain(p: GearRefProfileExtended): CurveChain {
  return new CurveChain(
    p.rdCurve,
    p.toothCurve,
    p.raCurve,
    p.toothCurveMirror,
    p.roConnector2.clone().reverse(),
    p.roCurve.clone().reverse(),
    p.roConnector0,
  );
}

/** Port of trim_reference_profile. */
export function trimReferenceProfile(
  toothCurve: Curve,
  refCurves: GearRefCircles,
  fillet: FilletParam,
  pitchAngle: number,
): GearRefProfile {
  const rACurve = refCurves.rACurve;
  const rDCurve = refCurves.rDCurve;

  if (!(fillet.tipFillet > 0)) {
    const raGuess = -pitchAngle / 8 / rACurve.length;
    let solTip = findCurveIntersect(toothCurve, rACurve, [0.9, raGuess]);
    if (!solTip.success) {
      solTip = findCurveIntersect(
        toothCurve, rACurve, [0.9, raGuess], IntersectMethod.MINDISTANCE,
      );
    }
    const solcheck = vnorm(vsub(toothCurve.at(solTip.x[0]), rACurve.at(solTip.x[1])));
    if ((solTip.success || solcheck < 1e-5) && toothCurve.at(solTip.x[0])[1] < 0) {
      toothCurve.setEndOn(solTip.x[0]);
    } else {
      const solMid = findCurvePlaneIntersect(toothCurve, { planeNormal: UP, guess: 1 });
      toothCurve.setEndOn(solMid.x[0]);
    }
  }

  if (!(fillet.rootFillet > 0)) {
    const rdGuess = -pitchAngle / 2 / rDCurve.length;
    let solRoot = findCurveIntersect(toothCurve, rDCurve, [0.3, rdGuess]);
    let solcheck = vnorm(vsub(toothCurve.at(solRoot.x[0]), rDCurve.at(solRoot.x[1])));
    if (!solRoot.success) {
      const solRoot2 = findCurveIntersect(
        toothCurve, rDCurve, [0, rdGuess], IntersectMethod.MINDISTANCE,
      );
      const solcheck2 = vnorm(vsub(toothCurve.at(solRoot2.x[0]), rDCurve.at(solRoot2.x[1])));
      if (solRoot2.success || solcheck2 < 1e-5) {
        solcheck = solcheck2;
        solRoot = solRoot2;
      }
    }
    const angleCheck = Math.atan2(
      toothCurve.at(solRoot.x[0])[1],
      toothCurve.at(solRoot.x[0])[0],
    );
    if ((solRoot.success || solcheck < 1e-5) && angleCheck > -pitchAngle / 2) {
      toothCurve.setStartOn(solRoot.x[0]);
    } else {
      const planeNorm = rotateVectorZ(UP, -pitchAngle / 2);
      const guess = solRoot.success || solcheck < 1e-5 ? solRoot.x[0] : 0.1;
      let solMid2 = findCurvePlaneIntersect(toothCurve, { planeNormal: planeNorm, guess });
      const solBot = minimizeND(
        (x) => {
          const p = toothCurve.at(x[0]);
          return Math.hypot(p[0], p[1]);
        },
        [solMid2.x[0]],
      );
      if (solBot.x[0] > solMid2.x[0]) solMid2 = solBot;
      toothCurve.setStartOn(solMid2.x[0]);
    }
  }

  const toothMirror = new MirroredCurve(toothCurve, UP);
  toothMirror.reverse();
  const toothRotate = new RotatedCurve(toothMirror, -pitchAngle, OUT);

  const pa1 = toothCurve.at(1);
  const pa2 = toothMirror.at(0);
  const centerA = vscale(OUT, (pa1[2] + pa2[2]) / 2);
  let raCurveOut: ArcCurve;
  if (vnorm(vsub(pa1, pa2)) > 1e-10) {
    raCurveOut = ArcCurve.from2PointCenter(pa1, pa2, centerA);
  } else {
    raCurveOut = new ArcCurve({
      radius: vnorm(vsub(refCurves.rACurve.at(0), refCurves.rACurve.center)),
      angle: 0,
      center: centerA,
      yaw: Math.atan2(pa1[1], pa1[0]),
      active: false,
    });
  }

  const pd1 = toothCurve.at(0);
  const pd2 = toothRotate.at(1);
  const centerD = vscale(OUT, (pd1[2] + pd2[2]) / 2);
  let rdCurveOut: ArcCurve;
  if (vnorm(vsub(pd1, pd2)) > 1e-10) {
    rdCurveOut = ArcCurve.from2PointCenter(pd2, pd1, centerD);
  } else {
    rdCurveOut = new ArcCurve({
      radius: vnorm(vsub(refCurves.rDCurve.at(0), refCurves.rDCurve.center)),
      angle: 0,
      center: centerD,
      yaw: Math.atan2(pd1[1], pd1[0]),
      active: false,
    });
  }

  const profile = new CurveChain(rdCurveOut, toothCurve, raCurveOut, toothMirror);
  const start = profile.at(0);
  const end = profile.at(1);
  const angle0 = Math.atan2(start[1], start[0]);
  const angle1 = Math.atan2(end[1], end[0]);

  const roCurveOut = new ArcCurve({
    radius: refCurves.rOCurve.radius,
    center: refCurves.rOCurve.center,
    angle: angle1 - angle0,
    yaw: angle0,
  });

  return {
    raCurve: raCurveOut,
    rdCurve: rdCurveOut,
    roCurve: roCurveOut,
    toothCurve,
    toothCurveMirror: toothMirror,
    pitchAngle,
    transform: new GearTransform(),
  };
}

export interface ProfileRecipeAtZ {
  toothGenerator: ToothGenerator;
  cone: ConicDataParams;
  limits: ToothLimitParam;
  pitchAngle: number;
  transform: GearTransform;
  fillet: FilletParam;
}

/** Port of generate_reference_profile. */
export function generateReferenceProfile(input: ProfileRecipeAtZ): GearRefProfile {
  const conic = new GearPolarTransform(input.cone.coneAngle, input.cone.baseRadius);
  const refCurves = generateReferenceCircles(
    input.toothGenerator.params.pitchRadius, input.limits, conic,
  );
  let toothCurve = input.toothGenerator.generateToothCurve();

  if (input.fillet.tipReduction > 0) {
    const rAh = applyTipReduction(
      toothCurve,
      conic.polarTransform(refCurves.rACurve.at(0))[0],
      conic.polarTransform(refCurves.rDCurve.at(0))[0],
      input.fillet.tipReduction,
      conic,
    );
    const pa = conic.inversePolarTransform([rAh, 0, 0]);
    refCurves.rACurve = ArcCurve.fromPointCenterAngle(pa, vscale(OUT, pa[2]), 2 * PI);
  }
  if (input.fillet.tipFillet > 0) {
    const chain = toothCurve instanceof CurveChain ? toothCurve : new CurveChain(toothCurve);
    toothCurve = applyFillet(chain, input.pitchAngle, refCurves.rACurve, input.fillet.tipFillet, -1);
  }
  if (input.fillet.rootFillet > 0) {
    const chain = toothCurve instanceof CurveChain ? toothCurve : new CurveChain(toothCurve);
    toothCurve = applyFillet(chain, input.pitchAngle, refCurves.rDCurve, input.fillet.rootFillet, 1);
  }
  const profile = trimReferenceProfile(toothCurve, refCurves, input.fillet, input.pitchAngle);
  profile.transform = input.transform.clone();
  return profile;
}

/** Port of generate_profile_extensions (conic branch). */
export function generateProfileExtensions(
  profile: GearRefProfile,
  coneData: ConicDataParams,
): GearRefProfileExtended {
  if (coneData.coneAngle === 0) throw new Error('bevel-only port: coneAngle must be != 0');
  const coneCenter = vscale(OUT, coneData.baseRadius / Math.tan(conicGamma(coneData)));

  const toothStartPoint = profile.toothCurve.at(0);
  const toothStartPlaneNormal = vcross(OUT, vnormalize(toothStartPoint));
  const solRoMid = findCurvePlaneIntersect(profile.roCurve, {
    planeNormal: toothStartPlaneNormal,
    guess: 0.5,
  });
  const roMidpoint = profile.roCurve.at(solRoMid.x[0]);
  const roCurveTooth = profile.roCurve.clone();
  roCurveTooth.setStartOn(solRoMid.x[0]);
  const roCurveDedendum = profile.roCurve.clone();
  roCurveDedendum.setEndOn(solRoMid.x[0]);

  const rdConnector = ArcCurve.from2PointCenter(
    profile.toothCurve.at(0),
    profile.toothCurveMirror.at(1),
    profile.rdCurve.center,
  );

  const raConnector = ArcCurve.from2PointCenter(
    rotateVectorZ(profile.raCurve.at(1), -profile.pitchAngle),
    profile.raCurve.at(0),
    profile.raCurve.center,
  );

  const profileStart = profileChain(profile).at(0);
  const roConnector0 = ArcCurve.from2PointCenter(
    profile.roCurve.at(0), profileStart, coneCenter,
  );
  const roConnector1 = ArcCurve.from2PointCenter(
    roMidpoint, profile.toothCurve.at(0), coneCenter,
  );
  const roConnector2 = ArcCurve.from2PointCenter(
    profile.roCurve.at(1), profile.toothCurveMirror.at(1), coneCenter,
  );
  const toothCenterline = ArcCurve.from2PointCenter(
    rdConnector.at(0.5), profile.raCurve.at(0.5), coneCenter,
  );

  return {
    ...profile,
    roConnector0,
    roConnector1,
    roConnector2,
    rdConnector,
    raConnector,
    roCurveTooth,
    roCurveDedendum,
    toothCenterline,
  };
}

export interface GearToothParamData {
  numTeeth: number;
  numCutoutTeeth: number;
  insideTeeth: boolean;
}

export function pitchAngleOf(p: GearToothParamData): number {
  return (2 * PI) / p.numTeeth;
}

export function numTeethAct(p: GearToothParamData): number {
  return Math.floor(p.numTeeth - p.numCutoutTeeth);
}

/** Port of the Gear manager class (conic/bevel only). */
export class Gear {
  zVals: number[];
  module: number;
  toothParam: GearToothParamData;
  shapeRecipe: (z: number) => ProfileRecipeAtZ;
  transform: GearTransform;

  constructor(opts: {
    zVals: number[];
    module: number;
    toothParam: GearToothParamData;
    shapeRecipe: (z: number) => ProfileRecipeAtZ;
    transform: GearTransform;
  }) {
    this.zVals = [...opts.zVals];
    this.module = opts.module;
    this.toothParam = opts.toothParam;
    this.shapeRecipe = opts.shapeRecipe;
    this.transform = opts.transform;
  }

  get pitchAngle(): number {
    return pitchAngleOf(this.toothParam);
  }

  curveGenAtZ(z: number): GearRefProfile {
    return generateReferenceProfile(this.shapeRecipe(z));
  }

  /** Cone of the recipe at z with the combined (global * local) transform applied. */
  coneAtZ(z: number): { cone: ConicDataParams; transform: GearTransform } {
    const recipe = this.shapeRecipe(z);
    const trf = this.transform.mul(recipe.transform);
    return { cone: recipe.cone, transform: trf };
  }

  /** Sphere through the profile at z: [center, R] in final space. */
  sphereDataAtZ(z: number): { center: Vec3; R: number } {
    const { cone, transform } = this.coneAtZ(z);
    const centerUntransformed = vscale(
      OUT, cone.baseRadius / Math.tan(conicGamma(cone)),
    );
    // ConicData.center applies the plain (non-angle) transform; the cone is
    // rotationally symmetric about Z so the gear angle does not matter.
    const center = new Transform({
      center: transform.center,
      orientation: transform.orientation,
      scale: transform.scale,
    }).apply(centerUntransformed);
    const R = Math.abs(
      (cone.baseRadius / Math.sin(conicGamma(cone))) * transform.scale,
    );
    return { center, R };
  }
}
