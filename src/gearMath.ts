/**
 * Port of py_gearworks.gearmath mesh-placement helpers (bevel branch).
 */
import { OUT, RIGHT, Vec3 } from './defs.js';
import { ConicDataParams, conicGamma, conicSphericalRadius } from './gearTeeth.js';
import { GearTransform, Transform } from './transform.js';
import {
  fromRotvec, Mat3, mclone, mmul, mtranspose, mvmul,
  vadd, vcross, vdot, vnorm, vnormalize, vscale, vsub,
} from './vec.js';

export function angleOfVectorInXy(v: Vec3): number {
  return Math.atan2(v[1], v[0]);
}

/** Python-style modulo into [0, 1). */
function mod1(x: number): number {
  return ((x % 1) + 1) % 1;
}

/** Cone parameters bound to a gear transform (py_gearworks ConicData). */
export interface ConicData extends ConicDataParams {
  transform: Transform;
}

/** Spherical center (tip) of the cone, in final (transformed) space. */
export function conicCenter(c: ConicData): Vec3 {
  return c.transform.apply(vscale(OUT, c.baseRadius / Math.tan(conicGamma(c))));
}

/** Port of py_gearworks.gearmath.calc_bevel_gear_placement_vector. */
export function calcBevelGearPlacementVector(
  targetDirNorm: Vec3,
  coneData1: ConicData,
  coneData2: ConicData,
  insideRing1 = false,
  insideRing2 = false,
  offset = 0,
): Vec3 {
  const gamma1 = conicGamma(coneData1);
  const gamma2 = conicGamma(coneData2);
  const R = conicSphericalRadius(coneData1) * Math.abs(coneData1.transform.scale);

  let angleRef: number;
  if (insideRing2) {
    angleRef = gamma2 - gamma1 + offset / R;
  } else if (insideRing1) {
    angleRef = gamma2 - gamma1 - offset / R;
  } else {
    angleRef = gamma1 + gamma2 + offset / R;
  }
  const z2 = coneData2.transform.zAxis;
  const rotAx = vnormalize(vcross(targetDirNorm, z2));
  const rot1 = fromRotvec(vscale(rotAx, angleRef));
  const centerH1 = R * Math.cos(gamma1);
  const diffVector = mvmul(rot1, vscale(z2, -centerH1));
  return vadd(conicCenter(coneData2), diffVector);
}

/** Port of py_gearworks.gearmath.calc_mesh_orientation. */
export function calcMeshOrientation(
  gear1ConeAngle: number,
  gear2ConeAngle: number,
  R: number,
  gear2: Transform,
  insideRing1 = false,
  insideRing2 = false,
  targetDir: Vec3 = RIGHT,
  offset = 0,
): Mat3 {
  const z2 = gear2.zAxis;
  let targetDirNorm = vsub(targetDir, vscale(z2, vdot(targetDir, z2)));
  const gamma1 = gear1ConeAngle / 2;
  const gamma2 = gear2ConeAngle / 2;
  if (vnorm(targetDirNorm) < 1e-12) {
    targetDirNorm = gear2.xAxis;
  } else {
    targetDirNorm = vnormalize(targetDirNorm);
  }
  if (gear1ConeAngle === 0 && gear2ConeAngle === 0) {
    return mclone(gear2.orientation);
  }
  let angleRef: number;
  if (insideRing2) {
    angleRef = gamma2 - gamma1 + offset / R;
  } else if (insideRing1) {
    angleRef = gamma2 - gamma1 - offset / R;
  } else {
    angleRef = gamma1 + gamma2 + offset / R;
  }
  const rotAx = vnormalize(vcross(targetDirNorm, z2));
  return mmul(fromRotvec(vscale(rotAx, angleRef)), gear2.orientation);
}

/** Port of py_gearworks.gearmath.calc_mesh_angle. */
export function calcMeshAngle(
  gearTransform1: GearTransform,
  gearTransform2: GearTransform,
  pitchAngle1: number,
  pitchAngle2: number,
  gear1InsideRing = false,
  gear2InsideRing = false,
): number {
  let centerDiffDir = vsub(gearTransform2.center, gearTransform1.center);
  if (vnorm(centerDiffDir) < 1e-12) {
    centerDiffDir = gearTransform2.xAxis;
  } else {
    centerDiffDir = vnormalize(centerDiffDir);
  }

  let contactDirGear1: Vec3;
  let contactDirOther: Vec3;
  let phaseOffset: number;
  let phaseSign: number;
  if (gear1InsideRing) {
    contactDirGear1 = centerDiffDir;
    contactDirOther = centerDiffDir;
    phaseOffset = 0;
    phaseSign = 1;
  } else if (gear2InsideRing) {
    contactDirGear1 = vscale(centerDiffDir, -1);
    contactDirOther = vscale(centerDiffDir, -1);
    phaseOffset = 0;
    phaseSign = 1;
  } else {
    contactDirGear1 = centerDiffDir;
    contactDirOther = vscale(centerDiffDir, -1);
    phaseOffset = 0.5;
    phaseSign = -1;
  }

  // v @ orientation (row-vector convention) == orientation^T @ v
  const angleOfOther = angleOfVectorInXy(
    mvmul(mtranspose(gearTransform2.orientation), contactDirOther),
  );
  const targetAngleGear1 = angleOfVectorInXy(
    mvmul(mtranspose(gearTransform1.orientation), contactDirGear1),
  );

  const phaseOfOther = mod1((gearTransform2.angle - angleOfOther) / pitchAngle2);

  return targetAngleGear1 + mod1(phaseSign * phaseOfOther + phaseOffset) * pitchAngle1;
}
