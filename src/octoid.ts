import { Vec3, PI } from './defs.js';
import { octoidFn, octoidContactFn } from './funcGen.js';

export interface OctoidParams {
  baseRad: number;
  sphereRad: number;
  alpha: number;
  angle?: number;
}

/** Spherical involute ("octoid") flank curve for bevel gears. */
export function octoid(t: number, params: OctoidParams): Vec3 {
  const { baseRad, sphereRad, alpha = (20 * PI) / 180, angle = 0 } = params;
  return octoidFn(t, baseRad, sphereRad, alpha, angle);
}

/** Contact curve of the octoid generation. */
export function octoidContact(t: number, params: OctoidParams): Vec3 {
  const { baseRad, sphereRad, alpha = (20 * PI) / 180, angle = 0 } = params;
  return octoidContactFn(t, baseRad, sphereRad, alpha, angle);
}
