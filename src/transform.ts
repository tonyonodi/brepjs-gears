/**
 * Port of py_gearworks.base_classes transform types.
 */
import { ORIGIN, Vec3 } from './defs.js';
import { fromEuler, IDENTITY3, Mat3, mclone, mmul, mvmul, vadd, vclone, vscale } from './vec.js';

export interface TransformDataLike {
  center: Vec3;
  orientation: Mat3;
  scale: number;
}

/** General 3D transform: rotate (orientation), scale, then translate. */
export class Transform implements TransformDataLike {
  center: Vec3;
  orientation: Mat3;
  scale: number;

  constructor({ center = ORIGIN, orientation = IDENTITY3, scale = 1 }: Partial<TransformDataLike> = {}) {
    this.center = vclone(center);
    this.orientation = mclone(orientation);
    this.scale = scale;
  }

  apply(p: Vec3): Vec3 {
    // points @ orientation.T * scale + center  (row-vector convention)
    return vadd(vscale(mvmul(this.orientation, p), this.scale), this.center);
  }

  get xAxis(): Vec3 {
    return [this.orientation[0][0], this.orientation[1][0], this.orientation[2][0]];
  }

  get zAxis(): Vec3 {
    return [this.orientation[0][2], this.orientation[1][2], this.orientation[2][2]];
  }

  mul(other: Transform): Transform {
    return new Transform({
      center: vadd(this.center, vscale(mvmul(this.orientation, other.center), this.scale)),
      orientation: mmul(this.orientation, other.orientation),
      scale: this.scale * other.scale,
    });
  }

  clone(): Transform {
    return new Transform(this);
  }
}

/** Gear transform: like Transform but with an extra rotation-progress angle about local Z. */
export class GearTransform extends Transform {
  angle: number;

  constructor(
    { center = ORIGIN, orientation = IDENTITY3, scale = 1, angle = 0 }:
      Partial<TransformDataLike & { angle: number }> = {},
  ) {
    super({ center, orientation, scale });
    this.angle = angle;
  }

  override apply(p: Vec3): Vec3 {
    // points @ rotZ(angle).T @ orientation.T * scale + center
    const rotated = mvmul(fromEuler('z', this.angle), p);
    return vadd(vscale(mvmul(this.orientation, rotated), this.scale), this.center);
  }

  override mul(other: Transform): GearTransform {
    const otherAngle = other instanceof GearTransform ? other.angle : 0;
    return new GearTransform({
      center: vadd(
        this.center,
        vscale(
          mvmul(mmul(this.orientation, fromEuler('z', this.angle)), other.center),
          this.scale,
        ),
      ),
      orientation: mmul(this.orientation, other.orientation),
      scale: this.scale * other.scale,
      angle: this.angle + otherAngle,
    });
  }

  override clone(): GearTransform {
    return new GearTransform(this);
  }
}
