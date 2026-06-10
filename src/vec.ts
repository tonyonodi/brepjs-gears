import { Vec3 } from './defs.js';

export type Mat3 = [Vec3, Vec3, Vec3]; // row-major

export const IDENTITY3: Mat3 = [
  [1, 0, 0],
  [0, 1, 0],
  [0, 0, 1],
];

export function vadd(a: Vec3, b: Vec3): Vec3 {
  return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
}

export function vsub(a: Vec3, b: Vec3): Vec3 {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}

export function vscale(a: Vec3, s: number): Vec3 {
  return [a[0] * s, a[1] * s, a[2] * s];
}

export function vmul(a: Vec3, b: Vec3): Vec3 {
  return [a[0] * b[0], a[1] * b[1], a[2] * b[2]];
}

export function vdot(a: Vec3, b: Vec3): number {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

export function vcross(a: Vec3, b: Vec3): Vec3 {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];
}

export function vnorm(a: Vec3): number {
  return Math.hypot(a[0], a[1], a[2]);
}

export function vnormalize(a: Vec3): Vec3 {
  const n = vnorm(a);
  return [a[0] / n, a[1] / n, a[2] / n];
}

export function vclone(a: Vec3): Vec3 {
  return [a[0], a[1], a[2]];
}

/** Matrix * vector (column-vector convention). */
export function mvmul(m: Mat3, v: Vec3): Vec3 {
  return [
    m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
    m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
    m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
  ];
}

export function mmul(a: Mat3, b: Mat3): Mat3 {
  const out: number[][] = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
  ];
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      out[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
    }
  }
  return out as Mat3;
}

export function mtranspose(m: Mat3): Mat3 {
  return [
    [m[0][0], m[1][0], m[2][0]],
    [m[0][1], m[1][1], m[2][1]],
    [m[0][2], m[1][2], m[2][2]],
  ];
}

export function mclone(m: Mat3): Mat3 {
  return [vclone(m[0]), vclone(m[1]), vclone(m[2])];
}

/** Rotation matrix about X axis. */
export function rotX(a: number): Mat3 {
  const c = Math.cos(a);
  const s = Math.sin(a);
  return [
    [1, 0, 0],
    [0, c, -s],
    [0, s, c],
  ];
}

/** Rotation matrix about Y axis. */
export function rotY(a: number): Mat3 {
  const c = Math.cos(a);
  const s = Math.sin(a);
  return [
    [c, 0, s],
    [0, 1, 0],
    [-s, 0, c],
  ];
}

/** Rotation matrix about Z axis. */
export function rotZ(a: number): Mat3 {
  const c = Math.cos(a);
  const s = Math.sin(a);
  return [
    [c, -s, 0],
    [s, c, 0],
    [0, 0, 1],
  ];
}

const AXIS_ROT: Record<string, (a: number) => Mat3> = { x: rotX, y: rotY, z: rotZ };

/**
 * Equivalent of scipy Rotation.from_euler with lowercase (extrinsic) axes.
 * Rotations are applied in sequence about fixed axes:
 * R = R_n @ ... @ R_2 @ R_1 (first axis in the string is applied first).
 */
export function fromEuler(seq: string, angles: number | number[]): Mat3 {
  const arr = typeof angles === 'number' ? [angles] : angles;
  if (seq.length !== arr.length) throw new Error('euler seq/angles mismatch');
  let m = IDENTITY3;
  for (let i = 0; i < seq.length; i++) {
    m = mmul(AXIS_ROT[seq[i]](arr[i]), m);
  }
  return m;
}

/** Rotation matrix from a rotation vector (axis * angle). */
export function fromRotvec(rv: Vec3): Mat3 {
  const angle = vnorm(rv);
  if (angle < 1e-300) return IDENTITY3;
  const [x, y, z] = vscale(rv, 1 / angle);
  const c = Math.cos(angle);
  const s = Math.sin(angle);
  const C = 1 - c;
  return [
    [c + x * x * C, x * y * C - z * s, x * z * C + y * s],
    [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
    [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
  ];
}

/** Rotate vector v around axis by angle (axis assumed normalized). */
export function rotateAround(v: Vec3, axis: Vec3, angle: number): Vec3 {
  return mvmul(fromRotvec(vscale(axis, angle)), v);
}
