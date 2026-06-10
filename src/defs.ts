export const PI = Math.PI;

/** Numerical differentiation global 'small step' (matches py_gearworks defs.DELTA). */
export const DELTA = 1e-6;
export const DEG2RAD = PI / 180;
export const RAD2DEG = 180 / PI;

export type Vec3 = [number, number, number];

export const ORIGIN: Vec3 = [0, 0, 0];
export const RIGHT: Vec3 = [1, 0, 0];
export const LEFT: Vec3 = [-1, 0, 0];
export const UP: Vec3 = [0, 1, 0];
export const DOWN: Vec3 = [0, -1, 0];
export const OUT: Vec3 = [0, 0, 1];
export const IN: Vec3 = [0, 0, -1];
