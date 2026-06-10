import { readFileSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const FIXTURES = join(dirname(fileURLToPath(import.meta.url)), '..', 'fixtures');

export function loadJson<T = any>(name: string): T {
  return JSON.parse(readFileSync(join(FIXTURES, name), 'utf8')) as T;
}

export function loadStepBlob(name: string): Blob {
  const buf = readFileSync(join(FIXTURES, name));
  return new Blob([buf]);
}

let kernelReady: Promise<void> | null = null;

/** Initialize the brepjs occt-wasm kernel once. */
export function initKernel(): Promise<void> {
  if (!kernelReady) {
    kernelReady = (async () => {
      const { OcctKernel } = await import('occt-wasm');
      const { registerKernel, OcctWasmAdapter } = await import('brepjs');
      const kernel = await OcctKernel.init();
      registerKernel('occt-wasm', OcctWasmAdapter.fromKernel(kernel));
    })();
  }
  return kernelReady;
}

export function dist3(a: readonly number[], b: readonly number[]): number {
  return Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

export interface GearFixture {
  name: string;
  params: {
    number_of_teeth: number;
    module: number;
    height: number;
    cone_angle: number;
    helix_angle: number;
  };
  z_vals: number[];
  module: number;
  cone_angle: number;
  pitch_angle: number;
  transform: { center: number[]; angle: number; scale: number; orientation: number[][] };
  sphere_data: { z: number; center: number[]; R: number }[];
  slices: {
    z: number;
    local_transform: { center: number[]; angle: number; scale: number };
    recipe_limits: { h_a: number; h_d: number; h_o: number };
    ra_radius: number;
    rd_radius: number;
    ro_radius: number;
    tooth_curve: { t: number[]; points: number[][] };
    profile_transformed?: { t: number[]; points: number[][] };
  }[];
  n_points_vert: number;
  n_z_tweens: number;
  solid_stats: {
    volume: number;
    area: number;
    bbox_min: number[];
    bbox_max: number[];
    center_of_mass: number[];
  };
}
