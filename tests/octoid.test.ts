import { describe, it, expect } from 'vitest';
import { octoid } from '../src/index.js';
import { loadJson, dist3 } from './helpers.js';

interface OctoidCase {
  base_rad: number;
  sphere_rad: number;
  alpha: number;
  angle: number;
  t: number[];
  points: number[][];
}

const cases = loadJson<OctoidCase[]>('octoid_fixtures.json');

describe('octoid (spherical involute) curve', () => {
  for (const c of cases) {
    it(`matches py_gearworks for base_rad=${c.base_rad}, sphere_rad=${c.sphere_rad}, alpha=${c.alpha.toFixed(3)}, angle=${c.angle}`, () => {
      for (let i = 0; i < c.t.length; i++) {
        const p = octoid(c.t[i], {
          baseRad: c.base_rad,
          sphereRad: c.sphere_rad,
          alpha: c.alpha,
          angle: c.angle,
        });
        expect(dist3(p, c.points[i]), `t=${c.t[i]}`).toBeLessThan(1e-9);
      }
    });
  }
});
