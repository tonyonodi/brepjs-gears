import { describe, it, expect } from 'vitest';
import { BevelGear } from '../src/index.js';
import { loadJson, dist3, GearFixture } from './helpers.js';

const FIXTURE_NAMES = [
  'bevel_straight_z16.json',
  'bevel_spiral_z8.json',
  'bevel_spiral_z21.json',
];

function gearFromFixture(f: GearFixture): BevelGear {
  return new BevelGear({
    numberOfTeeth: f.params.number_of_teeth,
    module: f.params.module,
    height: f.params.height,
    coneAngle: f.params.cone_angle,
    helixAngle: f.params.helix_angle,
  });
}

describe.each(FIXTURE_NAMES)('profile generation – %s', (name) => {
  const fixture = loadJson<GearFixture>(name);

  it('produces matching profile slices', () => {
    const gear = gearFromFixture(fixture);

    for (const slice of fixture.slices) {
      const sample = gear.profileAtZ(slice.z);

      expect(sample.localTransform.scale).toBeCloseTo(slice.local_transform.scale, 9);
      expect(sample.localTransform.angle).toBeCloseTo(slice.local_transform.angle, 9);
      expect(dist3(sample.localTransform.center, slice.local_transform.center)).toBeLessThan(1e-9);

      expect(sample.raRadius).toBeCloseTo(slice.ra_radius, 7);
      expect(sample.rdRadius).toBeCloseTo(slice.rd_radius, 7);
      expect(sample.roRadius).toBeCloseTo(slice.ro_radius, 7);

      // tooth curve chain, sampled at the same parameters
      for (let i = 0; i < slice.tooth_curve.t.length; i++) {
        const p = sample.toothCurvePoint(slice.tooth_curve.t[i]);
        expect(
          dist3(p, slice.tooth_curve.points[i]),
          `z=${slice.z} tooth t=${slice.tooth_curve.t[i]}`,
        ).toBeLessThan(1e-6);
      }

      // full closed one-pitch profile with transform applied
      if (slice.profile_transformed) {
        for (let i = 0; i < slice.profile_transformed.t.length; i++) {
          const p = sample.profilePointTransformed(slice.profile_transformed.t[i]);
          expect(
            dist3(p, slice.profile_transformed.points[i]),
            `z=${slice.z} profile t=${slice.profile_transformed.t[i]}`,
          ).toBeLessThan(1e-6);
        }
      }
    }
  });
});
