import { describe, it, expect, beforeAll } from 'vitest';
import { BevelGear } from '../src/index.js';
import { loadJson, loadStepBlob, dist3, initKernel, GearFixture } from './helpers.js';

const FIXTURE_NAMES = [
  'bevel_straight_z16',
  'bevel_spiral_z8',
  'bevel_spiral_z21',
];

beforeAll(async () => {
  await initKernel();
});

// Comparison tolerances.
//
// Notes on what "equal" can mean here:
// - OCCT mass properties on these B-reps carry noise: the py_gearworks
//   reference itself reports a center-of-mass y of ~0.008 for a perfectly
//   16-fold-symmetric gear. Volume noise is a few 0.1%.
// - py_gearworks' splitter keeps the tiny self-intersection wedge of the
//   undercut trochoid loop at the tooth root; this port trims it, which is a
//   ~0.1% volume difference localized below the working flank.
// - Direct boolean subtraction of two nearly-coincident complex solids is
//   numerically unreliable in OCCT, so the strong geometric checks are
//   per-plane section areas and boundary-sample distances instead.
const VOLUME_RTOL = 5e-3;
const BBOX_ATOL = 1e-2;
const COM_ATOL = 2e-2;
const SLAB_RTOL = 1e-2;
// Boundary samples on the undercut trochoid (which contains the loop wedge —
// see note above) deviate up to ~6e-3; the working flank stays within 2e-3.
const BOUNDARY_ATOL = 8e-3;

describe.each(FIXTURE_NAMES)('bevel gear solid – %s', (name) => {
  const fixture = loadJson<GearFixture>(`${name}.json`);

  it('builds a valid solid matching the py_gearworks reference', async () => {
    const {
      isOk, unwrap, importSTEP, measureVolume, measureVolumeProps,
      getBounds, getSolids, isShape3D, isValid, box, intersect,
      vertex, measureDistance,
    } = await import('brepjs');

    const gear = new BevelGear({
      numberOfTeeth: fixture.params.number_of_teeth,
      module: fixture.params.module,
      height: fixture.params.height,
      coneAngle: fixture.params.cone_angle,
      helixAngle: fixture.params.helix_angle,
    });

    const solid = gear.buildSolid() as any;
    expect(isShape3D(solid)).toBe(true);
    expect(isValid(solid)).toBe(true);

    const ref = fixture.solid_stats;

    const vol = unwrap(measureVolume(solid));
    expect(Math.abs(vol - ref.volume) / ref.volume, 'relative volume error').toBeLessThan(VOLUME_RTOL);

    const bb = getBounds(solid) as any;
    const bbMin = [bb.xMin, bb.yMin, bb.zMin];
    const bbMax = [bb.xMax, bb.yMax, bb.zMax];
    for (let i = 0; i < 3; i++) {
      expect(Math.abs(bbMin[i] - ref.bbox_min[i]), `bbox min[${i}]`).toBeLessThan(BBOX_ATOL);
      expect(Math.abs(bbMax[i] - ref.bbox_max[i]), `bbox max[${i}]`).toBeLessThan(BBOX_ATOL);
    }

    const props = unwrap(measureVolumeProps(solid)) as any;
    const com: number[] = Array.isArray(props.centerOfMass)
      ? props.centerOfMass
      : [props.centerOfMass.x, props.centerOfMass.y, props.centerOfMass.z];
    expect(dist3(com, ref.center_of_mass), 'center of mass').toBeLessThan(COM_ATOL);

    // load reference STEP
    const stepRes = await importSTEP(loadStepBlob(`${name}.step`));
    expect(isOk(stepRes)).toBe(true);
    const refSolids = getSolids(unwrap(stepRes) as any);
    expect(refSolids.length).toBe(1);
    const refSolid = refSolids[0];

    // z-distribution of material: volume inside horizontal slabs must match.
    // (Plane sections are fragile on grazing cuts; solid ∩ box booleans are
    // robust and compare the same property integrated over each slab.)
    const zMin = ref.bbox_min[2];
    const zMax = ref.bbox_max[2];
    const rMax = Math.max(...ref.bbox_max.slice(0, 2).map(Math.abs)) * 2.5;
    const slabFractions = [0, 0.25, 0.5, 0.75, 1];
    for (let s = 0; s < slabFractions.length - 1; s++) {
      const lo = zMin + (zMax - zMin) * slabFractions[s];
      const hi = zMin + (zMax - zMin) * slabFractions[s + 1];
      const slab = box(rMax, rMax, hi - lo, { at: [0, 0, (lo + hi) / 2], centered: true });
      const vMine = unwrap(measureVolume(unwrap(intersect(solid, slab as any)) as any));
      const vRef = unwrap(measureVolume(unwrap(intersect(refSolid as any, slab as any)) as any));
      expect(
        Math.abs(vMine - vRef),
        `slab volume z in [${lo.toFixed(2)}, ${hi.toFixed(2)}] (mine=${vMine}, ref=${vRef})`,
      ).toBeLessThan(SLAB_RTOL * Math.max(vRef, ref.volume * 0.02));
    }

    // boundary proximity: sampled points of the analytic tooth profile must
    // lie on or inside the reference solid (within tolerance). Gears are
    // built at center=origin, angle=0, so final space = module * profile.
    const module = fixture.params.module;
    for (const zr of [0.25, 0.75]) {
      const z = fixture.z_vals[0] + (fixture.z_vals[1] - fixture.z_vals[0]) * zr;
      const sample = gear.profileAtZ(z);
      for (let i = 0; i <= 12; i++) {
        const p = sample.profilePointTransformed(i / 12);
        const pFinal: [number, number, number] = [p[0] * module, p[1] * module, p[2] * module];
        const d = unwrap(measureDistance(refSolid as any, vertex(pFinal) as any));
        expect(d, `boundary point t=${i / 24} z=${z}`).toBeLessThan(BOUNDARY_ATOL);
      }
    }
  });
});
