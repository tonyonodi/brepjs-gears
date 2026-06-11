import { describe, it, expect } from 'vitest';
import { BevelGear, conicCenter, Vec3 } from '../src/index.js';
import { Mat3, mmul, mvmul, rotZ } from '../src/vec.js';
import { loadJson, dist3, initKernel } from './helpers.js';

interface MeshToFixture {
  name: string;
  params: {
    z1: number;
    z2: number;
    cone_angle_1: number;
    cone_angle_2: number;
    module: number;
    target_dir: number[];
    other_center: number[];
    other_angle: number;
    other_orientation: number[][];
  };
  expected: {
    center: number[];
    orientation: number[][];
    angle: number;
  };
}

const FIXTURES = loadJson<MeshToFixture[]>('meshto_fixtures.json');

// meshTo is pure transform math; the reference values come straight from the
// py_gearworks gearmath functions (fixtures/gen_meshto_fixtures.py).
const ATOL = 1e-9;

function buildPair(params: MeshToFixture['params']): { pinion: BevelGear; wheel: BevelGear } {
  const pinion = new BevelGear({
    numberOfTeeth: params.z1,
    coneAngle: params.cone_angle_1,
    module: params.module,
  });
  const wheel = new BevelGear({
    numberOfTeeth: params.z2,
    coneAngle: params.cone_angle_2,
    module: params.module,
    center: params.other_center as Vec3,
    angle: params.other_angle,
  });
  wheel.gearcore.transform.orientation = params.other_orientation as Mat3;
  return { pinion, wheel };
}

describe.each(FIXTURES)('meshTo – $name', ({ params, expected }) => {
  it('matches the py_gearworks placement', () => {
    const { pinion, wheel } = buildPair(params);
    pinion.meshTo(wheel, params.target_dir as Vec3);

    const trf = pinion.gearcore.transform;
    expect(dist3(trf.center, expected.center)).toBeLessThan(ATOL);
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        expect(Math.abs(trf.orientation[i][j] - expected.orientation[i][j])).toBeLessThan(ATOL);
      }
    }
    expect(Math.abs(trf.angle - expected.angle)).toBeLessThan(ATOL);
  });

  it('makes the pitch cone apexes coincide', () => {
    const { pinion, wheel } = buildPair(params);
    pinion.meshTo(wheel, params.target_dir as Vec3);
    expect(dist3(conicCenter(pinion.coneData), conicCenter(wheel.coneData))).toBeLessThan(1e-9);
  });
});

describe('meshTo – backlash and angle bias', () => {
  const params = FIXTURES[0].params;

  it('angleBias shifts the mesh angle within the backlash', () => {
    const { pinion, wheel } = buildPair(params);
    pinion.meshTo(wheel, params.target_dir as Vec3, 0.5, 0);
    const centered = pinion.gearcore.transform.angle;

    pinion.meshTo(wheel, params.target_dir as Vec3, 0.5, 1);
    const biased = pinion.gearcore.transform.angle;

    const backlashAct = params.module * 0.5;
    const expectedShift =
      backlashAct / 2 / pinion.pitchRadius / Math.cos(pinion.params.pressureAngle);
    expect(Math.abs(biased - centered - expectedShift)).toBeLessThan(1e-12);
  });

  it('builds the meshed pinion solid in place (kernel)', async () => {
    await initKernel();
    const { measureVolume, getBounds, unwrap } = await import('brepjs');
    const mkPinion = () => new BevelGear({
      numberOfTeeth: params.z1,
      coneAngle: params.cone_angle_1,
      module: params.module,
      height: 4,
    });

    const local = mkPinion().buildSolid();
    const meshed = mkPinion();
    const wheel = new BevelGear({
      numberOfTeeth: params.z2,
      coneAngle: params.cone_angle_2,
      module: params.module,
      height: 4,
    });
    meshed.meshTo(wheel, params.target_dir as Vec3);
    const placed = meshed.buildSolid();

    // rigid placement: volume is preserved
    const vLocal = unwrap(measureVolume(local as never)) as number;
    const vPlaced = unwrap(measureVolume(placed as never)) as number;
    expect(Math.abs(vPlaced - vLocal) / vLocal).toBeLessThan(1e-6);

    // the bounding-box center follows the gear transform (the gear is nearly
    // rotationally symmetric about its axis, so AABB distortion is small)
    const cb = (s: unknown): Vec3 => {
      const b = getBounds(s as never) as unknown as {
        xMin: number; xMax: number; yMin: number; yMax: number; zMin: number; zMax: number;
      };
      return [(b.xMin + b.xMax) / 2, (b.yMin + b.yMax) / 2, (b.zMin + b.zMax) / 2];
    };
    const trf = meshed.gearcore.transform;
    const rot = mmul(trf.orientation, rotZ(trf.angle));
    const expected = mvmul(rot, cb(local)).map((v, i) => v + trf.center[i]) as Vec3;
    expect(dist3(cb(placed), expected)).toBeLessThan(0.5);
  }, 240000);

  it('defaults backlash to the sum of both gears\' backlash parameters', () => {
    const mk = (backlash: number) => {
      const pair = buildPair(params);
      const pinion = new BevelGear({
        numberOfTeeth: params.z1,
        coneAngle: params.cone_angle_1,
        module: params.module,
        backlash,
      });
      return { pinion, wheel: pair.wheel };
    };

    const a = mk(0.2);
    a.pinion.meshTo(a.wheel, params.target_dir as Vec3, undefined, 1);

    const b = mk(0);
    b.pinion.meshTo(b.wheel, params.target_dir as Vec3, 0.2, 1);

    expect(Math.abs(a.pinion.gearcore.transform.angle - b.pinion.gearcore.transform.angle))
      .toBeLessThan(1e-12);
  });
});
