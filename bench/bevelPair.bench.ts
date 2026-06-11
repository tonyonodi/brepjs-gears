/**
 * Benchmark: generating a meshed pair of spiral (helixed) bevel gears.
 *
 * Solid construction is OCCT-bound and takes tens of seconds per gear, so each
 * bench runs a fixed (small) number of full iterations instead of tinybench's
 * default time-budget sampling. Bump BENCH_ITERATIONS for tighter statistics:
 *
 *   BENCH_ITERATIONS=3 npm run bench
 *
 * The bench callbacks must be `async`: tinybench detects async-ness of a plain
 * function by calling it (once in warmup, once in run), which here would cost
 * two unmeasured ~40s gear builds per bench.
 */
import { bench, describe } from 'vitest';
import { BevelGear, BevelGearParams } from '../src/index.js';
import { initKernel } from '../tests/helpers.js';

await initKernel();

const ITERATIONS = Number(process.env.BENCH_ITERATIONS ?? 1);
const OPTS = {
  iterations: ITERATIONS,
  warmupIterations: 0,
  warmupTime: 0,
  time: 0,
  throws: true,
};

// the z8/z21 spiral pair on perpendicular axes from examples/bevel-pair.mts
const z1 = 8;
const z2 = 21;
const gamma = Math.atan2(z1, z2);
const beta = Math.PI / 6;

const pinionParams: BevelGearParams = {
  numberOfTeeth: z1,
  module: 2,
  height: 10,
  coneAngle: 2 * gamma,
  helixAngle: beta,
};
const wheelParams: BevelGearParams = {
  numberOfTeeth: z2,
  module: 2,
  height: 10,
  coneAngle: 2 * (Math.PI / 2 - gamma),
  helixAngle: -beta,
};

describe('spiral bevel pair (z8/z21, m2, beta=30deg)', () => {
  bench('full pair: construct, meshTo, build both solids', async () => {
    const pinion = new BevelGear(pinionParams);
    const wheel = new BevelGear(wheelParams);
    pinion.meshTo(wheel);
    pinion.buildSolid();
    wheel.buildSolid();
  }, OPTS);

  bench('pinion (z8) solid only', async () => {
    new BevelGear(pinionParams).buildSolid();
  }, OPTS);

  bench('wheel (z21) solid only', async () => {
    new BevelGear(wheelParams).buildSolid();
  }, OPTS);

  bench('profile math only: construct + meshTo + reference profiles', async () => {
    const pinion = new BevelGear(pinionParams);
    const wheel = new BevelGear(wheelParams);
    pinion.meshTo(wheel);
    for (const gear of [pinion, wheel]) {
      const [z0, z1g] = [gear.gearcore.zVals[0], gear.gearcore.zVals[1]];
      for (let i = 0; i < 4; i++) {
        gear.profileAtZ(z0 + ((z1g - z0) * i) / 3);
      }
    }
  }, OPTS);
});
