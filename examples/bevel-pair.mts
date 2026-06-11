/**
 * Generate a pair of spiral bevel gears and export them as STEP files.
 * Run with: npx tsx examples/bevel-pair.mts
 */
import { OcctKernel } from 'occt-wasm';
import { registerKernel, OcctWasmAdapter, exportSTEP, unwrap } from 'brepjs';
import { writeFileSync } from 'node:fs';
import { BevelGear } from '../src/index.js';

const kernel = await OcctKernel.init();
registerKernel('occt-wasm', OcctWasmAdapter.fromKernel(kernel));

const z1 = 8;
const z2 = 21;
const gamma = Math.atan2(z1, z2);
const beta = Math.PI / 6;

const pinion = new BevelGear({
  numberOfTeeth: z1,
  module: 2,
  height: 10,
  coneAngle: 2 * gamma,
  helixAngle: beta,
});
const wheel = new BevelGear({
  numberOfTeeth: z2,
  module: 2,
  height: 10,
  coneAngle: 2 * (Math.PI / 2 - gamma),
  helixAngle: -beta,
});

// place the pinion against the wheel, so the exported pair meshes
pinion.meshTo(wheel);

for (const [name, gear] of [['pinion', pinion], ['wheel', wheel]] as const) {
  console.time(name);
  const solid = gear.buildSolid();
  const blob = unwrap(exportSTEP(solid as never));
  writeFileSync(`${name}.step`, Buffer.from(await blob.arrayBuffer()));
  console.timeEnd(name);
}
