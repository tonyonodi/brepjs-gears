/**
 * B-rep construction of the bevel gear solid using brepjs.
 *
 * Port of py_gearworks.conv_build123d.GearBuilder (cone_angle != 0 path), with
 * one architectural change: instead of fitting rational Bezier patches to the
 * profile slices (scipy optimization + Face.make_bezier_surface), the
 * analytically-exact profile slices are lofted with OCCT ThruSections.
 * For straight bevels (2 sections) the ruled loft is geometrically identical;
 * for spiral bevels the loft interpolates the same sections the Bezier fit
 * approximates.
 */
import {
  box, cone, cut, fuse, getFaces, getSolids, getBounds, interpolateCurve, intersect,
  isOk, loft, measureDistance, measureVolume, rotate, split,
  sphere, translate, unwrap, vertex, wire,
} from 'brepjs';
import type { AnyShape, Edge, Face, Result, Shape3D, Solid, Wire } from 'brepjs';
import { RAD2DEG, Vec3 } from './defs.js';
import type { BevelGear } from './bevelGear.js';
import { Curve, CurveChain, findCurveIntersect } from './curve.js';
import { mmul, rotZ, toRotvec, vnorm, vscale, vsub } from './vec.js';
import { GearRefProfile, numTeethAct, profileChain } from './core.js';
import { GearTransform } from './transform.js';

const SIDE_SURFACE_EXTENSION_RATIO = 0.01;
const OVERSAMPLING_RATIO = 3;
const POINTS_PER_CURVE = 33;

function expect<T>(r: Result<T, unknown>, what: string): T {
  if (!isOk(r)) {
    throw new Error(`${what} failed: ${JSON.stringify((r as { error?: unknown }).error)}`);
  }
  return unwrap(r) as T;
}

/** Sample an active curve and build a BSpline edge through the points. */
function curveToEdge(
  curve: Curve,
  transform: GearTransform,
  { mirrorY = false, reversed = false }: { mirrorY?: boolean; reversed?: boolean } = {},
): Edge {
  const pts: Vec3[] = [];
  for (let i = 0; i < POINTS_PER_CURVE; i++) {
    let t = i / (POINTS_PER_CURVE - 1);
    if (reversed) t = 1 - t;
    const p = curve.at(t);
    const q: Vec3 = mirrorY ? [p[0], -p[1], p[2]] : p;
    pts.push(transform.apply(q));
  }
  return expect<Edge>(interpolateCurve(pts), 'interpolateCurve');
}

/**
 * The undercut trochoid kept by py_gearworks' trim still contains its loop, so
 * the undercut and flank curves cross once more above the chain junction —
 * making the one-pitch profile self-intersecting. py_gearworks feeds those
 * self-intersecting surfaces to the splitter, which discards the loop wedge.
 * OCCT's splitter chokes on it here instead, so we clip both curves at that
 * secondary crossing before building the wire (same resulting solid).
 */
function deCrossToothCurves(toothCurves: Curve[]): Curve[] {
  if (toothCurves.length < 2) return toothCurves;
  const [undercut, flank] = [toothCurves[0], toothCurves[1]];
  let best: { su: number; sf: number } | null = null;
  for (let i = 0; i <= 6; i++) {
    for (let j = 0; j <= 6; j++) {
      const sol = findCurveIntersect(undercut, flank, [0.05 + (0.9 * i) / 6, 0.05 + (0.9 * j) / 6]);
      const su = sol.x[0];
      const sf = sol.x[1];
      const d = vnorm(vsub(undercut.at(su), flank.at(sf)));
      if (d > 1e-9) continue;
      if (su < 0.01 || su > 0.999 || sf < 0.001 || sf > 0.99) continue; // junction or out of range
      if (su > 0.999 - 1e-6 && sf < 1e-3) continue;
      if (!best || sf > best.sf) best = { su, sf };
    }
  }
  if (!best) return toothCurves;
  const undercut2 = undercut.clone();
  const flank2 = flank.clone();
  undercut2.setEndOn(best.su);
  flank2.setStartOn(best.sf);
  return [undercut2, flank2, ...toothCurves.slice(2)];
}

/**
 * One-pitch profile wire (rd -> tooth -> ra -> mirrored tooth) at a given z.
 * The mirrored tooth is built per sub-curve (in reverse order) so that the
 * kink between undercut and flank stays an exact wire vertex instead of being
 * smoothed over by a single interpolated edge.
 */
function profileWireAtZ(profile: GearRefProfile, globalTrf: GearTransform): Wire {
  // compose the global gear transform (module scale, angle, center) with the
  // per-slice local transform so the solid is built directly in final space
  const trf = globalTrf.mul(profile.transform);
  const toothCurves = deCrossToothCurves(
    (profile.toothCurve instanceof CurveChain
      ? profile.toothCurve.getCurves()
      : [profile.toothCurve]
    ).filter((c) => c.active),
  );

  const edges: Edge[] = [];
  if (profile.rdCurve.active) edges.push(curveToEdge(profile.rdCurve, trf));
  for (const c of toothCurves) edges.push(curveToEdge(c, trf));
  if (profile.raCurve.active) edges.push(curveToEdge(profile.raCurve, trf));
  for (let i = toothCurves.length - 1; i >= 0; i--) {
    edges.push(curveToEdge(toothCurves[i], trf, { mirrorY: true, reversed: true }));
  }
  return expect<Wire>(wire(edges), 'wire assembly');
}

/** Port of GearBuilder.gen_ref_solid: the blank that gets cut by tooth surfaces. */
function genRefSolid(gear: BevelGear): Shape3D {
  const core = gear.gearcore;
  const z0 = core.zVals[0];
  const z1 = core.zVals[core.zVals.length - 1];

  const profile0 = core.curveGenAtZ(z0);
  const profile1 = core.curveGenAtZ(z1);

  // blank is built in final space: sphere data uses the global gear transform
  const { center: center0, R: R0abs } = core.sphereDataAtZ(z0);
  const { center: center1, R: R1abs } = core.sphereDataAtZ(z1);
  const R0 = Math.abs(R0abs);
  const R1 = Math.abs(R1abs);

  const sph1 = translate(sphere(R0), center0);
  const sph2 = translate(sphere(R1), center1);

  const union = expect<Shape3D>(fuse(sph1, sph2), 'sphere fuse');
  const lens = expect<Shape3D>(intersect(sph1, sph2), 'sphere intersect');
  let refSolid: Shape3D = expect<Shape3D>(cut(union, lens), 'sphere cut');

  const trf0 = core.transform.mul(profile0.transform);
  const trf1 = core.transform.mul(profile1.transform);
  const cO0 = trf0.apply(profile0.roCurve.center);
  const cO1 = trf1.apply(profile1.roCurve.center);
  const hO = cO1[2] - cO0[2];
  const rO0 = profile0.roCurve.radius * trf0.scale;
  const rO1 = profile1.roCurve.radius * trf1.scale;
  const rOCone = cone(rO0, rO1, hO, { at: [0, 0, cO0[2]] });

  // if the spherical symmetric difference yields several solids, keep the one
  // closest to the gear body (mirror of build123d sort_by_distance)
  const pieces = getSolids(refSolid as AnyShape);
  if (pieces.length > 1) {
    const target = vertex([
      (cO0[0] + cO1[0]) / 2, (cO0[1] + cO1[1]) / 2, (cO0[2] + cO1[2]) / 2,
    ]);
    let best: Solid = pieces[0];
    let bestD = Infinity;
    for (const p of pieces) {
      const d = expect<number>(
        measureDistance(p as never, target as never) as never, 'measureDistance',
      );
      if (d < bestD) {
        bestD = d;
        best = p;
      }
    }
    refSolid = best as Shape3D;
  }

  // outside-teeth path: fuse with the outer cone, then keep only z >= cO0.z
  refSolid = expect<Shape3D>(fuse(refSolid, rOCone, { simplify: true }), 'blank fuse');
  // box used as a half-space: removes everything below the plane z = cO0.z
  const L = 8 * Math.max(R0, R1);
  const belowBox = box(L, L, L, { at: [0, 0, cO0[2] - L / 2], centered: true });
  refSolid = expect<Shape3D>(cut(refSolid, belowBox), 'blank plane split');
  return refSolid;
}

/**
 * Move a solid built in axis-aligned (scale-only) space into its final pose:
 * rotate by orientation * rotZ(angle) about the origin, then translate.
 * Mirror of py_gearworks conv_build123d.apply_transform_part.
 */
function placeSolid(solid: Shape3D, trf: GearTransform): Shape3D {
  const rotvec = toRotvec(mmul(trf.orientation, rotZ(trf.angle)));
  const angle = vnorm(rotvec);
  let out = solid;
  if (angle > 1e-12) {
    out = rotate(out, angle * RAD2DEG, {
      at: [0, 0, 0],
      axis: vscale(rotvec, 1 / angle),
    }) as Shape3D;
  }
  if (vnorm(trf.center) > 0) {
    out = translate(out, trf.center) as Shape3D;
  }
  return out;
}

/** Build the gear solid in final (module-scaled, positioned) space. */
export function buildBevelGearSolid(gear: BevelGear): Shape3D {
  const core = gear.gearcore;
  // The construction assumes the gear axis is +Z through the origin (cones,
  // half-space cuts and the tooth pattern all use world Z), so build with the
  // module scale only and pose the finished solid afterwards — mirror of
  // py_gearworks GearBuilder, which builds a copy with an identity transform
  // and applies the real one to the resulting Part.
  const finalTrf = core.transform;
  core.transform = new GearTransform({ scale: finalTrf.scale });
  try {
    return placeSolid(buildAxisAlignedSolid(gear), finalTrf);
  } finally {
    core.transform = finalTrf;
  }
}

function buildAxisAlignedSolid(gear: BevelGear): Shape3D {
  const core = gear.gearcore;
  const z0 = core.zVals[0];
  const z1 = core.zVals[core.zVals.length - 1];
  const zdiff = z1 - z0;

  const nVert = gear.defaultNVert();
  const nZTweens = Math.max(2, Math.ceil((nVert - 2) * OVERSAMPLING_RATIO) + 2);

  // side surfaces use a slightly extended z range so the split tool fully
  // pierces the blank (mirror of side_surface_extension_ratio)
  const zLo = z0 - SIDE_SURFACE_EXTENSION_RATIO * zdiff;
  const zHi = z1 + SIDE_SURFACE_EXTENSION_RATIO * zdiff;

  const sliceWires: Wire[] = [];
  for (let i = 0; i < nZTweens; i++) {
    const z = zLo + ((zHi - zLo) * i) / (nZTweens - 1);
    sliceWires.push(profileWireAtZ(core.curveGenAtZ(z), core.transform));
  }

  // open section wires make ThruSections produce a shell of side faces
  const strip = expect<Shape3D>(
    loft(sliceWires, { ruled: nZTweens === 2 }),
    'tooth strip loft',
  );
  const stripFaces = getFaces(strip as AnyShape);

  const nTeeth = numTeethAct(core.toothParam);
  const pitchDeg = (core.pitchAngle * RAD2DEG);
  const allFaces: Face[] = [];
  for (let j = 0; j < nTeeth; j++) {
    for (const f of stripFaces) {
      allFaces.push(j === 0 ? f : rotate(f, pitchDeg * j, { at: [0, 0, 0], axis: [0, 0, 1] }));
    }
  }

  const blank = genRefSolid(gear);
  const splitRes = expect<AnyShape>(
    split(blank as never, allFaces as never[]) as never,
    'splitting blank with tooth surfaces',
  );
  const withVolumes = getSolids(splitRes).map((p) => ({
    p,
    vol: expect<number>(measureVolume(p as never) as never, 'volume'),
  }));
  const maxVol = Math.max(...withVolumes.map((w) => w.vol));
  // the splitter can leave tiny sliver solids along strip seams — drop them
  const pieces = withVolumes.filter((w) => w.vol > 1e-4 * maxVol).map((w) => w.p);
  if (pieces.length < 2) {
    throw new Error('Split of blank solid via gear surfaces produced a single piece');
  }
  // the gear body is the piece with the lowest bounding-box center
  // (mirror of build123d sort by bounding_box().center().Z)
  const withCenters = pieces.map((p) => {
    const b = getBounds(p as AnyShape) as { zMin: number; zMax: number };
    return { p, cz: (b.zMin + b.zMax) / 2 };
  });
  withCenters.sort((a, b) => a.cz - b.cz);
  // everything was already built in final space (transform composed in)
  return withCenters[0].p as Shape3D;
}
