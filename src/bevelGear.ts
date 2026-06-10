/**
 * Port of py_gearworks.wrapper.BevelGear (InvoluteGear.calc_params, bevel path).
 */
import { ORIGIN, OUT, PI, Vec3 } from './defs.js';
import { vscale } from './vec.js';
import {
  Gear, GearRefProfile, generateProfileExtensions, profileChain, ProfileRecipeAtZ,
} from './core.js';
import { GearTransform } from './transform.js';
import {
  OctoidTooth, OctoidUndercutTooth, ToothLimitParam,
} from './gearTeeth.js';
import { buildBevelGearSolid } from './builder.js';

/** Parameters for a bevel gear, mirroring py_gearworks.BevelGear. */
export interface BevelGearParams {
  numberOfTeeth: number;
  /** Full cone angle in radians (2 * pitch half-angle). Default PI/2. */
  coneAngle?: number;
  /** Helix (spiral) angle in radians. Default 0. */
  helixAngle?: number;
  /** Tooth surface height (along the cone surface), real units. Default 1. */
  height?: number;
  module?: number;
  enableUndercut?: boolean;
  rootFillet?: number;
  tipFillet?: number;
  tipTruncation?: number;
  profileShift?: number;
  addendumCoefficient?: number;
  dedendumCoefficient?: number;
  /** Pressure angle in radians. Default 20 degrees. */
  pressureAngle?: number;
  backlash?: number;
  crowning?: number;
  zAnchor?: number;
  center?: Vec3;
  angle?: number;
}

export interface ProfileSliceSample {
  z: number;
  localTransform: { center: Vec3; angle: number; scale: number };
  raRadius: number;
  rdRadius: number;
  roRadius: number;
  /** Sample the tooth curve chain at parameter t in [0, 1]. */
  toothCurvePoint(t: number): Vec3;
  /** Sample the one-pitch profile (with local transform applied) at t in [0, 1]. */
  profilePointTransformed(t: number): Vec3;
}

interface ResolvedParams extends Required<Omit<BevelGearParams, 'center' | 'angle'>> {
  center: Vec3;
  angle: number;
}

export class BevelGear {
  readonly params: ResolvedParams;
  readonly gearcore: Gear;

  constructor(params: BevelGearParams) {
    const p: ResolvedParams = {
      numberOfTeeth: params.numberOfTeeth,
      coneAngle: params.coneAngle ?? PI / 2,
      helixAngle: params.helixAngle ?? 0,
      height: params.height ?? 1.0,
      module: params.module ?? 1.0,
      enableUndercut: params.enableUndercut ?? true,
      rootFillet: params.rootFillet ?? 0.0,
      tipFillet: params.tipFillet ?? 0.0,
      tipTruncation: params.tipTruncation ?? 0.1,
      profileShift: params.profileShift ?? 0,
      addendumCoefficient: params.addendumCoefficient ?? 1.0,
      dedendumCoefficient: params.dedendumCoefficient ?? 1.2,
      pressureAngle: params.pressureAngle ?? (20 * PI) / 180,
      backlash: params.backlash ?? 0,
      crowning: params.crowning ?? 0,
      zAnchor: params.zAnchor ?? 0,
      center: params.center ?? ORIGIN,
      angle: params.angle ?? 0,
    };
    if (p.coneAngle === 0) {
      throw new Error('brepjs-gears BevelGear requires a non-zero cone angle');
    }
    this.params = p;
    this.gearcore = this.calcParams();
  }

  get beta(): number {
    return this.params.helixAngle;
  }

  get gamma(): number {
    return this.params.coneAngle / 2;
  }

  private calcParams(): Gear {
    const p = this.params;
    const rpRef = p.numberOfTeeth / 2;
    const pitchAngle = (2 * PI) / p.numberOfTeeth;
    const gamma = p.coneAngle / 2;

    const crowningFunc = (z: number, offset: number): number =>
      offset - ((z * 2) / p.height - 1) ** 2 * (p.crowning / rpRef) * 1e-3;

    const backlashAngleVal = p.backlash / 2 / (rpRef * Math.cos(p.pressureAngle));

    const toothAngle =
      pitchAngle / 4 +
      (p.profileShift * Math.tan(p.pressureAngle)) / rpRef -
      backlashAngleVal;

    const spiralCoeff = Math.tan(this.beta) / rpRef;
    const angleFunc = (z: number): number => z * spiralCoeff;

    const hD = p.dedendumCoefficient - p.profileShift;
    const hA = p.addendumCoefficient + p.profileShift;
    const hO = Math.max(hD + 0.5, 2);
    const limits: ToothLimitParam = { hD, hA, hO };
    // py_gearworks quirk kept: the non-undercut OctoidTooth is constructed
    // without ref_limits, so it uses the defaults for its lower-limit check.
    const defaultLimits: ToothLimitParam = { hA: 1, hD: 1.2, hO: 2 };

    const zH = p.height / p.module;
    const zVals = [-zH * p.zAnchor, zH * (1 - p.zAnchor)];

    const shapeRecipe = (z: number): ProfileRecipeAtZ => {
      const toothParams = {
        pressureAngle: p.pressureAngle,
        pitchRadius: rpRef,
        pitchIntersectAngle: crowningFunc(z, toothAngle),
        coneAngle: p.coneAngle,
        pitchAngle,
        refLimits: p.enableUndercut ? limits : defaultLimits,
      };
      const toothGenerator = p.enableUndercut
        ? new OctoidUndercutTooth(toothParams)
        : new OctoidTooth(toothParams);
      return {
        toothGenerator,
        limits,
        fillet: {
          rootFillet: p.rootFillet,
          tipFillet: p.tipFillet,
          tipReduction: p.tipTruncation,
        },
        cone: { coneAngle: p.coneAngle, baseRadius: rpRef },
        pitchAngle,
        transform: new GearTransform({
          scale: 1 - (z * 2 * Math.sin(gamma)) / p.numberOfTeeth,
          center: vscale(OUT, z * Math.cos(gamma)),
          angle: angleFunc(z),
        }),
      };
    };

    return new Gear({
      zVals,
      module: p.module,
      toothParam: {
        numTeeth: p.numberOfTeeth,
        numCutoutTeeth: 0,
        insideTeeth: false,
      },
      shapeRecipe,
      transform: new GearTransform({
        center: p.center,
        angle: p.angle,
        scale: p.module,
      }),
    });
  }

  /** Reference profile slice info at a given z (module-1 reference space). */
  profileAtZ(z: number): ProfileSliceSample {
    const profile: GearRefProfile = this.gearcore.curveGenAtZ(z);
    const chain = profileChain(profile);
    const trf = profile.transform;
    return {
      z,
      localTransform: {
        center: trf.center,
        angle: trf.angle,
        scale: trf.scale,
      },
      raRadius: profile.raCurve.radius,
      rdRadius: profile.rdCurve.radius,
      roRadius: profile.roCurve.radius,
      toothCurvePoint: (t: number) => profile.toothCurve.at(t),
      profilePointTransformed: (t: number) => trf.apply(chain.at(t)),
    };
  }

  /** Extended profile (with ring connectors) at z, used by the solid builder. */
  extendedProfileAtZ(z: number) {
    const recipe = this.gearcore.shapeRecipe(z);
    return generateProfileExtensions(this.gearcore.curveGenAtZ(z), recipe.cone);
  }

  /**
   * Number of vertical sections used for the solid construction
   * (port of InvoluteGear.build_part's n_vert heuristic).
   */
  defaultNVert(): number {
    const [z0, z1] = [this.gearcore.zVals[0], this.gearcore.zVals[this.gearcore.zVals.length - 1]];
    const angles: number[] = [];
    for (let i = 0; i < 20; i++) {
      const z = z0 + ((z1 - z0) * i) / 19;
      angles.push(this.gearcore.shapeRecipe(z).transform.angle);
    }
    const twist = Math.abs(Math.max(...angles) - Math.min(...angles));
    if (this.params.crowning === 0 && this.beta === 0) return 2;
    if (twist > PI / 6) return 3 + Math.floor(twist / (PI / 6));
    return 4;
  }

  /** Build the gear solid (brepjs Shape3D), in final (module-scaled) space. */
  buildSolid(): unknown {
    return buildBevelGearSolid(this);
  }
}
