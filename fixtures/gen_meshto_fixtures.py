"""Generate reference values for the brepjs-gears meshTo port.

Replicates the InvoluteGear.mesh_to bevel branch (wrapper.py:717-746) using the
py_gearworks gearmath functions directly, so it runs without build123d.
"""
import json
import sys

sys.path.insert(0, "/tmp/pgw_src")
import numpy as np
from py_gearworks.base_classes import ConicData, GearTransform
from py_gearworks.gearmath import (
    calc_bevel_gear_placement_vector,
    calc_mesh_angle,
    calc_mesh_orientation,
)
from scipy.spatial.transform import Rotation as scp_Rotation


def mesh_to_bevel(z1, z2, cone_angle_1, cone_angle_2, module, target_dir,
                  other_center, other_angle, other_orientation):
    """Returns (center, orientation, angle) for gear 1 meshed to gear 2."""
    t1 = GearTransform(scale=module)
    t2 = GearTransform(
        scale=module,
        center=np.array(other_center, dtype=float),
        angle=other_angle,
        orientation=np.array(other_orientation, dtype=float),
    )
    cone1 = ConicData(cone_angle=cone_angle_1, base_radius=z1 / 2, transform=t1)
    cone2 = ConicData(cone_angle=cone_angle_2, base_radius=z2 / 2, transform=t2)
    target_dir = np.array(target_dir, dtype=float)
    target_dir = target_dir / np.linalg.norm(target_dir)

    t1.center = calc_bevel_gear_placement_vector(
        target_dir, cone1, cone2, False, False, offset=0,
    )
    # wrapper passes gearcore.cone.R (module units); only used with offset=0
    R_module_units = abs((z1 / 2) / np.sin(cone_angle_1 / 2))
    t1.orientation = calc_mesh_orientation(
        cone_angle_1, cone_angle_2, R_module_units, t2,
        False, False, target_dir, offset=0,
    )
    t1.angle = calc_mesh_angle(
        t1, t2, 2 * np.pi / z1, 2 * np.pi / z2,
        gear1_inside_ring=False, gear2_inside_ring=False,
    )
    return t1


cases = []

# case 1: perpendicular pair, wheel at origin, default direction
g1 = np.arctan2(8, 21)
cases.append(dict(
    name="perpendicular z8/z21 m2, target RIGHT",
    z1=8, z2=21, cone_angle_1=2 * g1, cone_angle_2=2 * (np.pi / 2 - g1),
    module=2.0, target_dir=[1, 0, 0],
    other_center=[0, 0, 0], other_angle=0.0,
    other_orientation=np.eye(3).tolist(),
))

# case 2: same pair, target UP (y)
cases.append(dict(
    name="perpendicular z8/z21 m2, target UP",
    z1=8, z2=21, cone_angle_1=2 * g1, cone_angle_2=2 * (np.pi / 2 - g1),
    module=2.0, target_dir=[0, 1, 0],
    other_center=[0, 0, 0], other_angle=0.0,
    other_orientation=np.eye(3).tolist(),
))

# case 3: wheel pre-transformed (rotated about x by 30deg, offset, spun)
rot = scp_Rotation.from_euler("x", np.pi / 6).as_matrix()
cases.append(dict(
    name="z8/z21 m2, wheel transformed, oblique target",
    z1=8, z2=21, cone_angle_1=2 * g1, cone_angle_2=2 * (np.pi / 2 - g1),
    module=2.0, target_dir=[2, 1, 0.5],
    other_center=[3, -2, 5], other_angle=0.35,
    other_orientation=rot.tolist(),
))

# case 4: non-perpendicular shaft angle (75 deg), z10/z17, m1.5
sigma = np.deg2rad(75)
ratio = 17 / 10
g1b = np.arctan2(np.sin(sigma), ratio + np.cos(sigma))
cases.append(dict(
    name="75deg shaft z10/z17 m1.5, target RIGHT",
    z1=10, z2=17, cone_angle_1=2 * g1b, cone_angle_2=2 * (sigma - g1b),
    module=1.5, target_dir=[1, 0, 0],
    other_center=[0, 0, 0], other_angle=0.0,
    other_orientation=np.eye(3).tolist(),
))

out = []
for c in cases:
    t = mesh_to_bevel(
        c["z1"], c["z2"], c["cone_angle_1"], c["cone_angle_2"], c["module"],
        c["target_dir"], c["other_center"], c["other_angle"],
        c["other_orientation"],
    )
    out.append(dict(
        name=c["name"],
        params=c,
        expected=dict(
            center=t.center.tolist(),
            orientation=t.orientation.tolist(),
            angle=float(t.angle),
        ),
    ))

print(json.dumps(out, indent=2))
