from mating import load_step, get_faces, compute_fingerprint, face_relation, find_valid_correspondence
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
import numpy as np
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_OUT
from OCC.Core.TopoDS import topods
from OCC.Core.TopExp import TopExp_Explorer, topexp
from OCC.Core.BRepClass3d import BRepClass3d_SolidClassifier
from OCC.Core.gp import gp_Pnt
from OCC.Core.GeomAbs import (GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_Cone,
                              GeomAbs_Sphere, GeomAbs_Torus, GeomAbs_BezierSurface,
                              GeomAbs_BSplineSurface, GeomAbs_SurfaceOfRevolution,
                              GeomAbs_SurfaceOfExtrusion, GeomAbs_OffsetSurface,
                              GeomAbs_OtherSurface)
import itertools



TYPE_NAMES = {
    GeomAbs_Plane: "plane", GeomAbs_Cylinder: "cylinder", GeomAbs_Cone: "cone",
    GeomAbs_Sphere: "sphere", GeomAbs_Torus: "torus",
    GeomAbs_BezierSurface: "bezier", GeomAbs_BSplineSurface: "bspline",
    GeomAbs_SurfaceOfRevolution: "revolution",
    GeomAbs_SurfaceOfExtrusion: "extrusion",
    GeomAbs_OffsetSurface: "offset", GeomAbs_OtherSurface: "other",
}

def print_faces(path, solid_index=None, near=None, radius=15.0):
    for s_idx, solid in enumerate(load_step(path)):
        if solid_index is not None and s_idx != solid_index:
            continue
        faces = get_faces(solid)
        shown = 0
        print(f"--- Solid {s_idx} ({len(faces)} faces, near={near}, radius={radius}) ---")
        for f in faces:
            p = f["params"]
            if near is not None and np.linalg.norm(p["centroid"] - np.array(near)) > radius:
                continue
            shown += 1
            line = f"  Face[{f['index']:3d}] {f['type']:<10} area={p['area']:8.2f}"
            if "radius" in p:
                line += f" r={p['radius']:7.3f}"
            if "convex" in p:
                line += "  male  " if p["convex"] else "  female"
            if "axis_direction" in p:
                d = p["axis_direction"]
                line += f"  dir=({d[0]:+.3f},{d[1]:+.3f},{d[2]:+.3f})"
            c = p["centroid"]
            line += f"  c=({c[0]:7.2f},{c[1]:7.2f},{c[2]:7.2f})"
            print(line)
        print(f"--- shown {shown} of {len(faces)} ---")


def face_census(path):
    for s_idx, solid in enumerate(load_step(path)):
        counts, total = {}, 0
        exp = TopExp_Explorer(solid, TopAbs_FACE)
        while exp.More():
            name = TYPE_NAMES.get(BRepAdaptor_Surface(topods.Face(exp.Current())).GetType(), "?")
            counts[name] = counts.get(name, 0) + 1
            total += 1
            exp.Next()
        print(f"solid {s_idx}: {total} faces")
        print(f"   {counts}")


def check_outward(path, solid_index, face_ids, eps=1e-3):
    """Step a hair along each face's computed outward normal. If we land
    outside the solid, the normal was right."""
    solid = load_step(path)[solid_index]
    for f in get_faces(solid):
        if f["index"] not in face_ids:
            continue
        p = f["params"]
        if "axis_direction" not in p:
            continue
        pt = p["centroid"] + eps * p["axis_direction"]
        clf = BRepClass3d_SolidClassifier(solid, gp_Pnt(*pt), 1e-7)
        outside = clf.State() == TopAbs_OUT
        print(f"face {f['index']}: normal points {'OUT (correct)' if outside else 'INTO SOLID (flip is inverted)'}")


def why_not(part_path, sel_ids, assembly_path, cand_solid, cand_ids):
    """Feed the matcher the known-correct answer and print which test kills it."""
    sel_solid = load_step(part_path)[0]
    sel_faces = [f for f in get_faces(sel_solid) if f["index"] in sel_ids]
    cand_faces = [f for f in get_faces(load_step(assembly_path)[cand_solid])
                  if f["index"] in cand_ids]

    print(f"selection: {len(sel_faces)} faces, candidate: {len(cand_faces)} faces")
    for f in sel_faces + cand_faces:
        p = f["params"]
        print(f"  {f['type']:<9} area={p['area']:8.2f} "
              f"dir={np.round(p.get('axis_direction', [0,0,0]), 3)}")

    fp = compute_fingerprint(sel_faces)
    for (i, j), (sd, sdist) in fp.items():
        print(f"sel pair ({i},{j}): dot={sd:.4f} dist={sdist:.4f}")
    for i, j in itertools.combinations(range(len(cand_faces)), 2):
        cd, cdist = face_relation(cand_faces[i], cand_faces[j])
        print(f"cand pair ({i},{j}): dot={cd:.4f} dist={cdist:.4f}")

    mapping, fit = find_valid_correspondence(sel_faces, cand_faces, fp)
    print(f"result: mapping={mapping} fit={fit}")
    
def select_faces(path, solid_index=0, kind=None, area=None, radius=None,
                 near=None, within=None, direction=None, convex=None, verbose=True):
    out = []
    print(f"from{path}")
    for f in get_faces(load_step(path)[solid_index]):
        p = f["params"]
        if kind is not None and f["type"] != kind:
            continue
        if area is not None and not (area[0] <= p["area"] <= area[1]):
            continue
        if radius is not None and not (radius[0] <= p.get("radius", -1) <= radius[1]):
            continue
        if convex is not None and p.get("convex") != convex:
            continue
        if near is not None and np.linalg.norm(p["centroid"] - np.array(near)) > within:
            continue
        if direction is not None:
            d = p.get("axis_direction")
            if d is None or abs(np.dot(d, np.array(direction))) < 0.98:
                continue
        out.append(f["index"])
        if verbose:
            c = p["centroid"]
            print(f"  {f['index']:3d} {f['type']:<9} area={p['area']:8.2f} "
                  f"c=({c[0]:7.2f},{c[1]:7.2f},{c[2]:7.2f})")
    return out


# select_faces("key_shaft.step", kind="plane")
# select_faces("key_shaft.step", kind="cylinder")
# select_faces("hub.step", kind="cylinder", radius=(4.9, 5.1))
# select_faces("hub.step", kind="plane")
select_faces("stepfiles/abracket.STEP", kind="cylinder", verbose=True)

# why_not("key_shaft.step", [0, 2, 5], "hub.step", 0, [3, 4, 6])
# part_path = "stepfiles/rod.step"
# assembly_path = "stepfiles/abracket.step"
# selection_face_ids = [0]

# why_not(part_path, selection_face_ids,assembly_path,0,[0,1,8,15,19] )
# why_not(part_path, selection_face_ids, assembly_path, 19, [2, 4, 6, 9])
 
# why_not(part_path, selection_face_ids, assembly_path, 19, [2, 4, 6, 9])   
# why_not(part_path, selection_face_ids, assembly_path, 20, [2, 4, 6, 9])
# why_not(part_path, selection_face_ids, assembly_path, 21, [2, 4, 6, 9])