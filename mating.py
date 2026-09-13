from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.TopExp import TopExp_Explorer, topexp
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.GeomAbs import (GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_Cone,
                              GeomAbs_Sphere, GeomAbs_Torus, GeomAbs_BezierSurface,
                              GeomAbs_BSplineSurface, GeomAbs_SurfaceOfRevolution,
                              GeomAbs_SurfaceOfExtrusion, GeomAbs_OffsetSurface,
                              GeomAbs_OtherSurface)
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_SOLID, TopAbs_REVERSED, TopAbs_OUT
from OCC.Core.TopoDS import topods
from OCC.Core.TopTools import TopTools_IndexedDataMapOfShapeListOfShape, TopTools_ListIteratorOfListOfShape
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRepClass3d import BRepClass3d_SolidClassifier
from OCC.Core.gp import gp_Pnt
from OCC.Display.SimpleGui import init_display
from OCC.Core.Quantity import Quantity_Color, Quantity_TOC_RGB

try:
    from OCC.Core.BRepGProp import brepgprop
    _surface_props = brepgprop.SurfaceProperties
    _volume_props = brepgprop.VolumeProperties
except ImportError:
    from OCC.Core.BRepGProp import brepgprop_SurfaceProperties as _surface_props
    from OCC.Core.BRepGProp import brepgprop_VolumeProperties as _volume_props

import numpy as np
import itertools
from collections import deque


dot_product_tolerance = 1e-2  #supposedly step files can have minor deviations
distance_tolerance = 2.0 #using centroid distance so need more tolerance (face sizes can be different which would potentially impact centroid locations)
exact_tolerance = 0.5 #gaps below this = exact fit
max_clearance = 2.0 #max gap permissible
parallel_dot = 0.98 # again some slack bcs apparently step files can have minor deviations

max_hops = 2

TYPE_NAMES = {
    GeomAbs_Plane: "plane", GeomAbs_Cylinder: "cylinder", GeomAbs_Cone: "cone",
    GeomAbs_Sphere: "sphere", GeomAbs_Torus: "torus",
    GeomAbs_BezierSurface: "bezier", GeomAbs_BSplineSurface: "bspline",
    GeomAbs_SurfaceOfRevolution: "revolution",
    GeomAbs_SurfaceOfExtrusion: "extrusion",
    GeomAbs_OffsetSurface: "offset", GeomAbs_OtherSurface: "other",
}

SUPPORTED = {"plane", "cylinder", "cone"}


def load_step(path):  #loads the step file, loads the solids and returns them
    reader = STEPControl_Reader()
    status = reader.ReadFile(path)

    if status != IFSelect_RetDone:
        raise Exception("Error reading STEP file.")
    reader.TransferRoots() #converts STEP to OCCT
    shape = reader.OneShape() #wrapps multiple roots into a single 'shape'

    solids = []
    explorer = TopExp_Explorer(shape, TopAbs_SOLID)
    while explorer.More():
        solid = topods.Solid(explorer.Current())
        solids.append(solid)
        explorer.Next()
    return solids


def solid_volume(solid):  #used by find_all_mates() to skip the selection part's instance if it is included (however, risk of different parts with same volume or the same part repeated being excluded)
    props = GProp_GProps()  #empty container 
    _volume_props(solid, props)  #reads values and adds em to props
    return props.Mass()  #returns the volume (yeah even though it says mass)


def get_faces(solid): #for all faces, creates a dictionary with an internal index, its "face" data,  geometry (plane, cylinder etc), and geometric params
    faces = []
    explorer = TopExp_Explorer(solid, TopAbs_FACE)
    idx = 0
    while explorer.More():
        face = topods.Face(explorer.Current())
        adapter = BRepAdaptor_Surface(face)
        kind = TYPE_NAMES.get(adapter.GetType(), "unknown") #sets type as unknown if it isnt identified inside TYPE_NAMES

        props = GProp_GProps()
        _surface_props(face, props)
        c = props.CentreOfMass()

        params = {  #these are defined for all faces
            "area": props.Mass(), 
            "centroid": np.array([c.X(), c.Y(), c.Z()]),
            "u_span": adapter.LastUParameter() - adapter.FirstUParameter(), #for a cylinder, angle U param is 0 to 2pi around the cylinder, i.e. this tell us if the cylinder is closed or not
            "v_span": adapter.LastVParameter() - adapter.FirstVParameter(), #length/axis of cylinder; distance along slope from apex for cone
            
            # i think u span and v span for a plane are probably just the maximum horizontal and vertical distance 
        }
        flip = -1.0 if face.Orientation() == TopAbs_REVERSED else 1.0  #TopAbs_reversed stores whether or not the face is "reversed", this is apparently not accounted for by the geometric normal 

        if kind == "plane":
            ax = adapter.Plane().Axis()
            params["axis_location"] = np.array([ax.Location().X(), ax.Location().Y(), ax.Location().Z()])
            params["axis_direction"] = flip * np.array([ax.Direction().X(), ax.Direction().Y(), ax.Direction().Z()])
        elif kind == "cylinder":
            cyl = adapter.Cylinder()
            ax = cyl.Axis()
            params["radius"] = cyl.Radius()
            params["axis_location"] = np.array([ax.Location().X(), ax.Location().Y(), ax.Location().Z()])
            params["axis_direction"] = np.array([ax.Direction().X(), ax.Direction().Y(), ax.Direction().Z()]) #flip not required here because for a cylinder, it doesnt really matter where the axis points
            params["convex"] = flip > 0  #convex  True implies male (shaft), False implies Female (hole) 
        elif kind == "cone":
            cone = adapter.Cone()
            ax = cone.Axis()
            params["radius"] = cone.RefRadius()
            params["half_angle"] = cone.SemiAngle()
            params["axis_location"] = np.array([ax.Location().X(), ax.Location().Y(), ax.Location().Z()])
            params["axis_direction"] = np.array([ax.Direction().X(), ax.Direction().Y(), ax.Direction().Z()])
            params["convex"] = flip > 0

        faces.append({"index": idx, "face": face, "type": kind, "params": params})
        idx += 1
        explorer.Next()
    return faces  


def _map_extent(indexed_map):
    for attr in ("Extent", "Size", "Length"):  #different occ versions have different attributs for the max length basically, which is what this function helps us find
        if hasattr(indexed_map, attr):
            return getattr(indexed_map, attr)()
    raise AttributeError(
        "Could not find Extent()/Size()/Length() on this map object -- "
        "check your OCCT/pythonocc-core version."
    )


def adjacency(solid, faces): #returns a map which tells us which faces share an edge
    edge_to_face_map = TopTools_IndexedDataMapOfShapeListOfShape()  #basically an OCC data structure (kinda like a dictionary), in this case mapping edge to corresponding faces
    topexp.MapShapesAndAncestors(solid, TopAbs_EDGE, TopAbs_FACE, edge_to_face_map) #each edge has 2 faces, this creates a map between an edge and its 2 faces

    def face_index(shape): #lookup which basically returns the id of a face
        for f in faces:
            if f["face"].IsSame(shape):
                return f["index"]
        return None

    adjacency_list = {f["index"]: set() for f in faces}

    for i in range(1, _map_extent(edge_to_face_map) + 1 ):
        face_list = edge_to_face_map.FindFromIndex(i) #list of faces from edge i

        faces_for_edge = []
        face_iterator = TopTools_ListIteratorOfListOfShape(face_list)
        while face_iterator.More():  #i mean this loop is kinda unecessary, its only needed to convert a C++ linked list into a python list, face_list and faces_for_edge conceptually hold the same information
            faces_for_edge.append(face_iterator.Value()) 
            face_iterator.Next()

        indexes = [face_index(f) for f in faces_for_edge] #find the indexes of the above selected faces

        for a, b in itertools.combinations(indexes, 2): #for the 2 faces that neighbor this edge, in the adjacency list for a add b (i.e. b is a neighbor of a) and add b for a (i.e. a is a neighbor of b)
            if a is None or b is None:
                continue
            adjacency_list[a].add(b)
            adjacency_list[b].add(a)

    return adjacency_list


def neighborhood(start_idx, adjacency_list, max_hops):  #standard BFS search, returns all the faces within 2 hops of the current face
    visited = {start_idx : 0}
    queue = deque([start_idx]) 

    while queue:
        current = queue.popleft() #picks out the first element
        current_hops = visited[current]
        if current_hops >= max_hops:
            continue
        for neighbor in adjacency_list[current]:
            if neighbor not in visited:
                visited[neighbor] = current_hops + 1
                queue.append(neighbor)

    return set(visited.keys())


def face_relation(f1, f2):  #relationship between any 2 faces

    if f1["type"] == "plane" and f2["type"] == "plane": #if both faces are planes, take dot product of normals and return that along with distance between normals (if parallel) else return distances between centroids
        normal1 = f1["params"]["axis_direction"]
        normal2 = f2["params"]["axis_direction"]
        dot_product = np.dot(normal1, normal2)  #abs not taken as negative dot means planes are looking at each other vs +ve meaning planes look away from each other
        delta = f1["params"]["centroid"] - f2["params"]["centroid"]
        if abs(dot_product) > parallel_dot:
            distance = abs(np.dot(delta, normal1))    # parallel: true gap
        else:
            distance = np.linalg.norm(delta)          # otherwise: centroid distance
        return dot_product, distance

    elif f1["type"] in ("cylinder", "cone") and f2["type"] in ("cylinder", "cone"):
        axis1 = f1["params"]["axis_direction"]
        axis2 = f2["params"]["axis_direction"]
        loc1 = f1["params"]["axis_location"]
        loc2 = f2["params"]["axis_location"]
        dot_product = abs(np.dot(axis1, axis2))  #direction of cylinder/cone's axis not as important because its arbitrary

        cross = np.cross(axis1, axis2)
        cross_norm = np.linalg.norm(cross) #normalized, i.e. unit vector
        if cross_norm < 1e-9: #almost 0, i.e. parallel
            axes_gap = loc2 - loc1
            perpendicular_gap = axes_gap - np.dot(axes_gap, axis1) * axis1  #because location of points can also be arbitrary, this subtracts the component of the gap along any one axes from the gap, returning the perpendicular distance
            #also, here multiplied with axis1 because axes_gap is still a vector and the dot product will result in a scalar. we need a vector along axis one because that's what needs to be subtracted
            distance = np.linalg.norm(perpendicular_gap) #unit vector
        else:
            distance = abs(np.dot(loc2 - loc1, cross)) / cross_norm  #distance along cross product of axes basically

        return dot_product, distance

    else:
        if f1["type"] == "plane":
            plane = f1
            other = f2
        else:
            plane = f2
            other = f1

        normal = plane["params"]["axis_direction"]
        axis = other["params"]["axis_direction"]
        dot_product = abs(np.dot(normal, axis))
        delta = other["params"]["centroid"] - plane["params"]["centroid"]  #for one plane and the other a cylinder/cone, find distance between centroids along the normal and dot product of the axes
        distance = abs(np.dot(delta, normal))
        return dot_product, distance


def compute_fingerprint(face_group):  #runs on every pair of faces, returns i, j : dot, distance for each pair of faces

    fingerprint = {}
    n = len(face_group)
    for i, j in itertools.combinations(range(n), 2):
        fingerprint[(i, j)] = face_relation(face_group[i], face_group[j])

    return fingerprint


# def find_valid_correspondence(selection_faces, target_faces, selection_fingerprint):

    if len(selection_faces) != len(target_faces):  #if length of faces doesnt match, rejected
        return None, None

    def by_type(faces):  #returns ids of each type of face, i.e. ids for planes, ids for cylinders, ids for cones etc
        d = {}
        for i, f in enumerate(faces):
            d.setdefault(f["type"], []).append(i)
        return d

    selected_by_type = by_type(selection_faces)
    target_by_type = by_type(target_faces)

    if set(selected_by_type) != set(target_by_type):  #if selection has 3 planes and 2 cylinders and target has 2 planes and 3 cylinders, this will be rejected
        return None, None
    if any(len(selected_by_type[t]) != len(target_by_type[t]) for t in selected_by_type):
        return None, None

    type_order = list(selected_by_type)
    per_type_permutations = [list(itertools.permutations(target_by_type[t])) for t in type_order]

    for combo in itertools.product(*per_type_permutations):
        mapping = {}
        for t, perm in zip(type_order, combo):
            mapping.update(zip(selected_by_type[t], perm))

        clearances = []
        ok = True
        for sel_i, cand_i in mapping.items():
            sel_face, cand_face = selection_faces[sel_i], target_faces[cand_i]

            if sel_face["type"] in ("cylinder", "cone"):
                if sel_face["params"]["convex"] == cand_face["params"]["convex"]:  #i.e. shaft cannot mate with a shaft
                    ok = False
                    break
                if sel_face["type"] == "cone":
                    if abs(sel_face["params"]["half_angle"] - cand_face["params"]["half_angle"]) > dot_product_tolerance:  #i should probably change dot product tolerance to something else
                        ok = False
                        break
                sel_radius = sel_face["params"]["radius"]
                cand_radius = cand_face["params"]["radius"]
                if cand_face["params"]["convex"]: # convex True implies male (shaft), else female(hole)
                    gap = sel_radius - cand_radius
                else:
                    gap = cand_radius - sel_radius
                if gap < -exact_tolerance or gap > max_clearance:  #in case gap is -ve i.e. shaft might be slightly larger than the hole
                    ok = False
                    break
                clearances.append(gap)

        if not ok:
            continue

        for (i, j), (sel_dot, sel_dist) in selection_fingerprint.items():
            cand_i, cand_j = mapping[i], mapping[j]
            cand_dot, cand_dist = face_relation(target_faces[cand_i], target_faces[cand_j])  #face relation takes in any 2 faces and gives the dot and distance between them
            if abs(sel_dot - cand_dot) > dot_product_tolerance:
                ok = False
                break
            if abs(sel_dist - cand_dist) > distance_tolerance:
                ok = False
                break
            if abs(sel_dot) > parallel_dot:
                clearances.append(cand_dist - sel_dist)
        if not ok:
            continue

        if not clearances:
            fit_type = "unknown"
        elif max(abs(c) for c in clearances) <= exact_tolerance:
            fit_type = "exact"
        else:
            fit_type = "clearance"
        return mapping, fit_type

    return None, None

def find_valid_correspondence(selection_faces, target_faces, selection_fingerprint):  
    if len(selection_faces) != len(target_faces): #if length of faces doesnt match, rejected
        return None, None

    def by_type(faces):  #returns type of each face
        d = {}
        for f in faces:
            d.setdefault(f["type"], []).append(f)
        return d

    selected_by_type = by_type(selection_faces)
    target_by_type = by_type(target_faces)

    if set(selected_by_type) != set(target_by_type):
        return None, None
    if any(len(selected_by_type[t]) != len(target_by_type[t]) for t in selected_by_type):  #if selection has 3 planes and 2 cylinders and target has 2 planes and 3 cylinders, this will be rejected
        return None, None

    clearances = []
    for t in ("cylinder", "cone"):
        if t not in selected_by_type:
            continue
        
        sel_male = sum(1 for f in selected_by_type[t] if f["params"]["convex"])  #enforce that num(males) and num(females) need to be the same, otherwise a hole might mate with a hole, or a rod with a rod
        tgt_female = sum(1 for f in target_by_type[t] if not f["params"]["convex"])
        if sel_male != tgt_female:
            return None, None
        sel_female = sum(1 for f in selected_by_type[t] if not f["params"]["convex"])
        tgt_male = sum(1 for f in target_by_type[t] if f["params"]["convex"])
        if sel_female != tgt_male:
            return None, None
        
        if t not in selected_by_type:
            continue
        
        sel_curved = sorted(selected_by_type[t], key=lambda f: f["params"]["radius"])  
        tgt_curved = sorted(target_by_type[t], key=lambda f: f["params"]["radius"])
        
        for sf, tf in zip(sel_curved, tgt_curved):  #enforce that radii are within tolerance
            sr, tr = sf["params"]["radius"], tf["params"]["radius"]
            if tf["params"]["convex"]:
                gap = sr - tr
            else:
                gap = tr - sr
            if gap < -exact_tolerance or gap > max_clearance:
                return None, None
            clearances.append(gap)

    sel_sorted = sorted(selection_fingerprint.values()) #all dot-distance pairs from selected sorted
    cand_pairs = []
    for i, j in itertools.combinations(range(len(target_faces)), 2):
        cand_pairs.append(face_relation(target_faces[i], target_faces[j]))  #all target faces taken pairwise
    cand_sorted = sorted(cand_pairs)  #target faces sorted

    for (sd, ss), (cd, cs) in zip(sel_sorted, cand_sorted):
        if abs(sd - cd) > dot_product_tolerance:
            return None, None
        if abs(ss - cs) > distance_tolerance:
            return None, None
        if abs(sd) > parallel_dot:
            clearances.append(cs - ss)

    if not clearances:
        fit_type = "unknown"
    elif max(abs(c) for c in clearances) <= exact_tolerance:
        fit_type = "exact"
    else:
        fit_type = "clearance"

    return True, fit_type

def search_solid(selection_faces, sel_fingerprint, solid_faces, adjacency, max_hops):  #candidate generation, generates target faces
    type_counts = {}
    for f in selection_faces:
        type_counts[f["type"]] = type_counts.get(f["type"], 0) + 1
    relevant_types = set(type_counts)

    sel_areas = [f["params"]["area"] for f in selection_faces]
    area_hi = 8.0 * max(sel_areas)
    area_lo = 0.1*min(sel_areas)
    sel_radii = [f["params"]["radius"] for f in selection_faces
                 if "radius" in f["params"]]

    def viable(f):
        if f["type"] not in relevant_types:
            return False
        if f["type"] == "plane":
            return area_lo <= f["params"]["area"] <= area_hi  #to filter any plane based on the max and min areas in the selection of faces, so pretty rough filter but filter nevertheless
        r = f["params"].get("radius")
        if r is None:
            return False
        return any(abs(r - sr) <= max_clearance for sr in sel_radii) #any returns True as long as any of em are True

    tested_combinations = set()
    matches = []
    skipped = 0
    for start in solid_faces:
        if not viable(start):
            continue

        pool_idxs = neighborhood(start["index"], adjacency, max_hops)  #runs BFS
        pool = [f for f in solid_faces if f["index"] in pool_idxs and viable(f)]  #stores the viable faces inside pool

        by_type = {}
        for f in pool:
            by_type.setdefault(f["type"], []).append(f)  #groups faces by type, i.e. cylinder = 2, 1, 4; plane = 5,6
        if any(len(by_type.get(t, [])) < count for t, count in type_counts.items()):  #makes sure that the type counts that we just grouped are more than the ones in selection faces, i.e. if selection faces has 3 cylinders, pool should have 4
            continue

        type_order = list(type_counts)  #list of types 
        per_type_choices = [list(itertools.combinations(by_type[t], type_counts[t]))
                            for t in type_order]  #for each type, pick the required number of faces from the pool, each type will be a different list?

        total = 1
        for choices in per_type_choices:
            total *= len(choices)  #counts combinations in to be tested
        
        if total > 20000:
            print(f"[debug] start {start['index']}: {total} combos, skipping") #just to speed things up, if all possible choices exceed 20k, skip
            skipped+=1
            continue

        for combo_parts in itertools.product(*per_type_choices):  #generates all combinations of plane groups + cylinder groups
            combo = [f for part in combo_parts for f in part]  

            if not any(f["index"] == start["index"] for f in combo):  #if combination doesnt contain start, no point checking here because this combination will get generated again when starting from one of the other included faces
                continue

            key = frozenset(f["index"] for f in combo)  #frozenset so that the same faces arent computed again
            if key in tested_combinations:
                continue
            tested_combinations.add(key)

            # mapping, fit_type = find_valid_correspondence(
            #     selection_faces, combo, sel_fingerprint)
            # if mapping is not None:
            #     matches.append((combo, mapping, fit_type))
                
            result, fit_type = find_valid_correspondence(selection_faces, combo, sel_fingerprint)  #find matches
            if result is not None:
                matches.append((combo, None, fit_type))
    if skipped:
        print(f"[debug] skipped {skipped} start faces (combo cap)")

    return matches

def find_all_mates(selection_part_path, selection_face_ids, assembly_path, max_hops=max_hops):

    sel_solids = load_step(selection_part_path)  #load step file
    
    if not sel_solids:
        raise RuntimeError("No solids found in selection part file")
    elif len(sel_solids) > 1:
        raise RuntimeError(f"Selection file has {len(sel_solids)} solids, expected 1")
    
    sel_all_faces = get_faces(sel_solids[0])

    selection_faces = [f for f in sel_all_faces if f["index"] in selection_face_ids]
    if len(selection_faces) != len(selection_face_ids):
        raise RuntimeError("Some selected face IDs were not found")
    bad = [f["index"] for f in selection_faces if f["type"] not in SUPPORTED]
    if bad:
        raise RuntimeError(f"Unsupported surface type on selected faces {bad}")

    sel_fingerprint = compute_fingerprint(selection_faces)
    sel_volume = solid_volume(sel_solids[0])

    asm_solids = load_step(assembly_path) #load assy
    print(f"[debug] sel_all_faces count: {len(sel_all_faces)}")
    print(f"[debug] assembly has {len(asm_solids)} solid(s)")
    all_matches = []

    for part_index, solid in enumerate(asm_solids):
        if abs(solid_volume(solid) - sel_volume) < 1e-6 * max(1.0, sel_volume): #used by find_all_mates() to skip the selection part's instance if it is included (however, risk of different parts with same volume or the same part repeated being excluded)
            print(f"[debug] part_index {part_index}: SKIPPED (same volume as selection part)") 
            continue

        solid_faces = get_faces(solid)
        adj = adjacency(solid, solid_faces)
        matches = search_solid(selection_faces, sel_fingerprint, solid_faces, adj, max_hops)

        for combo, mapping, fit_type in matches:
            all_matches.append({
                "part_index": part_index,
                "face_ids": [f["index"] for f in combo],
                "fit_type": fit_type,
            })

    print(f"[debug] total matches found: {len(all_matches)}")
    for m in all_matches:
        print(f"  {m}")
    return all_matches





def matches_to_highlight_dict(matches):
    """Convert find_all_mates() output into {part_index: set(face_index)},
    merging face ids across multiple matches on the same part."""
    by_part = {}
    for m in matches:
        by_part.setdefault(m["part_index"], set()).update(m["face_ids"])
    return by_part


def render_highlighted(step_path, highlighted_by_solid, output_path,
                        default_rgb=(0.6, 0.6, 0.6), highlight_rgb=(0.0, 0.8, 0.2),
                        display=None, only_matched=False, transparency=0.75):


    if display is None:
        display, _, _, _ = init_display()
    display.EraseAll()

    default_color = Quantity_Color(*default_rgb, Quantity_TOC_RGB)
    highlight_color = Quantity_Color(*highlight_rgb, Quantity_TOC_RGB)

    solid_indices = set(highlighted_by_solid.keys()) if only_matched else None
    for s_idx, solid in enumerate(load_step(step_path)):
        if solid_indices is not None and s_idx not in solid_indices:
            continue
        highlighted_ids = highlighted_by_solid.get(s_idx, set())

        # base solid in default color, transparent so interior matches
        # (socket walls, bores) are not hidden behind the outer shell
        ais = display.DisplayShape(solid, color=default_color, update=False)
        if transparency:
            handle = ais[0] if isinstance(ais, (list, tuple)) else ais
            display.Context.SetTransparency(handle, transparency, False)

        # re-draw just the highlighted faces on top, in the highlight color
        if highlighted_ids:
            for f in get_faces(solid):
                if f["index"] in highlighted_ids:
                    display.DisplayColoredShape(f["face"], color=highlight_color, update=False)

    display.FitAll()
    display.View.Dump(output_path)
    return display

samples = {"rod": ["stepfiles/rod.step", "stepfiles/abracket.step", [0], "images/rod_selection.png", "images/rod_matches.png"],
           "keyshaft": ["stepfiles/key_shaft.step", "stepfiles/hub.step", [0, 2, 5], "images/keyshaft_selection.png", "images/keyshaft_matches.png"],
           "keyshaft_clear": ["stepfiles/key_shaft.step", "stepfiles/hub2clearance.step", [0, 2, 5], "images/keyshaft_selection_clearance.png", "images/keyshaft_matches_clearance.png"],
           "lathe" : ["stepfiles/only_key.step", "stepfiles/WM290V_RedrawStart.stp", [1, 3, 14, 16], "images/lathe_selection.png", "images/lathe_matches.png"],
           "onlychuck":["stepfiles/only_key.step", "stepfiles/onlychuck.step", [1, 3, 14, 16], "images/chuck_selection.png", "images/chuck_matches.png"]
               }


if __name__ == "__main__":

    test = samples["lathe"]

    matches = find_all_mates(
        selection_part_path=test[0],
        selection_face_ids=test[2],
        assembly_path=test[1],
        max_hops=max_hops,
    )

    highlight_dict = matches_to_highlight_dict(matches)
    print(f"matched {len(highlight_dict)} distinct part(s):")
    for part_idx, face_ids in highlight_dict.items():
        print(f"  part_index {part_idx}: faces {sorted(face_ids)}")

    from OCC.Display.SimpleGui import init_display
    display, _, _, _ = init_display()
    display = render_highlighted(test[0], {0: set(test[2])},
                                 test[3], display=display, transparency=0.0)
    display = render_highlighted(test[1], highlight_dict, test[4],
                                 display=display, only_matched=False)
    display.FitAll()
    display.Repaint()
    ok = display.View.Dump(test[4])
    print(f"DUMP : {ok}")



def pick_faces(path, solid_index=0):
    from OCC.Display.SimpleGui import init_display
    solid = load_step(path)[solid_index]
    faces = get_faces(solid)
    display, start_loop, _, _ = init_display()
    picked = []

    def on_select(shapes, x, y):
        for s in shapes:
            for f in faces:
                if f["face"].IsSame(s):
                    picked.append(f["index"])
                    print(f"picked {f['index']} ({f['type']}, area={f['params']['area']:.2f})")

    display.SetSelectionModeFace()
    display.register_select_callback(on_select)
    display.DisplayShape(solid, update=True)
    start_loop()
    print(f"selection_face_ids = {sorted(set(picked))}")
    return sorted(set(picked))
# pick_faces("stepfiles/only_key.step")

