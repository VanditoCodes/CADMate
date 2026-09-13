from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeBox
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCC.Core.gp import gp_Pnt, gp_Ax2, gp_Dir, gp_Trsf, gp_Vec
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.Core.STEPControl import STEPControl_Writer, STEPControl_AsIs

def moved(shape, dx, dy, dz):
    t = gp_Trsf()
    t.SetTranslation(gp_Vec(dx, dy, dz))
    return BRepBuilderAPI_Transform(shape, t, True).Shape()

def write_step(shapes, path):
    w = STEPControl_Writer()
    for s in shapes:
        w.Transfer(s, STEPControl_AsIs)
    w.Write(path)
    

shaft = BRepPrimAPI_MakeCylinder(10.0, 60.0).Shape()
key = moved(BRepPrimAPI_MakeBox(6.0, 3.0, 40.0).Shape(), -3.0, 10.0, 10.0)
key = moved(BRepPrimAPI_MakeBox(6.0, 4.0, 40.0).Shape(), -3.0, 8.0, 10.0)
keyed_shaft = BRepAlgoAPI_Fuse(shaft, key).Shape()

write_step([keyed_shaft], "keyed_shaft.step")

BORE_CLEAR = 0.1
KEY_CLEAR  = 0.2
KEY_DEPTH_EXTRA = 0.5

hub = moved(BRepPrimAPI_MakeBox(50.0, 50.0, 25.0).Shape(), -25.0, -25.0, 15.0)

bore = BRepPrimAPI_MakeCylinder(10.0 + BORE_CLEAR, 100.0).Shape()

slot = moved(BRepPrimAPI_MakeBox(6.0 + 2 * KEY_CLEAR,
                                 3.0 + KEY_DEPTH_EXTRA,
                                 100.0).Shape(),
             -(3.0 + KEY_CLEAR),
             8.0,
             0.0)

hub = BRepAlgoAPI_Cut(hub, bore).Shape()
hub = BRepAlgoAPI_Cut(hub, slot).Shape()
hub2 = moved(hub, 0.0, 0.0, 25.0)
write_step([shaft, key, hub, hub2], "keyed_assy.step")
write_step([hub], "hub.step")