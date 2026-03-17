# 3D Gmsh + FreeFEM++ Reference

## Overview

Rules and patterns for 3D FEM analysis using Gmsh (OCC kernel) for mesh generation and FreeFEM++ (`gmshload3`) for solving. These were learned from debugging real failures — follow them strictly.

---

## Rule 1: Use OCC 3D Primitives + BooleanFragments (Not Extrude)

**Problem:** `Extrude{} Layers{}` produces unstable surface IDs. BoundingBox queries pick up wrong faces → `gmshload3` fails with "Border element not in mesh".

**Solution:** Build geometry from OCC solid primitives and split with BooleanFragments.

```geo
SetFactory("OpenCASCADE");

// Create solid primitives
Box(1) = {0, 0, 0,       L, W, h1};       // Layer 1
Box(2) = {0, 0, h1,      L, W, h2};       // Layer 2
Box(3) = {0, 0, h1+h2,   L, W, h3};       // Layer 3

// Fragment to create conformal mesh at interfaces
BooleanFragments{ Volume{1}; Delete; }{ Volume{2,3}; Delete; }
```

**Why BooleanFragments?** It splits volumes at shared surfaces, producing conformal meshes with stable entity IDs.

---

## Rule 2: Never Assign Internal Interfaces to Physical Surface

**Problem:** Interior faces shared by two tetrahedra (e.g., layer interface at z=h1) crash `gmshload3` if defined as `Physical Surface`.

**Cause:** `gmshload3` treats Physical Surface elements as boundary elements. An internal face belongs to two volumes, so FreeFEM++ cannot classify it as a boundary → "Border element not in mesh".

**Solution:** Only define Physical Surface for **external** boundary faces. Omit all internal interface surfaces.

```geo
// CORRECT: only external boundaries
Physical Surface("bottom")  = {bottom_face};
Physical Surface("top")     = {top_face};
Physical Surface("sides")   = {side_faces[]};

// WRONG: do NOT include internal interfaces
// Physical Surface("interface_12") = {interface_face};  // CRASHES gmshload3
```

**To extract interface data**, use Python post-processing (see Rule 5).

---

## Rule 3: Classify Volumes by Centroid z-Coordinate

**Problem:** `Volume In BoundingBox` returns all volumes that *intersect* the bounding box, not those *contained* in it. Adjacent layers sharing a boundary get duplicated.

**Solution:** Compute each volume's centroid and classify based on that.

```geo
// Get all volumes after BooleanFragments
vols[] = Volume{:};

For i In {0:#vols[]-1}
    bb() = BoundingBox Volume{vols[i]};
    // bb = {xmin, ymin, zmin, xmax, ymax, zmax}
    zcenter = (bb(2) + bb(5)) / 2;

    If (zcenter < h1)
        layer1_vols[] += vols[i];
    ElseIf (zcenter < h1 + h2)
        layer2_vols[] += vols[i];
    Else
        layer3_vols[] += vols[i];
    EndIf
EndFor

Physical Volume("layer1") = {layer1_vols[]};
Physical Volume("layer2") = {layer2_vols[]};
Physical Volume("layer3") = {layer3_vols[]};
```

---

## Rule 4: Use `Th(i).x` for 3D Coordinate Access

**Problem:** `Wh(i).x` (FE space DOF coordinate) works in 2D but fails to compile with `mesh3`.

**Solution:** Always use mesh vertex coordinates directly.

```freefem
mesh3 Th = gmshload3("mesh.msh");
fespace Vh(Th, P1);
Vh u;  // solve ...

// Output node coordinates and solution
{
    ofstream f("results/solution.txt");
    for (int i = 0; i < Vh.ndof; i++)
        f << Th(i).x << " " << Th(i).y << " " << Th(i).z << " " << u[][i] << endl;
}
```

**Note:** For P1 elements, `ndof == nv` (number of vertices), so `Th(i).x` maps 1:1 to DOFs.

---

## Rule 5: Extract Interface Data via Python Post-Processing

**Problem:** Marking internal interfaces in the mesh (Physical Surface) causes gmshload3 errors (Rule 2).

**Solution:** Export full-field data from FreeFEM++, then filter in Python.

```python
import numpy as np

# Load full-field data exported from FreeFEM++
data = np.loadtxt("results/solution.txt")  # x, y, z, value
x, y, z, val = data[:, 0], data[:, 1], data[:, 2], data[:, 3]

# Extract interface at z ≈ target
z_interface = 0.008  # example: 8mm
tol = 1e-4
mask = np.abs(z - z_interface) < tol

x_if, y_if, val_if = x[mask], y[mask], val[mask]

# Visualize or analyze interface data
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

tri = Triangulation(x_if, y_if)
fig, ax = plt.subplots()
tc = ax.tricontourf(tri, val_if, levels=20, cmap="RdYlBu_r")
plt.colorbar(tc)
ax.set_title(f"Interface at z = {z_interface}")
fig.savefig("results/interface.png", dpi=150, bbox_inches="tight")
```

---

## Typical 3D Workflow Summary

1. **Gmsh `.geo`**: OCC primitives → BooleanFragments → centroid-based Physical Volume → external-only Physical Surface
2. **Generate mesh**: `gmsh -3 model.geo -o mesh.msh -format msh2`
3. **FreeFEM++ `.edp`**: `gmshload3("mesh.msh")` → solve → export with `Th(i).x/y/z`
4. **Python post-processing**: Load results → filter interfaces by z ≈ target → visualize
