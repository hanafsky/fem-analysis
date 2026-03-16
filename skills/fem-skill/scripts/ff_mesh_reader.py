#!/usr/bin/env python3
"""
FreeFEM++ .msh file reader → PyVista UnstructuredGrid converter.

FreeFem++ .msh format:
  Line 1: nv nt ne  (num_vertices, num_triangles, num_edges)
  Next nv lines: x y boundary_label
  Next nt lines: v1 v2 v3 region_label  (1-indexed)
  Next ne lines: v1 v2 boundary_label   (1-indexed)
"""
import numpy as np

try:
    import pyvista as pv
    HAS_PYVISTA = True
except ImportError:
    HAS_PYVISTA = False


def read_freefem_msh(filename: str) -> dict:
    """Read a FreeFem++ .msh file.
    
    Returns:
        dict with keys: points, triangles, edges, boundary_labels, region_labels
    """
    with open(filename) as f:
        nv, nt, ne = map(int, f.readline().split())
        
        points = np.zeros((nv, 3))
        boundary_labels_v = np.zeros(nv, dtype=int)
        for i in range(nv):
            parts = f.readline().split()
            points[i, 0] = float(parts[0])
            points[i, 1] = float(parts[1])
            boundary_labels_v[i] = int(parts[2])
        
        triangles = np.zeros((nt, 3), dtype=int)
        region_labels = np.zeros(nt, dtype=int)
        for i in range(nt):
            parts = f.readline().split()
            triangles[i] = [int(parts[0])-1, int(parts[1])-1, int(parts[2])-1]
            region_labels[i] = int(parts[3])
        
        edges = np.zeros((ne, 2), dtype=int)
        boundary_labels_e = np.zeros(ne, dtype=int)
        for i in range(ne):
            parts = f.readline().split()
            edges[i] = [int(parts[0])-1, int(parts[1])-1]
            boundary_labels_e[i] = int(parts[2])
    
    return {
        "points": points,
        "triangles": triangles,
        "edges": edges,
        "nv": nv, "nt": nt, "ne": ne,
        "vertex_boundary_labels": boundary_labels_v,
        "region_labels": region_labels,
        "edge_boundary_labels": boundary_labels_e,
    }


def to_pyvista_grid(mesh_data: dict, point_data: dict = None) -> "pv.UnstructuredGrid":
    """Convert FreeFem++ mesh to PyVista UnstructuredGrid.
    
    Args:
        mesh_data: Output of read_freefem_msh()
        point_data: dict of {name: np.array} for scalar fields on vertices
    
    Returns:
        PyVista UnstructuredGrid
    """
    if not HAS_PYVISTA:
        raise ImportError("PyVista is required for grid conversion")
    
    points = mesh_data["points"]
    triangles = mesh_data["triangles"]
    nt = mesh_data["nt"]
    
    cells = np.hstack([np.full((nt, 1), 3, dtype=int), triangles]).ravel()
    cell_types = np.full(nt, 5, dtype=np.uint8)  # VTK_TRIANGLE = 5
    
    grid = pv.UnstructuredGrid(cells, cell_types, points)
    
    if point_data:
        for name, data in point_data.items():
            grid.point_data[name] = data
    
    return grid


def make_deformed_grid(mesh_data: dict, u1: np.ndarray, u2: np.ndarray,
                       scale: float = None, point_data: dict = None) -> "pv.UnstructuredGrid":
    """Create a deformed mesh for structural visualization.
    
    Args:
        mesh_data: Output of read_freefem_msh()
        u1, u2: Displacement fields (x, y)
        scale: Deformation magnification. If None, auto-calculated.
        point_data: Additional scalar fields
    
    Returns:
        PyVista UnstructuredGrid with deformed coordinates
    """
    points = mesh_data["points"].copy()
    
    # Auto-scale: make max deformation ~10% of characteristic length
    if scale is None:
        L_char = max(np.ptp(points[:, 0]), np.ptp(points[:, 1]))
        max_disp = max(np.abs(u1).max(), np.abs(u2).max())
        if max_disp > 0:
            scale = 0.1 * L_char / max_disp
        else:
            scale = 1.0
    
    deformed = points.copy()
    deformed[:, 0] += u1 * scale
    deformed[:, 1] += u2 * scale
    
    triangles = mesh_data["triangles"]
    nt = mesh_data["nt"]
    cells = np.hstack([np.full((nt, 1), 3, dtype=int), triangles]).ravel()
    cell_types = np.full(nt, 5, dtype=np.uint8)
    
    grid = pv.UnstructuredGrid(cells, cell_types, deformed)
    
    if point_data:
        for name, data in point_data.items():
            grid.point_data[name] = data
    
    # Store scale factor as field data
    grid.field_data["deformation_scale"] = np.array([scale])
    
    return grid


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python ff_mesh_reader.py <mesh.msh> [solution.txt]")
        sys.exit(1)
    
    mesh_data = read_freefem_msh(sys.argv[1])
    print(f"Vertices: {mesh_data['nv']}")
    print(f"Triangles: {mesh_data['nt']}")
    print(f"Edges: {mesh_data['ne']}")
    print(f"X range: [{mesh_data['points'][:,0].min():.4f}, {mesh_data['points'][:,0].max():.4f}]")
    print(f"Y range: [{mesh_data['points'][:,1].min():.4f}, {mesh_data['points'][:,1].max():.4f}]")
    
    if len(sys.argv) >= 3:
        sol = np.loadtxt(sys.argv[2])
        print(f"Solution: {sol.shape[0]} values, range [{sol.min():.6e}, {sol.max():.6e}]")
