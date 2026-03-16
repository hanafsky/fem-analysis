#!/usr/bin/env python3
"""
FEM result visualization with PyVista (off-screen rendering).

Supports:
- Scalar field contour plots (temperature, Von Mises, etc.)
- Deformed mesh overlay (structural analysis)
- Multi-panel comparison views
- Automatic colormap and scale bar configuration
"""
import os
os.environ["PYVISTA_OFF_SCREEN"] = "true"

import numpy as np
import pyvista as pv
pv.OFF_SCREEN = True

from ff_mesh_reader import read_freefem_msh, to_pyvista_grid, make_deformed_grid


def plot_scalar_field(mesh_file: str, solution_file: str, output_png: str,
                      scalar_name: str = "Solution",
                      cmap: str = "viridis", show_edges: bool = True,
                      window_size: tuple = (1000, 700)):
    """Plot a scalar field on the mesh.
    
    Args:
        mesh_file: Path to FreeFem++ .msh file
        solution_file: Path to text file with nodal values
        output_png: Output image path
        scalar_name: Label for the scalar bar
        cmap: Matplotlib colormap name
        show_edges: Whether to show mesh edges
        window_size: Image dimensions (width, height)
    """
    mesh_data = read_freefem_msh(mesh_file)
    solution = np.loadtxt(solution_file)
    
    grid = to_pyvista_grid(mesh_data, {scalar_name: solution})
    
    plotter = pv.Plotter(off_screen=True, window_size=window_size)
    plotter.add_mesh(grid, scalars=scalar_name, cmap=cmap,
                     show_edges=show_edges, edge_color="gray", line_width=0.3)
    plotter.add_scalar_bar(title=scalar_name)
    plotter.view_xy()
    plotter.screenshot(output_png)
    
    print(f"Saved: {output_png}")
    print(f"  {scalar_name} range: [{solution.min():.6e}, {solution.max():.6e}]")
    return solution.min(), solution.max()


def plot_deformed_stress(mesh_file: str, u1_file: str, u2_file: str,
                         stress_file: str, output_png: str,
                         stress_name: str = "Von Mises [Pa]",
                         scale: float = None,
                         cmap: str = "jet",
                         window_size: tuple = (1200, 500)):
    """Plot deformed mesh with stress contour overlay.
    
    Args:
        mesh_file: Path to FreeFem++ .msh file
        u1_file, u2_file: Displacement component files
        stress_file: Stress field file
        output_png: Output image path
        stress_name: Stress field label
        scale: Deformation scale (None for auto)
        cmap: Colormap
        window_size: Image dimensions
    """
    mesh_data = read_freefem_msh(mesh_file)
    u1 = np.loadtxt(u1_file)
    u2 = np.loadtxt(u2_file)
    stress = np.loadtxt(stress_file)
    
    # Original mesh (wireframe reference)
    grid_orig = to_pyvista_grid(mesh_data)
    
    # Deformed mesh with stress
    grid_def = make_deformed_grid(mesh_data, u1, u2, scale=scale,
                                  point_data={stress_name: stress})
    
    actual_scale = grid_def.field_data["deformation_scale"][0]
    
    plotter = pv.Plotter(off_screen=True, window_size=window_size)
    plotter.add_mesh(grid_orig, style="wireframe", color="gray",
                     line_width=0.5, opacity=0.3)
    plotter.add_mesh(grid_def, scalars=stress_name, cmap=cmap, show_edges=False)
    plotter.add_scalar_bar(title=stress_name)
    plotter.add_text(f"Deformation x{actual_scale:.1f}", position="upper_left",
                     font_size=10, color="black")
    plotter.view_xy()
    plotter.screenshot(output_png)
    
    print(f"Saved: {output_png}")
    print(f"  Deformation scale: {actual_scale:.1f}x")
    print(f"  Max |u1|: {np.abs(u1).max():.6e}")
    print(f"  Max |u2|: {np.abs(u2).max():.6e}")
    print(f"  {stress_name} range: [{stress.min():.6e}, {stress.max():.6e}]")


def plot_multi_panel(mesh_file: str, data_files: dict, output_png: str,
                     cmaps: dict = None, window_size: tuple = (1600, 500)):
    """Create a multi-panel comparison plot.
    
    Args:
        mesh_file: Path to FreeFem++ .msh
        data_files: dict of {panel_name: solution_file_path}
        output_png: Output image path
        cmaps: dict of {panel_name: cmap_name} (optional)
        window_size: Image dimensions
    """
    mesh_data = read_freefem_msh(mesh_file)
    n_panels = len(data_files)
    
    if cmaps is None:
        cmaps = {name: "viridis" for name in data_files}
    
    plotter = pv.Plotter(off_screen=True, shape=(1, n_panels), window_size=window_size)
    
    for idx, (name, filepath) in enumerate(data_files.items()):
        plotter.subplot(0, idx)
        solution = np.loadtxt(filepath)
        grid = to_pyvista_grid(mesh_data, {name: solution})
        plotter.add_mesh(grid, scalars=name, cmap=cmaps.get(name, "viridis"),
                         show_edges=True, edge_color="gray", line_width=0.3)
        plotter.add_scalar_bar(title=name)
        plotter.add_text(name, position="upper_edge", font_size=10)
        plotter.view_xy()
    
    plotter.screenshot(output_png)
    print(f"Saved multi-panel: {output_png} ({n_panels} panels)")


def export_for_threejs(mesh_file: str, solution_file: str, 
                       output_json: str, scalar_name: str = "value"):
    """Export mesh + solution as JSON for Three.js React artifact.
    
    Args:
        mesh_file: FreeFem++ .msh
        solution_file: Nodal values
        output_json: Output JSON path
        scalar_name: Field name
    """
    import json
    
    mesh_data = read_freefem_msh(mesh_file)
    solution = np.loadtxt(solution_file)
    
    data = {
        "vertices": mesh_data["points"][:, :2].tolist(),
        "triangles": mesh_data["triangles"].tolist(),
        "scalars": {
            "name": scalar_name,
            "values": solution.tolist(),
            "min": float(solution.min()),
            "max": float(solution.max()),
        }
    }
    
    with open(output_json, "w") as f:
        json.dump(data, f)
    
    print(f"Exported: {output_json} ({mesh_data['nv']} vertices, {mesh_data['nt']} triangles)")


if __name__ == "__main__":
    import sys
    print("FEM Visualization Pipeline")
    print("Usage: Import functions or run with specific arguments")
    print("  plot_scalar_field(mesh, solution, output)")
    print("  plot_deformed_stress(mesh, u1, u2, stress, output)")
    print("  plot_multi_panel(mesh, {name: file, ...}, output)")
    print("  export_for_threejs(mesh, solution, output)")
