#!/usr/bin/env python3
"""Visualize the Cu-Steel thermal contact simulation with PyVista.

Two separate meshes are rendered side by side with a shared colormap so the
temperature jump at the interface (x = 0.05 m) is clearly visible. Also emits
a centerline (y = 0.015 m) profile plot that makes the Rc-induced discontinuity
explicit.

Inputs (from results/):
    mesh_cu.msh, mesh_st.msh
    T1_t{60,150,300}.txt, T2_t{60,150,300}.txt
    centerline.csv

Outputs (to results/):
    fig_t060.png, fig_t150.png, fig_t300.png      (2D fields)
    fig_centerline.png                             (T vs x line plots)
"""

import os
import sys
from pathlib import Path

# PyVista must be configured off-screen BEFORE import.
os.environ["PYVISTA_OFF_SCREEN"] = "true"

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyvista as pv

pv.OFF_SCREEN = True

# --- Paths -----------------------------------------------------------------
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = ROOT / "results"
sys.path.insert(0, str(ROOT / ".claude" / "skills" / "fem-analysis" / "scripts"))
from ff_mesh_reader import read_freefem_msh, to_pyvista_grid  # noqa: E402


# --- Helpers ---------------------------------------------------------------
def load_block(mesh_file: Path, sol_file: Path, name: str = "T"):
    mesh = read_freefem_msh(str(mesh_file))
    T = np.loadtxt(sol_file)
    assert T.shape[0] == mesh["nv"], (
        f"DOF mismatch: {sol_file.name} has {T.shape[0]} but mesh has {mesh['nv']} vertices"
    )
    grid = to_pyvista_grid(mesh, {name: T})
    return mesh, grid, T


def render_field(t_label: int, out_png: Path, vmin: float, vmax: float):
    cu_mesh = RESULTS / "mesh_cu.msh"
    st_mesh = RESULTS / "mesh_st.msh"
    cu_sol = RESULTS / f"T1_t{t_label}.txt"
    st_sol = RESULTS / f"T2_t{t_label}.txt"

    _, cu_grid, Tcu = load_block(cu_mesh, cu_sol, "T")
    _, st_grid, Tst = load_block(st_mesh, st_sol, "T")

    # Narrow window matching the 100x30 mm aspect ratio (~3.3:1), plus margin for
    # the color bar. Keeps the blocks large on the canvas.
    plotter = pv.Plotter(off_screen=True, window_size=(1500, 520))
    plotter.background_color = "white"

    scalar_bar_args = dict(
        title="T [K]",
        vertical=True,
        position_x=0.90,
        position_y=0.12,
        height=0.76,
        width=0.035,
        n_labels=7,
        title_font_size=18,
        label_font_size=16,
        fmt="%.1f",
        color="black",
    )

    plotter.add_mesh(
        cu_grid,
        scalars="T",
        cmap="coolwarm",
        clim=(vmin, vmax),
        show_edges=False,
        scalar_bar_args=scalar_bar_args,   # attach bar to this mapper
    )
    plotter.add_mesh(
        st_grid,
        scalars="T",
        cmap="coolwarm",
        clim=(vmin, vmax),
        show_edges=False,
        show_scalar_bar=False,             # same color mapping — no second bar
    )

    # Interface line at x = 0.05 m, y ∈ [0, 0.03]
    interface = pv.Line((0.05, 0.0, 0.0), (0.05, 0.03, 0.0), resolution=10)
    plotter.add_mesh(interface, color="black", line_width=3.0)

    plotter.add_text(
        f"t = {t_label} s    Cu (left)   |   Steel (right)",
        position="upper_edge",
        font_size=14,
        color="black",
    )
    plotter.add_text(
        f"Cu: [{Tcu.min():.2f}, {Tcu.max():.2f}] K     "
        f"Steel: [{Tst.min():.2f}, {Tst.max():.2f}] K",
        position="lower_edge",
        font_size=11,
        color="black",
    )

    plotter.view_xy()
    plotter.camera.parallel_projection = True
    plotter.reset_camera()
    plotter.camera.zoom(1.6)

    plotter.screenshot(str(out_png), transparent_background=False)
    plotter.close()
    print(f"  wrote {out_png.name}")


def render_centerline(out_png: Path):
    import csv

    centerline = RESULTS / "centerline.csv"
    data = {}   # time -> side -> list[(x,T)]
    with open(centerline, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            t = float(row["time"])
            data.setdefault(t, {"cu": [], "st": []})[row["side"]].append(
                (float(row["x"]), float(row["T"]))
            )

    wanted = [60.0, 150.0, 300.0]
    fig, ax = plt.subplots(figsize=(9, 5.5), dpi=130)

    colors = {60.0: "#d62728", 150.0: "#2ca02c", 300.0: "#1f77b4"}
    jump_rows = []
    for t in wanted:
        if t not in data:
            continue
        cu = sorted(data[t]["cu"])
        st = sorted(data[t]["st"])
        xcu = np.array([p[0] for p in cu]) * 1000  # mm
        Tcu = np.array([p[1] for p in cu])
        xst = np.array([p[0] for p in st]) * 1000
        Tst = np.array([p[1] for p in st])

        ax.plot(xcu, Tcu, "-", color=colors[t], lw=2, label=f"t = {int(t)} s  (Cu)")
        ax.plot(xst, Tst, "--", color=colors[t], lw=2, label=f"t = {int(t)} s  (Steel)")

        # Mark the jump visually with a short vertical bar
        ax.vlines(50.0, Tst[0], Tcu[-1], colors=colors[t], lw=1.8, alpha=0.8)
        ax.plot([50.0], [Tcu[-1]], "o", color=colors[t], ms=5)
        ax.plot([50.0], [Tst[0]], "o", color=colors[t], ms=5, mfc="white")

        jump_rows.append((int(t), Tcu[-1], Tst[0], Tcu[-1] - Tst[0]))

    ax.axvline(50.0, color="k", lw=1, ls=":", alpha=0.6, label="Contact x = 50 mm")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("T [K]")
    ax.set_title("Centerline temperature profile (y = 15 mm) — Rc = 1e-4 m²K/W")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower left", fontsize=8, ncol=2)

    # Jump table inset, positioned top-right so it never overlaps the curves
    # (all curves are high on the RHS of the plot near y ≈ 370 K).
    rows = "\n".join(
        f"t={t:>3d}s:  T_Cu={tcu:6.2f}  T_St={tst:6.2f}  ΔT={dT:6.3f} K"
        for t, tcu, tst, dT in jump_rows
    )
    ax.text(
        0.98, 0.98, rows,
        transform=ax.transAxes,
        fontsize=9,
        va="top",
        ha="right",
        family="monospace",
        bbox=dict(facecolor="white", edgecolor="gray", alpha=0.9, pad=6),
    )

    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)
    print(f"  wrote {out_png.name}")


# --- Main ------------------------------------------------------------------
def main():
    # Shared colormap range across all time snapshots
    all_T = []
    for t in (60, 150, 300):
        all_T.append(np.loadtxt(RESULTS / f"T1_t{t}.txt"))
        all_T.append(np.loadtxt(RESULTS / f"T2_t{t}.txt"))
    vmin = min(a.min() for a in all_T)
    vmax = max(a.max() for a in all_T)
    print(f"Shared colormap range: T ∈ [{vmin:.2f}, {vmax:.2f}] K")

    for t in (60, 150, 300):
        render_field(t, RESULTS / f"fig_t{t:03d}.png", vmin, vmax)

    render_centerline(RESULTS / "fig_centerline.png")


if __name__ == "__main__":
    main()
