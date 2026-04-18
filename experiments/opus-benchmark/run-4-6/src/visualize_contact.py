#!/usr/bin/env python3
"""Render Cu+Steel contact thermal snapshots with PyVista."""
import os
import sys
from pathlib import Path

os.environ["PYVISTA_OFF_SCREEN"] = "true"

import numpy as np
import pyvista as pv
pv.OFF_SCREEN = True

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / ".claude" / "skills" / "fem-analysis" / "scripts"))
from ff_mesh_reader import read_freefem_msh, to_pyvista_grid  # noqa: E402

RESULTS = REPO / "results"
SNAPS = [60, 150, 300]
VMIN, VMAX = 293.0, 373.0


def load_snapshot(tag: str):
    cu = read_freefem_msh(str(RESULTS / "mesh_cu.msh"))
    st = read_freefem_msh(str(RESULTS / "mesh_st.msh"))
    T_cu = np.loadtxt(RESULTS / f"T_cu_{tag}.txt")
    T_st = np.loadtxt(RESULTS / f"T_st_{tag}.txt")
    g_cu = to_pyvista_grid(cu, {"T": T_cu})
    g_st = to_pyvista_grid(st, {"T": T_st})
    return g_cu, g_st


def render_snapshots():
    # Render each snapshot as its own image, then composite vertically via matplotlib.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg

    tmp_pngs = []
    for t in SNAPS:
        tag = f"{t:03d}"
        g_cu, g_st = load_snapshot(tag)
        plotter = pv.Plotter(off_screen=True, window_size=(1600, 420))
        plotter.add_mesh(
            g_cu, scalars="T", cmap="inferno", clim=(VMIN, VMAX),
            show_edges=False,
            scalar_bar_args={
                "title": "T [K]",
                "n_labels": 6,
                "vertical": True,
                "position_x": 0.90,
                "position_y": 0.10,
                "width": 0.04,
                "height": 0.8,
                "label_font_size": 18,
                "title_font_size": 20,
                "color": "black",
            },
        )
        plotter.add_mesh(
            g_st, scalars="T", cmap="inferno", clim=(VMIN, VMAX),
            show_edges=False, show_scalar_bar=False,
        )
        interface = pv.Line((0.05, 0.0, 0.0), (0.05, 0.030, 0.0))
        plotter.add_mesh(interface, color="cyan", line_width=3)
        plotter.view_xy()
        plotter.camera.zoom(1.7)
        plotter.set_background("white")
        tmp = RESULTS / f"_snap_{tag}.png"
        plotter.screenshot(str(tmp))
        plotter.close()
        tmp_pngs.append((t, tmp))

    fig, axes = plt.subplots(3, 1, figsize=(14, 8.5))
    for ax, (t, png) in zip(axes, tmp_pngs):
        ax.imshow(mpimg.imread(png))
        ax.set_title(
            f"t = {t} s   (Cu left   |   contact surface (cyan)   |   Steel right)",
            fontsize=13,
        )
        ax.axis("off")
    fig.suptitle(
        "Transient heat conduction, Cu-Steel with contact Rc = 1×10⁻⁴ m²K/W",
        fontsize=15, y=0.995,
    )
    fig.tight_layout()
    out = RESULTS / "contact_thermal_snapshots.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    for _, p in tmp_pngs:
        p.unlink(missing_ok=True)
    print(f"Saved: {out}")


def render_interface_profile():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = {60: "#1f77b4", 150: "#ff7f0e", 300: "#2ca02c"}
    for t in SNAPS:
        data = np.loadtxt(RESULTS / f"iface_{t:03d}.txt")
        y, Tcu, Tst = data[:, 0] * 1000, data[:, 1], data[:, 2]
        ax.plot(y, Tcu, "-", color=colors[t], lw=2, label=f"Cu  side  t={t}s")
        ax.plot(y, Tst, "--", color=colors[t], lw=2, label=f"Steel side t={t}s")
    ax.set_xlabel("y [mm] along interface (x = 50 mm)")
    ax.set_ylabel("Temperature [K]")
    ax.set_title("Interface temperature (imperfect contact, Rc = 1e-4 m²K/W)")
    ax.grid(True, alpha=0.3)
    ax.legend(ncol=3, fontsize=9, loc="lower center")
    fig.tight_layout()
    out = RESULTS / "interface_profile.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    render_snapshots()
    render_interface_profile()
