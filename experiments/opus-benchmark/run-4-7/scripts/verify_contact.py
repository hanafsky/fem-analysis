#!/usr/bin/env python3
"""Verify the Rc contact condition:   q · Rc  ≟  ΔT  at the interface x = 0.05.

Because FreeFEM's buildmesh produces an **unstructured** triangulation, there
is no "column of vertices" to use for a one-sided FD. Instead we compute the
interface heat flux element-wise:

    For every boundary edge with the interface label, find its owning triangle,
    compute the constant P1 gradient of T on that triangle, and take
            q_n = -k · (∇T · n_out)
    Average over the interface, weighted by edge length.

This gives the physically correct outward heat flux for a P1 finite-element
solution.

Outputs  results/contact_report.md.
"""

import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = ROOT / "results"
sys.path.insert(0, str(ROOT / ".claude" / "skills" / "fem-analysis" / "scripts"))
from ff_mesh_reader import read_freefem_msh  # noqa: E402

RC_TARGET = 1.0e-4
K_CU = 400.0
K_ST = 50.0
X_INTERFACE = 0.05


def build_edge_triangle_map(tris: np.ndarray) -> dict:
    """Map (min_vidx, max_vidx) -> list of triangle indices containing that edge."""
    m = defaultdict(list)
    nt = len(tris)
    for ti in range(nt):
        for i in range(3):
            a = int(tris[ti, i])
            b = int(tris[ti, (i + 1) % 3])
            m[(min(a, b), max(a, b))].append(ti)
    return m


def triangle_grad(P: np.ndarray, Tvals: np.ndarray):
    """P1 gradient on a triangle. P shape (3,2), Tvals shape (3,)."""
    x0, x1, x2 = P[:, 0]
    y0, y1, y2 = P[:, 1]
    area2 = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)   # = 2·area (signed)
    dN0_dx = (y1 - y2) / area2
    dN1_dx = (y2 - y0) / area2
    dN2_dx = (y0 - y1) / area2
    dN0_dy = (x2 - x1) / area2
    dN1_dy = (x0 - x2) / area2
    dN2_dy = (x1 - x0) / area2
    dT_dx = dN0_dx * Tvals[0] + dN1_dx * Tvals[1] + dN2_dx * Tvals[2]
    dT_dy = dN0_dy * Tvals[0] + dN1_dy * Tvals[1] + dN2_dy * Tvals[2]
    return dT_dx, dT_dy, abs(area2) * 0.5


def interface_flux(mesh: dict, T: np.ndarray, k: float, iface_label: int, n_out: tuple):
    """Length-weighted-average outward heat flux `-k · ∇T · n_out` on the interface."""
    pts = mesh["points"][:, :2]
    tris = mesh["triangles"]
    edges = mesh["edges"]
    edge_labels = mesh["edge_boundary_labels"]

    e2t = build_edge_triangle_map(tris)

    fluxes = []
    lengths = []
    T_edge = []
    y_edge = []
    for ei in range(len(edges)):
        if int(edge_labels[ei]) != iface_label:
            continue
        v0, v1 = int(edges[ei, 0]), int(edges[ei, 1])
        key = (min(v0, v1), max(v0, v1))
        cand = e2t.get(key, [])
        if len(cand) != 1:
            # Boundary edges belong to exactly one triangle; skip otherwise.
            continue
        ti = cand[0]
        tri = tris[ti]
        dT_dx, dT_dy, _ = triangle_grad(pts[tri], T[tri])
        q = -k * (dT_dx * n_out[0] + dT_dy * n_out[1])

        p0 = pts[v0]
        p1 = pts[v1]
        L = float(np.linalg.norm(p1 - p0))
        fluxes.append(q)
        lengths.append(L)
        T_edge.append(0.5 * (T[v0] + T[v1]))
        y_edge.append(0.5 * (p0[1] + p1[1]))

    fluxes = np.asarray(fluxes)
    lengths = np.asarray(lengths)
    T_edge = np.asarray(T_edge)
    y_edge = np.asarray(y_edge)

    q_weighted = float(np.sum(fluxes * lengths) / np.sum(lengths))
    T_weighted = float(np.sum(T_edge * lengths) / np.sum(lengths))

    # Interface nodal T (average of nodes that live exactly on x = X_INTERFACE)
    iface_node_mask = np.isclose(pts[:, 0], X_INTERFACE, atol=1e-9)
    T_node_mean = float(T[iface_node_mask].mean())
    T_node_std = float(T[iface_node_mask].std())

    return {
        "q_avg": q_weighted,
        "q_min": float(fluxes.min()),
        "q_max": float(fluxes.max()),
        "q_std": float(fluxes.std()),
        "num_edges": int(len(fluxes)),
        "T_edge_avg": T_weighted,
        "T_node_mean": T_node_mean,
        "T_node_std": T_node_std,
        "y_edge": y_edge,
        "T_edge": T_edge,
        "fluxes": fluxes,
    }


def analyse_snapshot(t: int) -> dict:
    mesh_cu = read_freefem_msh(str(RESULTS / "mesh_cu.msh"))
    mesh_st = read_freefem_msh(str(RESULTS / "mesh_st.msh"))
    T1 = np.loadtxt(RESULTS / f"T1_t{t}.txt")
    T2 = np.loadtxt(RESULTS / f"T2_t{t}.txt")

    # Cu: interface is the RIGHT side of Th1 → label 2, outward n = +x̂
    cu = interface_flux(mesh_cu, T1, K_CU, iface_label=2, n_out=(+1.0, 0.0))
    # Steel: interface is the LEFT side of Th2 → label 4, outward n = -x̂
    st = interface_flux(mesh_st, T2, K_ST, iface_label=4, n_out=(-1.0, 0.0))

    # Signed interface flux, positive if heat flows Cu → Steel
    q_cu_out = cu["q_avg"]        # outward from Cu (positive = leaving Cu)
    q_st_in = -st["q_avg"]        # into Steel (positive = entering)
    q_iface = 0.5 * (q_cu_out + q_st_in)   # best-estimate signed interface flux
    flux_continuity_err = abs(q_cu_out - q_st_in) / max(abs(q_iface), 1e-12) * 100

    dT = cu["T_node_mean"] - st["T_node_mean"]
    Rc_eff = dT / q_iface if abs(q_iface) > 1e-12 else float("nan")
    rel_err = abs(Rc_eff - RC_TARGET) / RC_TARGET * 100 if not np.isnan(Rc_eff) else float("nan")

    return {
        "t": t,
        "T_cu": cu["T_node_mean"],
        "T_st": st["T_node_mean"],
        "dT": dT,
        "q_cu_out": q_cu_out,
        "q_st_in": q_st_in,
        "q_iface": q_iface,
        "Rc_eff": Rc_eff,
        "rel_err_pct": rel_err,
        "flux_continuity_err_pct": flux_continuity_err,
        "num_iface_edges_cu": cu["num_edges"],
        "num_iface_edges_st": st["num_edges"],
    }


def main():
    lines = [
        "# Thermal Contact Resistance Verification",
        "",
        f"- Target Rc = **{RC_TARGET:.3e} m²K/W**",
        f"- λ_Cu = {K_CU} W/mK   λ_St = {K_ST} W/mK",
        "- Interface at x = 0.05 m, height 0.03 m (conformal meshes)",
        "- Flux method: length-weighted average of element-wise P1 gradients on",
        "  each boundary edge along the interface.",
        "",
        "## Per-snapshot check",
        "",
        "| t [s] | T_Cu⁻ [K] | T_St⁺ [K] | ΔT [K] | q_Cu_out [W/m²] | q_St_in [W/m²] | Rc_eff [m²K/W] | &#124;err&#124; [%] | flux-cont. err [%] |",
        "|------:|----------:|----------:|-------:|----------------:|---------------:|---------------:|-------------------:|-------------------:|",
    ]

    results = []
    for t in (60, 150, 300):
        r = analyse_snapshot(t)
        results.append(r)
        lines.append(
            f"| {r['t']} "
            f"| {r['T_cu']:.4f} "
            f"| {r['T_st']:.4f} "
            f"| {r['dT']:+.4f} "
            f"| {r['q_cu_out']:+.3f} "
            f"| {r['q_st_in']:+.3f} "
            f"| {r['Rc_eff']:.4e} "
            f"| {r['rel_err_pct']:.3f} "
            f"| {r['flux_continuity_err_pct']:.3f} |"
        )

    max_err = max(r["rel_err_pct"] for r in results)
    max_cont = max(r["flux_continuity_err_pct"] for r in results)

    lines += [
        "",
        "## Interpretation",
        "",
        "- `q_Cu_out` = `-λ_Cu · ∇T · n_out` on Cu interface (n_out = +x̂), averaged over",
        "  boundary edges (length-weighted).",
        "- `q_St_in`  = `+λ_St · ∇T · n_out_reversed` on Steel interface (n_out = -x̂), i.e.",
        "  the heat flux arriving into steel — should equal `q_Cu_out` by conservation.",
        "- `Rc_eff = (T_Cu⁻ − T_St⁺) / q_iface`, using the mean of the two fluxes.",
        "- `flux-continuity err = |q_Cu_out − q_St_in| / q_iface` — measures how well",
        "  discrete heat balance is satisfied at the interface (should be near 0).",
        "",
        "## Analytical steady-state cross-check",
        "",
        "1D thermal-resistance series (per unit area):",
        "- R_Cu = L/λ_Cu = 0.05/400 = 1.250 × 10⁻⁴",
        "- R_Rc =                   1.000 × 10⁻⁴",
        "- R_St = L/λ_St = 0.05/50  = 1.000 × 10⁻³",
        "- R_conv = 1/h              = 1.000 × 10⁻¹",
        "- R_tot ≈ 1.01350 × 10⁻¹ m²K/W",
        "- q_ss ≈ (373 − 293) / R_tot ≈ 789 W/m²",
        "- ΔT_contact_ss = q_ss · Rc ≈ 0.0789 K",
        "",
        "During the transient, interface flux is much higher than the eventual steady",
        "value because heat is still propagating into the cold steel, so ΔT_contact",
        "is larger and decays monotonically toward 0.08 K.",
        "",
        "## Summary",
        "",
        f"- Max |Rc_eff − Rc| / Rc  = **{max_err:.3f} %**",
        f"- Max flux-continuity err = **{max_cont:.3f} %**",
        "",
        (
            "**PASS**" if (max_err < 5.0 and max_cont < 5.0) else "**CHECK REQUIRED**"
        )
        + " — tolerance 5 % on both metrics.",
    ]

    out = RESULTS / "contact_report.md"
    out.write_text("\n".join(lines) + "\n")
    print(out.read_text())
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
