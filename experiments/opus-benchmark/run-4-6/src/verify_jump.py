#!/usr/bin/env python3
"""Verify that the interface temperature jump matches q·Rc.

Method:
1. Read interface samples (iface_###.txt) and the full temperature arrays.
2. Compute q_independent = λ_Cu · ∂T_Cu/∂x at x = W1 by finite difference
   using the two columns of Cu nodes closest to the interface.
3. Compare ΔT_measured with q_independent · Rc.
4. Write a markdown report.
"""
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / ".claude" / "skills" / "fem-analysis" / "scripts"))
from ff_mesh_reader import read_freefem_msh  # noqa: E402

RESULTS = REPO / "results"
LAMBDA_CU = 400.0
LAMBDA_ST = 50.0
Rc = 1.0e-4
W1 = 0.050
H = 0.030
SNAPS = [60, 150, 300]


def column_temperatures(mesh_file: Path, sol_file: Path, x_target: float, tol: float):
    mesh = read_freefem_msh(str(mesh_file))
    points = mesh["points"]
    T = np.loadtxt(sol_file)
    mask = np.abs(points[:, 0] - x_target) < tol
    if not mask.any():
        return None
    ys = points[mask, 1]
    Ts = T[mask]
    order = np.argsort(ys)
    return ys[order], Ts[order]


def flux_from_finite_difference(mesh_file, sol_file, x_interface, side):
    """Estimate heat flux at interface using nearest interior column."""
    mesh = read_freefem_msh(str(mesh_file))
    points = mesh["points"]
    T = np.loadtxt(sol_file)
    xs = points[:, 0]

    # Interface column (x ≈ x_interface)
    m_if = np.abs(xs - x_interface) < 1e-8
    if side == "cu":
        interior_candidates = xs[xs < x_interface - 1e-8]
        x_int = interior_candidates.max()
        m_in = np.abs(xs - x_int) < 1e-8
        dx = x_interface - x_int
        # Heat flux out of Cu in +x direction: q = -λ · dT/dx (convention)
        # With T_if > T_in typically, dT/dx > 0 at interface → q = -λ·+ < 0
        # We want magnitude and sign such that q > 0 means Cu → Steel:
        # use q = λ_Cu * (T_if - T_in)/dx (positive means heat flows right when T_if > T_in)
        T_if = T[m_if]; T_in = T[m_in]
        y_if = points[m_if, 1]; y_in = points[m_in, 1]
        # Sort by y, pair
        i_if = np.argsort(y_if); i_in = np.argsort(y_in)
        T_if, y_if = T_if[i_if], y_if[i_if]
        T_in, y_in = T_in[i_in], y_in[i_in]
        # Interpolate T_in onto y_if (they should match but safer)
        T_in_at_if = np.interp(y_if, y_in, T_in)
        # Actually the correct 1D flux from Cu into interface: heat moves from interior to surface,
        # so dT/dx at interface ≈ (T_if - T_in)/dx. But dT/dx<0 when T_if<T_in (cooler surface).
        # q (W/m^2) in +x direction from Cu body: q_Cu = λ_Cu * (T_in - T_if)/dx when interior hotter.
        # During warmup, T_in > T_if on Cu side (interior is hotter than interface because left side
        # is Dirichlet at 373 K and the interface is losing heat to steel).
        # So q_right (Cu -> Steel) = λ_Cu * (T_in - T_if)/dx  (positive when interior hotter).
        q = LAMBDA_CU * (T_in_at_if - T_if) / dx
        return y_if, q, T_if, T_in_at_if, dx
    elif side == "st":
        interior_candidates = xs[xs > x_interface + 1e-8]
        x_int = interior_candidates.min()
        m_in = np.abs(xs - x_int) < 1e-8
        dx = x_int - x_interface
        T_if = T[m_if]; T_in = T[m_in]
        y_if = points[m_if, 1]; y_in = points[m_in, 1]
        i_if = np.argsort(y_if); i_in = np.argsort(y_in)
        T_if, y_if = T_if[i_if], y_if[i_if]
        T_in, y_in = T_in[i_in], y_in[i_in]
        T_in_at_if = np.interp(y_if, y_in, T_in)
        # Steel side: flux into the interior in +x direction:
        # q_right = -λ_St * (T_in - T_if)/dx = λ_St * (T_if - T_in)/dx
        # During warmup, T_if > T_in on steel (interface fed by Cu), so q > 0.
        q = LAMBDA_ST * (T_if - T_in_at_if) / dx
        return y_if, q, T_if, T_in_at_if, dx


def analyze(t):
    tag = f"{t:03d}"
    iface = np.loadtxt(RESULTS / f"iface_{tag}.txt")
    y_face, Tcu_face, Tst_face = iface[:, 0], iface[:, 1], iface[:, 2]
    dT_measured = Tcu_face - Tst_face

    y_cu, q_cu, T_cu_if, T_cu_in, dx_cu = flux_from_finite_difference(
        RESULTS / "mesh_cu.msh", RESULTS / f"T_cu_{tag}.txt", W1, "cu"
    )
    y_st, q_st, T_st_if, T_st_in, dx_st = flux_from_finite_difference(
        RESULTS / "mesh_st.msh", RESULTS / f"T_st_{tag}.txt", W1, "st"
    )
    # Interpolate q onto iface y
    q_cu_i = np.interp(y_face, y_cu, q_cu)
    q_st_i = np.interp(y_face, y_st, q_st)
    q_avg = 0.5 * (q_cu_i + q_st_i)

    dT_predicted = q_avg * Rc
    rel_err = np.abs(dT_measured - dT_predicted) / np.maximum(np.abs(dT_measured), 1e-12)

    return {
        "t": t,
        "y": y_face,
        "Tcu": Tcu_face,
        "Tst": Tst_face,
        "dT_meas": dT_measured,
        "q_cu": q_cu_i,
        "q_st": q_st_i,
        "q_avg": q_avg,
        "dT_pred": dT_predicted,
        "rel_err": rel_err,
    }


def write_report(results):
    lines = []
    lines.append("# Contact-resistance verification report")
    lines.append("")
    lines.append(
        "Checks whether the observed temperature jump at the contact surface "
        "(x = 50 mm) equals the heat flux through the surface multiplied by the "
        "contact resistance Rc = 1×10⁻⁴ m²K/W."
    )
    lines.append("")
    lines.append("Flux is computed **independently** on each side by finite-difference "
                 "along x using the interface column and the nearest interior column:")
    lines.append("")
    lines.append("- Cu side:  `q_Cu = λ_Cu · (T_interior − T_interface) / Δx_Cu`")
    lines.append("- Steel side: `q_St = λ_St · (T_interface − T_interior) / Δx_St`")
    lines.append("")
    lines.append("`q_avg = (q_Cu + q_St)/2`.  Predicted jump: `ΔT_pred = q_avg · Rc`.")
    lines.append("")
    lines.append("## Summary (mid-interface y = 15 mm)")
    lines.append("")
    lines.append("| t [s] | T_Cu [K] | T_St [K] | ΔT_meas [K] | q_Cu [W/m²] | q_St [W/m²] | q_avg [W/m²] | ΔT_pred = q·Rc [K] | rel err [%] |")
    lines.append("|------:|---------:|---------:|------------:|------------:|------------:|-------------:|-------------------:|------------:|")
    worst = 0.0
    for r in results:
        mid = np.argmin(np.abs(r["y"] - H/2))
        line = "| {t} | {Tcu:.3f} | {Tst:.3f} | {dT:.4f} | {qcu:.2f} | {qst:.2f} | {qavg:.2f} | {dTp:.4f} | {err:.3f} |".format(
            t=r["t"], Tcu=r["Tcu"][mid], Tst=r["Tst"][mid],
            dT=r["dT_meas"][mid], qcu=r["q_cu"][mid], qst=r["q_st"][mid],
            qavg=r["q_avg"][mid], dTp=r["dT_pred"][mid], err=100*r["rel_err"][mid],
        )
        lines.append(line)
        worst = max(worst, float(r["rel_err"][mid]))
    lines.append("")
    pass_threshold = 0.05  # 5% — generous because flux FD is first-order
    status = "PASS" if worst < pass_threshold else "FAIL"
    lines.append(f"**Worst relative error at mid-interface: {100*worst:.3f}% → {status}**")
    lines.append(
        f"(Pass threshold = {100*pass_threshold:.1f}%. "
        "Flux is computed by first-order finite difference, so a few percent of "
        "error is expected from the FD discretization alone.)"
    )
    lines.append("")

    # Per-y table for t=300s (detailed)
    lines.append("## Per-y details at t = 300 s")
    lines.append("")
    lines.append("| y [mm] | T_Cu | T_St | ΔT_meas | q_Cu | q_St | ΔT_pred | rel err [%] |")
    lines.append("|-------:|-----:|-----:|--------:|-----:|-----:|--------:|------------:|")
    r300 = [x for x in results if x["t"] == 300][0]
    for i in range(0, len(r300["y"]), max(1, len(r300["y"])//8)):
        lines.append(
            "| {y:5.1f} | {Tcu:.3f} | {Tst:.3f} | {dT:.4f} | {qcu:.2f} | {qst:.2f} | {dTp:.4f} | {err:.3f} |".format(
                y=r300["y"][i]*1000, Tcu=r300["Tcu"][i], Tst=r300["Tst"][i],
                dT=r300["dT_meas"][i], qcu=r300["q_cu"][i], qst=r300["q_st"][i],
                dTp=r300["dT_pred"][i], err=100*r300["rel_err"][i],
            )
        )
    lines.append("")

    # Steady-state sanity check
    Rtot = W1/LAMBDA_CU + Rc + W1/LAMBDA_ST + 1.0/10.0
    q_steady = (373 - 293) / Rtot
    dT_steady = q_steady * Rc
    lines.append("## Steady-state analytical reference (1-D)")
    lines.append("")
    lines.append(
        "Treating the problem as 1-D series resistances "
        "`L_Cu/λ_Cu + Rc + L_St/λ_St + 1/h`:"
    )
    lines.append("")
    lines.append(f"- R_total = {Rtot:.5f} m²K/W")
    lines.append(f"- q_steady = (373 - 293) / R_total = **{q_steady:.2f} W/m²**")
    lines.append(f"- ΔT_interface at steady state = q · Rc = **{dT_steady:.4f} K**")
    lines.append(
        f"- Observed ΔT at t=300 s (mid-interface): **{r300['dT_meas'][len(r300['y'])//2]:.3f} K** "
        "(still warming up, so ΔT is larger than steady-state because transient heat "
        "flux into Cu is higher than the eventual steady flux)."
    )
    lines.append("")

    # Picard iteration info
    lines.append("## Notes")
    lines.append("")
    lines.append(
        "- The temperature is **represented on two disjoint meshes** (copper and steel) "
        "with independent P1 DOFs along the interface, so the jump is physical, not numerical."
    )
    lines.append(
        "- Per-step Picard iteration couples the two sides via the Robin-type interface "
        "flux `q = (T_Cu − T_St) / Rc`, which is mathematically equivalent to "
        "imposing `ΔT = q·Rc` directly. Independent verification uses flux from the "
        "interior gradient (finite difference), not the same Robin BC."
    )

    report = "\n".join(lines) + "\n"
    (RESULTS / "verification.md").write_text(report)
    print(f"Saved: {RESULTS / 'verification.md'}")
    print("")
    print("Summary (mid-interface):")
    for r in results:
        mid = np.argmin(np.abs(r["y"] - H/2))
        print(
            f"  t={r['t']:3d}s  ΔT_meas={r['dT_meas'][mid]:7.4f} K  "
            f"q_cu={r['q_cu'][mid]:8.2f} q_st={r['q_st'][mid]:8.2f} W/m² "
            f"ΔT_pred={r['dT_pred'][mid]:7.4f} K  rel_err={100*r['rel_err'][mid]:.3f}%"
        )


if __name__ == "__main__":
    results = [analyze(t) for t in SNAPS]
    write_report(results)
