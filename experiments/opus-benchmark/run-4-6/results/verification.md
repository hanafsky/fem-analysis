# Contact-resistance verification report

Checks whether the observed temperature jump at the contact surface (x = 50 mm) equals the heat flux through the surface multiplied by the contact resistance Rc = 1×10⁻⁴ m²K/W.

Flux is computed **independently** on each side by finite-difference along x using the interface column and the nearest interior column:

- Cu side:  `q_Cu = λ_Cu · (T_interior − T_interface) / Δx_Cu`
- Steel side: `q_St = λ_St · (T_interface − T_interior) / Δx_St`

`q_avg = (q_Cu + q_St)/2`.  Predicted jump: `ΔT_pred = q_avg · Rc`.

## Summary (mid-interface y = 15 mm)

| t [s] | T_Cu [K] | T_St [K] | ΔT_meas [K] | q_Cu [W/m²] | q_St [W/m²] | q_avg [W/m²] | ΔT_pred = q·Rc [K] | rel err [%] |
|------:|---------:|---------:|------------:|------------:|------------:|-------------:|-------------------:|------------:|
| 60 | 361.540 | 353.331 | 8.2085 | 82317.69 | 81701.22 | 82009.45 | 8.2009 | 0.092 |
| 150 | 367.979 | 364.211 | 3.7685 | 37750.04 | 37566.12 | 37658.08 | 3.7658 | 0.072 |
| 300 | 371.512 | 370.392 | 1.1205 | 11224.81 | 11171.43 | 11198.12 | 1.1198 | 0.061 |

**Worst relative error at mid-interface: 0.092% → PASS**
(Pass threshold = 5.0%. Flux is computed by first-order finite difference, so a few percent of error is expected from the FD discretization alone.)

## Per-y details at t = 300 s

| y [mm] | T_Cu | T_St | ΔT_meas | q_Cu | q_St | ΔT_pred | rel err [%] |
|-------:|-----:|-----:|--------:|-----:|-----:|--------:|------------:|
|   0.0 | 371.512 | 370.392 | 1.1204 | 11217.60 | 11181.35 | 1.1199 | 0.036 |
|   3.0 | 371.512 | 370.392 | 1.1204 | 11220.84 | 11174.69 | 1.1198 | 0.058 |
|   6.0 | 371.512 | 370.392 | 1.1204 | 11221.73 | 11174.61 | 1.1198 | 0.055 |
|   9.0 | 371.512 | 370.392 | 1.1204 | 11221.95 | 11175.12 | 1.1199 | 0.052 |
|  12.0 | 371.512 | 370.392 | 1.1204 | 11220.44 | 11176.16 | 1.1198 | 0.054 |
|  15.0 | 371.512 | 370.392 | 1.1205 | 11224.81 | 11171.43 | 1.1198 | 0.061 |
|  18.0 | 371.512 | 370.392 | 1.1204 | 11221.86 | 11177.06 | 1.1199 | 0.043 |
|  21.0 | 371.512 | 370.392 | 1.1205 | 11225.10 | 11172.45 | 1.1199 | 0.053 |
|  24.0 | 371.512 | 370.392 | 1.1205 | 11223.36 | 11170.80 | 1.1197 | 0.068 |
|  27.0 | 371.512 | 370.392 | 1.1204 | 11223.78 | 11170.44 | 1.1197 | 0.066 |
|  30.0 | 371.512 | 370.392 | 1.1203 | 11217.24 | 11180.16 | 1.1199 | 0.041 |

## Steady-state analytical reference (1-D)

Treating the problem as 1-D series resistances `L_Cu/λ_Cu + Rc + L_St/λ_St + 1/h`:

- R_total = 0.10123 m²K/W
- q_steady = (373 - 293) / R_total = **790.32 W/m²**
- ΔT_interface at steady state = q · Rc = **0.0790 K**
- Observed ΔT at t=300 s (mid-interface): **1.120 K** (still warming up, so ΔT is larger than steady-state because transient heat flux into Cu is higher than the eventual steady flux).

## Notes

- The temperature is **represented on two disjoint meshes** (copper and steel) with independent P1 DOFs along the interface, so the jump is physical, not numerical.
- Per-step Picard iteration couples the two sides via the Robin-type interface flux `q = (T_Cu − T_St) / Rc`, which is mathematically equivalent to imposing `ΔT = q·Rc` directly. Independent verification uses flux from the interior gradient (finite difference), not the same Robin BC.
