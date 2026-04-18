# Thermal Contact Resistance Verification

- Target Rc = **1.000e-04 m²K/W**
- λ_Cu = 400.0 W/mK   λ_St = 50.0 W/mK
- Interface at x = 0.05 m, height 0.03 m (conformal meshes)
- Flux method: length-weighted average of element-wise P1 gradients on
  each boundary edge along the interface.

## Per-snapshot check

| t [s] | T_Cu⁻ [K] | T_St⁺ [K] | ΔT [K] | q_Cu_out [W/m²] | q_St_in [W/m²] | Rc_eff [m²K/W] | &#124;err&#124; [%] | flux-cont. err [%] |
|------:|----------:|----------:|-------:|----------------:|---------------:|---------------:|-------------------:|-------------------:|
| 60 | 361.5400 | 353.3320 | +8.2080 | +82334.240 | +81546.123 | 1.0017e-04 | 0.170 | 0.962 |
| 150 | 367.9800 | 364.2111 | +3.7689 | +37611.679 | +37501.588 | 1.0035e-04 | 0.351 | 0.293 |
| 300 | 371.5120 | 370.3920 | +1.1200 | +11400.823 | +11152.900 | 9.9318e-05 | 0.682 | 2.199 |

## Interpretation

- `q_Cu_out` = `-λ_Cu · ∇T · n_out` on Cu interface (n_out = +x̂), averaged over
  boundary edges (length-weighted).
- `q_St_in`  = `+λ_St · ∇T · n_out_reversed` on Steel interface (n_out = -x̂), i.e.
  the heat flux arriving into steel — should equal `q_Cu_out` by conservation.
- `Rc_eff = (T_Cu⁻ − T_St⁺) / q_iface`, using the mean of the two fluxes.
- `flux-continuity err = |q_Cu_out − q_St_in| / q_iface` — measures how well
  discrete heat balance is satisfied at the interface (should be near 0).

## Analytical steady-state cross-check

1D thermal-resistance series (per unit area):
- R_Cu = L/λ_Cu = 0.05/400 = 1.250 × 10⁻⁴
- R_Rc =                   1.000 × 10⁻⁴
- R_St = L/λ_St = 0.05/50  = 1.000 × 10⁻³
- R_conv = 1/h              = 1.000 × 10⁻¹
- R_tot ≈ 1.01350 × 10⁻¹ m²K/W
- q_ss ≈ (373 − 293) / R_tot ≈ 789 W/m²
- ΔT_contact_ss = q_ss · Rc ≈ 0.0789 K

During the transient, interface flux is much higher than the eventual steady
value because heat is still propagating into the cold steel, so ΔT_contact
is larger and decays monotonically toward 0.08 K.

## Summary

- Max |Rc_eff − Rc| / Rc  = **0.682 %**
- Max flux-continuity err = **2.199 %**

**PASS** — tolerance 5 % on both metrics.
