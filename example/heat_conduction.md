# Heat Conduction Analysis

## Prompt

> Simulate steady-state heat conduction on a rectangular plate (1m x 0.5m).
> Left edge is held at 100 degrees C, right edge has convective cooling
> (h=10 W/m2K, ambient 20 degrees C), top and bottom are insulated.
> Material is steel (k=50 W/mK). Include a volumetric heat source of 10 kW/m3.

## Expected Behavior

1. Claude generates `heat.edp` with mixed boundary conditions (Dirichlet + Robin + Neumann)
2. Executes `FreeFem++ heat.edp`
3. Outputs mesh (`mesh.msh`) and temperature field (`temperature.txt`)
4. Visualizes temperature distribution using `coolwarm` or `hot` colormap
5. Self-evaluates: checks temperature gradient direction, BC satisfaction, physical plausibility
6. Presents the result with min/max temperature values

## Key Points

- Mixed BCs: Dirichlet (left), Robin/convection (right), Neumann/insulated (top, bottom)
- P1 elements are sufficient for thermal analysis
- Temperature should decrease from left to right
- Heat source raises the overall temperature above the simple conduction solution
- Colormap: `coolwarm` or `hot` is recommended for temperature fields

## Prompt Variations

**Transient analysis:**
> Run a transient heat conduction simulation. Start at 20 degrees C everywhere,
> then apply 100 degrees C on the left edge. Show temperature snapshots at t=0.1s, 0.5s, and 1.0s.

**With a heat sink:**
> Add a small rectangular region (0.4-0.6m x 0.2-0.3m) with high conductivity (k=400 W/mK)
> acting as a heat sink. Visualize the temperature field.

**Multi-material:**
> The left half is steel (k=50), the right half is aluminum (k=237).
> Show how the conductivity jump affects the temperature distribution.
