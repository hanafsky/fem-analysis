# Cantilever Beam Analysis

## Prompt

> Analyze a steel cantilever beam (1m x 0.1m) fixed on the left edge
> with a distributed load of 1 MPa on the right edge.
> Show the Von Mises stress and deformed shape.

## Expected Behavior

1. Claude generates `cantilever.edp` with plane stress formulation (E=210 GPa, nu=0.3)
2. Executes `FreeFem++ cantilever.edp`
3. Outputs mesh (`mesh.msh`), displacements (`u1.txt`, `u2.txt`), and stress (`vonmises.txt`)
4. Visualizes deformed mesh with Von Mises stress contour overlay using PyVista
5. Self-evaluates: checks deformation direction, stress concentration at fixed end, BC satisfaction
6. Presents the final image with key metrics (max displacement, max stress)

## Key Points

- Plane stress assumption (thin plate, sigma_zz = 0)
- P2 elements for displacement accuracy, projected to P1 for file output
- Auto-scaled deformation overlay (typically ~0.1 * L / max_displacement)
- Von Mises stress should peak near the fixed end (left edge)
- Analytical check: tip deflection delta = PL^3 / (3EI)

## Prompt Variations

**With a hole:**
> Same cantilever beam, but add a circular hole (radius 0.02m) at the center.
> Compare the stress concentration around the hole.

**Different material:**
> Analyze an aluminum cantilever beam (E=70 GPa, nu=0.33) with the same geometry.
> Compare the deflection with steel.

**With gravity:**
> Add self-weight (density 7800 kg/m^3) to the cantilever beam analysis.
