# Poisson Equation

## Prompt

> Solve the Poisson equation -laplacian(u) = 2*pi^2*sin(pi*x)*sin(pi*y)
> on a unit square with homogeneous Dirichlet boundary conditions.
> Compare the numerical solution with the exact solution u = sin(pi*x)*sin(pi*y).

## Expected Behavior

1. Claude generates `poisson.edp` with the given source term and BCs
2. Executes `FreeFem++ poisson.edp`
3. Outputs mesh (`mesh.msh`) and solution (`solution.txt`)
4. Computes and reports the L2 error against the exact solution
5. Visualizes the solution field using `viridis` colormap
6. Self-evaluates: checks symmetry, peak location at (0.5, 0.5), error magnitude

## Key Points

- Simplest PDE example — good for verifying the full pipeline works
- Exact solution is known, so numerical error can be quantified
- P1 elements on a 30x30 mesh should give L2 error ~ O(10^-3)
- Solution should peak at the center of the domain with value close to 1.0
- Homogeneous Dirichlet BCs on all four edges

## Prompt Variations

**L-shaped domain:**
> Solve the Poisson equation on an L-shaped domain.
> Use adaptive mesh refinement to capture the corner singularity.

**Neumann boundary:**
> Solve -laplacian(u) = 1 on a unit square with u=0 on the left edge
> and zero Neumann (natural) conditions on the other three edges.

**Visualization comparison:**
> Solve the Poisson equation on a unit square and show a side-by-side comparison
> of the numerical solution and the pointwise error.
