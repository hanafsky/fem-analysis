# Topology Optimization Reference (FreeFEM++)

## Overview

Topology optimization finds the optimal material distribution within a design domain
to minimize an objective (typically compliance) subject to a volume constraint.

Three main approaches implemented in FreeFEM++:
1. **SIMP** (Solid Isotropic Material with Penalization) — density-based
2. **Level-Set** method — boundary evolution via Hamilton-Jacobi
3. **Topological Derivative** — nucleation of holes

### Key Resources
- Allaire & Pantz, "Structural Optimization with FreeFEM++", SMO 32, 2006
- CMAP Toolbox: http://www.cmap.polytechnique.fr/~allaire/freefem_en.html
- 89-line code (Zhu et al., 2021): nonlinear SIMP in FreeFEM
- AutoFreeFEM (Allaire & Gfrerer, 2024): automatic shape derivative code generation
  - Repository: https://gitlab.tugraz.at/autofreefem/autofreefem

---

## Method 1: SIMP (Density-Based)

### Mathematical Formulation

Find density field θ(x) ∈ [0,1] to minimize:

```
min_θ  J(θ) = ∫_D σ(u) : ε(u) dx          (compliance)
s.t.   div(σ(u)) + f = 0  in D              (equilibrium)
       ∫_D θ dx = η·Vol(D)                   (volume constraint)
       E(x) = E_min + θ^p · (E_0 - E_min)   (SIMP interpolation)
```

where p ≥ 3 is the penalization exponent that drives θ toward 0 or 1.

### Template: Cantilever Compliance Minimization (SIMP)

```freefem
// =====================================================
// SIMP Topology Optimization — Cantilever Beam
// Based on CMAP Toolbox (Allaire et al.)
// =====================================================

// --- Design Parameters ---
real L = 2.0;          // Domain length
real H = 1.0;          // Domain height
real volfrac = 0.4;    // Target volume fraction (40%)
int niter = 80;        // Total optimization iterations
int niternp = 20;      // Non-penalized iterations (p=1)
real pmax = 4.0;       // Maximum SIMP exponent

// --- Material ---
real E0 = 1.0;         // Young's modulus (solid) [normalized]
real Emin = 1e-9;      // Young's modulus (void)
real nu = 0.3;         // Poisson's ratio
real mu0 = E0 / (2*(1+nu));
real lambda0 = E0*nu / ((1+nu)*(1-nu));  // plane stress

// --- Mesh ---
int nx = 80, ny = 40;
mesh Th = square(nx, ny, [L*x, H*y]);

// --- FE Spaces ---
fespace Vh(Th, [P2, P2]);   // Displacement (P2 for accuracy)
fespace Dh(Th, P1);          // Density field (P1)

// --- Initialize ---
Dh rho = volfrac;             // Uniform initial density
Dh rhoold, sensitivity;
real penal = 1.0;             // Start without penalization
real rfilter = 0.06;          // Helmholtz filter radius

// --- Optimization Loop ---
for (int iter = 0; iter < niter; iter++) {
    rhoold = rho;

    // Penalty ramp-up after non-penalized phase
    if (iter == niternp) penal = 2.0;
    if (iter > niternp && iter % 4 == 0 && penal < pmax)
        penal = min(penal + 0.5, pmax);

    // SIMP interpolated stiffness
    Dh Erho = Emin + rho^penal * (E0 - Emin);

    // Solve elasticity
    Vh [u1, u2], [v1, v2];
    solve Elasticity([u1, u2], [v1, v2]) =
        int2d(Th)(
            Erho * (
                lambda0*(dx(u1)+dy(u2))*(dx(v1)+dy(v2))
                + 2*mu0*(dx(u1)*dx(v1) + dy(u2)*dy(v2))
                + mu0*(dy(u1)+dx(u2))*(dy(v1)+dx(v2))
            )
        )
        - int1d(Th, 2)( (abs(y - H/2) < 0.05) * (-1.0) * v2 )
        + on(4, u1=0, u2=0);

    // Compliance
    real compliance = int2d(Th)(
        Erho * (lambda0*(dx(u1)+dy(u2))^2
               + 2*mu0*(dx(u1)^2 + dy(u2)^2)
               + mu0*(dy(u1)+dx(u2))^2)
    );

    // Sensitivity: dJ/dρ = -p·ρ^(p-1)·(E0-Emin)·ε:σ_unit
    Dh energyDens = lambda0*(dx(u1)+dy(u2))^2
                    + 2*mu0*(dx(u1)^2 + dy(u2)^2)
                    + mu0*(dy(u1)+dx(u2))^2;
    sensitivity = -penal * rho^(penal-1) * (E0-Emin) * energyDens;

    // Helmholtz PDE filter (avoids checkerboard)
    Dh sf, tf;
    solve SensFilter(sf, tf) =
        int2d(Th)(rfilter^2*(dx(sf)*dx(tf)+dy(sf)*dy(tf)) + sf*tf)
        - int2d(Th)(sensitivity * tf);
    sensitivity = sf;

    // Optimality Criteria (OC) update with bisection
    real l1 = 0, l2 = 1e9;
    Dh rhobackup = rho;
    real move = 0.2;
    while ((l2-l1)/(l1+l2) > 1e-3) {
        real lmid = 0.5*(l1+l2);
        for (int i = 0; i < Dh.ndof; i++) {
            real rhoNew = rhobackup[][i] * sqrt(-sensitivity[][i]/lmid);
            rhoNew = max(max(rhoNew, rhobackup[][i]-move), 1e-3);
            rhoNew = min(min(rhoNew, rhobackup[][i]+move), 1.0);
            rho[][i] = rhoNew;
        }
        real vol = int2d(Th)(rho) / (L*H);
        if (vol > volfrac) l1 = lmid; else l2 = lmid;
    }

    // Convergence monitoring
    real change = 0;
    for (int i = 0; i < Dh.ndof; i++)
        change = max(change, abs(rho[][i]-rhoold[][i]));
    real curvol = int2d(Th)(rho) / (L*H);

    if (iter % 10 == 0 || iter == niter-1)
        cout << "Iter " << iter << ": J=" << compliance
             << " vol=" << curvol << " p=" << penal
             << " Δρ=" << change << endl;
}

// --- Output ---
savemesh(Th, "topopt_mesh.msh");
{ ofstream f("topopt_rho.txt");
  for (int i=0; i<Dh.ndof; i++) f << rho[][i] << endl; }
cout << "ndof = " << Dh.ndof << endl;
```

### Variations

**Bridge / MBB beam** — change BCs and loads:
```freefem
// MBB beam: symmetry at left, roller at bottom-right corner
// + on(4, u1=0)    // symmetry: fix u1 only on left
// + on(1, u2=0)    // roller: fix u2 at bottom (or just one point)
// Load: top-left corner
// - int1d(Th, 3)( (x < 0.05) * (-1.0) * v2 )
```

**Multiple load cases** — solve for each, sum sensitivities:
```freefem
// for (int lc = 0; lc < nLoadCases; lc++) {
//     solve Elast_lc ... with load[lc]
//     sensitivity += weight[lc] * dJ_lc;
// }
```

**Self-weight (body force that depends on density)**:
```freefem
// Add: - int2d(Th)(rho * g * v2)   inside solve
// Additional sensitivity term required!
```

---

## Method 2: Level-Set Method

### Concept

Represent the shape Ω by a level-set function φ:
- φ(x) < 0  → inside Ω (solid)
- φ(x) > 0  → outside Ω (void)
- φ(x) = 0  → boundary ∂Ω

The shape evolves by solving a Hamilton-Jacobi equation:
```
∂φ/∂t + V·|∇φ| = 0
```
where V is the "velocity" derived from the shape derivative.

### Template: Level-Set Shape Optimization

```freefem
// Level-Set based compliance minimization
// Reference: Allaire, Jouve, Toader (2004)

real L = 2.0, H = 1.0;
real lagrange = 1.0;       // Lagrange multiplier for volume
real dt = 0.1;             // Pseudo-time step
int niter = 200;           // Optimization iterations
real E0 = 1.0, Emin = 1e-3, nu = 0.3;

// Mesh (fine, fixed)
int nx = 100, ny = 50;
mesh Th = square(nx, ny, [L*x, H*y]);

fespace Vh(Th, [P2, P2]);
fespace Wh(Th, P1);

// Initial level-set: perforated plate (multiple holes)
Wh phi;
phi = -cos(4*pi*x/L) * cos(4*pi*y/H) + 0.2;
// negative = solid, positive = void

for (int iter = 0; iter < niter; iter++) {
    // Ersatz material: use smoothed Heaviside of phi
    real eps = 0.01;
    Wh chi = 0.5*(1.0 - tanh(phi/eps));  // ≈1 in solid, ≈0 in void
    Wh Erho = Emin + chi * (E0 - Emin);

    real mu0 = E0/(2*(1+nu));
    real lambda0 = E0*nu/((1+nu)*(1-nu));

    // Solve elasticity with ersatz material
    Vh [u1, u2], [v1, v2];
    solve Elast([u1, u2], [v1, v2]) =
        int2d(Th)(
            Erho * (
                lambda0*(dx(u1)+dy(u2))*(dx(v1)+dy(v2))
                + 2*mu0*(dx(u1)*dx(v1)+dy(u2)*dy(v2))
                + mu0*(dy(u1)+dx(u2))*(dy(v1)+dx(v2))
            )
        )
        - int1d(Th, 2)( (abs(y-H/2)<0.05) * (-1.0) * v2 )
        + on(4, u1=0, u2=0);

    // Shape derivative = energy density - lagrange
    Wh shapeDeriv = -(
        lambda0*(dx(u1)+dy(u2))^2
        + 2*mu0*(dx(u1)^2+dy(u2)^2)
        + mu0*(dy(u1)+dx(u2))^2
    ) + lagrange;

    // Regularize velocity (solve auxiliary elliptic problem)
    Wh V, w;
    solve RegVelocity(V, w) =
        int2d(Th)(0.01*(dx(V)*dx(w)+dy(V)*dy(w)) + V*w)
        - int2d(Th)(shapeDeriv * w);

    // Update level-set: φ^{n+1} = φ^n + dt * V * |∇φ|
    // Simplified: direct descent
    phi = phi - dt * V;

    // Reinitialize level-set (periodically) to maintain signed distance
    if (iter % 5 == 0) {
        // Simple reinitialization by solving |∇φ| = 1
        // (or use adaptmesh for remeshing around φ=0)
    }

    // Monitor
    real vol = int2d(Th)(chi) / (L*H);
    real compliance = int2d(Th)(Erho * (
        lambda0*(dx(u1)+dy(u2))^2
        + 2*mu0*(dx(u1)^2+dy(u2)^2)
        + mu0*(dy(u1)+dx(u2))^2
    ));
    if (iter % 20 == 0)
        cout << "Iter " << iter << ": J=" << compliance
             << " vol=" << vol << endl;

    // Adjust Lagrange multiplier for volume
    real voltarget = 0.4;
    lagrange = lagrange + 0.5*(vol - voltarget);
}

savemesh(Th, "levelset_mesh.msh");
{ ofstream f("levelset_phi.txt");
  for (int i=0; i<Wh.ndof; i++) f << phi[][i] << endl; }
{ ofstream f("levelset_chi.txt");
  Wh chi = 0.5*(1.0 - tanh(phi/0.01));
  for (int i=0; i<Wh.ndof; i++) f << chi[][i] << endl; }
```

### Advanced: Level-Set with Mesh Evolution (mmg)

For higher quality boundaries, use the `mmg` remesher to adapt the mesh
to the zero level-set at each iteration. This requires `mmg` to be installed
and the `mshmet` plugin. See AutoFreeFEM for automated generation of this workflow.

---

## Method 3: Topological Derivative

### Concept

The topological derivative measures the sensitivity of the objective when
an infinitesimally small hole is created at point x. Used to decide where
to nucleate new holes.

For compliance with linear elasticity (plane stress):
```
DT(x) = π·( (λ+2μ)/(μ(λ+μ)) )·( 4μ·σ:ε + (λ-μ)·(tr σ)² / (λ+2μ) )
```

### Template: Topological Derivative Optimization

```freefem
// Topological derivative based optimization
// Reference: Luz Filho et al., SMO 66:74, 2023

// Same setup as SIMP, but the update rule differs:
// 1. Solve elasticity
// 2. Compute topological derivative DT(x) at all points
// 3. Evolve level-set: φ^{n+1} = φ^n + dt*(DT - lagrange)
// 4. Threshold: solid where φ < 0, void where φ > 0
// 5. Use adaptmesh for mesh refinement at φ=0

// Key FreeFEM feature: adaptmesh with iso-value tracking
// Th = adaptmesh(Th, phi, iso=0, nbvx=50000, err=0.01);
```

---

## Geometric Shape Optimization

### Concept

Unlike topology optimization, shape optimization moves the boundary
of an existing shape without creating/removing holes.

### Template: Cantilever Boundary Optimization

```freefem
// Geometric shape optimization
// The mesh boundary moves according to the shape gradient
// Uses movemesh to deform the domain

// Reference: cantilever_en.edp from CMAP Toolbox

real L = 2.0, H = 1.0;
real step = 0.02;     // Descent step size
int niter = 100;
real volmin = 0.3;    // Minimum volume fraction

// Initial mesh with boundary labels
border left(t=0, H)    { x=0; y=t;   label=1; }  // fixed
border top(t=0, L)      { x=t; y=H;   label=2; }  // free (optimizable)
border right(t=H, 0)    { x=L; y=t;   label=3; }  // loaded
border bottom(t=L, 0)   { x=t; y=0;   label=4; }  // free (optimizable)
mesh Th = buildmesh(left(20) + top(40) + right(20) + bottom(40));

for (int iter = 0; iter < niter; iter++) {
    // Solve elasticity on current mesh
    fespace Vh(Th, [P2, P2]);
    Vh [u1, u2], [v1, v2];
    // ... (same elasticity solve)

    // Compute shape gradient on boundary
    // Normal velocity = compliance_density - lagrange
    // Move boundary: new_boundary = old + step * velocity * normal

    // Use movemesh to deform
    // Th = movemesh(Th, [x + step*vx, y + step*vy]);

    // Remesh if quality degrades
    // Th = adaptmesh(Th, ...);
}
```

---

## AutoFreeFEM Integration

AutoFreeFEM (Python package) can automatically generate FreeFEM++ scripts
for shape optimization of arbitrary PDE-constrained problems.

### Installation
```bash
pip install autofreefem   # or: uv add autofreefem
```

### Usage Pattern
```python
from autofreefem import *

# Define fields and constants
u = VectorField("u", "P2", "Th", dirichlet={"4": [0, 0]})
v = VectorField("v", "P2", "Th")
E, nu = Constant("E", 1.0), Constant("nu", 0.3)

# Weak form (linear elasticity)
mu = E / (2*(1+nu))
lam = E*nu / ((1+nu)*(1-nu))
a = lam*div(u)*div(v) + 2*mu*inner(eps(u), eps(v))

# Objective (compliance)
J = integral(a.subs(v, u))

# Setup optimization → generates .edp file automatically
lag = Lagrangian(u, v, a, J)
lag.setup_optimization()
```

AutoFreeFEM automatically computes:
- Linearization (for nonlinear problems)
- Adjoint problem
- Shape derivative (via mathematical Lagrangian approach)
- Level-set mesh evolution code (using mmg)

---

## Visualization Notes

### For density-based (SIMP) results
```python
# In visualize.py, use binary colormap for final density
# cmap="binary" or threshold at 0.5:
# rho_binary = (rho > 0.5).astype(float)
```

### For level-set results
```python
# Visualize φ=0 iso-contour as boundary
# grid.contour([0], scalars="phi")  in PyVista
```

### Convergence plots
Track and plot: compliance vs iteration, volume vs iteration, max density change.
Export as CSV from FreeFEM++ and plot with matplotlib or in React artifact.

---

## Solver & Performance Notes

- SIMP on large meshes: use `solver=UMFPACK` or `solver=sparsesolver` for direct solvers
- For very large problems, FreeFEM++ supports PETSc/MPI parallelization
- Typical 2D problem (100×50 mesh, 80 iterations): ~30 seconds
- 3D problems require FreeFEM compiled with PETSc; see Allaire & Kelly's 3D toolbox
- IPOPT (`load "ff-Ipopt"`) available for constrained optimization with analytical gradients
- nlopt (`load "ff-NLopt"`) for derivative-free optimization methods
