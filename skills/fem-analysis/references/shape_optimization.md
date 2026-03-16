# Shape Optimization Reference (FreeFEM++)

## Overview

Shape optimization modifies the geometry boundary to optimize a functional,
without changing the topology (no new holes). Three sub-categories:

1. **Parametric (Sizing)**: Optimize a parameter field (e.g., plate thickness)
2. **Geometric**: Move the boundary of the shape
3. **Level-set based shape optimization**: Represent boundary implicitly, allows topology change

### Key Resources
- CMAP Toolbox: http://www.cmap.polytechnique.fr/~allaire/freefem_en.html
  - `plate.edp` — thickness optimization
  - `cantilever_en.edp` — geometric shape optimization
  - `levelset-cantilever.edp` — level-set method
- AutoFreeFEM: https://gitlab.tugraz.at/autofreefem/autofreefem
- FreeFEM docs: https://doc.freefem.org/documentation/algorithms-and-optimization.html
  - BFGS, NLCG, IPOPT, NLopt optimizers

---

## 1. Parametric Optimization (Thickness / Sizing)

### Template: Optimal Plate Thickness

```freefem
// Minimize compliance of an elastic plate with varying thickness h(x,y)
// h ∈ [hmin, hmax], subject to volume constraint ∫h = V_target

real L = 1.0, H = 1.0;
real hmin = 0.1, hmax = 1.0;
real hmean = 0.5;          // Target mean thickness
int niter = 50;
real step = 0.01;

mesh Th = square(40, 40, [L*x, H*y]);
fespace Vh(Th, [P2, P2]);
fespace Dh(Th, P1);

Dh h = hmean;  // Initial uniform thickness

real E = 1.0, nu = 0.3;
real mu = E/(2*(1+nu));
real lambda = E*nu/((1+nu)*(1-nu));

for (int iter = 0; iter < niter; iter++) {
    // Solve elasticity with thickness-dependent stiffness
    Vh [u1, u2], [v1, v2];
    solve Elast([u1, u2], [v1, v2]) =
        int2d(Th)(
            h * (
                lambda*(dx(u1)+dy(u2))*(dx(v1)+dy(v2))
                + 2*mu*(dx(u1)*dx(v1)+dy(u2)*dy(v2))
                + mu*(dy(u1)+dx(u2))*(dy(v1)+dx(v2))
            )
        )
        - int1d(Th, 2)(-1.0 * v2)
        + on(4, u1=0, u2=0);

    // Sensitivity: dJ/dh = -energy_density
    Dh sens = -(lambda*(dx(u1)+dy(u2))^2
                + 2*mu*(dx(u1)^2+dy(u2)^2)
                + mu*(dy(u1)+dx(u2))^2);

    // Project gradient to satisfy volume constraint
    real meanSens = int2d(Th)(sens) / (L*H);
    sens = sens - meanSens;

    // Update thickness
    h = h - step * sens;
    h = max(hmin, min(hmax, h));

    // Re-enforce volume
    real vol = int2d(Th)(h) / (L*H);
    h = h + (hmean - vol);
    h = max(hmin, min(hmax, h));
}

savemesh(Th, "thickness_mesh.msh");
{ ofstream f("thickness_h.txt");
  for (int i=0; i<Dh.ndof; i++) f << h[][i] << endl; }
```

---

## 2. Geometric Shape Optimization

### Concept

Move the boundary ∂Ω in the normal direction with velocity V_n derived
from the shape derivative:

```
dJ(Ω; V) = ∫_{∂Ω} V_n · g(u, p) ds
```

where g depends on the specific problem. For compliance:
```
g = σ(u):ε(u) - lagrange    (on the free boundary)
```

### Key FreeFEM++ Features

- `movemesh(Th, [x+dx, y+dy])` — deform mesh
- `adaptmesh(Th, ...)` — remesh when quality degrades
- `checkmovemesh(Th, [x+dx, y+dy])` — check validity before moving
- Boundary-only deformation: compute velocity field, extend to interior

### Template: Geometric Optimization with movemesh

```freefem
// Compliance minimization by moving the free boundary
real L = 2.0, H = 1.0;
int niter = 60;
real step = 0.005;
real E = 1.0, nu = 0.3;
real mu = E/(2*(1+nu));
real lambda = E*nu/((1+nu)*(1-nu));

// Initial mesh with labeled boundaries
border b1(t=0,L)    { x=t;   y=0;   label=1; }  // bottom (free)
border b2(t=0,H)    { x=L;   y=t;   label=2; }  // right (loaded)
border b3(t=L,0)    { x=t;   y=H;   label=3; }  // top (free)
border b4(t=H,0)    { x=0;   y=t;   label=4; }  // left (fixed)

int n = 30;
mesh Th = buildmesh(b1(n*2)+b2(n)+b3(n*2)+b4(n));

for (int iter = 0; iter < niter; iter++) {
    fespace Vh(Th, [P2, P2]);
    fespace Wh(Th, P1);

    Vh [u1, u2], [v1, v2];
    solve Elast([u1, u2], [v1, v2]) =
        int2d(Th)(
            lambda*(dx(u1)+dy(u2))*(dx(v1)+dy(v2))
            + 2*mu*(dx(u1)*dx(v1)+dy(u2)*dy(v2))
            + mu*(dy(u1)+dx(u2))*(dy(v1)+dx(v2))
        )
        - int1d(Th, 2)( (abs(y-H/2)<0.05)*(-1.0)*v2 )
        + on(4, u1=0, u2=0);

    // Shape gradient: velocity for boundary deformation
    // Extension-regularization: solve Laplacian to extend to interior
    Wh Vx, Vy, wx, wy;
    Wh energyDens = lambda*(dx(u1)+dy(u2))^2
                    + 2*mu*(dx(u1)^2+dy(u2)^2)
                    + mu*(dy(u1)+dx(u2))^2;

    // Solve extension problem: regularize shape gradient
    real alpha = 0.01;  // regularization parameter
    solve ExtendVel([Vx, Vy], [wx, wy]) =
        int2d(Th)(alpha*(dx(Vx)*dx(wx)+dy(Vx)*dy(wx)
                        +dx(Vy)*dx(wy)+dy(Vy)*dy(wy))
                  + Vx*wx + Vy*wy)
        + on(4, Vx=0, Vy=0)      // fixed boundary: no motion
        + on(2, Vx=0, Vy=0);     // loaded boundary: no motion

    // Deform mesh
    real minarea = checkmovemesh(Th, [x-step*Vx, y-step*Vy]);
    if (minarea > 0) {
        Th = movemesh(Th, [x-step*Vx, y-step*Vy]);
    }

    // Adapt mesh if quality is poor
    if (iter % 10 == 0)
        Th = adaptmesh(Th, [u1, u2], err=0.01, nbvx=10000);

    real compliance = int2d(Th)(energyDens);
    real vol = Th.measure;
    if (iter % 5 == 0)
        cout << "Iter " << iter << ": J=" << compliance
             << " Vol=" << vol << endl;
}

savemesh(Th, "shape_opt_mesh.msh");
```

### Important: Mesh Quality

Geometric shape optimization moves the mesh, which can cause:
- Tangled elements (negative area) → check with `checkmovemesh`
- Very stretched elements → use `adaptmesh` periodically
- Volume loss → add volume constraint via Lagrange multiplier

---

## 3. FreeFEM++ Optimization Solvers

### Built-in Optimizers

```freefem
// BFGS (quasi-Newton) — for few parameters
real[int] x0(3);  // Initial guess
BFGS(J, dJ, x0, eps=1e-6, nbiter=100);
// J: func real → real (objective)
// dJ: func real[int] → real[int] (gradient)

// NLCG (Nonlinear Conjugate Gradient) — for FE functions
NLCG(J, dJ, x0, eps=1e-6, nbiter=100);
// Better for large-scale (density as FE function)
```

### IPOPT (Interior Point Optimizer)

```freefem
load "ff-Ipopt"
// For constrained optimization with analytical Jacobians
// IPOPT(J, dJ, constraints, Jacobian, ...)
// See FreeFEM doc: algorithms-and-optimization.html
```

### NLopt

```freefem
load "ff-NLopt"
// Derivative-free methods: COBYLA, BOBYQA, Nelder-Mead
// Gradient-based: MMA, SLSQP, LBFGS
```

---

## Combining with AI Agent Workflow

### Workflow for AI-Driven Optimization

```
User: "Optimize the shape of a bracket with these loads and constraints"
    ↓
AI: Choose method (SIMP for topology, Level-Set for shape+topology)
    ↓
AI: Generate .edp file with appropriate template
    ↓
AI: Run FreeFEM++ → capture iteration history
    ↓
AI: Visualize intermediate + final density/shape with PyVista
    ↓
AI: Review result (check symmetry, load paths, manufacturing feasibility)
    ↓
AI: Generate convergence plot + final design summary
    ↓
AI: Optionally export to React 3D viewer for interactive exploration
```

### AI Self-Evaluation Checklist for Optimization Results

1. Did the optimizer converge? (check final Δρ < threshold)
2. Is the volume constraint satisfied? (within 1% of target)
3. Is the design binary? (for SIMP: check % of intermediate densities)
4. Are load paths visible? (continuous material from load to support)
5. Is the design symmetric if the problem is symmetric?
6. Are there floating islands of material? (disconnected regions = bad)
7. Is mesh-independence achieved? (re-run on finer mesh to verify)
