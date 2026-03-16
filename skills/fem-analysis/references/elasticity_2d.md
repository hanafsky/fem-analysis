# 2D Elasticity Reference (FreeFEM++)

## Plane Stress / Plane Strain

### Governing Equations

Plane stress (thin plate): σ_zz = 0
Plane strain (long body): ε_zz = 0

Lamé parameters:
- Plane stress: λ = Eν/((1+ν)(1-ν)),  μ = E/(2(1+ν))
- Plane strain: λ = Eν/((1+ν)(1-2ν)), μ = E/(2(1+ν))

### Template: Cantilever Beam with Tip Load

```freefem
// === Parameters ===
real L = 1.0;     // Length [m]
real H = 0.1;     // Height [m]
real E = 210e9;   // Young's modulus [Pa]
real nu = 0.3;    // Poisson's ratio
real P = 1e6;     // Load [N/m]

// Lamé constants (plane stress)
real mu = E / (2*(1+nu));
real lambda = E*nu / ((1+nu)*(1-nu));

// === Geometry ===
border bottom(t=0, L) { x=t; y=0;   label=1; }
border right(t=0, H)  { x=L; y=t;   label=2; }
border top(t=L, 0)    { x=t; y=H;   label=3; }
border left(t=H, 0)   { x=0; y=t;   label=4; }

// === Mesh ===
int n = 40;
mesh Th = buildmesh(bottom(n*10) + right(n) + top(n*10) + left(n));

// === FE Space (P2 for accuracy) ===
fespace Vh(Th, P2);
Vh u1, u2, v1, v2;

// === Solve ===
solve Elasticity([u1, u2], [v1, v2]) =
    int2d(Th)(
        lambda*(dx(u1)+dy(u2))*(dx(v1)+dy(v2))
        + 2*mu*(dx(u1)*dx(v1) + dy(u2)*dy(v2))
        + mu*(dy(u1)+dx(u2))*(dy(v1)+dx(v2))
    )
    - int1d(Th, 2)(P/H * v2)    // Distributed load on right edge
    + on(4, u1=0, u2=0);         // Fixed left edge

// === Post-processing: Von Mises stress ===
fespace Wh(Th, P1);
Wh ux1 = u1, ux2 = u2;  // P2 → P1 projection
Wh sxx = lambda*(dx(u1)+dy(u2)) + 2*mu*dx(u1);
Wh syy = lambda*(dx(u1)+dy(u2)) + 2*mu*dy(u2);
Wh sxy = mu*(dy(u1)+dx(u2));
Wh sigmaVM = sqrt(sxx^2 + syy^2 - sxx*syy + 3*sxy^2);

// === Output ===
cout << "Max u1 = " << u1[].max << ", Min u1 = " << u1[].min << endl;
cout << "Max u2 = " << u2[].max << ", Min u2 = " << u2[].min << endl;
cout << "Max Von Mises = " << sigmaVM[].max << endl;
cout << "ndof_P1 = " << Wh.ndof << endl;

savemesh(Th, "mesh.msh");
{ ofstream f("u1.txt"); for(int i=0; i<Wh.ndof; i++) f << ux1[][i] << endl; }
{ ofstream f("u2.txt"); for(int i=0; i<Wh.ndof; i++) f << ux2[][i] << endl; }
{ ofstream f("vonmises.txt"); for(int i=0; i<Wh.ndof; i++) f << sigmaVM[][i] << endl; }
```

### Variations

**Plate with hole** — Add circular border:
```freefem
real cx = L/2, cy = H/2, R = H/4;
border hole(t=0, 2*pi) { x=cx+R*cos(t); y=cy+R*sin(t); label=5; }
mesh Th = buildmesh(bottom(n*10) + right(n) + top(n*10) + left(n) + hole(-n*3));
// Note: negative point count for hole = clockwise = interior hole
```

**Body force (gravity)**:
```freefem
real rho = 7800; real g = 9.81;
// Add to solve: - int2d(Th)(rho*g*v2)
```

**Pressure load on surface**:
```freefem
real pressure = 1e5;
// Add: - int1d(Th, 3)(-pressure * v2)  // Normal to top surface
```

### Common Materials

| Material | E [GPa] | ν |
|----------|---------|---|
| Steel | 210 | 0.3 |
| Aluminum | 70 | 0.33 |
| Copper | 120 | 0.34 |
| Titanium | 116 | 0.32 |
| Concrete | 30 | 0.2 |

### Analytical Verification

**Cantilever beam tip deflection** (Euler-Bernoulli):
δ = PL³/(3EI), where I = bH³/12 (b = unit thickness for plane stress)

**Stress at fixed end** (beam theory):
σ_max = M*c/I = P*L*(H/2)/I
