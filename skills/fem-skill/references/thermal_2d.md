# 2D Thermal Analysis Reference (FreeFEM++)

## Steady-State Heat Conduction

### Governing Equation
-∇·(k∇T) = Q  (Poisson equation with thermal conductivity k and heat source Q)

### Template: Heat Conduction with Mixed BCs

```freefem
// === Parameters ===
real Lx = 1.0, Ly = 0.5;    // Domain size [m]
real k = 50.0;                // Thermal conductivity [W/(m·K)]
real Q = 1e4;                 // Volumetric heat source [W/m³]
real Thot = 100.0;            // Prescribed temperature [°C]
real h = 10.0;                // Convection coefficient [W/(m²·K)]
real Tinf = 20.0;             // Ambient temperature [°C]

// === Geometry ===
border bottom(t=0, Lx)  { x=t;   y=0;   label=1; }
border right(t=0, Ly)   { x=Lx;  y=t;   label=2; }
border top(t=Lx, 0)     { x=t;   y=Ly;  label=3; }
border left(t=Ly, 0)    { x=0;   y=t;   label=4; }

int n = 30;
mesh Th = buildmesh(bottom(n*2) + right(n) + top(n*2) + left(n));

fespace Vh(Th, P1);
Vh T, v;

// === Solve ===
// Left: Dirichlet (T=Thot), Right: Convection (Robin), Top/Bottom: Insulated (Neumann=0)
solve HeatConduction(T, v) =
    int2d(Th)(k*(dx(T)*dx(v) + dy(T)*dy(v)))  // Diffusion
    + int1d(Th, 2)(h*T*v)                       // Convection (Robin BC)
    - int2d(Th)(Q*v)                             // Heat source
    - int1d(Th, 2)(h*Tinf*v)                     // Convection reference
    + on(4, T=Thot);                              // Dirichlet BC

// === Output ===
cout << "Max T = " << T[].max << endl;
cout << "Min T = " << T[].min << endl;
cout << "ndof = " << Vh.ndof << endl;

savemesh(Th, "mesh.msh");
{ ofstream f("temperature.txt"); for(int i=0; i<Vh.ndof; i++) f << T[][i] << endl; }
```

### Template: Transient Heat Conduction

```freefem
real rho = 7800;   // Density [kg/m³]
real cp = 500;     // Specific heat [J/(kg·K)]
real dt = 0.01;    // Time step [s]
int nsteps = 100;  // Number of time steps

fespace Vh(Th, P1);
Vh T, Told, v;
Told = Tinit;  // Initial condition

for (int n = 0; n < nsteps; n++) {
    solve HeatTransient(T, v) =
        int2d(Th)(rho*cp/dt * T*v)
        + int2d(Th)(k*(dx(T)*dx(v) + dy(T)*dy(v)))
        - int2d(Th)(rho*cp/dt * Told*v)
        - int2d(Th)(Q*v)
        + on(4, T=Thot);
    Told = T;
    
    // Save snapshots at intervals
    if (n % 10 == 0) {
        ofstream f("temperature_" + n + ".txt");
        for(int i=0; i<Vh.ndof; i++) f << T[][i] << endl;
    }
}
```

### Common Thermal Conductivities

| Material | k [W/(m·K)] | ρ [kg/m³] | cp [J/(kg·K)] |
|----------|-------------|-----------|----------------|
| Copper | 400 | 8960 | 385 |
| Aluminum | 237 | 2700 | 900 |
| Steel | 50 | 7800 | 500 |
| Concrete | 1.4 | 2300 | 880 |
| Air | 0.026 | 1.2 | 1005 |

### Recommended Colormap
Use `"coolwarm"` or `"hot"` for temperature fields.
