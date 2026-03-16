# 2D Stokes Flow Reference (FreeFEM++)

## Stokes Equations
-μΔu + ∇p = f
∇·u = 0

### Template: Lid-Driven Cavity

```freefem
real L = 1.0;
real mu = 1.0;   // Dynamic viscosity

int n = 30;
mesh Th = square(n, n, [L*x, L*y]);

// Taylor-Hood elements: P2 for velocity, P1 for pressure
fespace Vh(Th, [P2, P2, P1]);
Vh [u1, u2, p], [v1, v2, q];

solve Stokes([u1, u2, p], [v1, v2, q]) =
    int2d(Th)(
        mu*(dx(u1)*dx(v1) + dy(u1)*dy(v1)
          + dx(u2)*dx(v2) + dy(u2)*dy(v2))
        - p*(dx(v1) + dy(v2))
        - q*(dx(u1) + dy(u2))
        - 1e-10*p*q  // Pressure stabilization
    )
    + on(1, 2, 4, u1=0, u2=0)     // No-slip walls
    + on(3, u1=1, u2=0);           // Lid velocity

// Output on P1 mesh
fespace Wh(Th, P1);
Wh ux1=u1, ux2=u2, pp=p;
Wh speed = sqrt(u1^2 + u2^2);

cout << "Max speed = " << speed[].max << endl;
cout << "Max pressure = " << pp[].max << endl;

savemesh(Th, "mesh.msh");
{ ofstream f("u1.txt"); for(int i=0; i<Wh.ndof; i++) f << ux1[][i] << endl; }
{ ofstream f("u2.txt"); for(int i=0; i<Wh.ndof; i++) f << ux2[][i] << endl; }
{ ofstream f("pressure.txt"); for(int i=0; i<Wh.ndof; i++) f << pp[][i] << endl; }
{ ofstream f("speed.txt"); for(int i=0; i<Wh.ndof; i++) f << speed[][i] << endl; }
```

### Recommended Visualization
- Velocity magnitude: `cmap="plasma"`
- Pressure field: `cmap="RdBu_r"` (diverging)
- Streamlines: Use PyVista's `streamlines_from_source` after converting to vector field
