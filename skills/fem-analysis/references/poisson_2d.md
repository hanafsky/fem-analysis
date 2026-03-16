# 2D Poisson/Laplace Reference (FreeFEM++)

## Poisson Equation: -Δu = f

### Template: Basic Poisson on Square

```freefem
real L = 1.0;
int n = 30;
mesh Th = square(n, n, [L*x, L*y]);

fespace Vh(Th, P1);
Vh u, v;

func f = 2*pi^2*sin(pi*x)*sin(pi*y);  // Source term

solve Poisson(u, v) =
    int2d(Th)(dx(u)*dx(v) + dy(u)*dy(v))
    - int2d(Th)(f*v)
    + on(1, 2, 3, 4, u=0);  // Homogeneous Dirichlet on all boundaries

// Exact solution for verification: u_exact = sin(pi*x)*sin(pi*y)
Vh ue = sin(pi*x)*sin(pi*y);
Vh err = u - ue;
cout << "L2 error = " << sqrt(int2d(Th)(err^2)) << endl;
cout << "Max u = " << u[].max << endl;

savemesh(Th, "mesh.msh");
{ ofstream f("solution.txt"); for(int i=0; i<Vh.ndof; i++) f << u[][i] << endl; }
```

### Template: L-shaped Domain

```freefem
border b1(t=0,1)   { x=t;   y=0;   label=1; }
border b2(t=0,0.5) { x=1;   y=t;   label=2; }
border b3(t=1,0.5) { x=t;   y=0.5; label=3; }
border b4(t=0.5,1) { x=0.5; y=t;   label=4; }
border b5(t=0.5,0) { x=t;   y=1;   label=5; }
border b6(t=1,0)   { x=0;   y=t;   label=6; }

int n = 20;
mesh Th = buildmesh(b1(n) + b2(n/2) + b3(n/2) + b4(n/2) + b5(n/2) + b6(n));
```

### Adaptive Mesh Refinement

```freefem
// Solve → Adapt → Re-solve loop
for (int iter = 0; iter < 3; iter++) {
    solve Poisson(u, v) = int2d(Th)(dx(u)*dx(v)+dy(u)*dy(v)) - int2d(Th)(f*v) + on(1,2,3,4,u=0);
    Th = adaptmesh(Th, u, err=0.01, nbvx=10000);
    cout << "Iter " << iter << ": nv=" << Th.nv << ", nt=" << Th.nt << endl;
}
```
