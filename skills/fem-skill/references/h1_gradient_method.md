# H1 Gradient Method (Traction Method / 力法) — Azegami Method

## Overview

The H1 gradient method (H1勾配法) is a shape optimization technique developed by
Prof. Hideyuki Azegami (畔上秀行, Nagoya University). Also known as the
**Traction Method** (力法) or **Force Method**.

**Key idea**: Instead of directly moving boundary nodes along the shape gradient,
the shape gradient function is treated as a **surface traction (力)** and an auxiliary
elasticity-like BVP is solved over the entire domain to obtain a smooth velocity field.
This naturally produces smooth, mesh-independent shape updates without checkerboard
or oscillation artifacts.

### References
- Azegami, H. "Shape Optimization Problems" Morikita Publishing, 2016 (with sample code)
- OCW: https://ocw.nagoya-u.jp/courses/0799
- OPTISHAPE-TS (commercial solver by Quint Corp.) is based on this theory
- YouTube: https://youtu.be/F3MiA83hgeQ

### Advantages over Level-Set / SIMP
- Produces inherently smooth boundary updates (no regularization filter needed)
- Domain variation type: mesh moves with the shape → accurate boundary representation
- Natural extension to Newton method (H1 Newton法) for superlinear convergence
- Works for both shape optimization and topology optimization (θ method)
- Unified framework for elastic, Stokes, and multi-physics problems

---

## Mathematical Foundation

### Shape Optimization Problem (Algorithm 3.7.2 from the book)

```
min  f₀(Ω) = ∫_Γ_N p·u ds        (mean compliance)
s.t. f₁(Ω) = |Ω| - |Ω₀| = 0     (volume constraint)
     state equation: -div(σ(u)) = 0  in Ω
```

### Shape Gradient (形状勾配関数)

Two equivalent formulations:

**Domain integral type** (Eq. 9.11.19):
```
g₀(ψ) = 2·(σ∇u)ᵀ:∇ψ - (σ:ε)·div(ψ)
```

**Boundary integral type** (Eq. 9.11.43):
```
g₀(ψ) = -(σ:ε)·(n·ψ)    on free boundary Γ_free
```

Volume constraint gradient:
```
g₁(ψ) = div(ψ)           (domain integral type)
g₁(ψ) = n·ψ              (boundary integral type)
```

### H1 Gradient Method — Core Idea

Instead of updating the shape with φ = -g₀ (L² gradient, causes oscillations),
solve for φ in H¹ by the **Robin-type BVP**:

```
Find φ_{g0} ∈ H¹(Ω) such that:
  cₐ·( ∇φ_{g0}:∇ψ + c_Ω·φ_{g0}·ψ ) = -g₀(ψ)   ∀ψ ∈ H¹(Ω)
```

where:
- `cₐ` (typically 100): H1 inner product coefficient (controls step size)
- `c_Ω` (typically 0): mass term coefficient (0 for pure H1, >0 adds L2 regularization)

Similarly for the volume constraint:
```
Find φ_{g1} such that:
  cₐ·( ∇φ_{g1}:∇ψ + c_Ω·φ_{g1}·ψ ) = -g₁(ψ)   ∀ψ
```

### Constrained Update (Eq. 3.7.12)

Lagrange multiplier for volume constraint:
```
λ₁ = -(f₁ + <g₁, φ_{g0}>) / <g₁, φ_{g1}>
```

Combined velocity field:
```
φ_g = φ_{g0} + λ₁·φ_{g1}
```

Update shape: `Ω^{k+1} = (I + φ_g)(Ω^k)`

### Boundary Conditions for the H1 gradient problem

- **Dirichlet boundary (fixed)**: φ = 0 (cannot move)
- **Neumann boundary (loaded)**: φ = 0 (load application point fixed)
- **Free boundary**: no BC (this is where optimization happens)
- **Symmetry boundary**: φ_n = 0 (normal component zero)

---

## Algorithm (H1 Gradient Method)

```
Step 1: Initialize mesh Ω₀, set k=0
Step 2: Solve state equation → u^k, compute f₀^{init}, |Ω₀|
Step 3: Solve H1 gradient problems:
        H1gradientf0 → φ_{g0}    (objective gradient)
        H1gradientf1 → φ_{g1}    (constraint gradient)
Step 4: Compute Lagrange multiplier:
        λ₁ = -(|Ω^k|-|Ω₀| + <g₁,φ_{g0}>) / <g₁,φ_{g1}>
Step 5: Combined update:
        φ_g = φ_{g0} + λ₁·φ_{g1}
Step 6: Move mesh: Th = movemesh(Th, [x+φ_g1, y+φ_g2])
Step 7: Re-solve state, adapt mesh
Step 8: Check convergence: |f₀^{k} - f₀^{k-1}| / |f₀^{init}| < ε
```

---

## Template: 2D Cantilever — Domain Integral, H1 Gradient

Based on: `9.11.5_shape_elastic/2d-cantilever/domain_integral/grad/main.edp`

```freefem
// H1 Gradient Method for 2D Cantilever Compliance Minimization
// Based on Azegami's book, Algorithm 3.7.2

string save="plots";

// Material
real ey=15.;    // Young's modulus
real nup=0.35;  // Poisson's ratio
real lame=ey*nup/((1.+nup)*(1.-2.*nup));
real mu=ey/(2.*(1.+nup));

// Optimization parameters
int kmax=100;
real ca=100.;     // H1 coefficient (larger = smaller steps, more stable)
real cOmega=0.;   // Mass term (0 = pure H1 gradient)
real errelas=0.001;
real f0err=1.e-7;

// --- Macros for strain, stress, shape gradient ---
macro div(u) (dx(u#1)+dy(u#2))//
macro grad(u) [dx(u#1),dx(u#2),dy(u#1),dy(u#2)]//
macro prod(u,v) (u#1*v#1+u#2*v#2)//

macro E11(u) (dx(u#1))//
macro E12(u) ((dx(u#2)+dy(u#1))/2)//
macro E21(u) ((dx(u#2)+dy(u#1))/2)//
macro E22(u) (dy(u#2))//
macro E(u) [E11(u),E12(u),E21(u),E22(u)]//

macro S11(u) (lame*div(u)+2*mu*dx(u#1))//
macro S12(u) (2*mu*((dx(u#2)+dy(u#1))/2))//
macro S21(u) (2*mu*((dx(u#2)+dy(u#1))/2))//
macro S22(u) (lame*div(u)+2*mu*dy(u#2))//
macro S(u) [S11(u),S12(u),S21(u),S22(u)]//

// σ∇u product (needed for domain integral form)
macro Sgut11(u) (S11(u)*dx(u#1)+S12(u)*dx(u#2))//
macro Sgut12(u) (S11(u)*dy(u#1)+S12(u)*dy(u#2))//
macro Sgut21(u) (S21(u)*dx(u#1)+S22(u)*dx(u#2))//
macro Sgut22(u) (S21(u)*dy(u#1)+S22(u)*dy(u#2))//
macro Sgut(u) [Sgut11(u),Sgut12(u),Sgut21(u),Sgut22(u)]//

// Shape gradients (domain integral type)
macro g0(psi) (2*(Sgut(u)'*grad(psi))-(S(u)'*E(u))*(div(psi)))//
macro g1(psi) (div(psi))//

// --- Boundary labels ---
int Dirichlet=1, Neumann=2, Free=3, Middle=4;

// --- Geometry (half-beam with symmetry) ---
real width=10, height=1;
border b1(t=0,width)       {x=t;y=0;label=Middle;}
border b2(t=0,0.05*height) {x=width;y=t;label=Neumann;}
border b3(t=0.05*height,0.5*height) {x=width;y=t;label=Free;}
border b4(t=width,0)       {x=t;y=0.5*height;label=Free;}
border b5(t=0.5*height,0)  {x=0;y=t;label=Dirichlet;}

mesh Th=buildmesh(b1(100)+b2(10)+b3(10)+b4(100)+b5(10));

fespace Vh(Th,[P2,P2]), Eh(Th,P1);

Vh [u1,u2],[v1,v2];  // displacement + test
Vh [phi1,phi2]=[0,0]; // accumulated deformation
Vh [varphig1,varphig2],[varphi01,varphi02],[varphi11,varphi12],[psi1,psi2];
Vh [pn1,pn2]=[0,-1.]; // traction
Eh e; // strain energy density

// State equation
problem elasticity([u1,u2],[v1,v2])
  =int2d(Th)(S(u)'*E(v))
   -int1d(Th,Neumann)(pn1*v1+pn2*v2)
   +on(Dirichlet,u1=0,u2=0);

// H1 gradient for objective f0
problem H1gradientf0([varphi01,varphi02],[psi1,psi2])
  =int2d(Th)(ca*((grad(varphi0)'*grad(psi))+cOmega*prod(varphi0,psi)))
   +int2d(Th)(g0(psi))
   +on(Dirichlet,varphi01=0)
   +on(Neumann,varphi01=0,varphi02=0)
   +on(Middle,varphi02=0);

// H1 gradient for constraint f1 (volume)
problem H1gradientf1([varphi11,varphi12],[psi1,psi2])
  =int2d(Th)(ca*((grad(varphi1)'*grad(psi))+cOmega*prod(varphi1,psi)))
   +int2d(Th)(g1(psi))
   +on(Dirichlet,varphi11=0)
   +on(Neumann,varphi11=0,varphi12=0)
   +on(Middle,varphi12=0);

// --- Initialization ---
elasticity;
Th = adaptmesh(Th,[u1,u2],err=errelas);
real meas, measinit, f0, f0init, f0prev;
meas = int2d(Th)(1.);
elasticity;
measinit=meas;
f0=int1d(Th,Neumann)(prod(pn,u));
f0init=f0;

// --- Optimization loop ---
for(int k=1; k<kmax+1; k++){
  f0prev=f0;

  // Step 3-4: solve H1 gradient problems
  H1gradientf0;
  H1gradientf1;

  // Step 5: Lagrange multiplier
  real product0=int2d(Th)(g1(varphi0));
  real product1=int2d(Th)(g1(varphi1));
  real lambda1=-(meas-measinit+product0)/product1;

  // Step 6: combined velocity + move mesh
  [varphig1,varphig2]=[varphi01+lambda1*varphi11,varphi02+lambda1*varphi12];
  [phi1,phi2]=[phi1,phi2]+[varphig1,varphig2];
  Th=movemesh(Th,[x+varphig1,y+varphig2]);

  // Step 7: re-solve
  elasticity;
  Th=adaptmesh(Th,[u1,u2],err=errelas);
  meas=int2d(Th)(1.);
  elasticity;
  f0=int1d(Th,Neumann)(prod(pn,u));

  cout<<"k="<<k<<" f0/f0init="<<f0/f0init
      <<" meas/measinit="<<meas/measinit<<endl;

  // Step 8: convergence check
  if(abs(f0prev-f0)/f0init < f0err) break;
}

// Output
savemesh(Th, "result.msh");
e=0.5*(S(u)'*E(u));
{ ofstream f("result_energy.txt");
  for(int i=0;i<Eh.ndof;i++) f<<e[][i]<<endl; }
```

---

## Template: Boundary Integral Type

The difference from domain integral: the shape gradient is evaluated **only on the free boundary**.

```freefem
// Shape gradients (boundary integral type, Eq. 9.11.43, 9.11.52)
macro g0(psi) (-(S(u)'*E(u))*(N.x*psi#1+N.y*psi#2))//
macro g1(psi) (N.x*psi#1+N.y*psi#2)//

// H1 gradient for f0: RHS is boundary integral on Free
problem H1gradientf0([varphi01,varphi02],[psi1,psi2])
  =int2d(Th)(ca*(grad(varphi0)'*grad(psi)+cOmega*prod(varphi0,psi)))
   +int1d(Th,Free)(g0(psi))     // ← boundary integral, not int2d!
   +on(Dirichlet,varphi01=0)
   +on(Neumann,varphi01=0,varphi02=0)
   +on(Middle,varphi02=0);

// Lagrange multiplier: also boundary integral
product0 = int1d(Th,Free)(g1(varphi0));
product1 = int1d(Th,Free)(g1(varphi1));
```

**When to use which:**
- Domain integral: more stable, works well when mesh adapts aggressively
- Boundary integral: theoretically cleaner, but requires good boundary mesh quality

---

## H1 Newton Method (H1 Newton法)

For faster convergence near optimum, switch from gradient to Newton after kN iterations.
The Newton method uses the **shape Hessian** h₀ and h₁.

**Key additions for Newton:**

```freefem
// Additional macros for Hessian
macro SE(u,v) [S11(u)*E11(v)+S12(u)*E21(v),
               S11(u)*E12(v)+S12(u)*E22(v),
               S21(u)*E11(v)+S22(u)*E21(v),
               S21(u)*E12(v)+S22(u)*E22(v)]//

macro gugvt(u,v) [dx(u#1)*dx(v#1)+dx(u#2)*dx(v#2),
                  dx(u#1)*dy(v#1)+dx(u#2)*dy(v#2),
                  dy(u#1)*dx(v#1)+dy(u#2)*dx(v#2),
                  dy(u#1)*dy(v#1)+dy(u#2)*dy(v#2)]//

// Hessian integrands (Eq. 9.11.38, 9.11.39)
macro h0(psi,varphi) (
  (S(u)'*E(u))*(gradt(varphi)'*grad(psi)+div(varphi)*div(psi))
  +(Sgut(u)'*(gugvt(psi,varphi)+gugvt(varphi,psi)))
  -2.0*(SE(u,u)'*(grad(varphi)*div(psi)+grad(psi)*div(varphi)))
)//

macro h1(psi,varphi) (-(gradt(varphi)'*grad(psi))+div(varphi)*div(psi))//

// Newton problem for f0
real ch=0.01;   // Hessian coefficient
real cD1=100.;  // H1 regularization (Newton phase)
real cD0=0.;    // Mass term (Newton phase)

problem H1Newtonf0([varphi01,varphi02],[psi1,psi2])
  =int2d(Th)(ch*(h0(psi,varphi0)+lambda1*h1(psi,varphi0)))
   +int2d(Th)(cD1*(grad(varphi0)'*grad(psi))+cD0*prod(varphi0,psi))
   +int2d(Th)(g0(psi))
   +on(Dirichlet,varphi01=0)
   +on(Neumann,varphi01=0,varphi02=0)
   +on(Middle,varphi02=0);

// Typical strategy: gradient for k < kN, Newton for k >= kN
int kN=20;
for(k=1; k<kN; k++) { /* H1 gradient loop */ }
for(k=kN; k<kmax+1; k++) { /* H1 Newton loop */ }
```

---

## Available Example Problems

From the book's sample code (061461src.zip):

### Shape Optimization (Chapter 9.11.5)

| Problem | Geometry | Methods Available |
|---|---|---|
| **2d-cantilever** | Beam: 10×0.5, fixed left, load right | boundary_integral/grad, domain_integral/grad, domain_integral/Newton |
| **2d-L-shape** | L-shaped domain, fixed top-left, load right | domain_integral/grad, domain_integral/Newton |
| **2d-hole** | Square with circular hole, tension | boundary_integral/grad, domain_integral/grad, domain_integral/Newton |
| **2d-hook** | Hook shape, 2 Dirichlet zones + load | boundary_integral/grad, domain_integral/grad, domain_integral/Newton |
| **3d-cantilever** | 3D beam (requires Gmsh mesh) | domain_integral/grad, domain_integral/Newton |

### Shape Optimization — Stokes (Chapter 9.12.5)

| Problem | Description |
|---|---|
| **2d-iso-body** | Isolated body in Stokes flow, minimize drag |
| (domain_integral/grad, Newton) | Uses P2-P1 Taylor-Hood elements |

### Topology Optimization (Chapter 8.8.5, 8.9.5)

| Problem | Description |
|---|---|
| **8.8.5 topo_elastic / 2d-cantilever** | θ-type H1 gradient for density optimization |
| **8.9.5 topo_stokes / 2d-iso-body** | Stokes flow topology optimization |
| **8.9.5 topo_stokes / 2d-Y-tube** | Y-tube flow channel topology optimization |

---

## Topology Optimization — θ Method (H1 gradient of θ type)

A different formulation where the design variable is a density field θ(x) ∈ [0,1],
using the same H1 regularization concept.

```freefem
// Penalization: E(θ) = θ^α * E₀
real alpha=3.0;

// Density field and its derivative
fespace Wh(Th, P1);
Wh theta, phi, dphi;  // phi=smoothed density, dphi=φ'

// Shape gradient for θ-type (Eq. 8.8.14, 8.8.16)
macro g0(vartheta) (-alpha*phi^(alpha-1)*dphi*(S(u)'*E(u))*vartheta)//
macro g1(vartheta) (dphi*vartheta)//

// H1 gradient for θ
problem H1gradientf0(varphi0, psi)
  =int2d(Th)(ca*(gradtheta(varphi0)'*gradtheta(psi)) + cD*varphi0*psi)
   +int2d(Th)(g0(psi))
   + /* boundary conditions */;

// Update: theta = theta + vartheta_g
// Then project: theta = max(0, min(1, theta))
// phi = smoothed version of theta (convolution or additional BVP)
```

---

## Parameter Tuning Guide

| Parameter | Typical Value | Effect |
|---|---|---|
| `ca` | 10–200 | H1 coefficient. **Larger = smaller steps, more stable**. Start with 100 |
| `cOmega` | 0–1 | Mass term. 0 for most problems. >0 for L-shape to prevent rigid body motion |
| `cD1` (Newton) | 100 | H1 regularization in Newton phase |
| `ch` (Newton) | 0.01–0.1 | Hessian weight. Too large → divergence |
| `kN` | 10–30 | When to switch from gradient to Newton |
| `errelas` | 0.001–0.01 | Mesh adaptation error threshold |
| `f0err` | 1e-7 | Convergence threshold |

### Tips
- If mesh tangling occurs: increase `ca` or add mesh quality check with `checkmovemesh`
- If convergence is slow: decrease `ca` (more aggressive steps)
- For 3D problems: use `cOmega > 0` and coarser mesh initially
- L-shape needs `cOmega=1` (otherwise rigid body rotation occurs)
- Always run `adaptmesh` after `movemesh` to maintain mesh quality

---

## AI Agent Workflow

```
User: "Optimize the shape of this bracket to minimize compliance"
    ↓
AI: Read this reference → choose domain_integral/grad template
    ↓
AI: Map user's geometry to border definitions + boundary labels
    ↓
AI: Set material, loads, BCs in template
    ↓
AI: Generate .edp with H1 gradient method
    ↓
User: Run locally in FreeFEM++
    ↓
AI: Parse kf0.d for convergence, visualize final mesh with PyVista
    ↓
AI: Check: Is f0 decreasing? Is volume preserved? Is mesh valid?
```
