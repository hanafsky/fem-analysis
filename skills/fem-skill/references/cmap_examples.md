# CMAP Official Examples — Shape & Topology Optimization

Source: Allaire et al., CMAP Ecole Polytechnique (MAP 562 course)
URL: http://www.cmap.polytechnique.fr/~allaire/freefem_en.html

These are the three official, battle-tested FreeFEM++ scripts from the CMAP toolbox.
They serve as canonical reference implementations. Claude should adapt these
patterns rather than writing optimization code from scratch.

---

## Example 1: Geometric Shape Optimization (cantilever_en.edp)

**What it does**: Moves the free boundary of a cantilever to minimize compliance,
keeping the volume constant via a Lagrange multiplier.

**Key techniques**:
- Two-mesh strategy: fine mesh `Sh` for elasticity, coarse mesh `Th` for geometry
- Shape gradient = strain energy density on the free boundary
- `extension` problem: regularizes boundary gradient into interior velocity field
- `movemesh` with `checkmovemesh` for safe deformation
- `adaptmesh` at each iteration for mesh quality
- Adaptive step control: increase if gradient aligns, decrease/backstep if not

**Setup (cantilever):**
- Domain ~9×8 with fixed support on left (2 patches, label=1)
- Point load on right edge (label=3): g1=0, g2=-1
- Free boundary (label=2): top, bottom curves + circular holes as initial shape
- Material: E=15, nu=0.35

**Algorithm skeleton:**

```freefem
// 1. Solve elasticity on fine mesh Sh
elasticity;

// 2. Compute shape gradient = energy density on boundary (label=2)
energ = 2*mu*(dx(u)^2+dy(v)^2+((dx(v)+dy(u))^2)/2)
        + lambda*(dx(u)+dy(v))^2;

// 3. Update Lagrange multiplier (volume constraint)
muv = 0.5*muv + 0.5*∫(energ on ∂Ω_free)/∫(1 on ∂Ω_free)
      + pasmu*(volume-volume0)/volume0;
gradient = (energ - muv) * cut;   // cut=0 on fixed/loaded boundaries

// 4. Extend gradient into interior (regularization)
solve extension([eu,ev],[w,s]) =
    ∫(∇eu·∇w + ∇ev·∇s + eu·w + ev·s)
    - ∫_{∂Ω_free}((w·nx+s·ny) * gradient)
    + on(fixed+loaded, eu=0, ev=0);

// 5. Move coarse mesh, check for inverted triangles
minaire = checkmovemesh(Th, [x,y])/100;
while (checkmovemesh(Th, [x+pas*eui, y+pas*evi]) < minaire)
    pas = pas/2;
Th = movemesh(Th, [x+pas*eui, y+pas*evi]);

// 6. Adapt meshes
Sh = adaptmesh(Th, [u,v], [eu,ev], err=0.05);
Th = adaptmesh(Sh, 5, 0, 5, IsMetric=1, hmin=0.05);
```

**Output pattern for visualization:**
```freefem
savemesh(Sh, save+iter+".msh");
// Energy density for coloring:
ofstream file(save+iter+".bb");
file << energ[].n << "\n";
for (int j=0; j<energ[].n; j++) file << energ[][j] << endl;
```

---

## Example 2: Level-Set Method (levelset-cantilever.edp)

**What it does**: Represents shape via level-set φ (solid where φ<0),
advects φ along shape gradient velocity with reinitialization.

**Key techniques:**
- Ersatz material: `Xapprox = X + weak*(1-X)` (weak=0.0001 for void)
- Level-set advection via characteristics: `convect([-V*N1,-V*N2], T, phi)`
- Reinitialization to signed distance: solve |∇φ|=1 iteratively
- Velocity regularization: Helmholtz smoothing of shape gradient
- Line search: accept step only if objective decreases

**Setup:**
- Domain [-1,1]×[-0.5,0.5], Dirichlet on left, load on right center
- Initial shape: perforated plate `phi0 = -0.1 + sin(πkx·x)*sin(πky·(y-0.5))`
- Material: mu=8, lambda=1, weak phase = 0.0001

**Algorithm skeleton:**

```freefem
// 1. Characteristic function from level-set
X = (phi < 0);
Xapprox = X + weak*(1-X);

// 2. Solve elasticity with ersatz material
solve elasticity(...) = int2d(Th)(
    2*mu*Xapprox*(...) + lambda*Xapprox*(...)
) - int1d(Th,3)(g1*w+g2*s) + on(2,u=0,v=0);

// 3. Shape gradient = energy density - Lagrange multiplier
V = 2*mu*Xapprox*(dx(u)^2+dy(v)^2+((dy(u)+dx(v))^2)/2)
    + lambda*Xapprox*(dx(u)+dy(v))^2 - lag*X;

// 4. Regularize velocity (Helmholtz smoothing)
solve smoothing(Vreg,vv) =
    int2d(Th)(eps2*(dx(Vreg)*dx(vv)+dy(Vreg)*dy(vv)) + Vreg*vv)
    - int2d(Th)(V*vv);

// 5. Compute normal from level-set gradient
nabla = dx(phi)^2 + dy(phi)^2;
N1 = dx(phi)/sqrt(nabla+eps1^2);
N2 = dy(phi)/sqrt(nabla+eps1^2);

// 6. Advect level-set
phim = convect([-V*N1,-V*N2], T, phi);

// 7. Reinitialize to signed distance function
for (iterinit=1; iterinit<Ninit; iterinit++) {
    S = phim/sqrt(phim^2 + h^2*(dx(phim)^2+dy(phim)^2));
    M1 = dx(phim)/sqrt(dx(phim)^2+dy(phim)^2+eps1^2);
    M2 = dy(phim)/sqrt(dx(phim)^2+dy(phim)^2+eps1^2);
    phiinit = convect([-S*M1,-S*M2], TT, phim) + TT*S;
    phim = phiinit;
}

// 8. Accept/reject based on objective decrease
if (objectivem <= objective*(1+tol)) { accept } else { T=T/2; reject }
```

---

## Example 3: Homogenization / SIMP (cantilever.homog.struct.edp)

**What it does**: Optimizes material density θ using homogenization theory
(Hashin-Shtrikman bounds for rank-2 laminates), with composite penalization.

**Key techniques:**
- Hashin-Shtrikman optimal microstructure bounds
- Rank-2 laminate orientation from stress eigenvectors
- Full anisotropic elasticity tensor (6 independent components in 2D: A,B,C,D,E,F)
- Tensor stored as 3×3 symmetric matrix using Voigt notation
- Lagrange multiplier by bisection (dichotomy) for exact volume constraint
- Penalization: `theta = ((1-cos(π·θ))/2 + p·θ)/(p+1)` after niternp iterations
- Two-phase: strong material + weak phase (epsil=0.01)

**Setup:**
- Structured mesh on [-1,1]×[0,1], Dirichlet on left, load at center-right
- Y=1, nu=0.3, target volume fraction = 0.4
- 40 iterations without penalization + 20 with penalization

**Density update rule (analytical for compliance with HS bounds):**

```freefem
// Optimal density from stress eigenvalues
theta1 = sqrt(mu+kappa)/sqrt(4*mu*kappa) * (|vp1| + |vp2|);
theta = theta1 / sqrt(lagrange);
theta = min(1, theta);

// Lagrange multiplier by bisection to satisfy volume constraint
while (|1 - weight/weight0| > 1e-6) {
    lagrange = (lagmax+lagmin)/2;
    theta = min(1, theta1/sqrt(lagrange));
    weight = ∫θ / volume;
    if (weight < weight0) lagmin=lagrange; else lagmax=lagrange;
}
```

**Anisotropic tensor update (Hashin-Shtrikman rank-2 laminate):**

The homogenized tensor depends on:
- Stress eigenvectors (e1,e2): lamination directions
- Proportions m1,m2: derived from |vp2|/(|vp1|+|vp2|) and vice versa
- Uses double inversion formula for HS upper bound
- Adds weak phase: `A = epsil/(1-epsil)*(λ+2μ) + θ*A_laminate`

---

## Choosing the Right Method

| Problem | Recommended Method | Example File |
|---|---|---|
| Optimize shape of existing design (no new holes) | Geometric | cantilever_en.edp |
| Find optimal material layout from scratch | Homogenization/SIMP | cantilever.homog.struct.edp |
| Allow topology changes (new holes, merging) | Level-Set | levelset-cantilever.edp |
| Quick prototyping with simple OC update | SIMP (see topology_optimization.md) | custom template |
| Nonlinear / multi-physics optimization | AutoFreeFEM | https://gitlab.tugraz.at/autofreefem/autofreefem |

## Adapting to Custom Problems

To adapt these examples to a new problem, typically modify:
1. **Geometry**: border definitions and buildmesh
2. **Boundary conditions**: labels and `on(...)` clauses
3. **Loading**: `int1d(Th, label)(force*v)` terms
4. **Material**: E, nu (Lamé parameters)
5. **Volume constraint**: target fraction and Lagrange multiplier tuning
6. **Output**: add savemesh + ofstream for PyVista visualization pipeline
