---
name: fem-analysis
disable-model-invocation: true
description: |
  FEM (Finite Element Method) analysis skill. Build models, generate meshes, and run analyses with FreeFEM++, then visualize results with PyVista.
  Automates the full pipeline: natural language problem description -> FreeFEM++ code (.edp) generation -> solver execution -> result visualization -> AI evaluation.
  Use this skill when the user mentions: FEM, FEA, finite element method, structural analysis, stress analysis, thermal analysis, elasticity,
  cantilever, stress analysis, thermal analysis, Poisson equation, elasticity,
  FreeFEM, freefem++, mesh, mesh generation, Von Mises, deformation analysis, eigenvalue analysis,
  topology optimization, shape optimization, SIMP,
  level-set, compliance minimization, structural optimization,
  density method, optimality criteria, sensitivity filter, AutoFreeFEM,
  H1 gradient method, traction method, Azegami,
  domain variation, shape Hessian, mean compliance minimization,
  or any request involving solving PDEs on 2D/3D domains with finite elements.
  Also trigger when the user uploads geometry descriptions, asks about structural mechanics problems,
  wants to visualize stress/displacement/temperature fields, or wants to optimize a structure's shape or topology.
---

# FEM Analysis Skill (FreeFEM++ + PyVista)

## Overview

End-to-end FEM analysis pipeline (macOS / Linux / WSL):
1. **Model & Mesh**: Generate FreeFEM++ `.edp` code from natural language
2. **Solve**: Execute FreeFEM++ and capture results
3. **Visualize**: PyVista rendering for publication-quality images
4. **Evaluate**: AI reviews visualization (self-check loop)

## Platform Notes

### macOS

#### FreeFEM++ on macOS

FreeFEM++ is installed via Homebrew. The binary name may differ by version:

```bash
# Check which is available
which FreeFem++ || which freefem++
```

**Homebrew paths (Apple Silicon):** `/opt/homebrew/bin/FreeFem++`
**Homebrew paths (Intel Mac):** `/usr/local/bin/FreeFem++`

#### PyVista on macOS

macOS has native OpenGL support — no `xvfb` needed.

**For off-screen rendering** (Claude Code sessions, CI, SSH):
```python
import os
os.environ["PYVISTA_OFF_SCREEN"] = "true"  # BEFORE importing pyvista
import pyvista as pv
pv.OFF_SCREEN = True
```

**For interactive display** (direct terminal use):
PyVista auto-detects display on macOS — no special setup needed.

#### matplotlib on macOS

```python
# Off-screen (Claude Code / SSH):
import matplotlib
matplotlib.use("Agg")

# Interactive (iTerm2 / Terminal.app):
# auto-selects "MacOSX" backend — no config needed
```

### Linux / WSL

#### FreeFEM++ on Linux/WSL

```bash
sudo apt update && sudo apt install -y freefem++
```

Binary path: `/usr/bin/FreeFem++` (apt version)

If the apt version is too old or unavailable, download official binaries from:
https://github.com/FreeFem/FreeFem-sources/releases

#### PyVista on Linux/WSL (headless)

WSL and headless Linux require `xvfb` for off-screen rendering:

```bash
sudo apt install -y libgl1-mesa-glx xvfb
```

**For off-screen rendering:**
```python
import os
os.environ["PYVISTA_OFF_SCREEN"] = "true"  # BEFORE importing pyvista
import pyvista as pv
pv.OFF_SCREEN = True
```

If rendering still fails, wrap the command with `xvfb-run`:
```bash
xvfb-run python src/visualize.py
```

#### matplotlib on Linux/WSL

```python
import matplotlib
matplotlib.use("Agg")  # Always use Agg in headless environments
```

## Available Tools

- FreeFEM++ 4.x — PDE solver (official .dmg on macOS / `sudo apt install freefem++` on Linux)
- gmsh — Advanced mesh generation (`brew install gmsh` / `sudo apt install gmsh`)
- Python 3.10+ — via uv or system package manager
- PyVista + VTK — 3D visualization (native arm64 wheels on Apple Silicon)
- meshio — Mesh format conversion
- Node.js 18+ — React Viewer (via fnm or package manager)

## Workflow

### Step 1: Problem Understanding

Clarify with user if needed: geometry, physics type, materials, BCs, desired output.

### Step 2: Generate FreeFEM++ Code

**Read the appropriate reference first:**
- `references/elasticity_2d.md` — Plane stress/strain
- `references/thermal_2d.md` — Heat conduction
- `references/poisson_2d.md` — Poisson/Laplace
- `references/stokes_2d.md` — Stokes flow
- `references/topology_optimization.md` — SIMP, Level-Set, Topological Derivative
- `references/shape_optimization.md` — Parametric, Geometric, FreeFEM optimizers (BFGS/IPOPT/NLopt)
- `references/h1_gradient_method.md` — Azegami H1 gradient method (Traction Method): gradient/Newton/theta methods, domain/boundary integral formulations
- `references/cmap_examples.md` — CMAP official toolbox code examples (Geometric / Level-Set / Homogenization)
- `references/3d_gmsh_freefem.md` — 3D Gmsh + FreeFEM++ (gmshload3, OCC workflow, volume classification)

**Output directory**: Always create a `results/` subdirectory and write all outputs there to keep the workspace clean:
```freefem
system("mkdir -p results");
savemesh(Th, "results/mesh.msh");
{ ofstream f("results/solution.txt"); for(int i=0; i<Vh.ndof; i++) f << u[][i] << endl; }
```

**CRITICAL: iovtk plugin is NOT available on this system.**
Always use text-based output as shown above.

**P2→P1 projection**: When solving with P2, project before saving:
```freefem
fespace Wh(Th, P1);
Wh ux = u;  // P2 → P1 projection
```

### Step 3: Execute

```bash
mkdir -p results
FreeFem++ problem.edp 2>&1
```

If not in PATH:
```bash
/opt/homebrew/bin/FreeFem++ problem.edp 2>&1    # macOS Apple Silicon
/usr/local/bin/FreeFem++ problem.edp 2>&1        # macOS Intel
/usr/bin/FreeFem++ problem.edp 2>&1              # Linux/WSL (apt)
```

### Step 4: Visualize

Use `scripts/ff_mesh_reader.py` for mesh parsing. Choose the renderer based on the task:

**matplotlib (preferred for 2D)**:
- Scalar field contour plots (temperature, stress, velocity, etc.)
- Streamline plots (Stokes flow)
- Better layout control, cleaner labels, publication-quality output
- Use `matplotlib.tri.Triangulation` for triangle mesh plotting

**PyVista (use when needed)**:
- Deformed mesh overlay (original wireframe + deformed shape) via `scripts/visualize.py`
- 3D visualization or interactive rotation
- Multi-panel comparison views
- Set `os.environ["PYVISTA_OFF_SCREEN"] = "true"` BEFORE importing pyvista

Read input from and save images to `results/` directory.

### Step 5: AI Self-Evaluation Loop

1. `view` the screenshot
2. Check: deformation direction, stress concentration locations, BC satisfaction, mesh adequacy
3. If issues → fix .edp → re-run → re-visualize
4. Present final result

### Step 6: Interactive 3D (Optional)

Export mesh+solution as JSON → React artifact with Three.js/Plotly for rotation/zoom.

## Common Gotchas

1. **iovtk unavailable** → always use text output + Python conversion
2. **2D visualization**: Prefer matplotlib over PyVista — better layout, cleaner scalar bars, publication-quality output. PyVista's off-screen rendering tends to produce small meshes with cut-off labels in 2D.
3. **Set PYVISTA_OFF_SCREEN** before import (only when using PyVista)
4. Auto-scale deformation: `scale = 0.1 * L / max_displacement`
4. FreeFem .msh is NOT gmsh format — use custom parser in scripts/
5. Boundary labels: 1=bottom, 2=right, 3=top, 4=left (rectangle convention)
6. **Homebrew binary name**: Check `FreeFem++` vs `freefem++` (case varies by version)
7. **Apple Silicon (arm64)**: VTK/PyVista wheels are native — no Rosetta needed
8. **macOS Gatekeeper**: First run may need `xattr -d com.apple.quarantine $(which FreeFem++)`
9.  **3D Gmsh: use OCC solids, not Extrude** — Box/Cylinder + BooleanFragments for stable surface IDs. Extrude causes gmshload3 "Border element not in mesh" errors.
10. **3D: never Physical Surface internal interfaces** — Interior faces shared by two tets crash gmshload3. Extract interface data via Python post-processing (z ≈ target filter).
11. **3D Volume classification: use centroid z, not BoundingBox** — `Volume In BoundingBox` returns intersecting (not contained) volumes, causing duplicates.
12. **3D FreeFEM++ coordinates: use `Th(i).x`** — `Wh(i).x` is 2D-only. For mesh3, always use `Th(i).x/y/z`.

## Local Dev Setup

See `references/project_setup.md` for full scaffold instructions (macOS with Homebrew / Linux with apt).
