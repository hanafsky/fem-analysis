# FEM Analysis Skill

A [Claude Code custom skill](https://docs.anthropic.com/en/docs/claude-code/skills) for end-to-end Finite Element Method (FEM) analysis using FreeFEM++ and PyVista.

## What it does

Given a natural language problem description, this skill automates:

1. **Model & Mesh** — Generates FreeFEM++ `.edp` code
2. **Solve** — Executes FreeFEM++ and captures results
3. **Visualize** — Renders publication-quality images with PyVista
4. **Evaluate** — AI reviews the visualization and iterates if needed

## Supported analysis types

- 2D elasticity (plane stress/strain)
- Thermal conduction
- Poisson/Laplace equations
- Stokes flow
- Topology optimization (SIMP, Level-Set)
- Shape optimization (H1 gradient method, parametric)

## Directory structure

```
fem-analysis/
├── skills/
│   └── fem-analysis/
│       ├── SKILL.md              # Skill definition (Claude Code)
│       ├── references/           # FEM theory & code references
│       │   ├── elasticity_2d.md
│       │   ├── thermal_2d.md
│       │   ├── poisson_2d.md
│       │   ├── stokes_2d.md
│       │   ├── topology_optimization.md
│       │   ├── shape_optimization.md
│       │   ├── h1_gradient_method.md
│       │   ├── cmap_examples.md
│       │   └── project_setup.md
│       ├── scripts/              # Python utilities
│       │   ├── ff_mesh_reader.py # FreeFem++ mesh parser
│       │   └── visualize.py      # PyVista rendering
│       └── templates/
│           └── FEMViewer.tsx     # React 3D viewer component
├── example/                      # Usage examples
├── README.md
└── README_JA.md
```

## Platform support

- **macOS** (Apple Silicon / Intel) — FreeFEM++ via official .dmg installer
- **Linux / WSL** — FreeFEM++ via apt

See `skills/fem-analysis/references/project_setup.md` for detailed setup instructions.

## Installation

### 1. Install the skill

Claude Code discovers skills from two locations. Choose one:

**Per-project** (available only in that project):
```bash
# From your project root
mkdir -p .claude/skills
cp -r /path/to/fem-analysis/skills/fem-analysis .claude/skills/
```

**Personal / global** (available in all your projects):
```bash
mkdir -p ~/.claude/skills
cp -r /path/to/fem-analysis/skills/fem-analysis ~/.claude/skills/
```

After copying, verify the structure:
```
.claude/skills/fem-analysis/
├── SKILL.md
├── references/
├── scripts/
└── templates/
```

### 2. Install prerequisites

The skill generates and runs FEM code, so the following tools must be installed on your system.

**macOS:**
```bash
# FreeFEM++ — download .dmg from official releases:
#   https://github.com/FreeFem/FreeFem-sources/releases
# Then install:
sudo cp -rf /Volumes/<mounted-dmg>/FreeFem++.app /Applications/
sudo xattr -rc /Applications/FreeFem++.app
# Add to PATH (adjust version number as needed):
#   bash/zsh: echo 'export PATH="/Applications/FreeFem++.app/Contents/ff-4.15.1/bin:$PATH"' >> ~/.zprofile
#   fish:     fish_add_path /Applications/FreeFem++.app/Contents/ff-4.15.1/bin

# gmsh (optional, for complex geometry)
brew install gmsh

# Python environment (via uv)
curl -LsSf https://astral.sh/uv/install.sh | sh
mkdir fem-workbench && cd fem-workbench
uv init
uv add pyvista meshio numpy vtk matplotlib scipy
```

**Linux / WSL (apt):**
```bash
sudo apt update && sudo apt install -y freefem++ gmsh
# If apt version is too old or unavailable, download from:
#   https://github.com/FreeFem/FreeFem-sources/releases
sudo apt install -y libgl1-mesa-glx xvfb   # for headless rendering

# Python environment (via uv)
curl -LsSf https://astral.sh/uv/install.sh | sh
mkdir fem-workbench && cd fem-workbench
uv init
uv add pyvista meshio numpy vtk matplotlib scipy
```

### 3. Verify

```bash
FreeFem++ --version
python -c "import pyvista; print(pyvista.__version__)"
```

For the full project scaffold (React viewer, vite config, etc.), see [`skills/fem-analysis/references/project_setup.md`](skills/fem-analysis/references/project_setup.md).

## License

MIT
