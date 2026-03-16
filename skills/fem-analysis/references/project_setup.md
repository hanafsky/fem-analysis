# Local Development Setup

## Project Structure

```
fem-workbench/
├── pyproject.toml          # Python deps (uv)
├── src/
│   ├── ff_mesh_reader.py   # FreeFem++ mesh parser
│   ├── visualize.py        # PyVista visualization
│   └── run_analysis.py     # CLI: edp generation → solve → visualize
├── problems/               # FreeFem++ .edp files (generated or manual)
│   └── cantilever.edp
├── results/                # Output meshes, solutions, images
├── viewer/                 # React + Vite interactive 3D viewer
│   ├── package.json
│   ├── vite.config.ts
│   ├── index.html
│   └── src/
│       ├── App.tsx
│       ├── FEMViewer.tsx   # Three.js mesh + scalar visualization
│       └── main.tsx
└── README.md
```

## Step 1: Homebrew Essentials

```bash
# Homebrew itself (skip if already installed)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# FreeFEM++
brew install freefem

# gmsh (optional, for complex 3D geometry)
brew install gmsh

# Verify
FreeFem++ --version
# If "command not found":
#   ls /opt/homebrew/bin/FreeFem*    # Apple Silicon
#   ls /usr/local/bin/FreeFem*       # Intel
```

**macOS Gatekeeper note:** First run may be blocked.
```bash
xattr -d com.apple.quarantine $(which FreeFem++)
```

## Step 2: Python Environment (uv)

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
source $HOME/.local/bin/env

mkdir fem-workbench && cd fem-workbench
uv init
uv add pyvista meshio numpy vtk matplotlib
uv add --dev pytest ipython
```

**pyproject.toml:**
```toml
[project]
name = "fem-workbench"
version = "0.1.0"
description = "FEM analysis workbench with FreeFEM++ and PyVista"
requires-python = ">=3.10"
dependencies = [
    "pyvista>=0.44",
    "meshio>=5.3",
    "numpy>=1.24",
    "vtk>=9.3",
    "matplotlib>=3.7",
]

[project.scripts]
fem-run = "src.run_analysis:main"
```

**Apple Silicon note:** PyVista/VTK provide native arm64 wheels — no Rosetta needed.

## Step 3: Node.js (fnm recommended)

```bash
brew install fnm
echo 'eval "$(fnm env --use-on-cd --shell zsh)"' >> ~/.zshrc
source ~/.zshrc
fnm install 22
fnm use 22
```

## Step 4: React Viewer (Vite + Three.js)

```bash
cd fem-workbench
npm create vite@latest viewer -- --template react-ts
cd viewer
npm install
npm install three @react-three/fiber @react-three/drei
npm install -D @types/three tailwindcss @tailwindcss/vite
```

**vite.config.ts:**
```typescript
import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

export default defineConfig({
  plugins: [react()],
  server: {
    port: 5173,
    proxy: {
      '/results': {
        target: 'http://localhost:3001',
        changeOrigin: true,
      }
    }
  }
})
```

## Step 5: Workflow

```bash
# 1. Solve
cd results && FreeFem++ ../problems/cantilever.edp 2>&1 | tee solve.log && cd ..

# 2. Visualize (off-screen for Claude Code)
uv run python -c "
import os; os.environ['PYVISTA_OFF_SCREEN']='true'
import sys; sys.path.insert(0, 'src')
from visualize import plot_scalar_field
plot_scalar_field('results/mesh.msh', 'results/vonmises.txt',
                  'results/vonmises.png', 'Von Mises')
"

# 3. Export for React viewer
uv run python -c "
import sys; sys.path.insert(0, 'src')
from visualize import export_for_threejs
export_for_threejs('results/mesh.msh', 'results/vonmises.txt',
                   'viewer/public/result.json', 'Von Mises')
"

# 4. Launch viewer
cd viewer && npm run dev
# Open http://localhost:5173
```

## FEMViewer Component

Loads exported JSON and renders:
- Triangle mesh with vertex coloring
- Colormap legend
- Orbit controls for rotation/zoom
- Toggle between scalar fields
- Deformation animation (if displacement data available)

## macOS Troubleshooting

| Issue | Fix |
|---|---|
| FreeFem++ not found | `echo 'export PATH="/opt/homebrew/bin:$PATH"' >> ~/.zshrc && source ~/.zshrc` |
| Gatekeeper blocks FreeFem++ | `xattr -d com.apple.quarantine $(which FreeFem++)` |
| VTK import error (Apple Silicon) | `uv run pip install --force-reinstall vtk` |
| Corporate proxy | Set http_proxy/https_proxy in ~/.zshrc + npm config set proxy |

## Alternative: WSL/Linux Setup

If developing on WSL (Windows Subsystem for Linux) or native Linux instead of macOS, replace the Homebrew steps with apt-based installation.

### Step 1: System Dependencies (Ubuntu/WSL)

```bash
# FreeFEM++
sudo apt update && sudo apt install -y freefem++

# gmsh (optional, for complex geometry)
sudo apt install -y gmsh

# For PyVista off-screen rendering in WSL (headless)
sudo apt install -y libgl1-mesa-glx xvfb
```

### Step 2: Python Environment (uv)

Same as macOS — uv works identically on Linux:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
source $HOME/.local/bin/env

mkdir fem-workbench && cd fem-workbench
uv init
uv add pyvista meshio numpy vtk matplotlib
uv add --dev pytest ipython
```

### Step 3: Node.js

```bash
# Using fnm
curl -fsSL https://fnm.vercel.app/install | bash
source ~/.bashrc
fnm install 22
fnm use 22
```

### Step 4: Workflow

Same as macOS, but use `xvfb-run` if PyVista rendering fails:

```bash
# If off-screen rendering fails:
xvfb-run uv run python src/visualize.py
```

### WSL/Linux Troubleshooting

| Issue | Fix |
|---|---|
| FreeFem++ not found | `sudo apt install -y freefem++` then verify with `which FreeFem++` |
| PyVista blank/crash in WSL | `sudo apt install -y xvfb` then `xvfb-run python ...` |
| VTK import error | `uv run pip install --force-reinstall vtk` |
| No display (headless) | Always set `PYVISTA_OFF_SCREEN=true` before importing pyvista |
