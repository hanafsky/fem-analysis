# FEM 解析スキル

[Claude Code カスタムスキル](https://docs.anthropic.com/en/docs/claude-code/skills)。FreeFEM++ と PyVista を使った有限要素法（FEM）解析のエンドツーエンド自動化ツールです。

## 機能

自然言語の問題記述から、以下を自動化します：

1. **モデル＆メッシュ** — FreeFEM++ の `.edp` コードを生成
2. **求解** — FreeFEM++ を実行し結果を取得
3. **可視化** — PyVista で論文品質の画像をレンダリング
4. **評価** — AI が可視化結果をレビューし、必要に応じて修正を繰り返す

## 対応解析タイプ

- 2D 弾性解析（平面応力/平面ひずみ）
- 熱伝導解析
- ポアソン/ラプラス方程式
- ストークス流れ
- トポロジー最適化（SIMP, Level-Set）
- 形状最適化（畔上H1勾配法, パラメトリック）

## ディレクトリ構成

```
fem-analysis/
├── skills/
│   └── fem-analysis/
│       ├── SKILL.md              # スキル定義（Claude Code）
│       ├── references/           # FEM 理論・コードリファレンス
│       │   ├── elasticity_2d.md
│       │   ├── thermal_2d.md
│       │   ├── poisson_2d.md
│       │   ├── stokes_2d.md
│       │   ├── topology_optimization.md
│       │   ├── shape_optimization.md
│       │   ├── h1_gradient_method.md
│       │   ├── cmap_examples.md
│       │   └── project_setup.md
│       ├── scripts/              # Python ユーティリティ
│       │   ├── ff_mesh_reader.py # FreeFem++ メッシュパーサー
│       │   └── visualize.py      # PyVista レンダリング
│       └── templates/
│           └── FEMViewer.tsx     # React 3D ビューアコンポーネント
├── example/                      # 使用例
├── README.md                     # English
└── README_JA.md                  # 日本語（このファイル）
```

## 対応プラットフォーム

- **macOS**（Apple Silicon / Intel）— 公式 .dmg インストーラーで FreeFEM++ をインストール
- **Linux / WSL** — apt で FreeFEM++ をインストール

詳細なセットアップ手順は `skills/fem-analysis/references/project_setup.md` を参照してください。

## インストール

### 1. スキルの配置

Claude Code は以下の 2 箇所からスキルを自動検出します。用途に応じて選んでください：

**プロジェクト単位**（そのプロジェクトでのみ有効）：
```bash
# プロジェクトルートで実行
mkdir -p .claude/skills
cp -r /path/to/fem-analysis/skills/fem-analysis .claude/skills/
```

**個人 / グローバル**（すべてのプロジェクトで有効）：
```bash
mkdir -p ~/.claude/skills
cp -r /path/to/fem-analysis/skills/fem-analysis ~/.claude/skills/
```

コピー後、以下の構成になっていることを確認：
```
.claude/skills/fem-analysis/
├── SKILL.md
├── references/
├── scripts/
└── templates/
```

### 2. 前提ツールのインストール

スキルは FEM コードの生成・実行を行うため、以下のツールがシステムに必要です。

**macOS：**
```bash
# FreeFEM++ — 公式リリースから .dmg をダウンロード：
#   https://github.com/FreeFem/FreeFem-sources/releases
# インストール：
sudo cp -rf /Volumes/<マウントされたdmg>/FreeFem++.app /Applications/
sudo xattr -rc /Applications/FreeFem++.app
# PATH に追加（バージョン番号は適宜変更）：
#   bash/zsh: echo 'export PATH="/Applications/FreeFem++.app/Contents/ff-4.15.1/bin:$PATH"' >> ~/.zprofile
#   fish:     fish_add_path /Applications/FreeFem++.app/Contents/ff-4.15.1/bin

# gmsh（任意、複雑な形状向け）
brew install gmsh

# Python 環境（uv 経由）
curl -LsSf https://astral.sh/uv/install.sh | sh
mkdir fem-workbench && cd fem-workbench
uv init
uv add pyvista meshio numpy vtk matplotlib scipy
```

**Linux / WSL（apt）：**
```bash
sudo apt update && sudo apt install -y freefem++ gmsh
# apt版が古い/存在しない場合は公式リリースからダウンロード：
#   https://github.com/FreeFem/FreeFem-sources/releases
sudo apt install -y libgl1-mesa-glx xvfb   # ヘッドレス環境でのレンダリング用

# Python 環境（uv 経由）
curl -LsSf https://astral.sh/uv/install.sh | sh
mkdir fem-workbench && cd fem-workbench
uv init
uv add pyvista meshio numpy vtk matplotlib scipy
```

### 3. 動作確認

```bash
FreeFem++ --version
python -c "import pyvista; print(pyvista.__version__)"
```

プロジェクトの完全なセットアップ（React ビューア、vite 設定など）は [`skills/fem-analysis/references/project_setup.md`](skills/fem-analysis/references/project_setup.md) を参照してください。

## ライセンス

MIT
