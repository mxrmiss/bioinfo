# tool_readme.txt
# Swiss-Prot Annotation Pipeline - Tool Requirements & Checks

本流程由 Snakemake 调度（可选），核心计算使用 DIAMOND。
下面列出“必须/推荐/可选”工具、用途，以及对应的存在性验证命令（可直接复制执行）。

======================================================================
1) 必须工具（Required）
======================================================================

[1.1] Python 3
- 用途：运行 scripts/*.py
- 建议：Python >= 3.9（你脚本使用了现代 typing/Pathlib）
- 验证命令：
python3 --version
python3 -c "import sys; print(sys.version)"

[1.2] PyYAML（Python 包）
- 用途：读取 config.yaml
- 验证命令：
python3 -c "import yaml; print('PyYAML OK', yaml.__version__)"

[1.3] DIAMOND
- 用途：建库（diamond makedb）与比对（diamond blastp）
- 验证命令：
diamond version
diamond --help | head -n 5


======================================================================
2) 推荐工具（Recommended）
======================================================================

[2.1] aria2c
- 用途：高速下载 Swiss-Prot（支持 metalink / 分块 / 断点续传）
- 注：若没有 aria2c，01_db.py 会回退到 wget；再无则用 Python 下载（更慢）
- 验证命令：
aria2c -v

[2.2] wget
- 用途：aria2c 不存在时的下载回退方案
- 验证命令：
wget --version | head -n 2


======================================================================
3) 可选工具（Optional，但强烈建议用于全流程）
======================================================================

[3.1] Snakemake
- 用途：一键按依赖关系跑完整流程（db→align→filter→annot→reports）
- 注：不用 snakemake 也能按顺序直接运行各脚本
- 验证命令：
snakemake --version
snakemake --help | head -n 5


======================================================================
4) 环境快速自检（建议复制执行）
======================================================================

# 4.1 检查 python + pyyaml
python3 -c "import sys,yaml; print('Python OK', sys.version.split()[0]); print('PyYAML OK', yaml.__version__)"

# 4.2 检查 diamond
diamond version

# 4.3 检查下载工具（至少有其一：aria2c 或 wget）
command -v aria2c >/dev/null 2>&1 && aria2c -v | head -n 1
command -v wget  >/dev/null 2>&1 && wget --version | head -n 1

# 4.4 检查 snakemake（可选）
command -v snakemake >/dev/null 2>&1 && snakemake --version


======================================================================
5) 常见问题提示（极简）
======================================================================

- 运行 scripts/01_db.py 报 “import yaml failed”
  说明缺 PyYAML：请先安装 pyyaml（在你自己的环境里装即可）

- 运行比对时报 diamond 不存在
  说明 DIAMOND 未安装或未在 PATH 中：先确保 diamond version 能运行

- 下载非常慢
  建议安装 aria2c，并在 config.yaml 中保持 db.use_metalink: true

