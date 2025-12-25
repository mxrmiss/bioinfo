#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
10_mcmctree_run_and_publish.py —— 两阶段调度（verbatim），配置优先，缺配也能独立跑

特性：
  1) 优先读取 config.yaml:mcmctree（仅用于定位路径与日志回显，不改模板内容）。
  2) 无配置文件或无该小节时，使用脚本顶部默认参数也能独立运行。
  3) 进入 work_dir 后，按模板原样直跑：Stage1→复制 out.BV→in.BV→Stage2。
  4) 实时把子进程输出打印到屏幕并写入日志；失败时打包最小排查包。
  5) Stage3 从 FigTree.tre 提取超时钟树，对 Newick 进行 CAFE 友好清洗后，
     写入 config.inputs.ultrametric_tree 指定路径。
  6) 若 work_dir 内已经存在 FigTree.tre，则跳过 Stage1/Stage2，仅执行 Stage3。
  7) 所有日志文件写入 config.paths.logs_dir 指定路径下的 10_mcmctree 子目录，
     若配置缺失，则默认使用 ./logs/10_mcmctree。

依赖：
  - Python 标准库 + PyYAML（pip/conda 安装 pyyaml 即可）
  - 外部程序：mcmctree（PAML 4.10.9）
使用：
  - 直接运行，无命令行参数；需要时在脚本顶部修改默认参数。
"""

from __future__ import annotations
import os, sys, shutil, tarfile, time, datetime, subprocess, re
from pathlib import Path
from typing import Optional

# ===================== 顶部默认参数（无配置文件时生效） =====================
# 默认工作目录（包含 mcmctree.ctl、mcmctree2.ctl、concat.clean.phy、species_calib.nobl.trees）
DEFAULT_WORK_DIR: Path = Path("results/06_cafe/mcmctree")

# 模板文件名（位于 work_dir 内；脚本不会改动）
STAGE1_CTL_NAME: str = "mcmctree.ctl"
STAGE2_CTL_NAME: str = "mcmctree2.ctl"

# 阶段1生成/阶段2需要的 BV 文件名（固定）
OUT_BV_NAME: str = "out.BV"
IN_BV_NAME: str = "in.BV"

# 当 in.BV 已存在时是否覆盖
OVERWRITE_INBV: bool = True

# 是否写 .done 哨兵
WRITE_SENTINELS: bool = True

# 配置文件路径（会优先使用此文件；若不存在则走默认参数）
CONFIG_PATH: Path = Path("config.yaml")

# ===================== 下面内容一般不需要改动 =====================

def ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def ts_compact() -> str:
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def which_exec(name: str) -> Optional[str]:
    return shutil.which(name)

def run_streaming(cmd: list[str], cwd: Path, log_file: Path, prefix: str) -> int:
    """实时流式运行外部命令：同时写屏幕与日志。"""
    with open(log_file, "a", encoding="utf-8") as lf:
        lf.write(f"[{ts()}] [CMD] (cwd={cwd}) $ {' '.join(cmd)}\n")
        lf.flush()

        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            line = line.rstrip("\n")
            print(f"[{prefix}] {line}")
            lf.write(f"[{ts()}] [{prefix}] {line}\n")
        proc.wait()
        lf.write(f"[{ts()}] [EXIT] code={proc.returncode}\n")
        lf.flush()
        return int(proc.returncode)

def summarize_file(p: Path) -> str:
    if not p.exists():
        return f"{p} [MISSING]"
    size = p.stat().st_size
    mtime = datetime.datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S")
    return f"{p} [size={size} bytes, mtime={mtime}]"

def copy_outBV_to_inBV(work_dir: Path, overwrite: bool, logf: Path) -> None:
    """将 Stage1 产生的 out.BV 复制为 in.BV，供 Stage2 使用。"""
    src = work_dir / OUT_BV_NAME
    dst = work_dir / IN_BV_NAME
    with open(logf, "a", encoding="utf-8") as lf:
        lf.write(f"[{ts()}] [INFO] Preparing BV files: {summarize_file(src)}\n")
    if not src.exists():
        raise FileNotFoundError(f"未发现 {src}，Stage1 可能失败，无法继续 Stage2。")
    if dst.exists() and not overwrite:
        with open(logf, "a", encoding="utf-8") as lf:
            lf.write(f"[{ts()}] [INFO] {dst} 已存在，按策略不覆盖。\n")
        print(f"[middle] {dst} 已存在（保留原样）。")
        return
    shutil.copy2(src, dst)
    with open(logf, "a", encoding="utf-8") as lf:
        lf.write(f"[{ts()}] [INFO] Copied {src.name} -> {dst.name}\n")
    print(f"[middle] 已复制 {src.name} -> {dst.name}")

def clean_newick_for_cafe(newick: str) -> str:
    """
    将 Newick 字符串清洗为 CAFE 友好格式：
      - 去除括号、逗号、冒号周围多余空格；
      - 冒号后紧跟分支长度，不保留 ' : 0.123' 这种形式；
      - 压缩连续空白为单个空格，并最终尽可能减少无意义空格。
    """
    s = newick.strip()
    # 把各种换行/制表等白字符统一成空格，方便后续处理
    s = re.sub(r"\s+", " ", s)

    # 去掉括号内外多余空格
    s = re.sub(r"\(\s+", "(", s)
    s = re.sub(r"\s+\)", ")", s)

    # 逗号左右不留空格
    s = re.sub(r"\s*,\s*", ",", s)

    # 冒号前不保留空格
    s = re.sub(r"\s+:", ":", s)

    # 冒号后紧跟数字（分支长度）
    s = re.sub(r":\s*([0-9.eE+-]+)", r":\1", s)

    s = s.strip()
    if not s.endswith(";"):
        s = s + ";"
    return s

def extract_ultrametric_newick_from_figtree(figtree_path: Path) -> str:
    """从 FigTree.tre 中提取一棵去除注释并清洗格式的超时钟树（Newick 字符串）。"""
    if not figtree_path.exists():
        raise FileNotFoundError(f"未发现 FigTree 文件：{figtree_path}")
    txt = figtree_path.read_text(encoding="utf-8")

    # 匹配第一条 tree 语句
    m = re.search(r"tree\s+\S+\s*=\s*(.+?);", txt, flags=re.IGNORECASE | re.DOTALL)
    if not m:
        raise ValueError("在 FigTree.tre 中未发现 tree 语句")
    s = m.group(1)

    # 去掉 [&...] 注释块（例如 &95%={L,U} 等）
    s = re.sub(r"\[&[^\]]*\]", "", s)

    # CAFE 友好格式清洗
    s = clean_newick_for_cafe(s)
    return s

def publish_ultrametric_tree(work_dir: Path, cfg: Optional[dict], log_path: Path) -> Path:
    """
    Stage3：发布供 CAFE 使用的超时钟树。
      1) 在 work_dir 下读取 FigTree.tre。
      2) 从中提取去注释、已清洗格式的 ultrametric 树。
      3) 写入 config.inputs.ultrametric_tree 指定路径。
    """
    figtree = work_dir / "FigTree.tre"
    with open(log_path, "a", encoding="utf-8") as lf:
        lf.write(f"[{ts()}] [INFO] Stage3: 准备从 {figtree} 发布超时钟树\n")

    if cfg is None or not isinstance(cfg, dict):
        raise RuntimeError(
            "Stage3 需要 config.yaml 中的 inputs.ultrametric_tree 设置，"
            "但当前未能正确读取配置。"
        )

    inputs = cfg.get("inputs") if isinstance(cfg, dict) else None
    if not isinstance(inputs, dict) or not inputs.get("ultrametric_tree"):
        raise RuntimeError("config.yaml: inputs.ultrametric_tree 未设置，无法发布超时钟树。")

    out_rel = str(inputs["ultrametric_tree"])
    out_path = Path(out_rel)

    ensure_dir(out_path.parent)
    newick = extract_ultrametric_newick_from_figtree(figtree)
    # 写出为单行 Newick，末尾补换行符
    out_path.write_text(newick + "\n", encoding="utf-8")

    with open(log_path, "a", encoding="utf-8") as lf:
        lf.write(f"[{ts()}] [INFO] Stage3: ultrametric tree written to {out_path}\n")

    print(f"[Stage3] 已从 {figtree} 发布超时钟树到 {out_path}")
    return out_path

def pack_debug_bundle(work_dir: Path, debug_dir: Path, tb_prefix: str, extra_files: list[Path], log_path: Path) -> Path:
    """打包最小排查包，便于后续排错。"""
    ensure_dir(debug_dir)

    ts_str = ts_compact()
    tar_path = debug_dir / f"{tb_prefix}_{ts_str}.tar.gz"

    # 默认加入的文件：ctl、种子、BV、主日志
    candidates = [
        work_dir / STAGE1_CTL_NAME,
        work_dir / STAGE2_CTL_NAME,
        work_dir / "SeedUsed",
        work_dir / OUT_BV_NAME,
        work_dir / IN_BV_NAME,
        work_dir / "out.txt",
        work_dir / "mcmc.txt",
        log_path,
    ]
    for x in extra_files:
        candidates.append(x)

    files_to_add = [p for p in candidates if p.exists()]

    with tarfile.open(tar_path, "w:gz") as tf:
        for p in files_to_add:
            arcname = p.relative_to(work_dir.parent) if work_dir in p.parents else p.name
            tf.add(p, arcname=str(arcname))

    return tar_path

def load_config(path: Path) -> Optional[dict]:
    """加载 YAML 配置；若文件不存在或解析失败则返回 None。"""
    if not path.exists():
        return None
    try:
        import yaml  # 延迟导入，以免环境无此库时影响其它脚本

        with open(path, "r", encoding="utf-8") as f:
            cfg = yaml.safe_load(f)
        if not isinstance(cfg, dict):
            print(f"[WARN] 配置文件 {path} 解析结果不是字典，将忽略。")
            return None
        return cfg
    except Exception as e:
        print(f"[WARN] 读取配置 {path} 失败：{e}；将使用默认参数。")
        return None

def resolve_workdir_from_config(cfg: dict) -> Optional[Path]:
    """
    从 config.yaml 中解析 mcmctree.work_dir；若未设置则返回 None。

    预期配置结构（示例）：
      mcmctree:
        work_dir: "results/06_cafe/mcmctree"
    """
    mc = cfg.get("mcmctree")
    if not isinstance(mc, dict):
        return None
    wd = mc.get("work_dir")
    if not wd:
        return None
    return Path(str(wd))

def write_sentinel(work_dir: Path, name: str) -> None:
    """在 work_dir 下写一个 .<name>.done 文件作为完成标记。"""
    sent = work_dir / f".{name}.done"
    sent.write_text(f"{ts()}  done by 10_mcmctree_run_and_publish.py\n", encoding="utf-8")

def main() -> int:
    # 0) 读取配置（若存在，以配置为主；否则使用默认）
    cfg = load_config(CONFIG_PATH)
    if cfg is None:
        work_dir = DEFAULT_WORK_DIR
        print(f"[init] 未发现或未能解析 config.yaml，使用默认工作目录：{work_dir}")
        cfg_seq = "concat.clean.phy"     # 仅用于日志回显
        cfg_tree = "species_calib.nobl.trees"
        paths_cfg = {}
        # 默认模板目录（与当前项目根目录下的 templates/ 对应）
        templates_dir = Path("templates")
    else:
        wd = resolve_workdir_from_config(cfg)
        work_dir = wd if wd is not None else DEFAULT_WORK_DIR
        mc = cfg.get("mcmctree", {}) if isinstance(cfg, dict) else {}
        cfg_seq = mc.get("seqfile", "concat.clean.phy")
        cfg_tree = mc.get("treefile", "species_calib.nobl.trees")
        # 从 config.yaml:mcmctree.templates 读取模板目录，默认仍为 "templates"
        templates_dir = Path(mc.get("templates", "templates"))
        paths_cfg = cfg.get("paths", {}) if isinstance(cfg, dict) else {}
        print(f"[init] 使用配置文件：{CONFIG_PATH}")
        print(f"[init] mcmctree.work_dir = {work_dir}")
        print(f"[init] mcmctree.seqfile  = {cfg_seq}")
        print(f"[init] mcmctree.treefile = {cfg_tree}")

    # 1) 日志/调试目录
    #    日志根目录优先使用 config.paths.logs_dir，否则默认为 ./logs
    logs_root = Path(paths_cfg.get("logs_dir", "logs"))
    LOG_DIR = logs_root / "10_mcmctree"
    DEBUG_DIR = work_dir / "_debug"
    DEBUG_TARBALL_PREFIX = "mcmctree_debug"

    ensure_dir(work_dir)
    ensure_dir(LOG_DIR)
    ensure_dir(DEBUG_DIR)
    log_path = LOG_DIR / f"run_{ts_compact()}.log"

    with open(log_path, "a", encoding="utf-8") as lf:
        lf.write(f"[{ts()}] [INIT] work_dir={work_dir}\n")
        lf.write(f"[{ts()}] [INIT] logs_root={logs_root}\n")

    figtree_path = work_dir / "FigTree.tre"

    # ========= 快速路径：已有 FigTree.tre，仅发布超时树 =========
    if figtree_path.exists():
        print(f"[init] 检测到已有 FigTree.tre：{figtree_path}")
        print("[init] 跳过 Stage1/Stage2，仅执行 Stage3 发布超时树。")
        with open(log_path, "a", encoding="utf-8") as lf:
            lf.write(f"[{ts()}] [INIT] FigTree.tre 已存在，跳过 Stage1/Stage2，仅执行 Stage3。\n")

        try:
            publish_ultrametric_tree(work_dir, cfg, log_path)
        except Exception as e:
            print(f"[ERR] Stage3 发布超时钟树失败：{e}")
            with open(log_path, "a", encoding="utf-8") as lf:
                lf.write(f"[{ts()}] [ERR] Stage3 failed (fast path): {e}\n")
            tb = pack_debug_bundle(work_dir, DEBUG_DIR, DEBUG_TARBALL_PREFIX, extra_files=[figtree_path], log_path=log_path)
            print(f"[ERR] 已生成排查包：{tb}")
            return 5

        # 成功收尾（快速路径也给出摘要）
        print("[OK] 快速路径：仅发布超时树完成。产物摘要：")
        for f in ["out.txt", "mcmc.txt", OUT_BV_NAME, IN_BV_NAME]:
            print("     - " + summarize_file(work_dir / f))
        if WRITE_SENTINELS:
            write_sentinel(work_dir, "mcmctree_stage1")
            write_sentinel(work_dir, "mcmctree_stage2")
        print(f"[OK] 日志文件：{log_path}")
        return 0

    # ========= 正常两阶段路径：Stage1 + Stage2 + Stage3 =========

    # 2) 查找 mcmctree 可执行文件
    mctree = which_exec("mcmctree")
    if not mctree:
        print("[ERR] 未在 PATH 中找到 mcmctree 可执行文件。")
        with open(log_path, "a", encoding="utf-8") as lf:
            lf.write(f"[{ts()}] [ERR] mcmctree not found in PATH\n")
        return 2
    print(f"[init] 使用 mcmctree 可执行文件：{mctree}")

    # 2.5) 若工作目录中缺少 ctl 模板，则从 templates_dir 自动拷贝
    for name in (STAGE1_CTL_NAME, STAGE2_CTL_NAME):
        src = templates_dir / name
        dst = work_dir / name
        if not dst.exists():
            if src.exists():
                shutil.copy2(src, dst)
                print(f"[init] 已从模板目录复制 {src} -> {dst}")
            else:
                print(f"[WARN] 期望的模板文件不存在：{src}")

    # 3) 基本输入文件检查：模板 & 关键输入
    required = [
        work_dir / STAGE1_CTL_NAME,
        work_dir / STAGE2_CTL_NAME,
        work_dir / "concat.clean.phy",
        work_dir / "species_calib.nobl.trees",
    ]
    missing = [p for p in required if not p.exists()]
    if missing:
        print("[ERR] 下列必需文件缺失：")
        for p in missing:
            print(f"      - {p}")
        with open(log_path, "a", encoding="utf-8") as lf:
            lf.write(f"[{ts()}] [ERR] Missing required files:\n")
            for p in missing:
                lf.write(f"  - {p}\n")
        tb = pack_debug_bundle(work_dir, DEBUG_DIR, DEBUG_TARBALL_PREFIX, extra_files=[], log_path=log_path)
        print(f"[ERR] 已生成排查包：{tb}")
        return 3

    # 4) Stage1
    print(f"[Stage1] 开始：{STAGE1_CTL_NAME}")
    t0 = time.time()
    code1 = run_streaming([mctree, STAGE1_CTL_NAME], work_dir, log_path, prefix="S1")
    t1 = time.time()
    print(f"[Stage1] 结束，退出码={code1}，耗时={t1 - t0:.1f}s")
    if code1 != 0:
        tb = pack_debug_bundle(work_dir, DEBUG_DIR, DEBUG_TARBALL_PREFIX, extra_files=[], log_path=log_path)
        print(f"[ERR] Stage1 失败，已生成排查包：{tb}")
        return code1

    # 5) out.BV -> in.BV
    try:
        copy_outBV_to_inBV(work_dir, OVERWRITE_INBV, logf=log_path)
    except Exception as e:
        print(f"[ERR] 中间步骤失败：{e}")
        tb = pack_debug_bundle(work_dir, DEBUG_DIR, DEBUG_TARBALL_PREFIX, extra_files=[], log_path=log_path)
        print(f"[ERR] 已生成排查包：{tb}")
        return 4

    # 6) Stage2
    print(f"[Stage2] 开始：{STAGE2_CTL_NAME}")
    t2 = time.time()
    code2 = run_streaming([mctree, STAGE2_CTL_NAME], work_dir, log_path, prefix="S2")
    t3 = time.time()
    print(f"[Stage2] 结束，退出码={code2}，耗时={t3 - t2:.1f}s")
    if code2 != 0:
        tb = pack_debug_bundle(work_dir, DEBUG_DIR, DEBUG_TARBALL_PREFIX, extra_files=[], log_path=log_path)
        print(f"[ERR] Stage2 失败，已生成排查包：{tb}")
        return code2

    # Stage3: 从 FigTree.tre 发布 ultrametric 超时钟树（供 CAFE 使用）
    try:
        publish_ultrametric_tree(work_dir, cfg, log_path)
    except Exception as e:
        print(f"[ERR] Stage3 发布超时钟树失败：{e}")
        with open(log_path, "a", encoding="utf-8") as lf:
            lf.write(f"[{ts()}] [ERR] Stage3 failed: {e}\n")
        tb = pack_debug_bundle(work_dir, DEBUG_DIR, DEBUG_TARBALL_PREFIX, extra_files=[], log_path=log_path)
        print(f"[ERR] 已生成排查包：{tb}")
        return 5

    # 7) 成功收尾
    print("[OK] 两阶段运行完成。产物摘要：")
    for f in ["out.txt", "mcmc.txt", OUT_BV_NAME, IN_BV_NAME]:
        print("     - " + summarize_file(work_dir / f))
    if WRITE_SENTINELS:
        write_sentinel(work_dir, "mcmctree_stage1")
        write_sentinel(work_dir, "mcmctree_stage2")
    print(f"[OK] 日志文件：{log_path}")
    return 0

if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n[INTERRUPTED] 用户中断。")
        sys.exit(130)

