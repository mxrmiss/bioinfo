#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
12_cafe5_run_models.py —— CAFE5 教程增强版（两阶段 / 极端家族修剪 / 失败率统计）

功能总览（与 config.yaml.cafe5 对应）：
  1) 读取 config.yaml 中的 cafe5 & paths 配置。
  2) 从 results/06_cafe/input/ 读取：
       - family.tsv            （第一列为 Desc = OG ID）
       - utree_for_cafe.nwk    （超时树）
  3) 可选“教程增强”：
       - two_stage_large.enable:
           按 copy_threshold 拆分为 primary 与 large 两个 family 表。
       - autofix_rounds:
           每个集合做若干轮“极端家族自动删除”：
             * 解析 CAFE5 日志中的
                 “Families with largest size differentials:”
               区块，拿到 OG 列表；
             * 从 family 表中删掉这些 OG，生成 family.autofix*.tsv；
             * 直到 CAFE5 不再提示 “You may want to try removing...”
               或达到 max_rounds。
       - 解析高失败率家族：
             “The following families had failure rates >20% of the time:”
         输出 flags/high_fail_ogs.list。
  4) 针对每个模型（默认只有 "global"）创建目录：
       results/06_cafe/models/<model>/
         ├─ family.primary.tsv
         ├─ family.large.tsv
         ├─ primary_global/
         ├─ large/                  （若 large 非空才跑）
         ├─ flags/high_fail_ogs.list
         ├─ sentinels/
         └─ run.log                 （主流程日志）
  5) 全部模型成功后，在 cafe_run_dir 下写入 .cafe.done 哨兵。

注意：
  - 不接受命令行参数，所有参数集中在 config.yaml 顶部。
  - 本脚本只修复路径与 CAFE5 调用逻辑，不对 error_model/multi_lambda 做复杂操作，
    但会安全地忽略这些配置（即它们存在也不会报错）。
"""

from __future__ import annotations

import os
import sys
import io
import re
import shutil
import subprocess
import logging
import datetime
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional

import yaml


# ======================== 通用工具 ========================

DEFAULT_CONFIG = "config.yaml"


def _expand_publish_placeholders(obj, publish_dir: str):
    """把 <publish_dir> 占位符替换成真实路径。"""
    if isinstance(obj, str):
        return obj.replace("<publish_dir>", publish_dir)
    if isinstance(obj, list):
        return [_expand_publish_placeholders(x, publish_dir) for x in obj]
    if isinstance(obj, dict):
        return {k: _expand_publish_placeholders(v, publish_dir) for k, v in obj.items()}
    return obj


def load_config(config_path: str = DEFAULT_CONFIG) -> Dict[str, Any]:
    """读取 YAML 配置，并展开 publish_dir 占位符。"""
    p = Path(config_path)
    if not p.exists():
        raise FileNotFoundError(f"[ERR] 未找到配置文件：{p}")
    cfg = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    pub = cfg.get("publish_dir")
    if pub:
        cfg["inputs"] = _expand_publish_placeholders(cfg.get("inputs", {}), str(pub))
    return cfg


def ensure_dir(p: Path) -> Path:
    """确保目录存在。"""
    p.mkdir(parents=True, exist_ok=True)
    return p


def rm_tree(p: Path):
    """安全删除整个目录树。"""
    if p.is_dir():
        shutil.rmtree(p)


def get_logger(name: str, logfile: Path, level: int = logging.INFO) -> logging.Logger:
    """同时打日志到文件 + 终端。"""
    ensure_dir(logfile.parent)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.handlers.clear()

    fmt = logging.Formatter(
        "[%(asctime)s] %(levelname)s - %(message)s",
        "%Y-%m-%d %H:%M:%S",
    )

    fh = logging.FileHandler(logfile, encoding="utf-8")
    fh.setFormatter(fmt)
    fh.setLevel(level)
    logger.addHandler(fh)

    sh = logging.StreamHandler(stream=sys.stdout)
    sh.setFormatter(fmt)
    sh.setLevel(level)
    logger.addHandler(sh)

    # 让 print 也立即 flush
    class _Flush(io.TextIOBase):
        def __init__(self, s):
            self.s = s

        def write(self, x):
            self.s.write(x)
            self.s.flush()
            return len(x)

    sys.stdout = _Flush(sys.stdout)
    sys.stderr = _Flush(sys.stderr)
    return logger


def banner(log: logging.Logger, text: str):
    bar = "=" * max(10, len(text) + 2)
    log.info(bar)
    log.info(f" {text} ")
    log.info(bar)


# ======================== family.tsv 处理 ========================

def split_family_primary_large(
    family_tsv: Path,
    out_primary: Path,
    out_large: Path,
    copy_threshold: int,
    log: logging.Logger,
) -> Tuple[int, int]:
    """
    按“单物种最大拷贝数阈值”拆分 family.tsv 为 primary 与 large 两个集合。

    约定：
      - 第一列为 Desc = OG ID（字符串），从第二列开始都是各物种拷贝数。
      - copy_threshold: 若某行任一物种基因数 >= 阈值，则划入 large 集合。
    """
    txt = family_tsv.read_text(encoding="utf-8").rstrip("\n").splitlines()
    if not txt:
        raise RuntimeError(f"[ERR] family.tsv 为空：{family_tsv}")

    header = txt[0]
    primary_lines = [header]
    large_lines = [header]

    for ln in txt[1:]:
        if not ln.strip():
            continue
        cols = ln.rstrip("\n").split("\t")
        if len(cols) <= 1:
            primary_lines.append(ln)
            continue
        counts = []
        for c in cols[1:]:
            c = c.strip()
            if c in ("", ".", "NA", "NaN", "-", "?"):
                continue
            try:
                counts.append(int(c))
            except ValueError:
                # 非数字一律忽略，不参与 max
                continue
        max_cnt = max(counts) if counts else 0
        if max_cnt >= copy_threshold:
            large_lines.append(ln)
        else:
            primary_lines.append(ln)

    out_primary.write_text("\n".join(primary_lines) + "\n", encoding="utf-8")
    out_large.write_text("\n".join(large_lines) + "\n", encoding="utf-8")

    n_primary = max(0, len(primary_lines) - 1)
    n_large = max(0, len(large_lines) - 1)
    log.info(
        f"[SPLIT] 阈值 {copy_threshold}: primary={n_primary} 行, large={n_large} 行"
    )
    return n_primary, n_large


# ======================== CAFE5 调用与解析 ========================

EXTREME_BLOCK_RE = re.compile(
    r"Families with largest size differentials:\R"
    r"((?:.*?OG\d+:\s*\d+\R?)+)",
    flags=re.M,
)

OG_KV_RE = re.compile(r"(OG\d+)\s*:\s*(\d+)")

FAIL_BLOCK_RE = re.compile(
    r"The following families had failure rates >20% of the time:\R"
    r"((?:.*?OG\d+\s+had\s+\d+\s+failures\R?)+)",
    flags=re.M,
)

FAIL_LINE_RE = re.compile(r"(OG\d+)\s+had\s+(\d+)\s+failures")


def run_cafe_once(
    cafe_bin: str,
    family_tsv: Path,
    tree_nwk: Path,
    threads: int,
    gamma_k: int,
    pvalue: float,
    cwd: Path,
    log: logging.Logger,
    label: str,
) -> str:
    """
    跑一次 CAFE5，返回 stdout 文本，用于后续解析极端家族 / 高失败率等。

    关键：family_tsv 与 tree_nwk 一律使用绝对路径，避免 cwd 更改导致找不到文件。
    """
    ensure_dir(cwd)
    family_abs = family_tsv.resolve()
    tree_abs = tree_nwk.resolve()

    cmd = [
        cafe_bin,
        "-i",
        str(family_abs),
        "-t",
        str(tree_abs),
        "-c",
        str(threads),
        "-k",
        str(gamma_k),
        "-P",
        str(pvalue),
    ]
    log.info(f"[CMD][{label}] {' '.join(cmd)} (cwd={cwd})")

    proc = subprocess.run(
        cmd,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    out = proc.stdout or ""
    for ln in out.splitlines():
        log.info(f"[{label}] {ln}")

    if proc.returncode != 0:
        raise RuntimeError(f"[ERR] CAFE5 运行失败（{label}），退出码 {proc.returncode}")

    return out


def parse_extreme_ogs(stdout_text: str) -> List[str]:
    """
    从 CAFE5 输出中解析“极端家族” OG 列表（取最后一个块）。
    """
    blocks = EXTREME_BLOCK_RE.findall(stdout_text)
    if not blocks:
        return []
    last_blk = blocks[-1]
    ogs: List[str] = []
    for m in OG_KV_RE.finditer(last_blk):
        ogs.append(m.group(1))
    return ogs


def drop_ogs_from_family(
    in_tsv: Path,
    out_tsv: Path,
    removed_list_tsv: Path,
    og_ids: List[str],
    log: logging.Logger,
) -> int:
    """
    从 family.tsv 中删除指定 OG（按第一列 Desc 精确匹配），并写出新表与删除列表。
    返回删除的行数。
    """
    if not og_ids:
        return 0
    og_set = set(og_ids)

    txt = in_tsv.read_text(encoding="utf-8").rstrip("\n").splitlines()
    if not txt:
        raise RuntimeError(f"[ERR] family 文件为空：{in_tsv}")

    header = txt[0]
    kept = [header]
    removed = [header]

    n_removed = 0
    for ln in txt[1:]:
        if not ln.strip():
            continue
        cols = ln.split("\t")
        desc = cols[0].strip() if cols else ""
        if desc in og_set:
            removed.append(ln)
            n_removed += 1
        else:
            kept.append(ln)

    out_tsv.write_text("\n".join(kept) + "\n", encoding="utf-8")
    removed_list_tsv.write_text("\n".join(removed) + "\n", encoding="utf-8")
    log.info(
        f"[FILTER] {in_tsv.name}: 解析到 {len(og_ids)} 个极端 OG，实际删除 {n_removed} 行，写入 {out_tsv.name}"
    )
    return n_removed


def parse_high_fail_ogs(stdout_text: str) -> List[str]:
    """
    解析“失败率 >20%”的 OG 列表（最后一个块）。
    """
    blocks = FAIL_BLOCK_RE.findall(stdout_text)
    if not blocks:
        return []
    last_blk = blocks[-1]
    ogs: List[str] = []
    for m in FAIL_LINE_RE.finditer(last_blk):
        ogs.append(m.group(1))
    return ogs


def ensure_gamma_results(out_dir: Path):
    """
    确认 Gamma_* 结果文件存在，否则视为 CAFE5 异常结束。
    """
    gamma = out_dir / "Gamma_results.txt"
    if not gamma.is_file():
        raise RuntimeError(f"[ERR] 未找到 Gamma_results.txt：{gamma}")


def run_with_autofix(
    cafe_bin: str,
    family_tsv: Path,
    tree_nwk: Path,
    threads: int,
    gamma_k: int,
    pvalue: float,
    cwd: Path,
    log: logging.Logger,
    label_prefix: str,
    max_autofix_rounds: int,
    model_dir: Path,
) -> Tuple[str, List[str]]:
    """
    针对一个 family 表执行多轮 CAFE5 + 自动极端家族剔除。

    返回：
      - 最后一轮 CAFE5 的 stdout
      - 高失败率 OG 列表（来自最后一轮）
    """
    current_tsv = family_tsv
    stdout_last = ""
    high_fail_ogs: List[str] = []

    for rnd in range(1, max_autofix_rounds + 1):
        label = f"{label_prefix} ROUND-{rnd}"
        log.info(f"[{label}] ---- CAFE5 输出开始 ----")
        stdout_last = run_cafe_once(
            cafe_bin=cafe_bin,
            family_tsv=current_tsv,
            tree_nwk=tree_nwk,
            threads=threads,
            gamma_k=gamma_k,
            pvalue=pvalue,
            cwd=cwd,
            log=log,
            label=label,
        )
        log.info(f"[{label}] ---- CAFE5 输出结束 ----")

        # 尝试解析极端家族
        extreme_ogs = parse_extreme_ogs(stdout_last)
        need_suggest = (
            "You may want to try removing the top few families" in stdout_last
            or "Failed to initialize any reasonable values" in stdout_last
        )

        if need_suggest and extreme_ogs and rnd < max_autofix_rounds:
            # 进行下一轮 autofix
            new_tsv = model_dir / f"family.autofix{rnd}.tsv"
            removed_tsv = model_dir / f"autofix_removed_round{rnd}.tsv"
            n_removed = drop_ogs_from_family(
                in_tsv=current_tsv,
                out_tsv=new_tsv,
                removed_list_tsv=removed_tsv,
                og_ids=extreme_ogs,
                log=log,
            )
            if n_removed == 0:
                # 解析到了 OG，但 family 表中一个都删不掉，说明列格式有问题
                raise RuntimeError(
                    "[ERR] 解析到极端 OG，但在 family.tsv 中未命中任何行，请检查首列是否为 OG ID/分隔符是否为 TAB。"
                )
            log.info(
                f"[CLEAN] 第 {rnd} 轮自动修正完成，进入下一轮。"
            )
            current_tsv = new_tsv
            continue
        else:
            if need_suggest and rnd == max_autofix_rounds:
                log.warning(
                    f"[WARN] 已达到 max_autofix_rounds={max_autofix_rounds}，"
                    f"仍然检测到 CAFE5 建议删除极端家族，请人工检查 family.tsv。"
                )
            break

    # 确认结果文件存在
    ensure_gamma_results(cwd)

    # 最后一轮解析高失败率 OG
    high_fail_ogs = parse_high_fail_ogs(stdout_last)
    if high_fail_ogs:
        log.info(f"[QC] 检测到高失败率家族 {len(high_fail_ogs)} 个。")

    return stdout_last, high_fail_ogs


# ======================== 主逻辑 ========================

def main():
    cfg = load_config()
    paths = cfg["paths"]
    cafe_cfg = cfg.get("cafe5", {}) or {}

    cafe_enable = bool(cafe_cfg.get("enable", True))
    threads = int(cafe_cfg.get("threads", 4))
    gamma_k = int(cafe_cfg.get("gamma_k", 3))
    pvalue = float(cafe_cfg.get("pvalue", 0.05))
    models = cafe_cfg.get("models", ["global"]) or ["global"]
    max_autofix_rounds = int(cafe_cfg.get("autofix_rounds", 1))

    two_stage_cfg = cafe_cfg.get("two_stage_large", {}) or {}
    two_stage_enable = bool(two_stage_cfg.get("enable", False))
    copy_threshold = int(two_stage_cfg.get("copy_threshold", 100))

    # 其它增强（先读出来，虽然本脚本暂时不深入使用）
    error_cfg = cafe_cfg.get("error_model", {}) or {}
    multi_cfg = cafe_cfg.get("multi_lambda", {}) or {}

    project_root = Path(".").resolve()
    cafe_run_dir = Path(paths["cafe_run_dir"]).resolve()
    logs_dir = Path(paths["logs_dir"]).resolve()

    # 输入路径
    input_dir = cafe_run_dir / "input"
    family_tsv = input_dir / "family.tsv"
    tree_nwk = input_dir / "utree_for_cafe.nwk"

    if not cafe_enable:
        # 若关闭，则写 .cafe.done 后退出
        ensure_dir(cafe_run_dir)
        (cafe_run_dir / ".cafe.done").write_text(
            "cafe5 disabled in config\n", encoding="utf-8"
        )
        return

    if not family_tsv.is_file():
        raise FileNotFoundError(f"[ERR] family.tsv 不存在：{family_tsv}")
    if not tree_nwk.is_file():
        raise FileNotFoundError(f"[ERR] utree_for_cafe.nwk 不存在：{tree_nwk}")

    # 日志
    logger = get_logger("aphylo.12", logs_dir / "12_cafe5_run_models.log")
    banner(logger, "APhylo 12 — CAFE5（教程增强版）")

    logger.info(f"[IN] family: {family_tsv}")
    logger.info(f"[IN] tree:   {tree_nwk}")
    logger.info(f"[CFG] threads = {threads}, gamma_k = {gamma_k}, pvalue = {pvalue}")
    logger.info(f"[CFG] models = {models}")
    logger.info(f"[CFG] max_autofix_rounds = {max_autofix_rounds}")
    logger.info(
        f"[CFG] two_stage_large.enable = {two_stage_enable}, "
        f"copy_threshold = {copy_threshold}"
    )

    cafe_bin = shutil.which("cafe5") or shutil.which("cafexp")
    if not cafe_bin:
        raise RuntimeError("[ERR] 未在 PATH 中找到 cafe5/cafexp 可执行文件")

    models_root = cafe_run_dir / "models"
    ensure_dir(models_root)

    # 每个模型单独跑
    for model in models:
        model_dir = models_root / model
        logger.info("")
        logger.info(f"================ [MODEL {model}] 开始 ================")

        # 若已有旧目录，重命名为备份，避免内容混在一起
        if model_dir.exists():
            ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            backup = model_dir.with_name(f"{model_dir.name}.bak_{ts}")
            logger.warning(f"[MODEL {model}] 检测到旧目录，重命名为 {backup}")
            model_dir.rename(backup)
        ensure_dir(model_dir)

        # 预先准备 primary / large family 表
        primary_tsv = model_dir / "family.primary.tsv"
        large_tsv = model_dir / "family.large.tsv"

        if two_stage_enable:
            n_primary, n_large = split_family_primary_large(
                family_tsv=family_tsv,
                out_primary=primary_tsv,
                out_large=large_tsv,
                copy_threshold=copy_threshold,
                log=logger,
            )
        else:
            # 不拆分：全部写入 primary，large 只保留表头
            txt = family_tsv.read_text(encoding="utf-8")
            primary_tsv.write_text(txt, encoding="utf-8")
            header = txt.splitlines()[0] if txt else ""
            large_tsv.write_text((header + "\n") if header else "", encoding="utf-8")
            logger.info("[SPLIT] two_stage_large.disable: 全部家族归入 primary。")

        # 清理 / 创建子目录
        primary_dir = ensure_dir(model_dir / "primary_global")
        large_dir = ensure_dir(model_dir / "large")
        flags_dir = ensure_dir(model_dir / "flags")
        sentinels_dir = ensure_dir(model_dir / "sentinels")

        # PRIMARY 阶段
        stdout_primary, high_fail_primary = run_with_autofix(
            cafe_bin=cafe_bin,
            family_tsv=primary_tsv,
            tree_nwk=tree_nwk,
            threads=threads,
            gamma_k=gamma_k,
            pvalue=pvalue,
            cwd=primary_dir,
            log=logger,
            label_prefix=f"MODEL={model} PRIMARY-GLOBAL",
            max_autofix_rounds=max_autofix_rounds,
            model_dir=model_dir,
        )

        # LARGE 阶段（只有在 two_stage 模式 + large 非空时才跑）
        high_fail_large: List[str] = []
        n_large_data = sum(1 for _ in large_tsv.read_text(encoding="utf-8").splitlines()) - 1
        if two_stage_enable and n_large_data > 0:
            logger.info(
                f"[MODEL {model}] 检测到 large 集合 {n_large_data} 行，启动 LARGE 阶段。"
            )
            stdout_large, high_fail_large = run_with_autofix(
                cafe_bin=cafe_bin,
                family_tsv=large_tsv,
                tree_nwk=tree_nwk,
                threads=threads,
                gamma_k=gamma_k,
                pvalue=pvalue,
                cwd=large_dir,
                log=logger,
                label_prefix=f"MODEL={model} LARGE-GLOBAL",
                max_autofix_rounds=max_autofix_rounds,
                model_dir=model_dir,
            )
        else:
            logger.info(
                f"[MODEL {model}] 未检测到 large 数据或 two_stage_large.disable，跳过 LARGE 阶段。"
            )

        # 写出高失败率 OG 汇总
        high_fail_all = sorted(set(high_fail_primary) | set(high_fail_large))
        hf_file = flags_dir / "high_fail_ogs.list"
        if high_fail_all:
            hf_file.write_text("\n".join(high_fail_all) + "\n", encoding="utf-8")
            logger.info(
                f"[MODEL {model}] 写出高失败率家族列表：{hf_file} （{len(high_fail_all)} 个 OG）"
            )
        else:
            hf_file.write_text("", encoding="utf-8")
            logger.info(f"[MODEL {model}] 未检测到高失败率家族，写入空文件：{hf_file}")

        # 模型级别 sentinel
        (sentinels_dir / f"{model}.done").write_text(
            f"{model} finished at {datetime.datetime.now()}\n",
            encoding="utf-8",
        )
        logger.info(f"[MODEL {model}] 完成。")

    # 全局 sentinel
    (cafe_run_dir / ".cafe.done").write_text(
        f"CAFE5 finished at {datetime.datetime.now()}\n",
        encoding="utf-8",
    )
    logger.info("")
    logger.info("[GLOBAL] 收工，全部模型完成，写出 .cafe.done。")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n")
        sys.exit(2)