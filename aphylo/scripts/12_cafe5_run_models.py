#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Run CAFE5 (primary + optional large) driven by config.yaml."""

from __future__ import annotations
import sys
import io
import re
import shutil
import subprocess
import logging
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional

import yaml

# ===================== 基础工具 =====================

CONFIG_PATH = "config.yaml"


def load_config(path: str = CONFIG_PATH) -> Dict[str, Any]:
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"[ERR] 未找到配置文件: {p}")
    cfg = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    return cfg


def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def get_logger(logfile: Path) -> logging.Logger:
    ensure_dir(logfile.parent)
    logger = logging.getLogger("aphylo.cafe5")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter(
        "[%(asctime)s] %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S"
    )

    fh = logging.FileHandler(logfile, encoding="utf-8")
    fh.setFormatter(fmt)
    fh.setLevel(logging.INFO)
    logger.addHandler(fh)

    sh = logging.StreamHandler(stream=sys.stdout)
    sh.setFormatter(fmt)
    sh.setLevel(logging.INFO)
    logger.addHandler(sh)

    # 让 print 也立即刷新，方便看进度
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


# ===================== family.tsv 处理 =====================

def split_primary_large(
    family_src: Path,
    model_dir: Path,
    threshold: int,
    primary_tag: str,
    large_tag: str,
    logger: logging.Logger,
) -> Tuple[Path, Optional[Path]]:
    """按最大拷贝数阈值拆分 primary / large 集。首列 Desc，第二列 Orthogroup。"""

    text = family_src.read_text(encoding="utf-8").rstrip("\n").splitlines()
    if not text:
        raise RuntimeError(f"[ERR] family 文件为空: {family_src}")

    header = text[0]
    primary_lines = [header]
    large_lines: List[str] = [header]

    for ln in text[1:]:
        if not ln.strip():
            continue
        parts = ln.split("\t")
        # 从第三列开始都是物种拷贝数
        counts = []
        for x in parts[2:]:
            x = x.strip()
            if x == "":
                continue
            # 粗暴转 int，不可转就算 0
            try:
                counts.append(int(x))
            except ValueError:
                counts.append(0)
        is_large = any(c >= threshold for c in counts)
        if is_large:
            large_lines.append(ln)
        else:
            primary_lines.append(ln)

    primary_file = model_dir / f"family.{primary_tag}.tsv"
    primary_file.write_text("\n".join(primary_lines) + "\n", encoding="utf-8")
    logger.info(f"[SPLIT] primary 行数: {len(primary_lines)-1}")

    large_file: Optional[Path] = None
    if len(large_lines) > 1:
        large_file = model_dir / f"family.{large_tag}.tsv"
        large_file.write_text("\n".join(large_lines) + "\n", encoding="utf-8")
        logger.info(f"[SPLIT] large 行数: {len(large_lines)-1}")
    else:
        logger.info("[SPLIT] 无 large 家族，后续 large 阶段跳过")

    return primary_file, large_file


def extract_extreme_ogs_from_runlog(runlog: Path) -> List[str]:
    """从 run.log 最后一个 'Families with largest size differentials:' 块中抽取 OG 列表。"""

    if not runlog.is_file():
        return []

    lines = runlog.read_text(encoding="utf-8", errors="ignore").splitlines()
    idxs = [
        i
        for i, ln in enumerate(lines)
        if "Families with largest size differentials:" in ln
    ]
    if not idxs:
        return []

    start = idxs[-1] + 1
    ogs: List[str] = []
    for ln in lines[start:]:
        t = ln.strip()
        if not t:
            break
        if "You may want to try removing" in t:
            break
        if "Failed to initialize any reasonable values" in t:
            break
        m = re.search(r"(OG\d+)", t)
        if m:
            ogs.append(m.group(1))

    return ogs


def filter_family_remove_ogs(
    family_in: Path, family_out: Path, ogs_to_remove: List[str]
) -> List[str]:
    """从 family.tsv 中按 Orthogroup 列剔除指定 OG，返回实际剔除的 OG 列表。"""

    rm_set = set(ogs_to_remove)
    if not rm_set:
        return []

    lines = family_in.read_text(encoding="utf-8").rstrip("\n").splitlines()
    if not lines:
        return []

    header = lines[0]
    kept = [header]
    removed: List[str] = []

    for ln in lines[1:]:
        if not ln.strip():
            continue
        parts = ln.split("\t")
        if len(parts) < 2:
            kept.append(ln)
            continue
        og = parts[1].strip()
        if og in rm_set:
            removed.append(og)
        else:
            kept.append(ln)

    family_out.write_text("\n".join(kept) + "\n", encoding="utf-8")
    return sorted(set(removed))


# ===================== 运行 CAFE5 =====================

def run_cafe5_once(
    cafe5_bin: str,
    family_file: Path,
    tree_file: Path,
    threads: int,
    gamma_k: int,
    pvalue: float,
    work_dir: Path,
    logger: logging.Logger,
    tag: str,
) -> Path:
    """在指定 work_dir 中跑一轮 CAFE5，输出 run.log 路径。"""

    ensure_dir(work_dir)
    runlog = work_dir / "run.log"

    cmd = [
        cafe5_bin,
        "-i",
        str(family_file.resolve()),
        "-t",
        str(tree_file.resolve()),
        "-c",
        str(threads),
    ]
    if gamma_k and gamma_k > 0:
        cmd += ["-k", str(gamma_k)]
    if pvalue is not None:
        cmd += ["-P", str(pvalue)]

    logger.info(f"[{tag}] CMD: {' '.join(cmd)}")

    with runlog.open("w", encoding="utf-8") as fh:
        proc = subprocess.Popen(
            cmd,
            cwd=str(work_dir),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            fh.write(line)
            fh.flush()
            # 原样流到屏幕
            print(line, end="")
        ret = proc.wait()

    if ret != 0:
        raise RuntimeError(f"[ERR] CAFE5 运行失败（{tag}），退出码 {ret}")

    logger.info(f"[{tag}] CAFE5 结束")
    return runlog


# ===================== 主流程 =====================

def main() -> None:
    cfg = load_config()

    cafe_cfg: Dict[str, Any] = cfg.get("cafe5") or {}
    paths_cfg: Dict[str, Any] = cfg.get("paths") or {}
    two_stage_cfg: Dict[str, Any] = cfg.get("two_stage_large") or {}

    models = cafe_cfg.get("models")
    if not isinstance(models, list) or not models:
        raise SystemExit("[ERR] config.cafe5.models 必须是非空列表，如 ['global']")

    threads = int(cafe_cfg.get("threads", 4))
    gamma_k = int(cafe_cfg.get("gamma_k", 3))
    pvalue = float(cafe_cfg.get("pvalue", 0.05))
    max_autofix_rounds = int(cafe_cfg.get("max_autofix_rounds", 0))

    two_stage_enable = bool(two_stage_cfg.get("enable", False))
    copy_threshold = int(two_stage_cfg.get("copy_threshold", 100))
    primary_tag = str(two_stage_cfg.get("primary_tag", "primary"))
    large_tag = str(two_stage_cfg.get("large_tag", "large"))

    inputs_cfg: Dict[str, Any] = cfg.get("inputs") or {}
    family_src = Path(inputs_cfg["family_tsv"])
    tree_file = Path(inputs_cfg["ultrametric_tree"])

    cafe_root = Path(paths_cfg["cafe_run_dir"])
    models_root = cafe_root / "models"
    logs_dir = Path(paths_cfg["logs_dir"])
    log_file = logs_dir / "12_cafe5_run_models.log"

    # 每次运行前清理旧 models 和旧日志
    if models_root.exists():
        shutil.rmtree(models_root)
    ensure_dir(models_root)
    if log_file.exists():
        log_file.unlink()

    logger = get_logger(log_file)
    logger.info("===================================")
    logger.info(" APhylo 12 — CAFE5 教程增强版（精简注释）")
    logger.info("===================================")
    logger.info(f"[IN] family: {family_src}")
    logger.info(f"[IN] tree:   {tree_file}")
    logger.info(f"[CFG] models: {models}")
    logger.info(f"[CFG] max_autofix_rounds: {max_autofix_rounds}")
    logger.info(f"[CFG] two_stage_large.enable: {two_stage_enable}")

    cafe5_bin = str(cfg.get("binaries", {}).get("cafe5", "cafe5"))

    if not family_src.is_file():
        raise SystemExit(f"[ERR] family_tsv 不存在: {family_src}")
    if not tree_file.is_file():
        raise SystemExit(f"[ERR] ultrametric_tree 不存在: {tree_file}")

    # 目前 models 通常只有 "global"，还是按循环写，方便之后扩展
    for model_name in models:
        model_dir = ensure_dir(models_root / model_name)
        logger.info(f"[MODEL {model_name}] 工作目录: {model_dir}")

        # 拷贝原始 family，留一份原始记录
        family_raw = model_dir / "family.raw.tsv"
        shutil.copy2(family_src, family_raw)

        # 拆 primary / large
        if two_stage_enable:
            primary_family, large_family = split_primary_large(
                family_raw,
                model_dir,
                threshold=copy_threshold,
                primary_tag=primary_tag,
                large_tag=large_tag,
                logger=logger,
            )
        else:
            primary_family = model_dir / f"family.{primary_tag}.tsv"
            shutil.copy2(family_raw, primary_family)
            large_family = None
            logger.info("[SPLIT] two_stage_large.disable: 全部家族走 primary")

        # ---------- primary: 自动极端修剪 + 运行 ----------
        current_family = primary_family
        primary_dir = ensure_dir(model_dir / f"{primary_tag}_global")

        for round_id in range(1, max_autofix_rounds + 1):
            tag = f"MODEL={model_name} PRIMARY-GLOBAL ROUND-{round_id}"
            logger.info(f"[{tag}] 开始 CAFE5 运行")
            runlog = run_cafe5_once(
                cafe5_bin=cafe5_bin,
                family_file=current_family,
                tree_file=tree_file,
                threads=threads,
                gamma_k=gamma_k,
                pvalue=pvalue,
                work_dir=primary_dir,
                logger=logger,
                tag=tag,
            )

            extreme_ogs = extract_extreme_ogs_from_runlog(runlog)
            if not extreme_ogs:
                logger.info(f"[{tag}] 未检测到极端 OG，终止自动修正循环")
                break

            next_family = model_dir / f"family.autofix{round_id+1}.tsv"
            removed = filter_family_remove_ogs(
                family_in=current_family,
                family_out=next_family,
                ogs_to_remove=extreme_ogs,
            )
            if not removed:
                logger.info(
                    f"[{tag}] 建议剔除的 OG 在 family 中未命中，终止自动修正循环"
                )
                break

            removed_file = model_dir / f"autofix_removed_round{round_id}.tsv"
            removed_file.write_text(
                "Orthogroup\n" + "\n".join(removed) + "\n", encoding="utf-8"
            )
            logger.info(
                f"[{tag}] 实际剔除 {len(removed)} 个 OG，写出: {removed_file.name}"
            )

            current_family = next_family

        # 若 max_autofix_rounds 为 0，则上面的循环不会进，直接跑了一次 primary

        # ---------- large: 只跑一轮，不做极端修剪 ----------
        if large_family is not None and large_family.is_file():
            large_dir = ensure_dir(model_dir / large_tag)
            tag_l = f"MODEL={model_name} LARGE"
            logger.info(f"[{tag_l}] 开始 CAFE5 运行（不做极端修剪）")
            run_cafe5_once(
                cafe5_bin=cafe5_bin,
                family_file=large_family,
                tree_file=tree_file,
                threads=threads,
                gamma_k=gamma_k,
                pvalue=pvalue,
                work_dir=large_dir,
                logger=logger,
                tag=tag_l,
            )
        else:
            logger.info(f"[MODEL {model_name}] 无 large 集，跳过 large 阶段")

    # 运行完成，写 sentinel
    done_flag = cafe_root / ".cafe.done"
    done_flag.write_text("ok\n", encoding="utf-8")
    logger.info(f"[DONE] 所有 CAFE5 模型运行完成，写出 {done_flag}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n")
        sys.exit(2)