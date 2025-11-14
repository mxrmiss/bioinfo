#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
12_cafe5_run_models.py —— CAFE5 运行脚本（教程增强版）
===================================================
功能概述：
  1) 读取 config.yaml 中 cafes 配置与路径。
  2) 拆分 family.tsv 为：
       - family.primary.tsv  —— 主要家族
       - family.large.tsv    —— 大拷贝家族（可选）
  3) 对 primary 阶段执行“极端 OG 自动剔除”多轮迭代：
       - 解析 CAFE5 日志中的 “Families with largest size differentials”
       - 从 family.primary.tsv 中删除对应 Orthogroup（第一列为 Desc，第二列为 OG）
       - 最多 max_autofix_rounds 轮（配置中设置）
  4) large 阶段只运行一次 CAFE5，不再做极端 OG 剔除。
  5) 解析日志中 “had XX failures” 行，统计高失败率 OG，
     写入 models/<model>/flags/high_fail_ogs.list。
  6) 生成 .cafe.done 哨兵文件，供 13_cafe5_aggregate.py 检查。

说明：
  - 所有参数统一从 config.yaml 读取；本脚本不使用命令行参数。
  - CAFE5 需要在当前环境中可直接执行（cafe5 命令可用）。
"""

from __future__ import annotations

import sys
import io
import os
import re
import logging
import subprocess
import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple

import yaml

# ======================== 基础工具 ========================

DEFAULT_CONFIG = "config.yaml"


def _expand_publish_placeholders(obj, publish_dir: str):
    """把 <publish_dir> 占位符展开（与 13 脚本保持一致）"""
    if isinstance(obj, str):
        return obj.replace("<publish_dir>", publish_dir)
    if isinstance(obj, list):
        return [_expand_publish_placeholders(x, publish_dir) for x in obj]
    if isinstance(obj, dict):
        return {k: _expand_publish_placeholders(v, publish_dir) for k, v in obj.items()}
    return obj


def load_config(path: str = DEFAULT_CONFIG) -> Dict[str, Any]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"[ERR] 未找到配置文件：{p}")
    cfg = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    pub = cfg.get("publish_dir")
    if pub:
        cfg["inputs"] = _expand_publish_placeholders(cfg.get("inputs", {}), str(pub))
        cfg["cafes"] = _expand_publish_placeholders(cfg.get("cafes", {}), str(pub))
    return cfg


def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def get_logger(logfile: Path) -> logging.Logger:
    ensure_dir(logfile.parent)
    logger = logging.getLogger("aphylo.12")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter(
        "[%(asctime)s] %(levelname)s - %(message)s",
        "%Y-%m-%d %H:%M:%S",
    )

    fh = logging.FileHandler(logfile, encoding="utf-8")
    fh.setFormatter(fmt)
    fh.setLevel(logging.INFO)
    logger.addHandler(fh)

    sh = logging.StreamHandler(stream=sys.stdout)
    sh.setFormatter(fmt)
    sh.setLevel(logging.INFO)
    logger.addHandler(sh)

    # 替换 stdout/stderr，保证流式输出立即刷新
    class _Flush(io.TextIOBase):
        def __init__(self, s): self.s = s
        def write(self, x):
            self.s.write(x)
            self.s.flush()
            return len(x)

    sys.stdout = _Flush(sys.stdout)
    sys.stderr = _Flush(sys.stderr)

    return logger


def banner(log: logging.Logger, text: str):
    bar = "=" * max(20, len(text) + 4)
    log.info(bar)
    log.info(f"  {text}  ")
    log.info(bar)


# ======================== family.tsv 处理 ========================

def split_family_by_copy(
    family_in: Path,
    model_dir: Path,
    threshold: int,
    logger: logging.Logger,
) -> Tuple[Path, Optional[Path]]:
    """
    按 copy 阈值把 family.tsv 拆成 primary / large 两份.
    family.tsv 格式：
      Desc  Orthogroup  sp1  sp2  ...
    """
    primary_out = model_dir / "family.primary.tsv"
    large_out = model_dir / "family.large.tsv"

    ensure_dir(model_dir)

    hdr: Optional[str] = None
    primary_lines: List[str] = []
    large_lines: List[str] = []

    with family_in.open(encoding="utf-8") as f:
        for i, line in enumerate(f):
            if i == 0:
                hdr = line.rstrip("\n")
                primary_lines.append(hdr)
                large_lines.append(hdr)
                continue
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            # 第一列 Desc，第二列 Orthogroup，从第三列起是拷贝数
            try:
                counts = [int(x) for x in cols[2:]]
            except ValueError:
                # 若存在非整数，直接丢进 primary，避免误删
                primary_lines.append(line.rstrip("\n"))
                continue
            max_cnt = max(counts) if counts else 0
            if max_cnt >= threshold:
                large_lines.append(line.rstrip("\n"))
            else:
                primary_lines.append(line.rstrip("\n"))

    primary_out.write_text("\n".join(primary_lines) + "\n", encoding="utf-8")

    if len(large_lines) > 1:
        large_out.write_text("\n".join(large_lines) + "\n", encoding="utf-8")
        logger.info(
            f"[SPLIT] 阈值 {threshold}: primary≈{len(primary_lines)-1} 行, "
            f"large≈{len(large_lines)-1} 行"
        )
        return primary_out, large_out
    else:
        if large_out.exists():
            large_out.unlink()
        logger.info(
            f"[SPLIT] 阈值 {threshold}: 未发现 large 家族，全部作为 primary 处理"
        )
        return primary_out, None


# ======================== 日志解析：极端 OG & 高失败率 ========================

EXTREME_BLOCK_RE = re.compile(
    r"Families with largest size differentials:\s*\n"
    r"((?:.*OG\d+:\s*\d+\s*\n)+)",
    re.MULTILINE,
)

OG_LINE_RE = re.compile(r"(OG\d+)\s*:\s*\d+")

FAIL_LINE_RE = re.compile(r"\b(OG\d+)\b.*\bhad\s+(\d+)\s+failures\b", re.IGNORECASE)


def extract_extreme_ogs_from_log(runlog: Path) -> List[str]:
    """从 run.log 中提取最后一个“极端家族”块里的 OG 列表"""
    if not runlog.is_file():
        return []
    txt = runlog.read_text(encoding="utf-8", errors="ignore")
    blocks = EXTREME_BLOCK_RE.findall(txt)
    if not blocks:
        return []
    last = blocks[-1]
    ogs = OG_LINE_RE.findall(last)
    # 去重保持顺序
    seen = set()
    res = []
    for og in ogs:
        if og not in seen:
            seen.add(og)
            res.append(og)
    return res


def remove_ogs_from_family(
    family_path: Path,
    ogs: List[str],
    out_path: Path,
    removed_out: Path,
) -> int:
    """
    从 family_path 中删除给定 OG（第二列 Orthogroup），写到 out_path，
    被删除行写到 removed_out。
    返回删除的行数。
    """
    og_set = set(ogs)
    removed: List[str] = []
    kept: List[str] = []

    with family_path.open(encoding="utf-8") as f:
        for i, line in enumerate(f):
            if i == 0:
                kept.append(line.rstrip("\n"))
                continue
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            og = cols[1] if len(cols) > 1 else ""
            if og in og_set:
                removed.append(line.rstrip("\n"))
            else:
                kept.append(line.rstrip("\n"))

    out_path.write_text("\n".join(kept) + "\n", encoding="utf-8")
    if removed:
        removed_out.write_text("\n".join(removed) + "\n", encoding="utf-8")
    return len(removed)


def extract_high_fail_ogs(runlog: Path) -> List[str]:
    """从 run.log 中抓取 'OGxxxxx had NN failures' 行，返回 OG 列表"""
    if not runlog.is_file():
        return []
    txt = runlog.read_text(encoding="utf-8", errors="ignore")
    ogs: List[str] = []
    seen = set()
    for og, _cnt in FAIL_LINE_RE.findall(txt):
        if og not in seen:
            seen.add(og)
            ogs.append(og)
    return ogs


# ======================== 运行 CAFE5 的封装 ========================

def run_cafe5(
    logger: logging.Logger,
    tag: str,
    cmd: List[str],
    workdir: Path,
    runlog: Path,
) -> int:
    """
    以流式方式运行 CAFE5:
      - stdout/stderr 合并
      - 每行加上 tag 写入 runlog & logger.info
    """
    ensure_dir(workdir)
    logger.info(f"[{tag}] CMD: {' '.join(cmd)} (cwd={workdir})")

    with runlog.open("a", encoding="utf-8") as lf:
        lf.write(
            f"========== {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} "
            f"{tag} START ==========\n"
        )
        lf.flush()
        proc = subprocess.Popen(
            cmd,
            cwd=str(workdir),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            line = line.rstrip("\n")
            pref = f"[{tag}] {line}"
            logger.info(pref)
            lf.write(pref + "\n")
        proc.wait()
        lf.write(
            f"========== {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} "
            f"{tag} END (RC={proc.returncode}) ==========\n"
        )
        lf.flush()
    if proc.returncode != 0:
        logger.error(f"[{tag}] CAFE5 退出码 {proc.returncode}")
    else:
        logger.info(f"[{tag}] CAFE5 运行完成")
    return proc.returncode


def move_previous_outputs(run_dir: Path, round_idx: int, logger: logging.Logger):
    """
    把上一轮 primary_global 输出文件移到 _round{idx}_tmp 子目录，仅保留 run.log。
    """
    tmp = ensure_dir(run_dir / f"_round{round_idx}_tmp")
    moved = 0
    for p in run_dir.iterdir():
        if not p.is_file():
            continue
        if p.name == "run.log":
            continue
        dest = tmp / p.name
        p.replace(dest)
        moved += 1
    if moved:
        logger.info(f"[PRIMARY-GLOBAL] 第 {round_idx} 轮产物已移入 {tmp}")


# ======================== 主流程 ========================

def main():
    cfg = load_config()
    paths_cfg = cfg.get("paths", {}) or {}
    cafes_cfg = cfg.get("cafes", {}) or {}

    cafe_run_dir = Path(paths_cfg.get("cafe_run_dir", "results/06_cafe"))
    logs_dir = Path(paths_cfg.get("logs_dir", "logs"))
    models_cfg = cafes_cfg.get("models", {}) or {}
    if not models_cfg:
        raise SystemExit("[ERR] config.cafes.models 为空，至少需要一个模型（如 global）。")

    two_stage_cfg = cafes_cfg.get("two_stage_large", {}) or {}
    ts_enable = bool(two_stage_cfg.get("enable", False))
    ts_threshold = int(two_stage_cfg.get("copy_threshold", 100))
    max_rounds = int(two_stage_cfg.get("max_autofix_rounds", 3))

    logger = get_logger(logs_dir / "12_cafe5_run_models.log")
    banner(logger, "APhylo 12 — CAFE5（教程增强版）")

    logger.info(f"[CFG] 运行目录: {cafe_run_dir}")
    ensure_dir(cafe_run_dir)

    models_root = ensure_dir(cafe_run_dir / "models")

    all_high_fail: List[str] = []

    # -------- 单模型循环（目前一般只有 global） --------
    for model_name, mcfg in models_cfg.items():
        logger.info(f"[MODEL {model_name}] ===============================")

        tree_path = Path(mcfg["tree"])
        family_in = Path(mcfg["family"])
        threads = int(mcfg.get("threads", 4))
        extra = str(mcfg.get("cmd_extra", "") or "").strip()

        if not tree_path.is_file():
            raise SystemExit(f"[ERR] 模型 {model_name} 的树不存在：{tree_path}")
        if not family_in.is_file():
            raise SystemExit(f"[ERR] 模型 {model_name} 的 family.tsv 不存在：{family_in}")

        logger.info(f"[MODEL {model_name}] [IN] family: {family_in}")
        logger.info(f"[MODEL {model_name}] [IN] tree:   {tree_path}")
        logger.info(f"[MODEL {model_name}] [CFG] max_autofix_rounds = {max_rounds}")

        model_dir = ensure_dir(models_root / model_name)
        flags_dir = ensure_dir(model_dir / "flags")

        # ---------- 拆分 primary / large ----------
        if ts_enable:
            primary_family, large_family = split_family_by_copy(
                family_in, model_dir, ts_threshold, logger
            )
        else:
            primary_family = model_dir / "family.primary.tsv"
            large_family = None
            # 直接拷贝原 family
            primary_family.write_text(
                family_in.read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            logger.info(f"[SPLIT] two_stage_large.disable=true：全部家族作为 primary 处理")

        # ---------- PRIMARY 阶段：多轮极端 OG 剔除 ----------
        primary_dir = ensure_dir(model_dir / "primary_global")
        primary_runlog = primary_dir / "run.log"
        high_fail_this_model: List[str] = []

        current_family = primary_family

        for r in range(1, max_rounds + 1):
            tag = f"MODEL={model_name} PRIMARY-GLOBAL ROUND={r}"
            logger.info(f"[{tag}] CAFE5 输出目录: {primary_dir}")

            cmd = [
                "cafe5",
                "-i",
                str(current_family),
                "-t",
                str(tree_path),
                "-c",
                str(threads),
                "-k",
                "3",
                "-P",
                "0.05",
            ]
            if extra:
                cmd.extend(extra.split())

            rc = run_cafe5(logger, tag, cmd, primary_dir, primary_runlog)
            if rc != 0:
                raise SystemExit(f"[ERR] {tag} 运行失败，退出。")

            # 解析极端 OG
            extreme_ogs = extract_extreme_ogs_from_log(primary_runlog)
            logger.info(f"[{tag}] 检测到极端 OG 数量：{len(extreme_ogs)}")

            if not extreme_ogs:
                logger.info(f"[{tag}] 未检测到极端 OG，停止自动修正。")
                break

            if r == max_rounds:
                logger.info(
                    f"[{tag}] 已达最大轮数({max_rounds})，不再继续剔除极端 OG。"
                )
                break

            # 输出被删除的家族列表
            removed_out = model_dir / f"autofix_removed_round{r}.tsv"
            new_family = model_dir / "family.primary.tmp.tsv"
            removed_cnt = remove_ogs_from_family(
                current_family,
                extreme_ogs,
                new_family,
                removed_out,
            )
            logger.info(
                f"[{tag}] 从 {current_family.name} 中剔除 {removed_cnt} 个极端家族，"
                f"写入 {removed_out.name}"
            )
            # 轮次产物归档
            move_previous_outputs(primary_dir, r, logger)
            # 替换 family.primary.tsv
            new_family.replace(primary_family)
            current_family = primary_family

        # 收集 high-fail OG（primary）
        hf_primary = extract_high_fail_ogs(primary_runlog)
        logger.info(
            f"[MODEL {model_name}] PRIMARY-GLOBAL 检测到 high-fail 家族：{len(hf_primary)}"
        )
        high_fail_this_model.extend(hf_primary)

        # ---------- LARGE 阶段：只跑一次，不做剔除 ----------
        if large_family and large_family.is_file():
            large_dir = ensure_dir(model_dir / "large")
            large_runlog = large_dir / "run.log"
            tag = f"MODEL={model_name} LARGE"

            logger.info(f"[{tag}] CAFE5 输出目录: {large_dir}")
            cmd = [
                "cafe5",
                "-i",
                str(large_family),
                "-t",
                str(tree_path),
                "-c",
                str(threads),
                "-k",
                "3",
                "-P",
                "0.05",
            ]
            if extra:
                cmd.extend(extra.split())

            rc = run_cafe5(logger, tag, cmd, large_dir, large_runlog)
            if rc != 0:
                raise SystemExit(f"[ERR] {tag} 运行失败，退出。")

            hf_large = extract_high_fail_ogs(large_runlog)
            logger.info(
                f"[MODEL {model_name}] LARGE 检测到 high-fail 家族：{len(hf_large)}"
            )
            high_fail_this_model.extend(hf_large)
        else:
            logger.info(
                f"[MODEL {model_name}] 未启用 large 阶段（无 large 家族或 two_stage_large.disable）"
            )

        # ---------- high-fail 汇总写入 flags ----------
        if high_fail_this_model:
            # 去重
            seen = set()
            uniq = []
            for og in high_fail_this_model:
                if og not in seen:
                    seen.add(og)
                    uniq.append(og)
            (flags_dir / "high_fail_ogs.list").write_text(
                "\n".join(uniq) + "\n", encoding="utf-8"
            )
            logger.info(
                f"[MODEL {model_name}] high_fail_ogs.list 写出 {len(uniq)} 条"
            )
            all_high_fail.extend(uniq)
        else:
            logger.info(f"[MODEL {model_name}] 未检测到 high-fail 家族。")

    # -------- 所有模型完成，写 .cafe.done --------
    done_flag = cafe_run_dir / ".cafe.done"
    done_flag.write_text(
        f"CAFE5 models finished at "
        f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n",
        encoding="utf-8",
    )
    logger.info("[DONE] 所有 CAFE5 模型运行完成，已写 .cafe.done")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n")
        sys.exit(2)