#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
12_cafe5_run_models.py —— 运行 CAFE5（教程增强版：两阶段 large / 自动剔除极端家族 / 标记高失败 OG）

核心特性：
  1) 支持 config.cafe5.models 中的多个模型（目前皇上一般只用 "global"）；
  2) 自动拆分 large 家族（copy 数 >= copy_threshold）：
       - primary 部分：family.primary.tsv  -> primary_global/ 目录
       - large   部分：family.large.tsv    -> large/ 目录
  3) 自动多轮剔除极端家族：
       - 检测日志中的 "Failed to initialize any reasonable values" +
         "Families with largest size differentials" 区块，
         每轮删除若干 OG，最多 max_autofix_rounds 轮；
       - 被删除的家族写入：
           flags/extreme_ogs_roundN.list
           autofix_removed_roundN.tsv
       - 每轮的测试运行放在 _roundN_tmp/ 目录，最终稳定结果放在 primary_global/。
  4) 记录高失败率 OG：
       - 解析 "The following families had failure rates >20% of the time" 区域，
         收集 OG 列表，写入 flags/high_fail_ogs.list（供 13 脚本使用）。
  5) 日志：
       - 主日志：paths.logs_dir/12_cafe5_run_models.log
       - 流式输出：CAFE5 的输出逐行同步到屏幕和日志（不会“憋大招”）。
  6) 输出路径已缩短：
       - CAFE5 调用统一使用 "-o ."，产物直接落在 primary_global/ 或 large/ 下，
         不再多出一层 results/ 子目录。

配置约定（config.yaml 中 cafe5 小节）：
  cafe5:
    enable: true
    threads: 30
    gamma_k: 3          # -k
    pvalue:  0.05       # -P
    models:
      - "global"

    # 自动修正轮数（可选，默认 0 表示不自动修正）
    max_autofix_rounds: 3

    # 两阶段：拆分 large 家族（可选）
    two_stage_large:
      enable: true
      copy_threshold: 100      # 某物种 copy >= 该值 → 认为是 large 家族
      keep_strategy: "split"   # 目前只实现 split；其它值等价于 split
      primary_tag: "primary"   # 仅用于日志标签
      large_tag:   "large"

其它配置（paths、binaries 部分）与之前脚本保持一致。
"""

from __future__ import annotations
import sys
import io
import logging
import subprocess
import re
import shutil
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional

import yaml

DEFAULT_CONFIG = "config.yaml"

# ====================== 通用工具：配置、路径、日志 ======================

def _expand_publish_placeholders(obj, publish_dir: str):
    if isinstance(obj, str):
        return obj.replace("<publish_dir>", publish_dir)
    if isinstance(obj, list):
        return [_expand_publish_placeholders(x, publish_dir) for x in obj]
    if isinstance(obj, dict):
        return {k: _expand_publish_placeholders(v, publish_dir) for k, v in obj.items()}
    return obj


def load_config(config_path: str = DEFAULT_CONFIG) -> Dict[str, Any]:
    """读取 config.yaml，并展开 <publish_dir> 占位符。"""
    p = Path(config_path)
    if not p.exists():
        raise FileNotFoundError(f"[ERR] 未找到配置文件：{p}")
    cfg = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    pub = cfg.get("publish_dir")
    if pub:
        cfg["inputs"] = _expand_publish_placeholders(cfg.get("inputs", {}), str(pub))
    return cfg


def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def need_dir(p: Path, what: str) -> Path:
    p = Path(p)
    if not p.is_dir():
        raise FileNotFoundError(f"[ERR] 缺少目录：{what} -> {p}")
    return p


def need_file(p: Path, what: str) -> Path:
    p = Path(p)
    if not p.is_file():
        raise FileNotFoundError(f"[ERR] 缺少文件：{what} -> {p}")
    return p


def clean_dir(p: Path):
    """清空目录（保留目录本身）。"""
    if p.is_dir():
        for child in p.iterdir():
            if child.is_dir():
                shutil.rmtree(child)
            else:
                try:
                    child.unlink()
                except FileNotFoundError:
                    pass
    p.mkdir(parents=True, exist_ok=True)


def get_logger(name: str, logfile: Path, level: int = logging.INFO) -> logging.Logger:
    """日志写文件 + 同步屏幕，并将 stdout/stderr 包装成自动 flush。"""
    ensure_dir(logfile.parent)
    lg = logging.getLogger(name)
    lg.setLevel(level)
    lg.handlers.clear()

    fmt = logging.Formatter(
        "[%(asctime)s] %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S"
    )
    fh = logging.FileHandler(logfile, encoding="utf-8")
    fh.setFormatter(fmt)
    fh.setLevel(level)
    sh = logging.StreamHandler(stream=sys.stdout)
    sh.setFormatter(fmt)
    sh.setLevel(level)
    lg.addHandler(fh)
    lg.addHandler(sh)

    class _Flush(io.TextIOBase):
        def __init__(self, s):
            self.s = s

        def write(self, x):
            self.s.write(x)
            self.s.flush()
            return len(x)

    sys.stdout = _Flush(sys.stdout)
    sys.stderr = _Flush(sys.stderr)
    return lg


def banner(logger: logging.Logger, text: str):
    bar = "=" * max(10, len(text) + 2)
    logger.info(bar)
    logger.info(f" {text} ")
    logger.info(bar)


def write_done(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    Path(path).touch()


# ====================== CAFE5 辅助：日志解析 ======================

def parse_extreme_ogs_from_text(text: str) -> List[str]:
    """
    从 CAFE5 输出中解析“极端家族”（Families with largest size differentials 块）。
    逻辑：取最后一个块，将其中的 OG0000012: 460 之类提取出来。
    """
    blocks = re.findall(
        r"Families with largest size differentials:\s*\n((?:.*OG\d+:\s*\d+\s*\n)+)",
        text,
        flags=re.I,
    )
    if not blocks:
        return []
    last = blocks[-1]
    ogs: List[str] = []
    for m in re.finditer(r"(OG\d+)\s*:\s*(\d+)", last):
        ogs.append(m.group(1))
    return ogs


def parse_high_fail_ogs_from_text(text: str) -> List[str]:
    """
    从输出中解析“失败率 >20%”的家族列表，返回 OG ID 列表。
    """
    lines = text.splitlines()
    res: List[str] = []
    flag = False
    for ln in lines:
        if "The following families had failure rates >20% of the time" in ln:
            flag = True
            continue
        if not flag:
            continue
        if not ln.strip():
            # 空行结束该块
            break
        m = re.search(r"(OG\d+)\s+had\s+(\d+)\s+failures", ln)
        if m:
            res.append(m.group(1))
    return res


def run_cafe_streaming(
    cmd: List[str],
    work_dir: Path,
    log: logging.Logger,
    model_tag: str,
    round_tag: str = "",
) -> str:
    """
    以流式方式运行 CAFE5：
      - 实时将 stdout 写入 logger（屏幕 + 日志）；
      - 将所有输出累积为字符串返回，用于后续正则解析。
    """
    tag = f"[MODEL={model_tag}{(' ' + round_tag) if round_tag else ''}]"
    log.info(f"{tag} CMD: " + " ".join(map(str, cmd)))

    proc = subprocess.Popen(
        cmd,
        cwd=str(work_dir),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True,
    )

    stdout_lines: List[str] = []
    log.info(f"{tag} ---- CAFE5 输出开始 ----")

    assert proc.stdout is not None
    for line in proc.stdout:
        line = line.rstrip("\n")
        stdout_lines.append(line)
        log.info(f"{tag} {line}")

    proc.wait()
    log.info(f"{tag} ---- CAFE5 输出结束 ----")

    text = "\n".join(stdout_lines)

    if proc.returncode != 0:
        raise RuntimeError(
            f"[ERR] CAFE5 退出码 {proc.returncode}：模型 {model_tag}；请查看日志。"
        )

    return text


# ====================== CAFE5 辅助：family 拆分与修正 ======================

def split_family_primary_large(
    family_tsv: Path,
    mdir: Path,
    copy_threshold: int,
    keep_strategy: str,
    log: logging.Logger,
    model: str,
) -> Tuple[Path, Optional[Path]]:
    """
    根据 copy_threshold 拆分 family.tsv 为 primary / large 两份，
    并写入模型目录下：
      - family.primary.tsv
      - family.large.tsv （若有 large 家族）
    返回 (primary_path, large_path_or_None)。
    """
    keep_strategy = (keep_strategy or "split").lower()
    tag = f"[MODEL {model}]"

    txt = family_tsv.read_text(encoding="utf-8").splitlines()
    if not txt:
        raise RuntimeError(f"[ERR] family.tsv 为空：{family_tsv}")

    header = txt[0]
    rows = txt[1:]

    primary_lines = [header]
    large_lines = [header]
    large_n = 0

    for ln in rows:
        if not ln.strip():
            continue
        parts = ln.rstrip("\n").split("\t")
        # 皇上确认：首列一定是 Desc（OG ID），后面都是拷贝数（整数或 0）
        counts = []
        for x in parts[1:]:
            x = x.strip()
            if not x:
                counts.append(0)
            else:
                try:
                    counts.append(int(x))
                except Exception:
                    counts.append(0)
        max_copy = max(counts) if counts else 0

        if max_copy >= copy_threshold:
            large_lines.append(ln.rstrip("\n"))
            large_n += 1
        else:
            primary_lines.append(ln.rstrip("\n"))

    primary_path = mdir / "family.primary.tsv"
    primary_path.write_text("\n".join(primary_lines) + "\n", encoding="utf-8")

    large_path: Optional[Path] = None
    if keep_strategy == "split" and large_n > 0:
        large_path = mdir / "family.large.tsv"
        large_path.write_text("\n".join(large_lines) + "\n", encoding="utf-8")

    log.info(
        f"{tag} 拆分完成：primary {len(primary_lines)-1} 行，large {large_n} 行；"
        f"阈值 copy_threshold={copy_threshold}"
    )
    return primary_path, large_path


def remove_ogs_from_family(
    src_family: Path,
    dst_family: Path,
    ogs_to_drop: List[str],
) -> int:
    """
    从 family TSV 中删除指定 OG（首列为 Desc），返回删除的行数。
    """
    og_set = set(ogs_to_drop)
    txt = src_family.read_text(encoding="utf-8").splitlines()
    if not txt:
        return 0
    header = txt[0]
    rows = txt[1:]

    kept = [header]
    removed = 0
    for ln in rows:
        if not ln.strip():
            continue
        parts = ln.rstrip("\n").split("\t")
        og = parts[0].strip()
        if og in og_set:
            removed += 1
            continue
        kept.append(ln.rstrip("\n"))

    dst_family.write_text("\n".join(kept) + "\n", encoding="utf-8")
    return removed


# ====================== 主流程：单模型运行（含多轮自动修正） ======================

def run_model_with_autofix(
    model: str,
    cafe_bin: str,
    tree_path: Path,
    family_orig: Path,
    out_root: Path,
    threads: int,
    gamma_k: Optional[int],
    pvalue: Optional[float],
    max_autofix_rounds: int,
    log: logging.Logger,
):
    """
    针对某一个 model（例如 "global"）：
      - 以 family_orig 为起点，最多进行 max_autofix_rounds 轮自动剔除极端家族；
      - 每轮运行 CAFE5，输出放到 _roundN_tmp/ 目录；
      - 最终稳定轮结果移入 primary_global/；
      - 解析 high-fail OG，写 flags/high_fail_ogs.list。
    """
    tag = f"[MODEL={model}]"
    mdir = ensure_dir(out_root / model)
    flags_dir = ensure_dir(mdir / "flags")
    sentinels = ensure_dir(mdir / "sentinels")
    primary_dir = ensure_dir(mdir / "primary_global")
    large_dir = ensure_dir(mdir / "large")

    # 开跑前：清空 primary_global / large / sentinels 目录（flags 保留历史标记）
    clean_dir(primary_dir)
    clean_dir(large_dir)
    clean_dir(sentinels)

    # 先根据 two_stage_large 把 family 拆成 primary / large
    # 注意：family_orig 是 paths.cafe_run_dir/input/family.tsv
    cfg = load_config()
    cafe_cfg = cfg.get("cafe5", {}) or {}
    two_stage = cafe_cfg.get("two_stage_large", {}) or {}
    ts_enable = bool(two_stage.get("enable", False))
    copy_threshold = int(two_stage.get("copy_threshold", 0) or 0)
    keep_strategy = str(two_stage.get("keep_strategy", "split") or "split")

    if ts_enable and copy_threshold > 0:
        family_primary, family_large = split_family_primary_large(
            family_orig,
            mdir,
            copy_threshold=copy_threshold,
            keep_strategy=keep_strategy,
            log=log,
            model=model,
        )
    else:
        # 不拆分，primary 就是原始 family，large 为空
        family_primary = family_orig
        (mdir / "family.primary.tsv").write_text(
            family_orig.read_text(encoding="utf-8"), encoding="utf-8"
        )
        family_large = None
        log.info(f"{tag} two_stage_large 关闭或阈值为 0，本轮不拆分 large 家族。")

    # ---------------- 第一部分：primary + 自动剔除极端 OG ----------------

    current_family = family_primary
    extreme_total: List[str] = []

    rounds = max(0, int(max_autofix_rounds or 0))
    if rounds == 0:
        rounds = 1  # 至少跑一轮（不做剔除）

    final_success = False
    last_stdout = ""

    for r in range(1, rounds + 1):
        round_tag = f"PRIMARY-GLOBAL ROUND-{r}"
        tmp_dir = ensure_dir(mdir / f"_round{r}_tmp")
        clean_dir(tmp_dir)

        cmd = [
            cafe_bin,
            "-i",
            str(current_family.resolve()),
            "-t",
            str(tree_path.resolve()),
            "-c",
            str(threads),
            "-o",
            ".",  # 结果直接写到 tmp_dir
        ]
        if gamma_k:
            cmd += ["-k", str(gamma_k)]
        if pvalue is not None:
            cmd += ["-P", str(pvalue)]

        log.info(f"{tag} [{round_tag}] CAFE5 启动准备")
        stdout_text = run_cafe_streaming(cmd, tmp_dir, log, model_tag=model, round_tag=round_tag)
        last_stdout = stdout_text

        # 检查是否还有 “Failed to initialize any reasonable values”
        if "Failed to initialize any reasonable values" not in stdout_text:
            log.info(f"{tag} [{round_tag}] 未检测到初始化失败提示，视为稳定轮。")
            # 将 tmp_dir 中的结果文件全部移入 primary_global/
            clean_dir(primary_dir)
            for p in tmp_dir.iterdir():
                target = primary_dir / p.name
                if target.exists():
                    if target.is_dir():
                        shutil.rmtree(target)
                    else:
                        target.unlink()
                shutil.move(str(p), str(target))
            final_success = True
            break

        if r >= rounds:
            # 已达最大修正轮数，仍然失败 → 报错
            log.error(
                f"{tag} [{round_tag}] 仍存在 'Failed to initialize any reasonable values'，"
                f"且已达到最大自动修正轮数（{rounds}）。"
            )
            raise RuntimeError(
                f"[ERR] 模型 {model} 自动修正用尽仍失败，请检查 family.tsv/树/物种匹配情况。"
            )

        # 解析极端 OG 并生成下一轮 family
        ogs = parse_extreme_ogs_from_text(stdout_text)
        if not ogs:
            log.error(
                f"{tag} [{round_tag}] 检测到初始化失败，但未能解析出任何极端家族 OG。"
            )
            raise RuntimeError(
                f"[ERR] 模型 {model} 无法从日志中提取极端 OG，请检查 CAFE5 输出。"
            )

        extreme_total.extend(ogs)
        (flags_dir / f"extreme_ogs_round{r}.list").write_text(
            "\n".join(ogs) + "\n", encoding="utf-8"
        )

        next_family = mdir / f"family.autofix{r}.tsv"
        removed = remove_ogs_from_family(current_family, next_family, ogs_to_drop=ogs)
        log.info(
            f"{tag} [FILTER] {next_family.name}: 从 {current_family.name} 删除 {removed} 行 OG，得到 {next_family.name}"
        )
        (mdir / f"autofix_removed_round{r}.tsv").write_text(
            "OG\n" + "\n".join(ogs) + "\n", encoding="utf-8"
        )

        current_family = next_family
        log.info(
            f"{tag} [CLEAN] 第 {r} 轮自动修正完成，进入下一轮。"
        )

    if not final_success:
        raise RuntimeError(
            f"[ERR] 模型 {model} 未能获得稳定 primary_global 结果，请检查日志。"
        )

    # ---------------- 第二部分：标记高失败家族（primary） ----------------

    high_fail_ogs = parse_high_fail_ogs_from_text(last_stdout)
    if high_fail_ogs:
        (flags_dir / "high_fail_ogs.list").write_text(
            "\n".join(high_fail_ogs) + "\n", encoding="utf-8"
        )
        log.info(
            f"{tag} 标记 high-fail 家族 {len(high_fail_ogs)} 个（failure rate >20%），"
            f"写入 {flags_dir/'high_fail_ogs.list'}"
        )
    else:
        log.info(f"{tag} 未检测到 high-fail 家族（failure rate >20% 区块缺失或为空）。")

    # ---------------- 第三部分：large 家族（若存在） ----------------

    if (mdir / "family.large.tsv").is_file():
        family_large = mdir / "family.large.tsv"
        if family_large.stat().st_size > 0:
            clean_dir(large_dir)
            round_tag = "LARGE"
            cmd = [
                cafe_bin,
                "-i",
                str(family_large.resolve()),
                "-t",
                str(tree_path.resolve()),
                "-c",
                str(threads),
                "-o",
                ".",  # 直接输出到 large_dir
            ]
            if gamma_k:
                cmd += ["-k", str(gamma_k)]
            if pvalue is not None:
                cmd += ["-P", str(pvalue)]

            log.info(f"{tag} [{round_tag}] 开始对 large 家族单独跑 CAFE5")
            _ = run_cafe_streaming(cmd, large_dir, log, model_tag=model, round_tag=round_tag)
            log.info(f"{tag} [{round_tag}] large 家族 CAFE5 完成。")
        else:
            log.info(f"{tag} large 家族文件为空，跳过第二阶段。")
    else:
        log.info(f"{tag} 未检测到 large 家族文件 family.large.tsv，跳过第二阶段。")

    # sentinel：标记该 model 已完成
    (sentinels / "primary_global.done").write_text("ok\n", encoding="utf-8")
    log.info(f"{tag} primary_global & large 全部完成。")


# ====================== 入口函数 ======================

def main():
    cfg = load_config()
    paths = cfg["paths"]
    cafe_cfg = cfg.get("cafe5", {}) or {}

    logs_dir = Path(paths["logs_dir"])
    log = get_logger("aphylo.12", logs_dir / "12_cafe5_run_models.log")
    banner(log, "APhylo 12 — CAFE5（教程增强版）")

    # 若关闭 CAFE5，直接写 done
    if not cafe_cfg.get("enable", True):
        log.info("cafe5.enable = false —— 跳过本步。")
        write_done(Path(paths["cafe_run_dir"]) / ".cafe.done")
        return

    threads = int(cafe_cfg.get("threads", 4) or 4)
    gamma_k = cafe_cfg.get("gamma_k", None)
    gamma_k = int(gamma_k) if gamma_k not in (None, "", "null") else None
    pvalue = cafe_cfg.get("pvalue", None)
    pvalue = float(pvalue) if pvalue not in (None, "", "null") else None
    models = cafe_cfg.get("models", ["global"]) or ["global"]
    max_autofix_rounds = int(cafe_cfg.get("max_autofix_rounds", 3) or 3)

    cafe_bin = cfg.get("binaries", {}).get("cafe5", "cafe5")

    # 输入：family.tsv + utree_for_cafe.nwk
    cafe_run_dir = Path(paths["cafe_run_dir"])
    input_dir = need_dir(cafe_run_dir / "input", "CAFE 输入目录")

    tsvs = sorted(input_dir.glob("*.tsv"))
    if not tsvs:
        raise FileNotFoundError("[ERR] CAFE 输入目录中未找到 family TSV（*.tsv）")
    family = need_file(tsvs[0], "CAFE family.tsv")

    nwks = sorted(input_dir.glob("*.nwk"))
    if not nwks:
        raise FileNotFoundError("[ERR] CAFE 输入目录中未找到树文件（*.nwk）")
    tree = need_file(nwks[0], "CAFE 超时钟树")

    log.info(f"[IN] family: {family.resolve()}")
    log.info(f"[IN] tree:   {tree.resolve()}")
    log.info(f"[CFG] max_autofix_rounds = {max_autofix_rounds}")

    models_root = ensure_dir(cafe_run_dir / "models")

    for model in models:
        run_model_with_autofix(
            model=str(model),
            cafe_bin=cafe_bin,
            tree_path=tree,
            family_orig=family,
            out_root=models_root,
            threads=threads,
            gamma_k=gamma_k,
            pvalue=pvalue,
            max_autofix_rounds=max_autofix_rounds,
            log=log,
        )

    # 全部模型完成，写全局 done
    write_done(cafe_run_dir / ".cafe.done")
    log.info("[DONE] 所有 CAFE5 模型运行完成。")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(str(e) + "\n")
        sys.exit(2)