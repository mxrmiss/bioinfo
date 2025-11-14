#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 12_cafe5_run_models.py —— 读取 config.yaml 按固定口径运行 CAFE5（两阶段可选、误差模型可选）

import os, sys, re, csv, shutil, subprocess, time, datetime
from pathlib import Path
from typing import List, Tuple, Dict, Any, Optional
import yaml
import logging

CFG_PATH = "config.yaml"

# ---------- 基础工具 ----------
def _expand_publish_placeholders(obj, publish_dir: str):
    if isinstance(obj, str): return obj.replace("<publish_dir>", publish_dir)
    if isinstance(obj, list): return [_expand_publish_placeholders(x, publish_dir) for x in obj]
    if isinstance(obj, dict): return {k:_expand_publish_placeholders(v, publish_dir) for k,v in obj.items()}
    return obj

def load_config(path: str) -> Dict[str, Any]:
    p = Path(path)
    cfg = yaml.safe_load(p.read_text(encoding="utf-8"))
    pub = cfg.get("publish_dir")
    if pub:
        cfg["inputs"] = _expand_publish_placeholders(cfg.get("inputs", {}), str(pub))
    return cfg

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p

def get_logger(logfile: Path) -> logging.Logger:
    ensure_dir(logfile.parent)
    lg = logging.getLogger("aphylo.cafe12")
    lg.setLevel(logging.INFO)
    lg.handlers.clear()
    fmt = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(logfile, encoding="utf-8")
    fh.setFormatter(fmt); fh.setLevel(logging.INFO)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt); sh.setLevel(logging.INFO)
    lg.addHandler(fh); lg.addHandler(sh)
    return lg

def run_stream(cmd: List[str], cwd: Path, tee_file: Path, logger: logging.Logger) -> int:
    ensure_dir(tee_file.parent)
    with tee_file.open("w", encoding="utf-8") as fout:
        proc = subprocess.Popen(cmd, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, bufsize=1, universal_newlines=True)
        for line in proc.stdout:
            s = line.rstrip("\n")
            logger.info(s)
            fout.write(s + "\n"); fout.flush()
        proc.wait()
        return proc.returncode

def move_results_up(dir_: Path):
    # CAFE5 默认在当前目录生成 results/ 子目录，这里将其内容提升一级
    res = dir_ / "results"
    if res.is_dir():
        for f in res.iterdir():
            target = dir_ / f.name
            if target.exists():
                if target.is_dir(): shutil.rmtree(target)
                else: target.unlink()
            if f.is_dir(): shutil.move(str(f), str(target))
            else: shutil.move(str(f), str(target))
        shutil.rmtree(res)

# ---------- 解析器 ----------
def parse_high_fail_ogs(text: str) -> List[str]:
    out = []
    flag = False
    for ln in text.splitlines():
        if "failure rates >20%" in ln:
            flag = True
            continue
        if flag:
            m = re.search(r'\b(OG\d+)\b.*?\bhad\b', ln)
            if m: out.append(m.group(1))
            elif ln.strip() == "" or ln.startswith("["):  # 读到空行或新段落即停止
                break
    return out

def parse_extreme_ogs(text: str, top_n: int = 20) -> List[Tuple[str,int]]:
    blks = re.findall(r"Families with largest size differentials:\s*([\s\S]*?)\n\s*(?:You may|Failed|----|$)", text)
    if not blks: return []
    last = blks[-1]
    out = []
    for og, n in re.findall(r"(OG\d+)\s*:\s*(\d+)", last):
        out.append((og, int(n)))
    return out[:top_n]

def split_primary_large(family_tsv: Path, out_primary: Path, out_large: Path, threshold: int,
                        keep_strategy: str, logger: logging.Logger) -> Tuple[int,int]:
    rows_p, rows_l = [], []
    with family_tsv.open("r", encoding="utf-8") as f:
        rdr = csv.reader(f, delimiter="\t")
        header = next(rdr)
        # 列定位：固定口径，首列 Desc，第二列 Orthogroup，其余物种列
        try:
            idx_desc = header.index("Desc")
            idx_og = header.index("Orthogroup")
        except ValueError:
            raise RuntimeError("family.tsv 表头应包含列：Desc 与 Orthogroup")
        species_idx = [i for i in range(len(header)) if i not in (idx_desc, idx_og)]
        rows_p.append(header); rows_l.append(header)
        for r in rdr:
            if not r: continue
            vals = []
            for i in species_idx:
                try:
                    vals.append(int(r[i]))
                except Exception:
                    vals.append(0)
            mx = max(vals) if vals else 0
            if mx >= threshold:
                rows_l.append(r)
                if keep_strategy == "split":  # split：large 进入 large，primary 仅保留非大拷贝
                    pass
                else:
                    rows_p.append(r)
            else:
                rows_p.append(r)
    with out_primary.open("w", encoding="utf-8", newline="") as fo:
        csv.writer(fo, delimiter="\t").writerows(rows_p)
    with out_large.open("w", encoding="utf-8", newline="") as fo:
        csv.writer(fo, delimiter="\t").writerows(rows_l)
    logger.info(f"[SPLIT] 阈值 {threshold}: primary={len(rows_p)-1} 行, large={len(rows_l)-1} 行")
    return len(rows_p)-1, len(rows_l)-1

def write_list(lst: List[str], path: Path):
    ensure_dir(path.parent)
    path.write_text("\n".join(lst) + ("\n" if lst else ""), encoding="utf-8")

# ---------- 核心流程 ----------
def run_one_model(cfg: Dict[str, Any], model_name: str, logger: logging.Logger):
    cafe = cfg["cafes"]
    errm = cfg.get("error_model", {}) or {}
    two_stage = cfg.get("two_stage_large", {}) or {}

    cafe_root = Path(cfg["paths"]["cafe_run_dir"])
    models_root = cafe_root / "models"
    model_dir = ensure_dir(models_root / model_name)
    flags_dir = ensure_dir(model_dir / "flags")

    # 输入路径（固定口径）
    family_tsv = Path(cfg["inputs"]["family_tsv"]).resolve()
    utree = Path(cfg["inputs"]["ultrametric_tree"]).resolve()
    if not family_tsv.is_file():
        raise FileNotFoundError(f"family_tsv 不存在：{family_tsv}")
    if not utree.is_file():
        raise FileNotFoundError(f"ultrametric_tree 不存在：{utree}")

    logger.info("===== CAFE5 开始 =====")
    logger.info(f"[IN] family: {family_tsv}")
    logger.info(f"[IN] tree  : {utree}")

    # 分集
    primary_tag = two_stage.get("primary_tag", "primary")
    large_tag   = two_stage.get("large_tag", "large")
    keep_strategy = two_stage.get("keep_strategy", "split")
    threshold = int(two_stage.get("copy_threshold", 100))
    do_two_stage = bool(two_stage.get("enable", True))

    fam_primary = model_dir / f"family.{primary_tag}.tsv"
    fam_large   = model_dir / f"family.{large_tag}.tsv"

    if do_two_stage:
        split_primary_large(family_tsv, fam_primary, fam_large, threshold, keep_strategy, logger)
    else:
        shutil.copy2(family_tsv, fam_primary)
        # large 空表保留表头
        with family_tsv.open("r", encoding="utf-8") as fi:
            header = fi.readline()
        fam_large.write_text(header, encoding="utf-8")
        logger.info("[SPLIT] 关闭两阶段：仅 primary 参与")

    # 公共 CAFE5 选项
    cafe_bin = cfg["binaries"]["cafe5"]
    threads  = int(cafe.get("threads", 1))
    gamma_k  = int(cafe.get("gamma_k", 3))
    pval     = cafe.get("pvalue", None)
    base_cmd = [cafe_bin, "-t", str(utree), "-c", str(threads), "-k", str(gamma_k)]
    if pval not in (None, "", 0, 0.0): base_cmd += ["-P", str(pval)]

    # 误差模型（按 apply_to 控制）
    em_mode = (errm.get("mode") or "off").lower()
    em_file = errm.get("file")
    em_apply = [x.lower() for x in (errm.get("apply_to") or [])]

    # ---------- 子任务执行器 ----------
    def run_subset(tag: str, family_file: Path, autofix: bool):
        sub_dir = ensure_dir(model_dir / (f"{tag}_global"))
        runlog = sub_dir / "run.log"
        # CAFE5 运行
        def _run_once(f: Path) -> str:
            cmd = base_cmd + ["-i", str(f)]
            # 误差模型开关
            if em_mode == "estimate" and tag in em_apply:
                cmd += ["-e"]
            elif em_mode == "use" and em_file and tag in em_apply:
                cmd += [f"-e{em_file}"]
            logger.info(f"[CMD] {tag.upper()} 开跑: " + " ".join(cmd))
            rc = run_stream(cmd, cwd=sub_dir, tee_file=runlog, logger=logger)
            if rc != 0:
                logger.error(f"[ERR] {tag} 运行返回码 {rc}")
            # 提升 results/*
            move_results_up(sub_dir)
            return runlog.read_text(encoding="utf-8", errors="ignore")

        # 首轮
        text = _run_once(family_file)

        # 记录高失败率家族
        high_fail = parse_high_fail_ogs(text)
        if high_fail:
            write_list(high_fail, flags_dir / "high_fail_ogs.list")
            logger.info(f"[{tag}] 失败率>20% 家族：{len(high_fail)}")

        # 自动剔除（仅 primary，large 不执行）
        if not autofix:
            return

        rounds = int(cafe.get("max_autofix_rounds", 0))
        if rounds <= 0: return
        src_file = family_file
        for r in range(1, rounds + 1):
            # 失败特征：初始化失败提示
            if "Failed to initialize any reasonable values" not in text and \
               "Families with largest size differentials" not in text:
                logger.info(f"[{tag}] 无需自动修正，停止。")
                break
            ex = parse_extreme_ogs(text, top_n=20)
            if not ex:
                logger.info(f"[{tag}] 未解析到极端家族，停止自动修正。")
                break
            # 生成剔除清单
            rm_list = {og for og, _ in ex}
            (model_dir / f"autofix_removed_round{r}.tsv").write_text(
                "\n".join([f"{og}\t{n}" for og, n in ex]) + "\n", encoding="utf-8"
            )
            # family 过滤
            dst = model_dir / f"family.autofix{r}.tsv"
            with src_file.open("r", encoding="utf-8") as fi, dst.open("w", encoding="utf-8") as fo:
                rdr = csv.reader(fi, delimiter="\t")
                w = csv.writer(fo, delimiter="\t")
                header = next(rdr); w.writerow(header)
                idx_og = header.index("Orthogroup")
                kept = 0
                for row in rdr:
                    if row and row[idx_og] not in rm_list:
                        w.writerow(row); kept += 1
            logger.info(f"[{tag}] 第 {r} 轮剔除 {len(rm_list)} 个家族；保留 {kept}。")
            # 重跑
            text = _run_once(dst)
            src_file = dst

    # primary（含自动剔除）
    run_subset(primary_tag, fam_primary, autofix=True)
    # large（不剔除）
    if do_two_stage:
        cnt_large = sum(1 for _ in open(fam_large, encoding="utf-8")) - 1
        if cnt_large > 0:
            run_subset(large_tag, fam_large, autofix=False)
            logger.info("[INFO] primary_global & large 全部完成。")
        else:
            logger.info("[INFO] large 为空，跳过。")

def main():
    cfg = load_config(CFG_PATH)

    # 固定路径（无兼容）
    cafe_root = Path(cfg["paths"]["cafe_run_dir"])
    models_root = cafe_root / "models"
    logs_root = Path(cfg["paths"]["logs_dir"])
    ensure_dir(cafe_root); ensure_dir(models_root.parent); ensure_dir(logs_root)

    # 每次运行清空 models 与本脚本日志
    if models_root.exists(): shutil.rmtree(models_root)
    for p in logs_root.glob("12_cafe5_run_models*.log"):
        p.unlink()

    logger = get_logger(logs_root / "12_cafe5_run_models.log")
    logger.info("========== APhylo 12 — CAFE5 ==========")
    logger.info(f"[DIR] models_root: {models_root.resolve()}")
    logger.info(f"[DIR] logs_root  : {logs_root.resolve()}")

    # 模型列表（固定口径：list）
    models = cfg["cafes"]["models"]
    if not isinstance(models, list) or not models:
        raise RuntimeError("[ERR] config.cafes.models 需为非空列表，例如 ['global']")

    for m in models:
        logger.info(f"----- MODEL {m} 开始 -----")
        run_one_model(cfg, m, logger)
        logger.info(f"----- MODEL {m} 结束 -----")

    # 写 done
    (cafe_root / ".cafe.done").write_text(datetime.datetime.now().isoformat()+"\n", encoding="utf-8")
    logger.info("[OK] 全部完成，已写 .cafe.done")

if __name__ == "__main__":
    main()