#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
01_db.py
================================================================================
[功能]
下载/校验/解压 UniProt Swiss-Prot reviewed FASTA，并构建 DIAMOND 数据库。

[皇上要求]
- 不接受命令行参数：参数集中在脚本顶部 + config.yaml
- 使用相对路径（以项目根目录为基准）
- 屏幕流式输出：下载需要真实可见的进度/summary
- 脚本可单独运行，也可被 Snakefile 调用
- 只下载必须文件：reldate.txt, README, RELEASE.metalink, uniprot_sprot.fasta.gz
- 默认走 HTTPS 直链；若 config.yaml 开关打开，则启用 aria2c metalink + select-file

[目录约定]
- 原始下载：data/db/raw/
- DIAMOND库：data/db/diamond/uniprot_sprot_{REL}.dmnd
- 留档：results/01_db/{REL}/
- 当前版本指针：results/01_db/CURRENT_REL.txt
- 日志：logs/

[依赖]
- diamond（必须）
- aria2c（推荐；无则回退 wget；再无则用 Python 下载）
================================================================================
"""

# ----------------------------- 参数区（皇上可改） -----------------------------

UNIPROT_BASE = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete"

URL_RELDATE  = f"{UNIPROT_BASE}/reldate.txt"
URL_README   = f"{UNIPROT_BASE}/README"
URL_METALINK = f"{UNIPROT_BASE}/RELEASE.metalink"
URL_SPROT_GZ = f"{UNIPROT_BASE}/uniprot_sprot.fasta.gz"

ARIA2C_CONN = 16
ARIA2C_SPLIT = 16
ARIA2C_CHUNK = "1M"
ARIA2C_CONTINUE = True

USE_WGET_FALLBACK = True

SPROT_FASTA_NAME = "uniprot_sprot.fasta"

# ----------------------------- import（含单跑修复） -----------------------------

import os
import re
import sys
import gzip
import shutil
import subprocess
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Optional, List, Tuple, Any, Dict

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts._common import (
    find_project_root,
    read_yaml_config,
    ensure_dir,
    md5sum_file,
)

# ----------------------------- 工具函数 -----------------------------

def which(cmd: str) -> Optional[str]:
    return shutil.which(cmd)

def stream_command(cmd: List[str], log_path: Path, cwd: Optional[Path] = None) -> int:
    """
    运行外部命令：stdout/stderr 合流，按“行”流式打印并写日志。
    适合以换行输出为主的程序（例如 diamond）。
    """
    ensure_dir(log_path.parent)
    with log_path.open("a", encoding="utf-8") as log:
        log.write("CMD: " + " ".join(cmd) + "\n")
        log.flush()

        p = subprocess.Popen(
            cmd,
            cwd=(cwd.as_posix() if cwd else None),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True,
        )
        assert p.stdout is not None
        for line in p.stdout:
            print(line.rstrip("\n"), flush=True)
            log.write(line)
            log.flush()
        return p.wait()

def download_small_file_python(url: str, out_path: Path) -> None:
    """
    用 Python 把 URL 当普通文件下载（不触发 aria2c Metalink 展开）。
    """
    ensure_dir(out_path.parent)
    if out_path.is_file() and out_path.stat().st_size > 0:
        print(f"[01_db] exists, skip download: {out_path}", flush=True)
        return

    print(f"[01_db] download (python): {url} -> {out_path}", flush=True)
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    if tmp.exists():
        tmp.unlink()

    with urllib.request.urlopen(url) as resp, tmp.open("wb") as f:
        shutil.copyfileobj(resp, f)

    tmp.replace(out_path)

    if not out_path.is_file() or out_path.stat().st_size == 0:
        raise RuntimeError(f"Downloaded file missing/empty: {out_path}")

def _cfg_get_db(cfg: Dict[str, Any]) -> Dict[str, Any]:
    d = cfg.get("db", {})
    if isinstance(d, dict):
        return d
    return {}

def download_big_sprot(
    url_https: str,
    metalink_path: Path,
    out_path: Path,
    log_path: Path,
    cfg_db: Dict[str, Any]
) -> None:
    """
    下载 uniprot_sprot.fasta.gz：

    默认（更稳）：HTTPS 直链
      aria2c ... -o uniprot_sprot.fasta.gz https://...

    可选（config 开关打开才启用）：metalink + select-file
      aria2c -M RELEASE.metalink --select-file=N ...

    方案 A：aria2c 不走 PIPE，继承终端输出（你能看到进度/summary）。
    """
    ensure_dir(out_path.parent)

    if out_path.is_file() and out_path.stat().st_size > 0:
        print(f"[01_db] exists, skip download: {out_path}", flush=True)
        return

    aria2 = which("aria2c")
    wget = which("wget") if USE_WGET_FALLBACK else None

    use_metalink = bool(cfg_db.get("use_metalink", False))
    select_file = int(cfg_db.get("metalink_select_file", 6))
    preferred_proto = str(cfg_db.get("metalink_preferred_protocol", "https")).strip().lower()
    if preferred_proto not in {"http", "https", "ftp"}:
        preferred_proto = "https"

    if aria2:
        ensure_dir(log_path.parent)

        base = [
            "aria2c",
            "--console-log-level=notice",
            "--log-level=notice",
            "-l", log_path.as_posix(),
            "-x", str(ARIA2C_CONN),
            "-s", str(ARIA2C_SPLIT),
            "-k", str(ARIA2C_CHUNK),
        ]
        if ARIA2C_CONTINUE:
            base.append("-c")

        if use_metalink:
            cmd = base + [
                "--metalink-preferred-protocol=" + preferred_proto,
                "-M", metalink_path.as_posix(),
                f"--select-file={select_file}",
            ]
            print(f"[01_db] download (aria2c metalink): select-file={select_file} preferred={preferred_proto}", flush=True)
            rc = subprocess.call(cmd, cwd=out_path.parent.as_posix())
            if rc != 0:
                raise RuntimeError(f"aria2c metalink download failed (exit={rc}). See {log_path}")

        else:
            cmd = base + ["-o", out_path.name, url_https]
            print(f"[01_db] download (aria2c https): {url_https}", flush=True)
            rc = subprocess.call(cmd, cwd=out_path.parent.as_posix())
            if rc != 0:
                raise RuntimeError(f"aria2c https download failed (exit={rc}). See {log_path}")

    elif wget:
        cmd = ["wget", "-O", out_path.as_posix(), url_https]
        print(f"[01_db] download (wget): {url_https}", flush=True)
        rc = stream_command(cmd, log_path)
        if rc != 0:
            raise RuntimeError(f"wget download failed (exit={rc}) for {url_https}")

    else:
        print("[01_db] aria2c/wget not found, fallback to python download (may be slow)", flush=True)
        download_small_file_python(url_https, out_path)

    if not out_path.is_file() or out_path.stat().st_size == 0:
        raise RuntimeError(f"Downloaded file missing/empty: {out_path}")

def extract_rel_and_date(reldate_path: Path) -> Tuple[str, str]:
    """
    从 reldate.txt 提取 REL 与 Swiss-Prot 日期。
    兼容格式：
      - UniProt Knowledgebase Release 2025_04 consists of:
      - UniProtKB/Swiss-Prot Release 2025_04 of 08-Oct-2025
      - 旧格式：Release 2025_04 ...
    """
    rel = ""
    sp_date = ""
    rel_from_sprot = ""

    lines = reldate_path.read_text(encoding="utf-8", errors="replace").splitlines()

    for line in lines:
        s = line.strip()

        if not rel:
            m = re.search(r"\bRelease\s+([0-9_]+)\b", s)
            if m:
                rel = m.group(1)

        if not rel_from_sprot:
            m = re.search(r"Swiss-Prot\s+Release\s+([0-9_]+)\b", s)
            if m:
                rel_from_sprot = m.group(1)

        if not sp_date:
            m = re.search(r"Swiss-Prot\s+Release\s+[0-9_]+\s+of\s+([0-9A-Za-z-]+)\b", s)
            if m:
                sp_date = m.group(1)

        if rel and sp_date:
            break

    if not rel:
        rel = rel_from_sprot

    if not rel:
        raise ValueError("Failed to parse REL from reldate.txt")

    if not sp_date:
        sp_date = "NA"

    return rel, sp_date

def parse_md5_from_metalink(metalink_path: Path, filename: str) -> str:
    """从 RELEASE.metalink 中解析指定文件的 md5，只做解析不做下载。"""
    tree = ET.parse(metalink_path.as_posix())
    root = tree.getroot()

    def local(tag: str) -> str:
        return tag.split("}", 1)[1] if "}" in tag else tag

    md5 = ""
    for node in root.iter():
        if local(node.tag) != "file":
            continue
        if node.attrib.get("name") != filename:
            continue
        for child in node.iter():
            if local(child.tag) == "hash":
                t = child.attrib.get("type", "").lower()
                if t == "md5":
                    md5 = (child.text or "").strip()
                    break
        if md5:
            break

    if not md5:
        raise ValueError(f"MD5 for {filename} not found in {metalink_path}")
    return md5

def md5sum_check(path: Path, expected_md5: str) -> None:
    got = md5sum_file(path)
    if got.lower() != expected_md5.lower():
        raise RuntimeError(f"MD5 mismatch for {path.name}: expected {expected_md5}, got {got}")
    print(f"[01_db] MD5 OK: {path.name}", flush=True)

def gunzip_to(src_gz: Path, dst: Path) -> None:
    if dst.is_file() and dst.stat().st_size > 0:
        print(f"[01_db] exists, skip gunzip: {dst}", flush=True)
        return
    ensure_dir(dst.parent)
    print(f"[01_db] gunzip: {src_gz} -> {dst}", flush=True)
    with gzip.open(src_gz.as_posix(), "rb") as fin, dst.open("wb") as fout:
        shutil.copyfileobj(fin, fout, length=1024 * 1024)

# ----------------------------- 主流程 -----------------------------

def main() -> None:
    root = find_project_root()
    os.chdir(root)

    cfg = read_yaml_config(root)
    cfg_db = _cfg_get_db(cfg)

    raw_dir = root / "data" / "db" / "raw"
    diamond_dir = root / "data" / "db" / "diamond"
    ensure_dir(raw_dir)
    ensure_dir(diamond_dir)
    ensure_dir(root / "results" / "01_db")
    ensure_dir(root / "logs")

    log_download_big = root / "logs" / "01_db.download_sprot.log"
    log_makedb = root / "logs" / "01_db.diamond_makedb.log"

    reldate = raw_dir / "reldate.txt"
    readme = raw_dir / "README"
    metalink = raw_dir / "RELEASE.metalink"
    sprot_gz = raw_dir / "uniprot_sprot.fasta.gz"
    sprot_fa = raw_dir / SPROT_FASTA_NAME

    print("[01_db] step: download metadata/files", flush=True)

    download_small_file_python(URL_RELDATE, reldate)
    download_small_file_python(URL_README, readme)
    download_small_file_python(URL_METALINK, metalink)

    download_big_sprot(
        url_https=URL_SPROT_GZ,
        metalink_path=metalink,
        out_path=sprot_gz,
        log_path=log_download_big,
        cfg_db=cfg_db
    )

    rel, sp_date = extract_rel_and_date(reldate)
    out = root / "results" / "01_db" / rel
    ensure_dir(out)

    print(f"[01_db] REL={rel} SwissProt_date={sp_date}", flush=True)

    print("[01_db] step: md5 check via RELEASE.metalink", flush=True)
    md5_expected = parse_md5_from_metalink(metalink, "uniprot_sprot.fasta.gz")
    md5sum_check(sprot_gz, md5_expected)

    print("[01_db] step: gunzip for diamond makedb", flush=True)
    gunzip_to(sprot_gz, sprot_fa)

    (out / "db_manifest.tsv").write_text(
        f"REL\t{rel}\n"
        f"SwissProt_date\t{sp_date}\n"
        f"DB_URL\t{URL_SPROT_GZ}\n",
        encoding="utf-8"
    )
    shutil.copy2(reldate.as_posix(), (out / "reldate.txt").as_posix())
    shutil.copy2(metalink.as_posix(), (out / "RELEASE.metalink").as_posix())
    shutil.copy2(readme.as_posix(), (out / "README").as_posix())

    (out / "md5_uniprot_sprot.fasta.gz.txt").write_text(
        f"{md5sum_file(sprot_gz)}  {sprot_gz.as_posix()}\n",
        encoding="utf-8"
    )

    dv = subprocess.run(["diamond", "version"], capture_output=True, text=True)
    (out / "diamond_version.txt").write_text((dv.stdout + dv.stderr).strip() + "\n", encoding="utf-8")

    db_prefix = diamond_dir / f"uniprot_sprot_{rel}"
    db_dmnd = db_prefix.with_suffix(".dmnd")

    if db_dmnd.is_file() and db_dmnd.stat().st_size > 0:
        print(f"[01_db] exists, skip makedb: {db_dmnd}", flush=True)
    else:
        print(f"[01_db] step: diamond makedb -> {db_dmnd}", flush=True)
        cmd = ["diamond", "makedb", "--in", sprot_fa.as_posix(), "--db", db_prefix.as_posix()]
        rc = stream_command(cmd, log_makedb)
        if rc != 0:
            raise RuntimeError(f"diamond makedb failed (exit={rc}). See {log_makedb}")
        if not db_dmnd.is_file() or db_dmnd.stat().st_size == 0:
            raise RuntimeError(f"DIAMOND db missing/empty after makedb: {db_dmnd}")

    (out / "diamond_db_fileinfo.txt").write_text(
        f"{db_dmnd.as_posix()}\t{db_dmnd.stat().st_size}\n",
        encoding="utf-8"
    )

    (root / "results" / "01_db" / "CURRENT_REL.txt").write_text(rel + "\n", encoding="utf-8")

    print(f"[01_db] DONE REL={rel} DB={db_dmnd}", flush=True)

if __name__ == "__main__":
    main()

