#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
提取并清洗 CDS + 蛋白质脚本（CDS-driven + 物种级并行 + 越界坐标过滤）

更新要点（最小必要改动）：
1) 并行：将“按物种处理”的主循环改为多进程池，物种之间并行；每个物种独占自己的 tmp 子目录。
2) 坐标越界修复：在 gffread 提取前，利用基因组序列长度过滤 GFF 中越界 feature（end>len 或 start<1），
   避免 gffread 的 `GffObj::getSpliced() improper genomic coordinate` 报错。
3) 仍走“CDS→翻译蛋白”的路线：最终蛋白由通过 QC 的 CDS 翻译得到，保证一一对应、可审计。
4) 其余流程、阈值、输出路径、清理策略保持不变。
"""

import sys, re, shutil, gzip, subprocess, datetime, csv, os
from pathlib import Path
from typing import List, Optional, Dict, Tuple
from multiprocessing import Pool, cpu_count

# ========================== 路径与阈值（可按需修改） ==========================
ANNOT_DIR = "data/annotation"          # GFF/GFF3 目录
GENOME_DIR = "data/genomic"            # 基因组 FASTA 目录
OUT_PEP_DIR = "data/proteomes"         # 输出蛋白目录
OUT_CDS_DIR = "data/cds"               # 输出CDS目录
TMP_DIR  = "data/.tmp_from_gff"        # 中间文件目录（每物种一个子目录）
LIST_DIR = "results/proteomes"         # 日志/列表目录
SPECIES_LIST = f"{LIST_DIR}/species.list"
MATCH_TSV = f"{LIST_DIR}/annot_genome_match.tsv"

SEQID_COVERAGE_MIN = 0.90              # GFF与FA序列名重合率阈值
RESUME_SKIP_PASS = True                # 已通过且有输出则跳过

# -------- 新增总开关：结束时是否删除 TMP 目录（默认 False=不删除） --------
CLEANUP_TMP_AT_END = False

# ========================== 依赖工具（AGAT现代脚本名） ==========================
AGAT_KEEP_LONGEST = "agat_sp_keep_longest_isoform.pl"
AGAT_FIX_DUPS     = "agat_sp_fix_features_locations_duplicated.pl"
AGAT_FIX_PHASES   = "agat_sp_fix_cds_phases.pl"   # 需要 -g -f
GFFREAD_BIN       = "gffread"
AGAT_EXTRA_ARGS   = ""                  # 预留额外参数（默认空）

# ========================== ID 规范化与QC参数 ==========================
STRIP_ID_PREFIXES = True
ID_PREFIX_REGEX   = r'^(rna\-|mrna\-|transcript:|cds:)'

PROT_MIN_LEN_AA   = 50                 # 蛋白最短长度（aa）
PROT_MAX_LEN_AA   = 5000               # 蛋白最长长度（aa）；设为 0 表示不启用
CDS_MIN_LEN_NT    = 150                # CDS 最短长度（nt，且需3的倍数）

# ========================== 并行参数 ==========================
# 默认使用 CPU 一半核心，至少 2；可用环境变量 N_WORKERS 覆盖
N_WORKERS = max(2, int(os.environ.get("N_WORKERS", max(1, (cpu_count() or 4)//2))))

# 限制一些底层库的线程数，避免过度抢核（对外部工具通常也安全）
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

# ========================== 扩展名与工具检测 ==========================
GFF_EXTS = [".gff.gz", ".gff3.gz", ".gff", ".gff3"]
FA_EXTS  = [".fa.gz", ".fna.gz", ".fasta.gz", ".fa", ".fna", ".fasta"]
TOOLS = ["awk", "gzip", "seqkit", GFFREAD_BIN, AGAT_KEEP_LONGEST, AGAT_FIX_DUPS, AGAT_FIX_PHASES]

MATCH_HEADER = ["species","gff_file","genome_file","cov_frac","status","pep_kept","cds_kept","note","timestamp"]

# ========================== 【新增】pep↔cds 映射输出参数（可关） ==========================
WRITE_PEP2CDS = True                                   # 是否写出蛋白-转录本映射
PEP2CDS_PER_SPECIES_DIR = "data/maps/pep2cds_by_species"  # 每物种分片
PEP2CDS_TSV = "data/maps/pep2cds.tsv"                  # 合并后的总表

# ========================== 运行辅助函数 ==========================
def need(cmd: str):
    """检查依赖命令是否存在；不存在则抛错"""
    if subprocess.call(f"command -v {cmd} >/dev/null 2>&1", shell=True) != 0:
        raise RuntimeError(f"[ERR] dependency not found: {cmd}")

def ensure_dir(p: Path): p.mkdir(parents=True, exist_ok=True)
def is_gz(path: Path) -> bool: return str(path).endswith(".gz")

def run(cmd: str, label: str = "CMD"):
    """执行 shell 命令；透传输出；失败则抛出异常（便于并行下单物种失败不影响全局）"""
    print(f"[{label}] {cmd}")
    res = subprocess.run(cmd, shell=True)
    if res.returncode != 0:
        raise RuntimeError(f"[ERR] command failed ({label}), exit code: {res.returncode}")

def bash_out(cmd: str) -> str:
    """获取命令stdout文本（stderr静默）"""
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    return p.stdout or ""

def _cat(path: Path) -> str:
    """根据是否压缩选择 cat 或 gzip -cd"""
    return f"gzip -cd '{path}'" if is_gz(path) else f"cat '{path}'"

def safe_int_last_line(s: str) -> int:
    s = (s or "").strip()
    if not s: return 0
    try:
        return int(s.splitlines()[-1].strip())
    except Exception:
        digits = re.findall(r"\d+", s)
        return int(digits[-1]) if digits else 0

def count_fasta_records(p: Path) -> int:
    return safe_int_last_line(bash_out(f"test -f '{p}' && grep -c '^>' '{p}' 2>/dev/null || echo 0"))

# ========================== 物种与路径 ==========================
def list_species_from_annotation(annot_dir: Path) -> List[str]:
    """从 annotation/ 中列出物种stem（去除各种中间后缀）"""
    stems_set = set()
    for f in sorted(annot_dir.iterdir()):
        if not f.is_file(): continue
        m = re.search(r'^(?P<stem>.+?)\.(gff3?|gff)(\.gz)?$', f.name)
        if not m: continue
        stem = m.group("stem")
        stem = re.sub(r'\.(agat_longest|keep_longest|agat.*|clean4agat|clean(?:ed)?|filtered.*|fixloc|fixphase)$', '', stem)
        stems_set.add(stem)
    return sorted(stems_set)

def find_by_stem(stem: str, directory: Path, exts: List[str]) -> Optional[Path]:
    """按 stem + 扩展名在目录中查找文件"""
    for e in exts:
        f = directory / f"{stem}{e}"
        if f.is_file(): return f
    return None

# ========================== 覆盖度检查 ==========================
def coverage_details(gff: Path, fa: Path) -> Tuple[float, int, int, int]:
    """
    计算 GFF 的 seqid 集合与 FASTA 的序列名集合的重合比例
    返回：覆盖率、GFF序列数、FASTA序列数、重叠数
    """
    env_prefix = "LC_ALL=C"
    gff_ids = bash_out(f"{env_prefix} " + _cat(gff) + r" | awk -F'\t' '$0!~/^#/ && NF>=1 {print $1}' | sort -u").splitlines()
    fa_ids  = bash_out(f"{env_prefix} " + _cat(fa)  + r" | grep '^>' | sed 's/^>//;s/ .*//' | sort -u").splitlines()
    if not gff_ids or not fa_ids:
        return 0.0, len(set(gff_ids)), len(set(fa_ids)), 0
    s_g = set(gff_ids); s_f = set(fa_ids)
    overlap = len(s_g & s_f)
    cov = overlap / len(s_g) if s_g else 0.0
    return cov, len(s_g), len(s_f), overlap

# ========================== 缓存读写（match表） ==========================
def load_match_table(path: Path) -> Dict[str, Dict[str, str]]:
    if not path.is_file(): return {}
    data: Dict[str, Dict[str, str]] = {}
    with path.open("r", encoding="utf-8") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for row in rdr:
            sp = row.get("species","").strip()
            if sp: data[sp] = row
    return data

def write_match_table(path: Path, rows: List[Dict[str, str]]):
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=MATCH_HEADER, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows: w.writerow(r)

# ========================== ID规范化与QC ==========================
def norm_ids_file(p: Path):
    """对一个只含 ID 的文件进行前缀清理（比如 rna-/mrna-/transcript: 等）"""
    if STRIP_ID_PREFIXES:
        run(f"awk '{{gsub(/^{ID_PREFIX_REGEX}/,\"\"); print}}' '{p}' > '{p}.norm' && mv '{p}.norm' '{p}'", "NORM")

def clean_cds_codon_safe(fna_in: Path, fna_out: Path, keep_ids_out: Path) -> int:
    """
    对 CDS 做“密码子安全性”筛选（CDS-driven关键步骤）：
    - 长度 ≥ CDS_MIN_LEN_NT 且 长度%3==0（直接对序列计算 length）
    - 翻译后不允许内部 '*'（允许末端 1 个 '*'）
    - 输出保留 ID 列表 keep.ids，并按该列表应用到 CDS（得到最终 codon-safe CDS）
    注：fx2tab 使用 -i -s（id, seq），避免列位次歧义
    """
    tmp_len = Path(str(fna_out) + ".len3.tmp")
    # 使用 -i -s：$1=id, $2=seq
    run(
        f"seqkit fx2tab -i -s '{fna_in}' | "
        f"awk '{{ if(length($2)>={CDS_MIN_LEN_NT} && length($2)%3==0) print \">\"$1\"\\n\"$2 }}' | "
        f"seqkit seq -w 0 > '{tmp_len}'",
        "QC:CDS_LEN_MOD3"
    )
    pep_tmp = Path(str(fna_out) + ".pep.tmp")
    run(f"seqkit translate -w 0 '{tmp_len}' > '{pep_tmp}'", "QC:TRANSLATE")

    # 统计内部终止：允许尾部1个'*'，内部不允许
    awk_stop = r"""{
      id=$1; pep=$2; t=pep; gsub(/\*/,"",t); sc=length(pep)-length(t);
      if (sc==0 || (sc==1 && substr(pep,length(pep),1)=="*")) print id;
    }"""
    run(f"seqkit fx2tab -i -s '{pep_tmp}' | awk '{awk_stop}' > '{keep_ids_out}'", "QC:KEEP_IDS")

    # 规范化 keep.ids（一次即可）
    norm_ids_file(keep_ids_out)

    # 应用 keep.ids 到 CDS（得到最终 CDS）
    n_ids = safe_int_last_line(bash_out(f"test -s '{keep_ids_out}' && wc -l < '{keep_ids_out}' || echo 0"))
    if n_ids > 0:
        run(f"seqkit grep -f '{keep_ids_out}' '{tmp_len}' | seqkit seq -w 0 > '{fna_out}'", "QC:APPLY_IDS")
    else:
        Path(fna_out).write_text("", encoding="utf-8")

    # 清理临时文件
    for p in [tmp_len, pep_tmp]:
        try: Path(p).unlink()
        except Exception: pass

    return count_fasta_records(fna_out)

def translate_cds_to_protein(final_cds: Path, out_pep: Path) -> int:
    """
    由“最终 CDS”直接翻译得到蛋白（CDS-driven 的产物）并做轻度 QC：
    - 去掉唯一的末端 '*'
    - 不允许内部 '*'
    - 非常见氨基酸 B/Z/J/U/O → 'X'
    - 过滤掉长度 < PROT_MIN_LEN_AA
    """
    pep_raw = Path(str(out_pep) + ".raw.tmp")
    run(f"seqkit translate -w 0 '{final_cds}' > '{pep_raw}'", "PROT:TRANSLATE_FROM_CDS")

    awk_script = r"""{
      id=$1; seq=toupper($2);
      gsub(/B|Z|J|U|O/,"X",seq);                # 非常见氨基酸 → X
      t=seq; gsub(/\*$/,"",t); seq=t;           # 去掉末端唯一 '*'
      t2=seq; gsub(/\*/,"",t2); sc=length(seq)-length(t2);
      if (sc>0) next;                           # 内部 '*' → 丢弃（理论上不会，因为前面已筛）
      if (length(seq)>0 && length(seq)<%MINLEN%) next;
      if (%MAXLEN%>0 && length(seq)>%MAXLEN%) next;
      print ">" id "\n" seq
    }""".replace("%MINLEN%", str(PROT_MIN_LEN_AA)).replace("%MAXLEN%", str(PROT_MAX_LEN_AA))

    pep_qc = Path(str(out_pep) + ".qc.tmp")
    run(f"seqkit fx2tab -i -s '{pep_raw}' | awk '{awk_script}' | seqkit seq -w 0 > '{pep_qc}'", "PROT:POST_QC")

    npep = count_fasta_records(pep_qc)
    if npep > 0:
        shutil.move(pep_qc, out_pep)
    else:
        Path(out_pep).write_text("", encoding="utf-8")
        try: pep_qc.unlink()
        except Exception: pass

    try: pep_raw.unlink()
    except Exception: pass

    return npep

# ========================== 越界坐标过滤（关键修复） ==========================
def filter_gff_by_contig_length(gff_in: Path, fa_plain: Path, gff_out: Path, report_out: Path) -> Tuple[int,int]:
    """
    生成基因组 contig 长度表，然后过滤 GFF 中任何超出长度范围的 feature 行：
    - 允许注释行原样通过
    - 仅对具有 9 列的 feature 行检查 $4(start)、$5(end)
    - 条件：1 <= start <= end <= len[seqid]
    - 被剔除的行计数写入报告
    返回 (kept, removed)
    """
    ensure_dir(gff_out.parent)
    contig_len = gff_out.parent / "contig.len.tsv"
    run(f"seqkit fx2tab -n -l '{fa_plain}' > '{contig_len}'", "LEN:CONTIG_TABLE")

    # 过滤并统计
    awk_filter = r"""BEGIN{FS=OFS="\t"}
NR==FNR {len[$1]=$2; next}
$0 ~ /^#/ {print; next}
NF<9 {next}
{
  seq=$1; s=$4; e=$5;
  if (seq in len && s>=1 && e>=s && e<=len[seq]) {print; kept++} else {removed++}
}
END{
  if (ENVIRON["__GFF_LEN_REPORT"]!="") {
    print "kept=" kept "\tremoved=" removed > ENVIRON["__GFF_LEN_REPORT"];
  }
}"""
    # 设置环境变量把统计写到报告文件
    env = os.environ.copy()
    env["__GFF_LEN_REPORT"] = str(report_out)
    print("[CLEAN:LENCHECK] filtering out-of-bound features by contig lengths ...")
    res = subprocess.run(
        f"awk '{awk_filter}' '{contig_len}' '{gff_in}' > '{gff_out}'",
        shell=True, env=env
    )
    if res.returncode != 0:
        raise RuntimeError("[ERR] GFF length-based filtering failed")

    # 读取报告数字
    kept, removed = 0, 0
    if report_out.is_file():
        txt = report_out.read_text().strip()
        m1 = re.search(r'kept=(\d+)', txt); m2 = re.search(r'removed=(\d+)', txt)
        kept = int(m1.group(1)) if m1 else 0
        removed = int(m2.group(1)) if m2 else 0
    return kept, removed

# 统计 GFF 中 CDS 行数（用于过滤前后对比 / 回退判断）
def _count_cds_lines(p: Path) -> int:
    return safe_int_last_line(bash_out(f"awk -F'\t' '$3==\"CDS\"{{c++}} END{{print c+0}}' '{p}'"))

# ========================== 【新增】写 pep↔cds 映射（每物种分片） ==========================
def _iter_fasta_ids(path: Path):
    """遍历 FASTA 的 ID（只取第一个空格前的 token）"""
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('>'):
                yield line[1:].strip().split()[0]

def _split_species_id(token: str, sp_fallback: str) -> Tuple[str, str]:
    """若 ID 形如 Species|GeneID，则拆出物种与基因；否则物种回退为当前物种"""
    if '|' in token:
        sp, pid = token.split('|', 1)
        return sp.strip(), pid.strip()
    return sp_fallback, token

# ========================== 【新增，仅新增】统一 ID 的“合同化清洗”函数 ==========================
def _sanitize_id_contract(x: str) -> str:
    """
    与下游契约一致的 ID 清洗规则：
    - 保留竖线 '|'（例如 gnl|WGS_AMQN|CAPTEDRAFT... 依然合法）
    - 仅允许 [A-Za-z0-9._:-|]，其它字符替换为 '_'
    - 去首尾空白
    """
    return re.sub(r"[^A-Za-z0-9._:\-\|]", "_", x.strip())

def write_pep2cds_for_species(sp: str, cds_fa: Path, pep_fa: Path, out_dir: Path) -> int:
    """
    根据“最终 CDS 与最终蛋白”顺序一一对应关系写出映射分片：
    species  protein_id  cds_id
    """
    ids_cds = list(_iter_fasta_ids(cds_fa))
    ids_pep = list(_iter_fasta_ids(pep_fa))

    # 建立基于 ID 的安全映射，而不是按行号硬配：
    # 思路：蛋白是从 CDS 翻译而来，ID 理应一致或仅差物种前缀；
    #      有些 CDS 可能被翻译后太短而被过滤掉，因此蛋白列表是CDS列表的子集；
    #      我们按“ID是否匹配”逐条输出，避免顺序错配。
    cds_clean_to_raw: Dict[str, str] = {}
    for raw_cds_id in ids_cds:
        sp_c, cds_clean = _split_species_id(raw_cds_id, sp)
        if cds_clean not in cds_clean_to_raw:
            cds_clean_to_raw[cds_clean] = raw_cds_id

    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / f"{sp}.tsv"
    newfile = not out_tsv.exists() or out_tsv.stat().st_size == 0

    written = 0
    with open(out_tsv, 'a', encoding='utf-8') as w:
        if newfile:
            w.write('species\tprotein_id\tcds_id\n')
        for pep_raw_id in ids_pep:
            sp_p, pep_clean = _split_species_id(pep_raw_id, sp)
            # 只在我们能找到同名CDS时写入，避免错配
            if pep_clean in cds_clean_to_raw:
                cds_raw_id = cds_clean_to_raw[pep_clean]
                sp_c2, cds_clean2 = _split_species_id(cds_raw_id, sp)

                # ----------【仅新增两行】在写出前进行“合同化清洗”----------
                pep_clean = _sanitize_id_contract(pep_clean)
                cds_clean2 = _sanitize_id_contract(cds_clean2)
                # -------------------------------------------------------

                w.write(f"{sp}\t{pep_clean}\t{cds_clean2}\n")
                written += 1

    return written

# ========================== 单物种处理（供并行调用） ==========================
def process_one_species(args) -> Dict[str, str]:
    """
    独立处理一个物种：AGAT修复 → 过滤 → 提取CDS → QC → 翻译蛋白
    返回写 match 表所需的行（字典），失败则抛异常（主进程捕获并汇总）
    """
    sp, gff_path, fa_path, cache_row = args
    sp_tmp = Path(TMP_DIR)/sp
    ensure_dir(sp_tmp)

    gff_plain = sp_tmp / "annot.gff3"
    fa_plain  = sp_tmp / "genome.fna"

    # 解压/复制到工作区
    if str(gff_path).endswith(".gz"):
        run(f"gzip -cd '{gff_path}' > '{gff_plain}'", f"PREP:GUNZIP_GFF ({sp})")
    else:
        run(f"cp -f '{gff_path}' '{gff_plain}'", f"PREP:CP_GFF ({sp})")

    if str(fa_path).endswith(".gz"):
        run(f"gzip -cd '{fa_path}' > '{fa_plain}'", f"PREP:GUNZIP_FA ({sp})")
    else:
        run(f"cp -f '{fa_path}' '{fa_plain}'", f"PREP:CP_FA ({sp})")

    # 0) 修复重复位置
    fixed1 = sp_tmp / "fixloc.gff3"
    run(f"{AGAT_FIX_DUPS} -g '{gff_plain}' -o '{fixed1}'", "AGAT:FIX_DUPS (1/3)")

    # 1) 修复 phase（需 genome）
    fixed2 = sp_tmp / "fixphase.gff3"
    run(f"{AGAT_FIX_PHASES} -g '{fixed1}' -f '{fa_plain}' -o '{fixed2}'", "AGAT:FIX_PHASE (2/3)")
    if (not fixed2.exists()) or fixed2.stat().st_size < 1024:
        print(f"[WARN] {sp}: fixphase 文件缺失或过小，回退使用 fixloc 结果。")
        fixed2 = fixed1

    # 2) 仅保留必要 feature —— 恢复 gene
    gff_clean = sp_tmp / "clean4agat.gff3"
    run(
        "awk -F'\\t' 'BEGIN{OFS=\"\\t\"} "
        "/^#/ {print; next} "
        "{t=$3; if (t==\"gene\"||t==\"mRNA\"||t==\"transcript\"||t==\"exon\"||t==\"CDS\"||"
        "t==\"five_prime_UTR\"||t==\"three_prime_UTR\"||t==\"start_codon\"||t==\"stop_codon\") print;}' "
        f"'{fixed2}' > '{gff_clean}'",
        "CLEAN:GFF_FILTER"
    )

    # 3) 每基因保留最长转录本
    gff_agat = sp_tmp / "agat_longest.gff3"
    run(f"{AGAT_KEEP_LONGEST} -g '{gff_clean}' -o '{gff_agat}' {AGAT_EXTRA_ARGS}", "AGAT:KEEP_LONGEST (3/3)")

    # --- 可选统计：keep_longest 前后 CDS 行数（便于排错）
    cds_before = _count_cds_lines(gff_clean)
    cds_after  = _count_cds_lines(gff_agat)
    print(f"[STAT] {sp}: CDS lines before/after keep_longest = {cds_before}/{cds_after}")

    # ★ 3.5) 基于 contig 长度的越界过滤（关键修复）
    gff_lenf = sp_tmp / "lenfilter.gff3"
    report   = sp_tmp / "lenfilter.report.txt"
    kept, removed = filter_gff_by_contig_length(gff_agat, fa_plain, gff_lenf, report)
    if removed > 0:
        print(f"[CLEAN:LENCHECK] {sp}: removed {removed} out-of-bound features; kept {kept}.")

    # —— 新增：若长度过滤后 CDS 为 0，则回退（避免 seqid 不匹配“团灭”）
    cds_after_lenf = _count_cds_lines(gff_lenf)
    print(f"[STAT] {sp}: CDS lines after len-filter = {cds_after_lenf}")
    if cds_after_lenf == 0:
        print(f"[WARN] {sp}: len-filter produced 0 CDS; fallback to keep_longest output.")
        gff_lenf = gff_agat

    # 4) 仅提取 CDS（不再提蛋白） —— 去掉 -V，保留 -E
    all_cds = sp_tmp / "all.cds.fna"
    run(f"{GFFREAD_BIN} -E '{gff_lenf}' -g '{fa_plain}' -x '{all_cds}'", "GFFREAD:EXTRACT_CDS")

    # 4.5) 断点统计：gffread 提取条数
    n_all_cds = count_fasta_records(all_cds)
    print(f"[STAT] {sp}: gffread CDS extracted = {n_all_cds}")
    if n_all_cds == 0:
        raise RuntimeError(f"[FATAL] {sp}: gffread produced 0 CDS (check gene/mRNA hierarchy).")

    # 5) 规范ID并写中间CDS（再做QC）
    out_cds = Path(OUT_CDS_DIR)  / f"{sp}.fna"
    ensure_dir(out_cds.parent)
    if STRIP_ID_PREFIXES:
        run(
            f"seqkit replace -p '^\\s*(\\S+).*' -r '$1' '{all_cds}' | "
            f"seqkit replace -p '{ID_PREFIX_REGEX}' -r '' | seqkit seq -w 0 > '{out_cds}'",
            "CDS:NORM_ID"
        )
    else:
        run(f"seqkit replace -p '^\\s*(\\S+).*' -r '$1' '{all_cds}' | seqkit seq -w 0 > '{out_cds}'", "CDS:NORM_ID")

    # 6) 密码子安全性 QC（keep.ids 产生于此，得到最终 CDS）
    keep_ids = sp_tmp / "keep.ids"
    n_cds = clean_cds_codon_safe(out_cds, out_cds, keep_ids)
    if n_cds == 0:
        raise RuntimeError(f"[FATAL] {sp}: codon-safe CDS count is 0.")

    # 7) 由“最终 CDS”直接翻译生成“最终蛋白”（CDS-driven）
    out_pep = Path(OUT_PEP_DIR) / f"{sp}.faa"
    ensure_dir(out_pep.parent)
    npep = translate_cds_to_protein(out_cds, out_pep)
    if npep == 0:
        raise RuntimeError(f"[FATAL] {sp}: translated protein count is 0 after QC.")

    # 【新增】写出 pep↔cds 映射分片（不影响既有输出）
    if WRITE_PEP2CDS:
        try:
            write_pep2cds_for_species(sp, out_cds, out_pep, Path(PEP2CDS_PER_SPECIES_DIR))
            print(f"[MAP] {sp}: pep2cds written.")
        except Exception as e:
            print(f"[WARN] {sp}: write pep2cds failed: {e}")

    # 8) 清理物种临时文件（不触碰 annotation/）
    for patt in ["*.fna","*.fa","*.fasta","*.cds.fna","*.protein.faa","*.gtf","*.tsv",
                 "keep*.ids","*.pep.tmp","*.len3.tmp","*.tmp","tmp.translate.faa",
                 "*.agat_longest.gff3","*.clean4agat.gff3","*.fixloc.gff3","*.fixphase.gff3",
                 "contig.len.tsv","lenfilter.report.txt","lenfilter.gff3"]:
        for f in sp_tmp.glob(patt):
            try: f.unlink()
            except Exception: pass

    # 返回给主进程写 match 表
    row = {
        "species": sp,
        "gff_file": cache_row.get("gff_file", str(gff_path)),
        "genome_file": cache_row.get("genome_file", str(fa_path)),
        "cov_frac": cache_row.get("cov_frac","0.0"),
        "status": "PASS_CODON_SAFE",
        "pep_kept": str(npep),
        "cds_kept": str(n_cds),
        "note": "",
        "timestamp": datetime.datetime.now().isoformat(timespec="seconds")
    }
    print(f"[DONE] {sp}: codon-safe CDS={n_cds}, PEP={npep}")
    return row

# ========================== 顶层 worker（为了解决 pickling 问题） ==========================
def _worker(argtuple):
    """
    进程池调用的顶层函数（必须在模块顶层，不能嵌在 main()）
    返回 ("OK", row) 或 ("ERR", {"species": sp, "error": msg})
    """
    try:
        return ("OK", process_one_species(argtuple))
    except Exception as e:
        sp = argtuple[0] if isinstance(argtuple, (list, tuple)) and argtuple else "?"
        return ("ERR", {"species": sp, "error": str(e)})

# ========================== 主流程 ==========================
def main():
    os.environ.setdefault("LC_ALL", "C")

    # 依赖检查 & 目录准备
    for t in TOOLS: need(t)
    for d in [OUT_PEP_DIR, OUT_CDS_DIR, TMP_DIR, LIST_DIR]: ensure_dir(Path(d))

    print("[STEP] 启动前清理临时目录 ...")
    if Path(TMP_DIR).exists():
        try:
            shutil.rmtree(TMP_DIR, ignore_errors=True)
            print(f"[CLEAN] removed {TMP_DIR}")
        except Exception as e:
            print(f"[WARN] failed to remove old tmp: {e}")

    annot = Path(ANNOT_DIR); genome = Path(GENOME_DIR)
    species = list_species_from_annotation(annot)
    if not species:
        print("[ERR] annotation 目录中未发现 GFF/GFF3。", file=sys.stderr); sys.exit(1)

    match_path = Path(MATCH_TSV)
    cache = load_match_table(match_path)
    rows_out: Dict[str, Dict[str, str]] = {sp: cache[sp] for sp in cache}

    print("[CHECK] 覆盖度检查（严格模式：任一失败→退出） ...")
    pairs = []
    for sp in species:
        gff = find_by_stem(sp, annot, GFF_EXTS)
        fa  = find_by_stem(sp, genome, FA_EXTS)
        if not gff or not fa:
            rows_out[sp] = {"species":sp,"gff_file":str(gff or ""), "genome_file":str(fa or ""),"cov_frac":"0.0","status":"FAIL_NO_FILE","pep_kept":"0","cds_kept":"0","note":"missing GFF or Genome","timestamp":datetime.datetime.now().isoformat(timespec="seconds")}
            print(f"[FATAL] {sp}: 缺少 GFF 或 Genome。", file=sys.stderr)
            write_match_table(match_path, [rows_out[k] for k in sorted(rows_out)])
            sys.exit(1)

        # 规范化工作区在每个 worker 中进行，这里只做覆盖度检查
        cov, gcnt, fcnt, ocnt = coverage_details(gff, fa)
        print(f"  - {sp}: coverage = {cov:.1%}")
        print(f"    · GFF seqid: {gcnt}  Genome seqid: {fcnt}  Overlap: {ocnt}")

        rows_out[sp] = {
            "species": sp,
            "gff_file": str(gff),
            "genome_file": str(fa),
            "cov_frac": f"{cov:.6f}",
            "status": ("PASS_CHECK" if cov>=SEQID_COVERAGE_MIN else "FAIL_COVERAGE"),
            "pep_kept": "0",
            "cds_kept": "0",
            "note": "",
            "timestamp": datetime.datetime.now().isoformat(timespec="seconds")
        }

        if cov < SEQID_COVERAGE_MIN:
            print(f"[FATAL] {sp}: GFF vs genome 覆盖率 {cov:.1%} < {SEQID_COVERAGE_MIN:.0%}。", file=sys.stderr)
            write_match_table(match_path, [rows_out[k] for k in sorted(rows_out)])
            sys.exit(4)

        out_pep = Path(OUT_PEP_DIR)/f"{sp}.faa"
        out_cds = Path(OUT_CDS_DIR)/f"{sp}.fna"
        if RESUME_SKIP_PASS and sp in cache and cache[sp].get("status","").startswith("PASS") \
           and out_pep.exists() and out_pep.stat().st_size>0 and out_cds.exists() and out_cds.stat().st_size>0:
            rows_out[sp] = {**cache[sp]}
            print(f"  - {sp}: 发现历史输出且状态PASS → 跳过提取")
            continue

        pairs.append((sp, gff, fa, rows_out[sp]))

    if not pairs and not any((Path(OUT_PEP_DIR)/f"{sp}.faa").is_file() for sp in species):
        print("[FATAL] 没有物种通过覆盖度检查，且无历史输出。", file=sys.stderr)
        write_match_table(match_path, [rows_out[k] for k in sorted(rows_out)])
        sys.exit(2)

    # ============ 并行跑物种 ============
    print(f"[STEP] 物种并行处理：workers={N_WORKERS}, tasks={len(pairs)}")
    results: List[Dict[str, str]] = []
    failures: List[Tuple[str, str]] = []

    with Pool(processes=N_WORKERS) as pool:
        for status, payload in pool.imap_unordered(_worker, pairs):
            if status == "OK":
                results.append(payload)
            else:
                failures.append((payload.get("species","?"), payload.get("error","")))

    # 将成功结果写入 rows_out
    for r in results:
        rows_out[r["species"]] = r
    # 报告失败但不影响其它物种产出
    if failures:
        print("\n[WARN] 以下物种处理失败（已跳过，详见错误信息）：")
        for sp, err in failures:
            print(f"  - {sp}: {err}")
            # 保留失败状态
            if sp in rows_out:
                rows_out[sp]["status"] = "FAIL_PROCESS"
                rows_out[sp]["note"] = err[:200]

    # species.list 基于已有 proteomes 生成
    ensure_dir(Path(LIST_DIR))
    existing = [p.stem for p in Path(OUT_PEP_DIR).glob("*.faa") if p.stat().st_size>0]
    Path(SPECIES_LIST).write_text("\n".join(sorted(set(existing))) + "\n", encoding="utf-8")

    # 最终清理与防污染（TMP 由开关控制）
    print("[STEP] 最终清理 & 污染防护 ...")
    try:
        if CLEANUP_TMP_AT_END:
            if Path(TMP_DIR).exists():
                shutil.rmtree(TMP_DIR, ignore_errors=True)
                print(f"[CLEAN] tmp dir {TMP_DIR} removed.")
        else:
            print(f"[KEEP] 按开关保留临时目录：{TMP_DIR}")

        # 保持原有注释目录“中间产物扫除”策略
        ann = Path(ANNOT_DIR)
        for bad in ann.glob("*"):
            if re.search(r'\.(clean4agat|agat_longest|keep_longest|fixloc|fixphase|agat.*)\.gff3(\.gz)?$', bad.name):
                try:
                    bad.unlink()
                    print(f"[CLEAN] removed stray middleware in annotation: {bad.name}")
                except Exception as e:
                    print(f"[WARN] cannot remove {bad}: {e}")
    except Exception as e:
        print(f"[WARN] final cleanup exception: {e}")

    # 写 match 表
    write_match_table(Path(MATCH_TSV), [rows_out[k] for k in sorted(rows_out)])

    # 【新增】合并每物种 pep↔cds 分片为全局映射（不影响既有产物）
    if WRITE_PEP2CDS:
        try:
            parts_dir = Path(PEP2CDS_PER_SPECIES_DIR)
            all_parts = sorted(parts_dir.glob('*.tsv'))
            if all_parts:
                Path(PEP2CDS_TSV).parent.mkdir(parents=True, exist_ok=True)
                with open(PEP2CDS_TSV, 'w', encoding='utf-8') as w:
                    w.write('species\tprotein_id\tcds_id\n')
                    for p in all_parts:
                        with open(p, 'r', encoding='utf-8', errors='ignore') as r:
                            header = True
                            for line in r:
                                if header and line.lower().startswith('species\tprotein_id\tcds_id'):
                                    header = False
                                    continue
                                header = False
                                w.write(line)
                print(f"[OUT] pep2cds: {PEP2CDS_TSV}")
            else:
                print("[WARN] pep2cds switch on but no per-species parts found.")
        except Exception as e:
            print(f"[WARN] merge pep2cds failed: {e}")

    print(f"[OUT] Proteomes: {OUT_PEP_DIR}/*.faa")
    print(f"[OUT] CDS      : {OUT_CDS_DIR}/*.fna")
    print(f"[OUT] Species  : {SPECIES_LIST}")
    print(f"[OUT] MatchTbl : {MATCH_TSV}")
    print("[ALL DONE] 🎉")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        # 主进程兜底报错（例如依赖缺失）
        print(str(e), file=sys.stderr)
        sys.exit(1)
