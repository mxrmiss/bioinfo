#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
零参数·一次学习·批量生成 samples.tsv / contrasts.tsv
- 仅对第一个样本询问 sample/group 名称，自动学习 group 规则应用到所有样本
- 智能识别 R1/R2（_R1/_2/.R1_001/-R2 等），支持多 lane（逗号合并）
- 默认 FASTQ 目录 = 脚本所在目录；交互简洁
"""

import os
import re
import sys
from collections import defaultdict
from itertools import combinations

# ========== 通用 ==========
EXT_RE = re.compile(r"\.(?:f(?:ast)?q)(?:\.(?:gz|bz2|xz|zst))?$", re.IGNORECASE)

def prompt_with_default(msg: str, default: str) -> str:
    s = input(f"{msg}（默认：{default}）：").strip()
    return s if s else default

# ========== 文件扫描与解析 ==========
READ_PATTERNS = [
    re.compile(r"^(?P<base>.+?)[._-]R(?P<read>[12])(?:[._-]?\d+)?\.(?:f(?:ast)?q)(?:\.(?:gz|bz2|xz|zst))?$", re.IGNORECASE),
    re.compile(r"^(?P<base>.+?)[._-](?P<read>[12])(?:[._-]?\d+)?\.(?:f(?:ast)?q)(?:\.(?:gz|bz2|xz|zst))?$", re.IGNORECASE),
]

def list_fastq_files(indir):
    if not os.path.isdir(indir):
        return []
    return sorted([f for f in os.listdir(indir) if EXT_RE.search(f or "")])

def split_read(fname):
    """
    稳健读端解析：去扩展名 -> 以 -/_/. 切分 -> 右向左找第一个 R1/R2 或 1/2
    base = 该标记左侧所有 token 用下划线连接
    """
    name = os.path.basename(fname)
    stem = re.sub(EXT_RE, "", name)
    tokens = re.split(r"[._-]+", stem)
    if not tokens:
        return None, None
    read_idx = None
    read = None
    for i in range(len(tokens)-1, -1, -1):
        t = tokens[i].upper()
        if re.fullmatch(r"R[12]", t):
            read = int(t[1]); read_idx = i; break
        if t in {"1","2"}:
            read = int(t); read_idx = i; break
    if read_idx is None:
        return None, None
    base_tokens = tokens[:read_idx]
    if not base_tokens:
        return None, None
    base = "_".join(base_tokens)
    return base, read

def find_pairs(indir):
    files = list_fastq_files(indir)
    if not files:
        print(f"× 在目录 {indir} 下未发现 FASTQ 文件。")
        return {}, None

    print(f"\n目录中的第一个文件是：{files[0]}")
    print(f"共找到 {len(files)} 个 FASTQ 文件。")

    pairs = defaultdict(lambda: {1: [], 2: []})
    unparsed = []
    for f in files:
        base, read = split_read(f)
        if base is None:
            unparsed.append(f); continue
        pairs[base][read].append(os.path.join(indir, f))

    cleaned = {}
    for base, d in pairs.items():
        if not d[1]:
            continue
        r1 = ",".join(d[1]); r2 = ",".join(d[2]) if d[2] else ""
        cleaned[base] = {1: r1, 2: r2}

    if unparsed:
        print(f"（提示）有 {len(unparsed)} 个文件未识别为 R1/R2，已忽略。")

    if not cleaned:
        print("× 未找到可用的 R1/R2 配对（至少需要 R1）。")
        return {}, None

    first_key = sorted(cleaned.keys())[0]
    return cleaned, first_key

# ========== “只问一次”的规则学习 ==========
def tokens_of_base(base: str):
    return re.split(r"[-_.]+", base)

def learn_group_rule_from_first(base_first: str, group_first: str):
    """
    从第一个样本的 (base_first, group_first) 学习 group 提取规则：
    规则优先级：
      1) group 等于 base 切分后的某个 token -> 记录其 token_index
      2) group 是第一个 token 的前缀（如 D0 匹配 D0、D1、D15 等）-> 记录 'prefix' 策略
      3) 回退：总是取 token[0]
    返回：callable f(base)->group
    """
    toks = tokens_of_base(base_first)
    # 1) 直接匹配某个 token
    for idx, tk in enumerate(toks):
        if group_first == tk:
            def f(base, i=idx):
                tt = tokens_of_base(base)
                return tt[i] if len(tt) > i else tt[0]
            return f

    # 2) 作为 token[0] 的前缀（常见 D0/D1/D8……）
    if toks:
        if toks[0].startswith(group_first):
            prefix = re.escape(re.match(r"^[A-Za-z]+", toks[0]).group(0)) if re.match(r"^[A-Za-z]+", toks[0]) else ""
            # 若首 token 是“字母+数字”，提取“相同字母 + 可变数字”
            if prefix:
                pattern = re.compile(rf"^{prefix}\d+", re.IGNORECASE)
                def f(base, pat=pattern):
                    tt = tokens_of_base(base)
                    m = pat.match(tt[0]) if tt else None
                    return m.group(0) if m else (tt[0] if tt else base)
                return f
            else:
                def f(base):
                    tt = tokens_of_base(base)
                    return tt[0] if tt else base
                return f

    # 3) 回退：取第一个 token
    def f(base):
        tt = tokens_of_base(base)
        return tt[0] if tt else base
    return f

# ========== 写表 ==========
def write_samples_tsv(rows, outpath):
    with open(outpath, "w", encoding="utf-8") as w:
        w.write("sample\tgroup\tfastq1\tfastq2\n")
        for r in rows:
            w.write(f"{r['sample']}\t{r['group']}\t{r['fastq1']}\t{r['fastq2']}\n")
    print(f"\n✔ 已生成样本信息表：{outpath}（{len(rows)} 条）")

# ========== contrasts（智能菜单） ==========
def _guess_baseline(groups):
    cand = {g.lower(): g for g in groups}
    for p in ["d0", "ck", "ctrl", "control", "con"]:
        if p in cand: return cand[p]
    for g in groups:
        if re.fullmatch(r"d0", g, re.IGNORECASE): return g
    return None

def _normalize_contrast(case, control): return f"{case}_vs_{control}"

def _mk_contrasts_from_baseline(groups, baseline):
    return [(_normalize_contrast(g, baseline), g, baseline) for g in sorted(groups) if g != baseline]

def _mk_contrasts_timeseries(groups):
    patt = re.compile(r"^d(\d+)$", re.IGNORECASE)
    d0, nums = None, []
    for g in groups:
        m = patt.match(g)
        if m:
            n = int(m.group(1))
            if n==0: d0 = g
            else: nums.append((n,g))
    if d0 and nums:
        nums.sort()
        return [(_normalize_contrast(g,d0), g, d0) for _,g in nums]
    return []

def _mk_contrasts_all_pairwise(groups):
    out=[]
    for a,b in combinations(sorted(groups),2):
        case, ctrl = (a,b) if a>b else (b,a)
        out.append((_normalize_contrast(case, ctrl), case, ctrl))
    return out

def _dedup_contrasts(cs):
    seen=set(); out=[]
    for c in cs:
        k=(c[1],c[2])
        if k in seen: continue
        seen.add(k); out.append(c)
    return out

def prompt_contrasts_smart(groups):
    groups=set(groups)
    if not groups:
        print("未检测到组别，跳过对比设计。"); return []
    baseline_guess = _guess_baseline(groups)

    print("\n=== 对比设计（contrasts） ===")
    print("可用组别：", ", ".join(sorted(groups)))
    if baseline_guess: print("基线组猜测：", baseline_guess)
    print("""请选择：
  1) 一键：所有组 vs 基线
  2) 一键：时间序列（D>0 vs D0）
  3) 一键：全对全（去重）
  4) 一句话自定义（示例：D1,D8,D15 vs D0）
  0) 跳过
""")
    choice = input("选择：").strip()
    contrasts=[]

    if choice=="1":
        base = prompt_with_default("请输入基线组", baseline_guess or "")
        if not base or base not in groups:
            print("× 基线组无效，已取消。"); return []
        contrasts = _mk_contrasts_from_baseline(groups, base)
    elif choice=="2":
        cs = _mk_contrasts_timeseries(groups)
        if not cs: print("× 未检测到 D0 或 D<数字>。"); return []
        contrasts = cs
    elif choice=="3":
        contrasts = _mk_contrasts_all_pairwise(groups)
    elif choice=="4":
        line = input("请输入一句（如 D1,D8,D15 vs D0）：").strip()
        m = re.match(r"^(.+?)\s+vs\s+(\S+)$", line, re.IGNORECASE)
        if not m: print("× 格式错误。"); return []
        lhs = [x.strip() for x in m.group(1).split(",") if x.strip()]
        base = m.group(2).strip()
        bad = [x for x in lhs+[base] if x not in groups]
        if bad: print("× 不存在的组：", ", ".join(bad)); return []
        contrasts = [(_normalize_contrast(a, base), a, base) for a in lhs if a!=base]
    else:
        print("已选择跳过对比设计。"); return []

    contrasts = _dedup_contrasts(contrasts)
    if contrasts:
        print("将生成：")
        for c in contrasts: print("  -", c[0])
    else:
        print("未生成任何对比。")
    return contrasts

def write_contrasts_tsv(contrasts, outpath):
    with open(outpath, "w", encoding="utf-8") as w:
        w.write("contrast\tcase\tcontrol\n")
        for c,a,b in contrasts:
            w.write(f"{c}\t{a}\t{b}\n")
    print(f"✔ 已生成对比设计表：{outpath}（{len(contrasts)} 条）")

# ========== 主流程 ==========
def main():
    print("=== RNA-seq 表格生成器（一次学习·全自动推断）===\n")

    # 默认目录 = 脚本所在目录
    script_dir = os.path.dirname(os.path.abspath(__file__))
    indir = prompt_with_default("请输入 FASTQ 所在目录", script_dir)

    files = list_fastq_files(indir)
    while not files:
        print("× 未发现 FASTQ 文件，请重新输入。")
        indir = prompt_with_default("请输入 FASTQ 所在目录", indir)
        files = list_fastq_files(indir)

    samples_out   = prompt_with_default("samples.tsv 输出路径", "samples.tsv")
    contrasts_out = prompt_with_default("contrasts.tsv 输出路径", "contrasts.tsv")

    pairs, first_key = find_pairs(indir)
    if not pairs:
        sys.exit(1)

    # —— 只问一次：第一个样本的 sample/group —— #
    print(f"\n第一个样本键：{first_key}")
    sample_first = prompt_with_default("请为第一个样本的 sample 名命名", first_key)
    # 默认 group：第一个 token
    default_group_first = tokens_of_base(first_key)[0] if tokens_of_base(first_key) else first_key
    group_first = prompt_with_default("请为第一个样本的 group 名命名", default_group_first)

    # 学习出 group 提取规则
    group_fn = learn_group_rule_from_first(first_key, group_first)

    # 按学到的规则，自动生成所有样本的行；sample 一律用样本键，首样本可重命名
    rows=[]; groups=set(); seen=set()
    for key in sorted(pairs.keys()):
        sample = sample_first if key == first_key else key
        # 防止首样本重名与后续键冲突
        base_name, i = sample, 2
        while sample in seen:
            sample = f"{base_name}_{i}"; i += 1
        group = group_fn(key)
        rows.append({"sample": sample, "group": group,
                     "fastq1": pairs[key][1], "fastq2": pairs[key][2]})
        seen.add(sample); groups.add(group)

    write_samples_tsv(rows, samples_out)

    # contrasts 智能菜单（可跳过）
    contrasts = prompt_contrasts_smart(groups)
    if contrasts:
        write_contrasts_tsv(contrasts, contrasts_out)

    print("\n✨ 完成！")

if __name__ == "__main__":
    main()
