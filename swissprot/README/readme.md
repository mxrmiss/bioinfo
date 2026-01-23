# Swiss-Prot Annotation Pipeline (DIAMOND + BestHit)

## 1. 目的
对每个物种的蛋白质序列（.faa）进行 UniProt Swiss-Prot（reviewed）同源注释：
- 用 DIAMOND blastp 比对到 Swiss-Prot
- 按阈值过滤命中
- 每个 query 蛋白选 1 条“最佳命中”(best hit)
- 解析 Swiss-Prot 标题字段，生成结构化注释表
- 汇总跨物种统计与长表，并把各物种注释表复制到 99_reports 目录

---

## 2. 输入与目录约定
项目根目录下需包含：
- config.yaml
- data/query/{species_id}.faa                # 每个物种一个蛋白 FASTA（header 为 protein_id）

脚本会自动创建/使用：
- data/db/raw/                               # Swiss-Prot 下载文件
- data/db/diamond/                           # DIAMOND 数据库（.dmnd）
- results/01_db/{REL}/                       # 下载与建库留档
- results/02_align/{REL}/{species}/          # DIAMOND 比对输出
- results/03_filter/{REL}/{species}/         # 过滤与 best hit
- results/04_annot/{REL}/{species}/          # 注释表 + QC
- results/99_reports/{REL}/                  # 汇总表 + 复制出的各物种注释表
- logs/                                      # 运行日志

---

## 3. 配置文件（config.yaml）
关键参数：
- species_prefixes: 物种列表（决定流程跑哪些物种）
- threads: DIAMOND 线程数
- diamond:
  - sensitivity: sensitive / very-sensitive / fast
  - evalue: DIAMOND blastp 的 evalue
  - max_target_seqs: 每条 query 最多输出多少命中（后续仍会 besthit 选 1）
- filter（03_filter.py 使用）：
  - EV_MAX: evalue 最大值
  - MIN_BITSCORE: bitscore 最小值
  - MIN_ALN_LEN: 对齐长度最小值（aa）
  - MIN_QCOV: query coverage 最小值
  - MIN_SCOV: subject coverage 最小值（当前默认 0，不限制）
- db（01_db.py 使用）：
  - use_metalink / metalink_select_file / metalink_preferred_protocol

---

## 4. 运行方式（推荐用 Snakemake）
在项目根目录运行：
- snakemake -j 1
（实际并行由 DIAMOND 的 threads 控制；species 逐个跑由规则控制）

也可单独运行脚本（按顺序）：
1) python scripts/01_db.py
2) python scripts/02_align.py        # 可用环境变量 SPECIES_ID 指定单个物种
3) python scripts/03_filter.py
4) python scripts/04_annot.py
5) python scripts/99_reports.py

---

## 5. 处理流程（每一步做什么）
### Step 01: 下载与建库（01_db.py）
- 下载：reldate.txt / README / RELEASE.metalink / uniprot_sprot.fasta.gz
- 用 RELEASE.metalink 提取 md5 并校验 uniprot_sprot.fasta.gz
- 解压出 uniprot_sprot.fasta
- diamond makedb 构建：data/db/diamond/uniprot_sprot_{REL}.dmnd
- 写入版本指针：results/01_db/CURRENT_REL.txt

### Step 02: 比对（02_align.py）
- DIAMOND blastp：query = data/query/{species}.faa vs Swiss-Prot dmnd
- 输出：results/02_align/{REL}/{species}/hits.diamond.tsv
- 同时记录：query 序列条数、md5、diamond 版本等

### Step 03: 过滤 + 每条 query 选 best hit（03_filter.py）
输入：hits.diamond.tsv（outfmt6，含 qlen/slen/stitle）
过滤阈值（默认）：
- evalue ≤ 1e-5
- bitscore ≥ 50
- alignment length ≥ 50 (aa)
- qcov ≥ 0.50  （qcov = qspan/qlen；qspan = |qend-qstart|+1；并 clamp 到[0,1]）
- scov ≥ 0.00  （scov = sspan/slen；当前不限制）
输出：
- hits.filtered.tsv：通过阈值的所有命中
- hits.besthit.tsv：每个 qseqid 仅保留 1 条最佳命中
  besthit 排序规则：evalue 越小越好 → bitscore 越大越好 → qcov 越大越好 → pident 越大越好
- summary.tsv：N_QUERY / N_FILTERED_HITS / N_BESTHIT / HIT_RATE

### Step 04: 解析 Swiss-Prot 标题生成注释表（04_annot.py）
输入：hits.besthit.tsv
从 sseqid/stitle 解析：
- accession / entry_name（如 sp|P12345|NAME）
- protein_name（OS= 之前的描述）
- organism（OS=...）
- gene_name（GN=...，若标题没有 GN 则为空）
输出：
- swissprot_annotation.tsv：结构化注释表（每个 query 一行）
- qc_missing_fields.tsv：统计 OS/GN 缺失率
- md5_swissprot_annotation.txt：注释表 md5

### Step 99: 汇总（99_reports.py）
输出到 results/99_reports/{REL}/：
- species_summary.tsv：每物种注释覆盖与缺失率汇总
- all_annotations.long.tsv：所有物种注释表拼成长表（首列 species_id）
- {species}.swissprot_annotation.tsv：把每个物种的最终注释表原样复制到 99 目录，便于集中取用

---

## 6. 单个物种注释表（swissprot_annotation.tsv）列含义与怎么看
位置：
- results/04_annot/{REL}/{species}/swissprot_annotation.tsv
以及（复制件）：
- results/99_reports/{REL}/{species}.swissprot_annotation.tsv

表头（每行对应 1 个 query 蛋白的 best hit 注释）：
1) protein_id
   - 你的查询蛋白 ID（来自 .faa 的 FASTA header，即 DIAMOND 的 qseqid）
   - 用它把注释回填到你的基因/蛋白列表里

2) swissprot_accession
   - Swiss-Prot 命中条目的 accession（去掉了可能的版本号后缀）
   - 用于唯一定位 Swiss-Prot 条目（稳定、可引用）

3) entry_name
   - Swiss-Prot 的 entry name（如某些条目会有 NAME）
   - 便于人工浏览与交叉核对（不是所有分析都必需）

4) protein_name
   - Swiss-Prot 标题中的蛋白名称描述（OS= 之前的部分）
   - 适合用于结果表“蛋白功能名称/注释描述”

5) organism
   - Swiss-Prot 条目的物种名（来自标题 OS=...）
   - 用于检查注释来源分布（例如是否大量命中到某些模式物种）

6) gene_name
   - Swiss-Prot 标题中的基因名（GN=...）
   - 可能为空：某些 Swiss-Prot 条目没有 GN 字段（属正常现象）

7) evalue
   - 该 best hit 的 e-value（越小越可信）
   - 主要用于质量把控与排序参考

8) bitscore
   - 该 best hit 的 bit score（越大越好）
   - 与 evalue 一起衡量匹配强度

9) pident
   - 对齐区段的百分比一致性（% identity）
   - 反映相似程度，但受对齐区段影响（需结合 coverage 看）

10) qcov
   - query coverage：命中覆盖 query 蛋白的比例（0-1）
   - 你当前阈值要求 qcov ≥ 0.5；越接近 1 通常越可信

11) scov
   - subject coverage：命中覆盖 Swiss-Prot 蛋白的比例（0-1）
   - 当前默认不设下限（MIN_SCOV=0），用于辅助判断是否是“只命中对方的一小段”

快速解读建议：
- 先看 qcov：接近 1 的通常注释最稳；接近 0.5 的属于“刚达标”，建议结合 pident/bitscore 判断
- 再看 evalue 与 bitscore：evalue 越小、bitscore 越大通常越可靠
- gene_name 为空不代表注释失败，只表示该 Swiss-Prot 标题未提供 GN

