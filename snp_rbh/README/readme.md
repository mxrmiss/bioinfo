# README.md

# snp_rbh — RBH-based divergence landscape (codon-level)

本项目用于比较两个物种的蛋白序列并识别 **RBH (reciprocal best hits)** 的一对一同源对，随后对每对同源基因做 **MAFFT 蛋白对齐 + PAL2NAL 回译 codon alignment**，统计 substitution/indel 差异，并在参考物种染色体上做 **chr 汇总 + 1Mb 滑窗汇总**，最后输出 3 张图（chr 柱状图 / 滑窗分面折线 / 滑窗热图）。

---

## 目录结构（必须保持）
snp_rbh/
├── data/
│   ├── annotation/      # GFF/GFF3（可 .gz）
│   ├── cds/             # CDS（当前脚本默认非 .gz）
│   ├── genomic/         # genome FASTA（可 .gz，用于 chr rename 自动推断 code）
│   └── proteomes/       # protein FASTA（当前脚本默认非 .gz）
├── scripts/             # 00~05 脚本 + pal2nal.pl
├── results/             # step0_chr ... step6_plots
└── logs/

你的当前数据布局示例：
- data/proteomes/Sinonovacula_constricta.faa
- data/proteomes/Sinonovacula_rivularis.faa
- data/cds/Sinonovacula_constricta.fna
- data/cds/Sinonovacula_rivularis.fna
- data/annotation/Sinonovacula_constricta.gff3.gz
- data/annotation/Sinonovacula_rivularis.gff.gz
- data/genomic/Sinonovacula_constricta.fa.gz
- data/genomic/Sinonovacula_rivularis.fasta.gz

---

## 流程概览（00 → 05）
00_chr_rename_map.py
  生成参考物种 chr 重命名映射 old_seqid -> chr01..（默认仅 >=10Mb 序列）
01_rbh_diamond.sh
  DIAMOND 双向 blastp（A→B、B→A）→ 覆盖度过滤 → 每 query 取 best hit → RBH pairs（第一列固定为参考物种）
02_make_pair_fastas.py
  从 RBH pairs 抽取两物种蛋白+CDS，生成每对一个小 fasta（.faa/.fna）
03_run_mafft_pal2nal.sh
  对每对蛋白做 MAFFT → 用 PAL2NAL 回译到 codon alignment（并行）
04_codon_stats_and_chr.py
  统计 substitution/indel/aligned_length，映射参考 GFF 坐标，输出 pair/chr/1Mb window 三张表
05_plot_density.R
  读取 step5_stats 画 3 张图（chr bar / window facet line / window heatmap），支持 METRIC_MODE 开关

---

## 快速运行（按顺序）
注意：所有脚本都使用相对路径并自动定位项目根目录；直接在 snp_rbh/ 下运行即可。

Step 0（可选但强烈建议）：chr 重命名映射
1) 编辑 scripts/00_chr_rename_map.py：
   - GENOME_FASTA = 参考物种 genome（例如 data/genomic/Sinonovacula_rivularis.fasta.gz）
   - CHROM_MIN_BP = 10_000_000（默认 10Mb）
2) 运行：
   python3 scripts/00_chr_rename_map.py
输出：
   results/step0_chr/<code>.chr_rename_map.tsv

Step 1：DIAMOND + RBH
1) 编辑 scripts/01_rbh_diamond.sh：
   - A_PROTEOME/B_PROTEOME：两物种蛋白
   - EVALUE/MIN_QCOV/MIN_SCOV：阈值
   - REFERENCE="A" 或 "B"：决定 rbh_pairs.tsv 第一列是谁（建议参考物种放第一列）
2) 运行：
   bash scripts/01_rbh_diamond.sh
输出（核心）：
   results/step2_rbh/rbh_pairs.tsv

Step 2：生成每对小 fasta（蛋白+CDS）
1) 编辑 scripts/02_make_pair_fastas.py：
   - PAIRS_TSV：results/step2_rbh/rbh_pairs.tsv
   - REF_PROTEOME/REF_CDS：必须对应 pairs 第一列的物种
   - OTH_PROTEOME/OTH_CDS：对应第二列物种
2) 运行：
   python3 scripts/02_make_pair_fastas.py
输出：
   results/step3_align_prot/pairs_faa/*.faa
   results/step4_codon/pairs_fna/*.fna

Step 3：MAFFT + PAL2NAL（并行）
1) 编辑 scripts/03_run_mafft_pal2nal.sh：
   - TOTAL_THREADS / JOBS / MAFFT_THREADS
   - PAL2NAL="scripts/pal2nal.pl"
2) 运行：
   bash scripts/03_run_mafft_pal2nal.sh
输出：
   results/step3_align_prot/mafft_aa/*.aa.aln
   results/step4_codon/codon_aln/*.codon.fa
日志：
   logs/mafft.<pair>.log
   logs/pal2nal.<pair>.log

Step 4：统计（pair / chr / 1Mb window）
1) 编辑 scripts/04_codon_stats_and_chr.py：
   - CODON_DIR：results/step4_codon/codon_aln
   - REF_GFF：参考物种 GFF（用于 transcript ID= 定位 chr/start/end）
   - REF_GENOME：参考 genome（用于自动推断 chr rename map 路径）
   - WINDOW_BP：默认 1_000_000
2) 运行：
   python3 scripts/04_codon_stats_and_chr.py
输出（核心，step5_stats 合同）：
   results/step5_stats/pair_stats.tsv
   results/step5_stats/chr_summary.tsv
   results/step5_stats/window_1Mb.tsv

Step 5：画图（3 张图）
1) 编辑 scripts/05_plot_density.R：
   - METRIC_MODE="pdistance" 或 "differences"
   - MIN_PAIRS_WINDOW：默认 0（不启用过滤）
   - OUTDIR：默认 results/step6_plots
2) 运行：
   Rscript scripts/05_plot_density.R
输出：
   results/step6_plots/density.<suffix>.chr_bar.png/pdf
   results/step6_plots/density.<suffix>.window_1Mb.facet_line.png/pdf
   results/step6_plots/density.<suffix>.window_1Mb.heatmap.png/pdf

---

## 关键口径说明（避免误解）
- 所有 “per Mb” 的密度分母都是 **aligned bases**（对齐且双方为 A/C/G/T 的位点数），不是染色体物理长度。
- pdistance 模式：严格 substitutions-only，展示为 (substitutions / aligned) * 1e6。
- differences 模式：substitutions + indel_bases，优先使用 step5_stats 的 density_per_Mb_aligned。
- 参考端 ref_id 自动判定：谁能在 REF_GFF 的 transcript(ID=) 里找到，谁就是参考端；pal2nal 输出序列顺序不影响统计。

