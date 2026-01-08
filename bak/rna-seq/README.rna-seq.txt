====================================================================================================
RNA-SEQ PIPELINE —— 超详细说明文档 README.txt
Version: stable | Format: pure text | Audience: researchers & pipeline users
Author: Li Ziqiang
====================================================================================================

本 README 说明当前 RNA-seq 项目中所有脚本的用途、输入依赖、配置字段与目录结构。

====================================================================================================
项目简介
====================================================================================================

本 RNA-seq 流水线假定已经具备：

  - 一套原始转录组测序数据（FASTQ，单端或双端）；
  - 对应物种的参考基因组及注释文件（genome.fa + GTF/GFF3）；
  - eggNOG-mapper 或本地功能注释结果（用于 GO/KEGG 注释）；
  - 完整的样本设计表与对比设计表（samples.tsv / contrasts.tsv）。

在此基础上，RNA-seq 流程负责完成：

  1) FASTQ 质控与清洗（fastp + cutadapt）；
  2) 参考准备与 Salmon 索引构建（gentrome + decoy-aware index）；
  3) Salmon 表达定量（转录本层）与 tximport 基因层聚合；
  4) 基于 DESeq2 的差异表达分析（DEG）；
  5) 基于 eggNOG/本地字典的 GO/KEGG 注释与 over-representation analysis（ORA）；
  6) 可选 GSEA 分析与统一的富集输入管理（背景集 / test 集）；
  7) 最终结果打包与 Methods skeleton/manifest 生成，面向投稿的 Source Data 发布。


====================================================================================================
一、整体流程与推荐运行顺序
====================================================================================================

推荐严格顺序如下（括号内为脚本文件名）：

1) 实验设计与对比定义（可选）
   - 00_make_sample_contrasts.py

2) 质控与表达定量
   - 01_qc_fastp_cutadapt.py
   - 02_refprep_salmon_index.py
   - 03_build_tx2gene_map.py
   - 04_salmon_quant.py
   - 05_matrix_aggregate.py
   - 05_tximport_aggregate.R

3) 差异表达模块
   - 06_deg_module.py
   - 06_f_deg.R

4) 注释与 GO/KEGG 字典
   - 07_prepare_emapper_annotations.py
   - go_local.R
   - kegg_local.sh

5) 富集分析 ORA/GSEA
   - 08_enrich_module.py
   - 08_g_enrich.R

6) 最终发布与装订
   - 09_publish_results.py


====================================================================================================
二、脚本逐个说明（功能 / 主要输入 / 特点 / 注意事项）
====================================================================================================

----------------------------------------
1) 00_make_sample_contrasts.py
----------------------------------------
功能：
  - 根据给定的 FASTQ 布局和分组信息，半自动构建并校验 samples.tsv 与 contrasts.tsv：
      * samples.tsv：每行一个 sample，包含 group、fastq 路径等；
      * contrasts.tsv：每行一个对比，定义 contrast/case/control。
  - 提供一次性体检：检查 fastq 是否存在、同一 sample 是否重复、多余/缺失分组等。

主要输入：
  - config.yaml：
      * data.samples_tsv
      * data.contrasts_tsv
  - 目录结构：
      * 原始 FASTQ 所在目录（由 config 中相关字段统一指定）。

特点：
  - 不参与后续统计运算，核心价值是“把实验设计变成机读 TSV 文件”；
  - 所有路径与 group 定义都写回 TSV，方便手工校对和版本管理。

注意事项：
  - 一旦 samples/contrasts 确认无误，后续无需频繁重跑本脚本；
  - 若实验设计发生变更（新增样本或对比），建议重新跑一次并更新 TSV。


----------------------------------------
2) 01_qc_fastp_cutadapt.py
----------------------------------------
功能：
  - 读取 samples.tsv 中的 FASTQ 信息，对每个样本执行：
      * fastp 质控过滤（Q 值、长度、N 碱基等）；
      * cutadapt 剪切接头（如配置）；
      * 输出 clean FASTQ 与 fastp JSON/HTML 报告。
  - 聚合生成样本级 QC 统计表与 outlier/reject 列表。

主要输入：
  - config.yaml：
      * data.samples_tsv
      * reference.genome_fa（可选，用于某些特殊过滤）
      * qc.thresholds.*（最小保留比例、最小 raw reads 等）
      * qc.outlier_rules.*（接头比例等判断规则）
      * resources.threads.qc
      * binaries.fastp / binaries.cutadapt
      * paths.qc_dir / paths.logs_dir
  - FASTQ 文件：按 samples_tsv 中路径定位。

特点：
  - 阈值全部来自 config.qc，便于跨项目统一标准；
  - 将“边缘样本”与“必剔样本”区分开来，防止误删。

注意事项：
  - 若 min_raw_reads / min_retention_rate 设置过高，会导致大量样本进入 rejects；
  - 建议先用小批样本试跑一次，确认 QC 报告符合预期，再批量跑全体。


----------------------------------------
3) 02_refprep_salmon_index.py
----------------------------------------
功能：
  - 根据参考基因组和 GTF/GFF3 注释，构建 Salmon 所需的 gentrome 和 decoy-aware 索引：
      * 使用 gffread 从 genome + GTF/GFF3 提取 transcripts.fa；
      * 构建 decoys.txt（染色体/contig 列表）；
      * 如有需要，生成 gentrome.fa（transcripts + decoys）；
      * 调用 salmon index 生成 ref/salmon_index/。

主要输入：
  - config.yaml：
      * reference.genome_fa
      * reference.annotation_gff / annotation_gtf
      * reference.emapper（仅记录，不直接在本脚本使用）
      * reference.salmon.index_dir / kmer_len / decoy_list / gentrome_fa / rebuild
      * resources.threads.salmon
      * binaries.salmon / binaries.gffread
      * paths.ref_dir / paths.logs_dir

特点：
  - 支持“已有索引则跳过”和“强制重建”等模式（rebuild: always/if_missing/never）；
  - decoy_list/gentrome_fa 可事先指定，也可由本脚本自动生成。

注意事项：
  - 一旦索引构建完毕，将被所有后续 quant 复用，一般无需重跑；
  - 参考基因组或注释文件变动时，必须重新运行本脚本并重建 index。


----------------------------------------
4) 03_build_tx2gene_map.py
----------------------------------------
功能：
  - 从 GTF/GFF3 和可选的 emapper 注释中抽取 transcript_id 与 gene_id 对应关系：
      * 生成 tx2gene.clean.tsv：两列（transcript_id, gene_id）；
      * 生成 tx2gene.blacklist.tsv：转录本数超阈值的基因（用于 QC）。

主要输入：
  - config.yaml：
      * reference.annotation_gff / annotation_gtf
      * annotations.*（ID 清洗策略、alias/strip prefix/suffix 等）
      * paths.maps_dir / paths.logs_dir

特点：
  - 所有 ID 规范化策略（大小写、前后缀等）集中在 annotations 段配置；
  - blacklists 有助于发现异常基因（如伪基因或注释错误）。

注意事项：
  - tx2gene.clean.tsv 是 tximport 聚合的唯一“ID 真相源”，后续请不要手工修改；
  - 若 emapper 等外部注释使用了不同的 ID 命名规则，必须在 annotations 中统一清洗。


----------------------------------------
5) 04_salmon_quant.py
----------------------------------------
功能：
  - 读取 samples.tsv、reference.salmon.index_dir，批量执行 Salmon quant：
      * 支持单端/双端；自动识别 clean/raw 作为输入；
      * 生成 quant.sf 及标准 Salmon 辅助文件。

主要输入：
  - config.yaml：
      * data.samples_tsv
      * reference.salmon.index_dir
      * salmon.*（library 类型、是否覆盖已有结果、额外参数等）
      * resources.threads.salmon
      * binaries.salmon
      * paths.quant_dir / paths.qc_dir / paths.logs_dir

特点：
  - 每个 sample 有独立的 quant 子目录，方便手工检查与删除；
  - 对已存在且完整的 quant 结果支持“跳过重跑”策略。

注意事项：
  - 文库类型建议统一用 -l A（自动推断），否则需要在 config.salmon 中显式指定；
  - 若清洗后 reads 极少，Salmon 可能对该样本给出警告或无法定量，应结合 QC 决定是否剔除。


----------------------------------------
6) 05_matrix_aggregate.py
----------------------------------------
功能：
  - 汇总所有 quant.sf，调用 tximport 生成基因层 counts/TPM 矩阵：
      * 构建设计好的 tximport_meta.tsv（sample 与 quant.sf 路径匹配）；
      * 调用 05_tximport_aggregate.R 完成聚合与矩阵统计。

主要输入：
  - config.yaml：
      * data.samples_tsv
      * data.contrasts_tsv
      * reference.maps（tx2gene.clean.tsv 的路径）
      * paths.quant_dir / paths.matrix_dir / paths.logs_dir
      * resources.threads.tximport
      * tximport.*（counts_from_abundance 等）

特点：
  - 将样本 QC 剔除信息（rejects.tsv）一并考虑，避免无效样本进入矩阵；
  - 所有路径与 meta 信息都写入 tximport_meta.tsv，便于复盘。

注意事项：
  - 运行前应确认 03 已生成 tx2gene.clean.tsv，且 quant.sf 数量与 samples.tsv 样本数一致；
  - 若 counts_from_abundance 调整策略改变，应重新运行本脚本及其后续 DEG/富集模块。


----------------------------------------
7) 05_tximport_aggregate.R
----------------------------------------
功能：
  - 在 05_matrix_aggregate.py 调用下执行 tximport：
      * 读取 tximport_meta.tsv 与 tx2gene.clean.tsv；
      * 生成基因层 counts 与 TPM 矩阵；
      * 计算 matrix_stats（样本库大小、分位数等）。

主要输入：
  - tximport_meta.tsv
  - tx2gene.clean.tsv
  - config.yaml 中的 tximport.* 与 paths.matrix_dir。

特点：
  - matrix_stats.tsv 直接为 Methods 提供一句话描述所需的数值；
  - 使用 tximport 官方推荐的 counts_from_abundance 策略，确保与主流文献一致。

注意事项：
  - 遇到某些 transcript 找不到 gene 的情况，会在日志中给出警告，应结合 annotations 策略排查；
  - 若后续 DEG 阶段发现“某些样本全 0 或极小”，优先检查本脚本输出及 tximport_meta.tsv 是否正确。


----------------------------------------
8) 06_deg_module.py
----------------------------------------
功能：
  - 根据 counts 矩阵与 contrasts.tsv，批量调度 06_f_deg.R 完成所有对比的 DESeq2 分析；
  - 统一管理 DEG 输出目录结构与日志。

主要输入：
  - config.yaml：
      * data.contrasts_tsv
      * paths.matrix_dir（counts）/ paths.deg_dir
      * resources.threads.deg
      * deg.*（lfc 阈值、FDR、是否使用 apeglm shrinkage 等）
      * logging.* / binaries.Rscript
  - gene_counts.tsv（05 步生成）。

特点：
  - 分离“对比设计”与“统计实现”：design 全在 contrasts.tsv 中定义；
  - 每个 contrast 产物都放在独立子目录，结构统一。

注意事项：
  - 若某 contrast 的组内样本数 < deg.min_samples_per_group，将被自动跳过并记录在日志中；
  - 若存在 batch 等协变量，需要在 06_f_deg.R 与 contrasts/design 中共同定义，保证模型一致。


----------------------------------------
9) 06_f_deg.R
----------------------------------------
功能：
  - 对一个给定 contrast 执行完整 DESeq2 流程：
      * 构建 DESeqDataSet（支持 batch 因子）；
      * 正规化与离散度估计；
      * 计算 log2FC、pvalue、padj；
      * 输出 DEG_all/DEG_up/DEG_down、varTopN、RLE 诊断等。

主要输入：
  - counts 矩阵（gene_counts.tsv 或子集）
  - contrast 信息（从 06_deg_module.py 传递的 design）
  - config.yaml 中的 deg.* 配置。

特点：
  - 所有阈值集中在 deg 段（lfc/fdr/diagnostics.varTopN 等），避免脚本内部硬编码；
  - 生成多种 QC 指标（varTopN.list / rle_range.tsv），方便后续绘图与批次效应评估。

注意事项：
  - independent_filter 一般建议开启，可提高 FDR 控制效率，但在极少基因项目中要慎用；
  - 如使用 apeglm 等 shrinkage 方法，需要保证相关 Bioc 包安装完整。


----------------------------------------
10) 07_prepare_emapper_annotations.py
----------------------------------------
功能：
  - 从 eggNOG-mapper 或本地功能注释结果中，构建 gene→GO / gene→KO / gene→pathway 三张长表；
  - 统计注释覆盖率（universe_coverage.tsv），作为富集分析 universe 的质量评估。

主要输入：
  - config.yaml：
      * reference.emapper
      * annotations.*（ID 清洗、允许的匹配模式等）
      * paths.annot_dir / paths.logs_dir

特点：
  - 所有 ID 清洗策略与 03 共用，确保 gene_id 在表达/DEG/注释三者之间完全一致；
  - 覆盖率报告帮助判断是否需要补充其他注释来源（如物种特异数据库）。

注意事项：
  - emapper TSV 中若 gene 列命名与流程预期不一致，需要在 annotations 段配置中调整列名映射；
  - 统一选择“以 gene 为单位”还是“以转录本为单位”，避免后续富集混淆。


----------------------------------------
11) go_local.R
----------------------------------------
功能：
  - 通过 GO.db 等 Bioconductor 包构建本地 GO 字典：
      * term2name.tsv：go_id → term_name；
      * obsolete_map.tsv：废弃 GOID 的替换与 consider 信息。

主要输入：
  - config.yaml：
      * go.dict.*（输出路径）
      * go.obsolete_policy / go.maxGS 等
      * paths.annot_dir / paths.logs_dir

特点：
  - 统一 GO term 名称来源为 GO.db，避免网上多个版本名称不一致的问题；
  - 提前处理 obsolete term，降低富集分析时的“找不到 term”错误率。

注意事项：
  - GO.db 等 Bioc 包版本应与 R 版本匹配，否则安装会失败；
  - obsolete_policy 若使用 “replace_or_consider”，需要在 08_g_enrich.R 中配合处理替代 ID。


----------------------------------------
12) kegg_local.sh
----------------------------------------
功能：
  - 从本地 KEGG TSV 或镜像版本构建：
      * term2gene.tsv：pathway → gene 映射；
      * term2name.tsv：pathway ID → 中文/英文名称。

主要输入：
  - config.yaml：
      * kegg.dict.*（输入与输出路径）
      * paths.annot_dir

特点：
  - 完全脱离线上 KEGG API，适应离线环境；
  - 允许根据项目需求过滤特定物种或特定层级的 KEGG 条目。

注意事项：
  - 需要预先准备好合法的 KEGG TSV 文件，并确保 ID 与 gene2ko/gene2pathway 一致；
  - 若未来更换 KEGG 源，只需更新本脚本与 config.kegg.dict，即可复用下游模块。


----------------------------------------
13) 08_enrich_module.py
----------------------------------------
功能：
  - 管理所有 ORA/GSEA 任务的输入与背景：
      * 从 DEG_up/DEG_down/PSG 等来源构建基因列表（inputs/*.list）；
      * 按对比/标签定义背景基因 universe（background/*.list）及 meta 信息；
      * 可选构建 GSEA 排名文件（gsea_ranks.tsv）。

主要输入：
  - config.yaml：
      * data.contrasts_tsv
      * paths.deg_dir / paths.annot_dir / paths.enrich_dir / paths.logs_dir
      * go.dict.* / kegg.dict.* / enrich.*（universe 策略、GSEA 选项等）
      * resources.threads.enrich
  - DEG 结果目录（06 输出）
  - 注释字典（gene2go / gene2pathway 等）

特点：
  - 把“什么基因集合需要做富集”这件事抽象为 tasks.tsv，易于记录和审计；
  - 背景集合与测试集合分离存放，有利于日后重跑 ORA/GSEA 而不依赖 DEG 目录结构。

注意事项：
  - universe 策略（例如“所有有表达的基因” vs “所有有注释的基因”）会影响富集结果，应在 config.enrich 中统一约定；
  - 若 GSEA 启用，需要确保 rank 文件与后续 08_g_enrich.R 读取逻辑匹配。


----------------------------------------
14) 08_g_enrich.R
----------------------------------------
功能：
  - 承接 08_enrich_module.py 的 tasks，完成：
      * GO ORA（BP/CC/MF 分本体或合并）；
      * KEGG ORA；
      * 可选 GO/KEGG GSEA；
      * 输出 by_term 与合并显著表。

主要输入：
  - 基因列表 inputs/*.list
  - 背景列表 background/*.list + meta
  - 注释字典：gene2go / gene2pathway / term2name / obsolete_map / KEGG term2gene 等
  - config.yaml 中 go.*, kegg.* 与 enrich.*。

特点：
  - 统一采用 hypergeometric / Fisher exact + BH/FDR；
  - 所有输出都是长表 TSV，方便后续作图脚本（气泡图、条形图等）直接消费。

注意事项：
  - 若输入集合太小或 universe 太小，可能出现“无有效条目”的情况，应先检查 inputs/background 是否合理；
  - GSEA 模块需要保证 rank 值的方向性定义与差异分析一致（例如按 log2FC 或 statistic 排序）。


----------------------------------------
15) 09_publish_results.py
----------------------------------------
功能：
  - 将流程中的关键 TSV 拷贝到统一的 Source Data 目录，并生成：
      * METHODS_README.txt：方法学说明 skeleton；
      * manifest.tsv：文件清单（含模块、用途等）；
      * publish_check.tsv：各模块产物完整性检查；
      * 可选针对每个 contrast 的 Excel 汇总文件。

主要输入：
  - config.yaml：
      * paths.matrix_dir / paths.deg_dir / paths.annot_dir / paths.enrich_dir / paths.publish_dir / paths.logs_dir
      * publish.*（是否生成 Excel、是否拷贝某些可选产物等）

特点：
  - 做“搬运工”而不是“再计算器”，所有内容来自前序模块；
  - 输出结构贴近投稿期刊要求，减少投稿前临时整理时间。

注意事项：
  - 新增下游脚本时，应在本脚本中同步添加需要打包的 TSV；
  - 建议在每次完整分析后跑一遍本脚本，并将 Source Data 目录纳入版本管理。


====================================================================================================
三、config.yaml 关键段与脚本之间的关系（简要）
====================================================================================================

1) project 与 data
   - project.name / species：仅用于日志与 Methods 描述，所有脚本可读取但不依赖其值做逻辑分支。
   - data.samples_tsv / contrasts_tsv：
       * 样本与对比的唯一真实来源；
       * 00/01/04/05/06/08/09 均直接或间接使用。

2) reference 与 salmon
   - reference.*：
       * genome_fa / annotation_gff/gtf：用于 02 提取转录本、03 解析基因、07 对齐 emapper；
       * emapper：07 读取注释的主入口。
   - reference.salmon.*：
       * index_dir / kmer_len / decoy_list / gentrome_fa / rebuild：主要由 02/04 消费。
   - 02_refprep_salmon_index.py 与 04_salmon_quant.py 是对这一段使用最重的两个脚本。

3) resources 与 binaries
   - resources.threads.*：控制 01/02/04/05/06/08 等所有“大头”任务的并行度；
   - resources.memory_gb.*：用于日志与可能的资源检查；
   - binaries.*：fastp/cutadapt/salmon/gffread/Rscript 等各脚本通用。

4) qc
   - 被 01_qc_fastp_cutadapt.py 完整消费：
       * thresholds.* 决定样本是否视为“成功质控”；
       * outlier_rules.* 决定哪些样本被标记为 outlier 而非直接剔除。

5) tximport 与 deg
   - tximport.*：
       * 05_matrix_aggregate.py 与 05_tximport_aggregate.R 使用；
       * 决定 counts_from_abundance 策略和矩阵统计行为。
   - deg.*：
       * 06_deg_module.py 与 06_f_deg.R 使用；
       * 决定 log2FC 阈值、FDR、诊断指标等。

6) annotations 与 emapper
   - annotations.*：
       * 03_build_tx2gene_map.py 与 07_prepare_emapper_annotations.py 共同使用；
       * 统一 ID 清洗、基因/转录本命名规则；
       * 防止“表达矩阵用一套 ID、注释用另一套 ID”的灾难。
   - reference.emapper：
       * 07 的主输入路径。

7) go 与 kegg
   - go.*：
       * go_local.R 与 08_g_enrich.R 使用；
       * 决定 GO 字典路径、OBSOLETE 处理策略、最大 gene set 大小等。
   - kegg.*：
       * kegg_local.sh 与 08_g_enrich.R 使用；
       * 决定 KEGG 字典路径、统计参数（p_adjust_method）等。

8) enrich 与 publish
   - enrich.*：
       * 08_enrich_module.py 与 08_g_enrich.R 使用；
       * 包括 ORA/GSEA 是否启用、universe 规则、输出位置等。
   - publish.*：
       * 09_publish_results.py 使用；
       * 决定是否生成 Excel、是否覆盖已有 Source Data 等行为。

9) paths 与 logging
   - paths.*：
       * 定义各模块结果目录（01_qc / 02_ref / 03_maps / 04_quant / 05_matrix / 06_deg / 07_annot / 08_enrich / 09_publish 及 logs 等）；
       * 所有脚本都是从这里读取路径，而不是内部硬编码。
   - logging.*：
       * 控制所有 Python 脚本的日志输出级别、是否带时间戳、日志文件位置。


====================================================================================================
四、项目目录与结果目录树状结构（简略）
====================================================================================================

以下为一个典型 RNA-seq 工程的目录结构示意，仅展示主要层级：

project/
├── config.yaml
├── README.txt
├── data/
│   ├── samples.tsv
│   └── contrasts.tsv
├── ref/
│   ├── genome.fa
│   ├── genome.gff3 / genome.gtf
│   ├── transcripts.fa
│   ├── decoys.txt
│   ├── gentrome.fa
│   └── salmon_index/
├── results/
│   ├── 01_qc/
│   │   ├── sample_qc.tsv
│   │   ├── summary.tsv
│   │   ├── outliers.tsv
│   │   ├── rejects.tsv
│   │   ├── raw_json/
│   │   ├── raw_html/
│   │   └── clean/
│   ├── 02_ref/
│   ├── 03_maps/
│   │   ├── tx2gene.clean.tsv
│   │   └── tx2gene.blacklist.tsv
│   ├── 04_quant/
│   │   └── <sample>/quant.sf
│   ├── 05_matrix/
│   │   ├── tximport_meta.tsv
│   │   ├── counts/
│   │   │   └── gene_counts.tsv
│   │   ├── tpms/
│   │   │   └── gene_tpm.tsv
│   │   └── matrix_stats.tsv
│   ├── 06_deg/
│   │   └── <contrast>/
│   │       ├── DEG_all.tsv
│   │       ├── DEG_up.tsv
│   │       ├── DEG_down.tsv
│   │       ├── design.txt
│   │       ├── varTopN.list
│   │       └── rle_range.tsv
│   ├── 07_annot/
│   │   ├── gene2go.tsv
│   │   ├── gene2ko.tsv
│   │   ├── gene2pathway.tsv
│   │   ├── universe_coverage.tsv
│   │   ├── go/
│   │   │   ├── term2name.tsv
│   │   │   └── obsolete_map.tsv
│   │   └── kegg/
│   │       ├── term2gene.tsv
│   │       └── term2name.tsv
│   ├── 08_enrich/
│   │   ├── tasks.tsv
│   │   ├── background/
│   │   ├── inputs/
│   │   ├── summary.tsv
│   │   └── <label>/
│   │       ├── GO_*_by_term_*.tsv
│   │       ├── GO_sig_*.tsv
│   │       ├── KEGG_by_term_*.tsv
│   │       └── gsea/
│   └── 09_publish/
│       └── source_data/
│           ├── matrix/
│           ├── deg/
│           ├── background/
│           ├── enrich/
│           ├── METHODS_README.txt
│           ├── manifest.tsv
│           └── publish_check.tsv
└── logs/


====================================================================================================
五、简要排查顺序（仅作提示）
====================================================================================================

1) 先看 01_qc：
   - sample_qc.tsv / summary.tsv 是否合理（总 reads、Q30、GC 分布）；
   - rejects.tsv 中是否出现出乎意料的大量样本。

2) 再看 02_ref 与 03_maps：
   - transcripts.fa 是否存在、长度合理；
   - tx2gene.clean.tsv 行数是否与转录本数大致匹配；
   - 是否存在大量 blacklisted 基因。

3) 检查 04_quant 与 05_matrix：
   - 任意抽查几个样本的 quant.sf 是否有非零 TPM；
   - gene_counts.tsv 中是否存在“全 0 列”或“极端大库”样本；
   - matrix_stats.tsv 中库大小分布是否正常。

4) 差异分析（06）：
   - 某个典型 contrast 的 DEG_all.tsv 中基因数是否与矩阵一致；
   - DEG_up/DEG_down 数量级是否符合生物学预期；
   - RLE 与 varTopN 诊断文件是否提示严重批次问题。

5) 注释与富集（07–08）：
   - universe_coverage.tsv 中 GO/KEGG 覆盖率是否过低；
   - tasks.tsv 是否覆盖了所有需要关注的对比与方向；
   - 若富集结果几乎为空，优先检查输入集合大小与 universe 定义。

6) 发布与装订（09）：
   - publish_check.tsv 中是否有模块被标记为“不完整”；
   - manifest.tsv 中的路径是否都能成功打开；
   - METHODS_README.txt 内容是否满足论文方法学撰写的基本需求。



产物附录

====================================================================================================================
文件/目录                              描述                                                         路径
                                                                                              来源脚本
====================================================================================================================

【00：实验设计与对比定义（可选）】
--------------------------------------------------------------------------------------------------------------------
samples.tsv                        样本设计表（sample、group、FASTQ 等）                                data/samples.tsv
contrasts.tsv                      对比设计表（contrast/case/control 定义）                            data/contrasts.tsv
                                                                                              00_make_sample_contrasts.py

【01：FASTQ 质控与清洗】
--------------------------------------------------------------------------------------------------------------------
sample_qc.tsv                      样本级 QC 表（reads 数、Q30、GC 等）                                 results/01_qc/sample_qc.tsv
summary.tsv                        整体 QC 汇总（与 sample_qc 结构相同）                                 results/01_qc/summary.tsv
outliers.tsv                       边缘可疑样本记录（待人工判断）                                      results/01_qc/outliers.tsv
rejects.tsv                        明确剔除样本列表（Fail-Fast）                                        results/01_qc/rejects.tsv
raw_json/                          fastp JSON 报告目录                                                results/01_qc/raw_json/
raw_html/                          fastp HTML 报告目录                                                results/01_qc/raw_html/
clean/                             清洗后的 FASTQ 目录（若启用 clean）                                   results/01_qc/clean/
                                                                                              01_qc_fastp_cutadapt.py

【02：参考准备与 Salmon 索引】
--------------------------------------------------------------------------------------------------------------------
transcripts.fa                     从基因组+GTF 提取的转录本序列                                          ref/transcripts.fa
decoys.txt                         decoy 染色体/contig ID 列表                                           ref/decoys.txt
gentrome.fa                        转录本+decoy 组合 gentrome 序列                                       ref/gentrome.fa
salmon_index/                      Salmon decoy-aware 索引目录                                           ref/salmon_index/
gentrome_decoy_summary.tsv         gentrome/decoy 摘要统计表                                    results/02_ref/gentrome_decoy_summary.tsv
                                                                                              02_refprep_salmon_index.py

【03：tx2gene 映射】
--------------------------------------------------------------------------------------------------------------------
tx2gene.clean.tsv                  转录本→基因映射表（transcript_id,gene_id）                           results/03_maps/tx2gene.clean.tsv
tx2gene.blacklist.tsv              转录本数超阈值基因黑名单                                             results/03_maps/tx2gene.blacklist.tsv
                                                                                              03_build_tx2gene_map.py

【04：Salmon 表达定量】
--------------------------------------------------------------------------------------------------------------------
04_quant/                          每样本 Salmon 输出总目录                                           results/04_quant/
<sample>/quant.sf                  Salmon 定量结果（TPM/counts/eff_length）                           results/04_quant/<sample>/quant.sf
04_salmon_quant.log                Salmon 定量模块运行日志                                            logs/04_salmon_quant.log
                                                                                              04_salmon_quant.py

【05：基因层表达矩阵（tximport 聚合）】
--------------------------------------------------------------------------------------------------------------------
tximport_meta.tsv                  tximport 元信息表（样本名与 quant.sf 路径）                        results/05_matrix/tximport_meta.tsv
gene_counts.tsv                    基因层 counts 矩阵（行=gene，列=sample）                           results/05_matrix/counts/gene_counts.tsv
gene_tpm.tsv                       基因层 TPM 矩阵                                                   results/05_matrix/tpms/gene_tpm.tsv
matrix_stats.tsv                   表达矩阵整体统计（基因数、库大小等）                              results/05_matrix/matrix_stats.tsv
                                                                                              05_matrix_aggregate.py + 05_tximport_aggregate.R

【06：差异表达分析（DEG）】
--------------------------------------------------------------------------------------------------------------------
06_deg/                            各对比 DEG 结果总目录                                              results/06_deg/
<contrast>/DEG_all.tsv             该对比所有基因 DESeq2 结果                                         results/06_deg/<contrast>/DEG_all.tsv
<contrast>/DEG_up.tsv              显著上调基因列表                                                   results/06_deg/<contrast>/DEG_up.tsv
<contrast>/DEG_down.tsv            显著下调基因列表                                                   results/06_deg/<contrast>/DEG_down.tsv
<contrast>/design.txt              DESeq2 设计式与样本分配记录                                       results/06_deg/<contrast>/design.txt
<contrast>/varTop100.list          最高变异基因 ID 列表（前 100）                                     results/06_deg/<contrast>/varTop100.list
<contrast>/rle_range.tsv           RLE 范围与正规化诊断指标                                           results/06_deg/<contrast>/rle_range.tsv
                                                                                              06_deg_module.py + 06_f_deg.R

【07：注释映射（GO/KEGG）】
--------------------------------------------------------------------------------------------------------------------
gene2go.tsv                        gene → GO ID 映射长表                                              results/07_annot/gene2go.tsv
gene2ko.tsv                        gene → KEGG KO 映射表                                             results/07_annot/gene2ko.tsv
gene2pathway.tsv                   gene → KEGG pathway ID 映射表                                     results/07_annot/gene2pathway.tsv
universe_coverage.tsv              注释覆盖率统计（GO/KO/pathway 覆盖情况）                          results/07_annot/universe_coverage.tsv
                                                                                              07_prepare_emapper_annotations.py

【GO 本地字典（go_local）】
--------------------------------------------------------------------------------------------------------------------
go/term2name.tsv                   GO term 英文名称字典（go_id,term_name）                            results/07_annot/go/term2name.tsv
go/obsolete_map.tsv                GO 废弃 ID 映射表（旧号→新号）                                     results/07_annot/go/obsolete_map.tsv
logs/go_local.log                  GO 本地字典构建日志                                                logs/go_local.log
                                                                                              go_local.R

【KEGG 本地字典（kegg_local）】
--------------------------------------------------------------------------------------------------------------------
kegg/term2gene.tsv                 KEGG pathway → gene 映射长表                                       results/07_annot/kegg/term2gene.tsv
kegg/term2name.tsv                 KEGG pathway ID → 名称字典                                         results/07_annot/kegg/term2name.tsv
                                                                                              kegg_local.sh

【08：富集分析 ORA/GSEA（输入与背景）】
--------------------------------------------------------------------------------------------------------------------
08_enrich/                         富集分析总目录                                                    results/08_enrich/
tasks.tsv                          ORA/GSEA 任务清单（label/test_set 等）                            results/08_enrich/tasks.tsv
background/<c>.list                每对比的背景基因宇宙列表                                          results/08_enrich/background/<contrast>.list
background/<c>.meta.tsv            每对比背景集合元信息                                              results/08_enrich/background/<contrast>.meta.tsv
inputs/*.list                      外部 tag 或 DEG 集合输入基因列表                                 results/08_enrich/inputs/*.list
<lbl>/gsea_ranks.tsv               GSEA 排名文件（若启用 GSEA）                                      results/08_enrich/<label>/gsea_ranks.tsv
                                                                                              08_enrich_module.py

【08：GO/KEGG ORA 结果与 GSEA 表】
--------------------------------------------------------------------------------------------------------------------
summary.tsv                        所有 label 的 ORA/GSEA 显著条目数汇总                             results/08_enrich/summary.tsv
<lbl>/GO_BP_by_term_*.tsv          GO-BP ORA 结果（all/up/down）                                     results/08_enrich/<label>/GO_BP_by_term_*.tsv
<lbl>/GO_CC_by_term_*.tsv          GO-CC ORA 结果（all/up/down）                                     results/08_enrich/<label>/GO_CC_by_term_*.tsv
<lbl>/GO_MF_by_term_*.tsv          GO-MF ORA 结果（all/up/down）                                     results/08_enrich/<label>/GO_MF_by_term_*.tsv
<lbl>/GO_sig_all.tsv               GO 三本体合并显著表（all）                                        results/08_enrich/<label>/GO_sig_all.tsv
<lbl>/GO_sig_up.tsv                GO 显著上调合并表                                                results/08_enrich/<label>/GO_sig_up.tsv
<lbl>/GO_sig_down.tsv              GO 显著下调合并表                                                results/08_enrich/<label>/GO_sig_down.tsv
<lbl>/KEGG_by_term_all.tsv         KEGG ORA 结果（all）                                              results/08_enrich/<label>/KEGG_by_term_all.tsv
<lbl>/KEGG_by_term_up.tsv          KEGG ORA 结果（up）                                               results/08_enrich/<label>/KEGG_by_term_up.tsv
<lbl>/KEGG_by_term_down.tsv        KEGG ORA 结果（down）                                             results/08_enrich/<label>/KEGG_by_term_down.tsv
<lbl>/gsea/GO_gsea.tsv             GO GSEA 结果（若启用 GSEA）                                       results/08_enrich/<label>/gsea/GO_gsea.tsv
<lbl>/gsea/KEGG_gsea.tsv           KEGG GSEA 结果                                                    results/08_enrich/<label>/gsea/KEGG_gsea.tsv
                                                                                              08_g_enrich.R

【09：最终发布与装订（Source Data + Excel）】
--------------------------------------------------------------------------------------------------------------------
source_data/                       面向投稿的汇总目录                                                results/09_publish/source_data/
source_data/matrix/                基因表达矩阵及统计副本                                            results/09_publish/source_data/matrix/
source_data/deg/<c>/               各对比 DEG 结果副本                                               results/09_publish/source_data/deg/<contrast>/
source_data/background/            背景基因 universe 及 meta 副本                                   results/09_publish/source_data/background/
source_data/enrich/<lbl>/          GO/KEGG ORA 与 GSEA 结果副本                                     results/09_publish/source_data/enrich/<label>/
METHODS_README.txt                 自动生成的方法学说明骨架                                          results/09_publish/source_data/METHODS_README.txt
manifest.tsv                       已打包文件清单（path/module/描述）                              results/09_publish/source_data/manifest.tsv
publish_check.tsv                  各模块产物缺失情况检查表                                         results/09_publish/source_data/publish_check.tsv
<contrast>.xlsx                    每对比 DEG+富集汇总 Excel（可选）                               results/09_publish/<contrast>.xlsx
                                                                                              09_publish_results.py
