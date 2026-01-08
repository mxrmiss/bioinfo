====================================================================================================
phylo —— 全链路超详细 README.txt
Version: stable | Format: pure text | Audience: researchers & pipeline users
Author: Li Ziqiang
====================================================================================================


01. 项目简介
----------------------------------------------------------------------------------------------------
phylo 是一个面向多物种系统发育基因组学（phylogenomics）的全自动流程。它从物种基因组
注释文件（GFF/GTF + Genome FASTA）出发，自动提取高质量 CDS 与蛋白序列，运行
OrthoFinder 建立正交家族，构建严格单拷贝（SCO）基因集合，对单基因对齐进行裁剪、
规范化、isoform 选择，拼接 AA 超矩阵，并使用 IQ-TREE 构建系统发育树。

最终，phylo 会将所有标准化产物以统一格式发布为 *aphylo-ready* 套件，直接用于
时间校准、正选择分析、CAFE5 家族扩张/收缩分析等更深层流程。

phylo 完整流程从 00 → 10，所有参数均在 config.yaml 中统一管理，不使用命令行参数。


02. 输入数据规范（Input Contracts）
----------------------------------------------------------------------------------------------------
所有输入数据必须满足以下约定，否则整个流水线可能会提前报错或产生不可复现的结果。

(1) Genome FASTA
    - 目录：data/genomic/
    - 要求：含 >chr/scaffold 的真实序列；seqid 必须与 GFF 内一致（至少 90% 覆盖度）
    - 支持：fa / fna / fasta 及其 .gz 压缩格式

(2) Annotation GFF/GTF
    - 目录：data/annotation/
    - 要求：结构完整的 GFF3/GTF，每条 feature 坐标必须在对应基因组序列范围内

(3) 物种命名
    - 最终所有 FASTA 表头必须统一为：Species|SeqID
    - species.list 将由脚本自动生成

(4) OrthoFinder 输入蛋白
    - 目录：data/proteomes/
    - 必须是规范化后的蛋白序列（脚本自动生成）

(5) CDS 数据
    - 目录：data/cds/
    - 与蛋白序列一一对应，由脚本自动生成

(6) pep2cds 全局映射表
    - data/maps/pep2cds.tsv
    - 由上游脚本自动生成，10_publish_aphylo_ready 的强依赖

(7) 配置文件
    - config.yaml
    - 所有参数必须写在此文件，不允许使用命令行参数


03. 脚本概述（Scripts Overview）
----------------------------------------------------------------------------------------------------
本章节说明所有脚本的用途、特点与注意事项。

----------------------------------------
[上游准备脚本：生成 PEP / CDS / pep2cds]
----------------------------------------

[A] extract_longest_from_gff.py
功能：
  - 以 GFF + Genome 为输入，按物种提取密码子安全（codon-safe）的 CDS 及其对应蛋白。
  - 在密码子安全校验后翻译蛋白，确保 PEP 与 CDS 一一对应。
  - 自动写出 pep2cds 映射：分物种版本 + 全局版本。
特点：
  - 物种并行处理；GFF 预修复；越界坐标过滤；内部“*”过滤；长度过滤。
  - 产生 annotation/genome 匹配表 annot_genome_match.tsv。
注意事项：
  - 覆盖度不达标直接终止；
  - 若 CDS 无法通过密码子安全校验会报错；
  - pep2cds.tsv 必须保留，后续发布硬依赖。

[B] rename_proteome_headers.py
功能：
  - 将所有蛋白 FASTA 表头改为 `Species|SeqID` 标准格式。
  - 可选备份原文件，保证安全可回退。
特点：
  - 自动根据文件名识别物种名；SeqID 清洗规则与 pep2cds 兼容。
注意事项：
  - 请在运行 OrthoFinder 前执行；
  - 不要手动修改表头，以免破坏映射一致性。


----------------------------------------
[主流程脚本：00 → 10]
----------------------------------------

[00] 00_of_run.py
功能：
  - 调用 OrthoFinder，生成全量 Orthogroups、MSA、GeneCount 等。
  - 写出 RESULTS_DIR.txt 作为锚点供后续查找。
注意事项：
  - 运行前必须保证 data/proteomes/*.faa 已标准化表头。

[01] 01_build_sco_keep_list.py
功能：
  - 构建严格 SCO 列表 ogs_selected.list。
特点：
  - 精确逐列扫描 OF 输出，不做推断。
注意事项：
  - ogs_selected.list 是后续所有严格步骤的核心白名单。

[02] 02_prepare_trimmed_alignment.py
功能：
  - 读取 OF 的 MSA，对每个 OG 进行裁剪，输出到 trim_norm 目录。
  - 可生成 colmask 列掩码（纯 0/1，长度=AA 列数）。
注意事项：
  - 裁剪后为空的文件会被忽略。

[03] 03_normalize_msas.py
功能：
  - 将每个裁剪 MSA 的表头统一替换为 Species|SeqID。
  - 写出完整映射表 header_map.tsv。
注意事项：
  - 对于异常表头会严格报错，请确保上游表头已规范。

[04] 04_collapse_isoforms.py
功能：
  - 同一 OG 内，每物种仅保留最长 isoform。
  - 生成 og_species_protein.tsv。
注意事项：
  - 必须在表头规范化后运行，否则会误判物种。

[05] 05_collapse_headers_to_species.py
功能：
  - 将表头从 Species|SeqID 转换为纯物种名（>Species）。
注意事项：
  - 此版本用于拼接超矩阵；是 phylo → aphylo 的接口。

[06] 06_filter_alignments_by_sco.py
功能：
  - 筛选可用于超矩阵的单基因对齐文件，写出 sco_filelist.txt。
注意事项：
  - filelist 使用绝对路径，确保复现稳定性。

[07] 07_concat_supermatrix.py
功能：
  - 基于 sco_filelist.txt 拼接 AA 超矩阵 supermatrix.faa。
  - 生成 partitions.txt（三列表：OG start end）。
特点：
  - 缺序物种在拼接时自动补全为 '-'。
注意事项：
  - 阈值（长度、最小物种数等）全部来自 config.yaml。

[08] 08_iqtree_concat.py
功能：
  - 使用 IQ-TREE 构建系统发育树 species_tree.nwk。
  - 自动注入分区。
注意事项：
  - partitions.txt 不存在时自动退化为单分区模式。

[09] 09_report_matrix.py
功能：
  - 从超矩阵统计每个物种的非 gap 数与占位比例，写出 matrix.tsv。
注意事项：
  - 若 supermatrix.faa 不存在，会报错。

[10] 10_publish_aphylo_ready.py
功能：
  - 将所有中间产物发布为 aphylo-ready 数据套件。
  - 检查 strict SCO 的物种集合是否严格等于物种树叶节点。
  - 构建 pep2cds_resolved.tsv、all_pep2cds_resolved.tsv。
  - 写出 manifest.json 与 QC_report.txt。
  - 严格验证 colmask、前缀白名单、CDS 匹配、ID 规范性。
注意事项：
  - STRICT_NORMALIZE_CDS=True 时，未能匹配的 CDS 会使发布失败。
  - 发布目录会在每次运行前被清空（若配置为 overwrite）。


04. 文件夹结构（Directory Layout）
----------------------------------------------------------------------------------------------------

项目运行结束后，目录结构应类似如下：

phylo/
├── config.yaml
├── data/
│   ├── annotation/
│   ├── genomic/
│   ├── proteomes/
│   ├── cds/
│   └── maps/
│       ├── pep2cds_by_species/
│       └── pep2cds.tsv
├── results/
│   ├── orthofinder/
│   ├── trim_norm/
│   ├── trim_norm_collapse/
│   ├── trim_norm_collapse_sco/
│   ├── supermatrix/
│   ├── trees/
│   ├── reports/
│   └── publish/
│       └── aphylo_ready/
│           ├── strict_sco_msa/
│           ├── colmask/
│           ├── ogs_selected.list
│           ├── sco_filelist.txt
│           ├── ogs_selected.published.list
│           ├── sco_filelist.published.txt
│           ├── family.tsv
│           ├── species_tree.nwk
│           ├── matrix.tsv
│           ├── og_species_protein.tsv
│           ├── pep2cds_resolved.tsv
│           ├── all_pep2cds_resolved.tsv
│           ├── excluded_reason.tsv
│           ├── manifest.json
│           └── QC_report.txt


05. config.yaml 说明（路径与参数总览）
----------------------------------------------------------------------------------------------------
主要字段（仅说明功能，不列 YAML）：

paths:
    input.proteome_dir         OrthoFinder 输入蛋白目录
    input.cds_dir              每物种 CDS 的输出路径
    reports_dir                各类报告（ogs_selected.list、sco_filelist.txt）
    msa_trim_dir               修剪后 MSA 存放
    normalize_dir              标准化后 MSA 存放
    collapse_dir               去同工型后 MSA 存放
    species_collapse_dir       物种层面 collapse 的 MSA（最终版本）
    supermatrix_dir            超矩阵输出目录（supermatrix.faa, partitions.txt）
    trees_dir                  IQ-TREE 输出目录
    publish_dir                发布为 APhylo 的目录

binaries:
    python                     指定 Python 可执行
    iqtree                     IQ-TREE 路径（如留空则自动查找）

species:
    alias_map                  物种名归一化映射，如：
                              Sinonovacula_constricta → Sconstricta

species_tree.concat:
    iqtree_threads             线程数
    iqtree_params              模型、bootstrap、alrt、extra 参数
    min_align_len              拼接前对 OG 对齐长度过滤
    min_taxa_per_locus         每个 OG 至少要有多少物种
    top_n_partitions           最长的若干个 OG 用于拼接（0 表示全部）

publish:
    include_ultrametric_tree   是否复制超时钟树
    db_prefix_allowlist        允许的数据库前缀，用于 pep2cds 检查

06. 可复现性（Reproducibility Contract）
----------------------------------------------------------------------------------------------------
要保证结果百分百复现，必须遵守以下规则：

1）所有参数都写在 config.yaml，不使用命令行参数  
2）所有 FASTA 表头均遵守 Species|SeqID  
3）RESULTS_DIR.txt 是 OrthoFinder 结果唯一入口，禁止手动修改  
4）不要手动改写 ogs_selected.list、sco_filelist.txt  
5）所有脚本都不应被修改，否则需重新发布 aphylo_ready  
6）发布脚本 10_publish_aphylo_ready.py 是最终入口，必须保持最新  
7）每次变动 PEP、CDS、pep2cds.tsv 都必须完整重跑 00–10  


07. 常见问题（FAQ）
----------------------------------------------------------------------------------------------------

Q1：为什么某些 OG 最终没有进入 strict_sco_msa？  
A1：因为严格要求“物种集合必须等于物种树叶节点”，缺任意一个物种都会被剔除。

Q2：为什么 colmask 会导致发布失败？  
A2：掩码长度必须精确等于 AA 列数，否则视为数据损坏。

Q3：为什么 pep2cds_unresolved.tsv 会出现大量条目？  
A3：多半是 CDS header 或 pep2cds.tsv 未统一；需检查 ID 前缀、alias_map 等。

Q4：是否可以只重跑 IQ-TREE？  
A4：可以，只要 supermatrix.faa 和 partitions.txt 没变。

Q5：OrthoFinder 结果能否复用？  
A5：可复用，只要蛋白表头与物种名不变，否则应重跑。



产物附录

====================================================================================================================
文件/目录                                 描述                                                         路径
                                                                                                       来源脚本
====================================================================================================================

【00：物种原始准备（蛋白/CDS/pep2cds）】
--------------------------------------------------------------------------------------------------------------------
data/proteomes/*.faa                      规范化后的物种蛋白序列（Species|SeqID 表头）                data/proteomes/
                                           用作 OrthoFinder 输入                                       extract_longest_from_gff.py
                                                                                                       rename_proteome_headers.py

data/cds/*.fna                             规范化 CDS，与 PEP 一一对应                                data/cds/
                                           后续 pep→cds 映射所需                                      extract_longest_from_gff.py

data/maps/pep2cds.tsv                      全局蛋白→CDS 对照表                                        data/maps/
                                           用于 strict_sco 发布与 codon 追踪                          extract_longest_from_gff.py

data/maps/pep2cds_by_species/*.tsv         每物种拆分的蛋白→CDS 映射                                  data/maps/pep2cds_by_species/
                                           发布阶段用于分物种解析                                     extract_longest_from_gff.py

results/proteomes/species.list             当前成功获得 PEP/CDS 的物种列表                            results/proteomes/
                                           用作 phylo 有效物种清单                                    extract_longest_from_gff.py

results/proteomes/annot_genome_match.tsv   GFF vs Genome 覆盖度/状态表                                results/proteomes/
                                           物种 QC 情况                                               extract_longest_from_gff.py


====================================================================================================================
【01：OrthoFinder 运行】
--------------------------------------------------------------------------------------------------------------------
results/orthofinder/Results_*/             OrthoFinder 主结果目录                                     results/orthofinder/
                                           包含 Orthogroups、MSA、GT、基因计数等                       00_of_run.py

results/orthofinder/RESULTS_DIR.txt        唯一锚点文件，指向有效 Results_*                           results/orthofinder/
                                           全后续流程依赖它定位 OF 输出                               00_of_run.py

Results_*/Orthogroups/Orthogroups.tsv      全家族矩阵                                                 results/orthofinder/Results_*/Orthogroups/
                                           后续严格 SCO 构建来源之一                                  OrthoFinder

Results_*/Orthogroups/Orthogroups.GeneCount.tsv  家族计数矩阵（给 CAFE）                             results/orthofinder/Results_*/Orthogroups/
                                                                                                       OrthoFinder

Results_*/MultipleSequenceAlignments/*.fa  单基因对齐（AA）                                          results/orthofinder/Results_*/MultipleSequenceAlignments/
                                           后续裁剪处理                                               OrthoFinder


====================================================================================================================
【02：严格 SCO 构建】
--------------------------------------------------------------------------------------------------------------------
results/reports/ogs_selected.list          严格 single-copy orthogroup 列表                            results/reports/
                                           每行一个 OG                                                01_build_sco_keep_list.py

results/reports/ogs_policy.json            SCO 选择策略说明                                           results/reports/
                                           记录使用了哪种来源（txt 或 column-scan）                    01_build_sco_keep_list.py

results/reports/.ogs_selected.done         严格 SCO 完成哨兵                                          results/reports/
                                           用于 Snakemake / 逻辑依赖                                   01_build_sco_keep_list.py


====================================================================================================================
【03：对齐裁剪（Trimmed MSA）】
--------------------------------------------------------------------------------------------------------------------
results/trim_norm/OG*.trim.faa             裁剪后的 MSA 文件（AA）                                     results/trim_norm/
                                           对应 OF MSA → trimal 处理                                  02_prepare_trimmed_alignment.py

results/colmask/OG*.colmask                对应每个裁剪 MSA 的掩码（0/1 串）                           results/colmask/
                                           长度=对齐列数                                              02_prepare_trimmed_alignment.py

results/trim_norm/.done                    裁剪完成哨兵                                               results/trim_norm/
                                           标记 OG*.trim.faa 均生成                                   02_prepare_trimmed_alignment.py


====================================================================================================================
【04：MSA 表头规范化（Species|SeqID）】
--------------------------------------------------------------------------------------------------------------------
results/trim_norm/OG*.trim.faa             表头已规范为 Species|SeqID                                 results/trim_norm/
                                           保证后续 isoform 折叠与拼接统一                             03_normalize_msas.py

results/reports/header_map.tsv             原表头→新表头的完整映射表                                  results/reports/
                                           记录 OG、旧 header、新 header、Species、SeqID               03_normalize_msas.py

results/trim_norm/.norm.done               表头规范化完成哨兵                                         results/trim_norm/
                                           后续折叠流程触发                                           03_normalize_msas.py


====================================================================================================================
【05：isoform 折叠（保留最长）】
--------------------------------------------------------------------------------------------------------------------
results/trim_norm_collapse/OG*.trim.faa    同一 OG、同一物种仅保留最长 isoform                        results/trim_norm_collapse/
                                           仍保留 Species|SeqID 表头                                  04_collapse_isoforms.py

results/reports/og_species_protein.tsv     (OG, Species, SeqID, length) 选择记录                      results/reports/
                                           供后续 strict 发布与 QC                                    04_collapse_isoforms.py

results/trim_norm_collapse/.done           isoform 折叠完成哨兵                                       results/trim_norm_collapse/
                                           下一阶段 species-only 使用                                 04_collapse_isoforms.py


====================================================================================================================
【06：物种-only 表头（>Species）】
--------------------------------------------------------------------------------------------------------------------
results/trim_norm_collapse_sco/OG*.trim.faa   表头仅物种名 >Species                                  results/trim_norm_collapse_sco/
                                              给超矩阵拼接、IQ-TREE、aphylo 用                       05_collapse_headers_to_species.py

results/trim_norm_collapse_sco/.done          species-only 转换完成哨兵                               results/trim_norm_collapse_sco/
                                              后续 filelist 专用                                      05_collapse_headers_to_species.py


====================================================================================================================
【07：严格 SCO filelist（绝对路径）】
--------------------------------------------------------------------------------------------------------------------
results/reports/sco_filelist.txt           每行一个绝对路径，对应 strict MSA                          results/reports/
                                           超矩阵拼接唯一入口                                         06_filter_alignments_by_sco.py

results/reports/.sco_filelist.done         filelist 完成哨兵                                         results/reports/
                                           供后续超矩阵流程检查                                       06_filter_alignments_by_sco.py


====================================================================================================================
【08：超矩阵拼接（AA Supermatrix）】
--------------------------------------------------------------------------------------------------------------------
results/supermatrix/supermatrix.faa       所有保留的 OG 拼接得到的 AA 超矩阵                         results/supermatrix/
                                          物种顺序统一、缺失位点填 '-'                                07_concat_supermatrix.py

results/supermatrix/partitions.txt        (OG, start, end) 三列表分区                                results/supermatrix/
                                          1-based 闭区间                                              07_concat_supermatrix.py

results/supermatrix/.supermatrix.done     超矩阵完成哨兵                                             results/supermatrix/
                                          IQ-TREE 需要它                                               07_concat_supermatrix.py


====================================================================================================================
【09：覆盖度报告】
--------------------------------------------------------------------------------------------------------------------
results/reports/matrix.tsv                覆盖度矩阵：total_sites / nongap_sites / occupancy%        results/reports/
                                          供 Methods QC 与发布                                       09_report_matrix.py

results/reports/.matrix.done              覆盖度完成哨兵                                             results/reports/
                                          标记 matrix.tsv 存在                                       09_report_matrix.py


====================================================================================================================
【10：IQ-TREE 构树】
--------------------------------------------------------------------------------------------------------------------
results/trees/species_tree.nwk            主物种树（AA 超矩阵构建）                                  results/trees/
                                          发布到 aphylo 的默认树                                     08_iqtree_concat.py

results/trees/species_tree.*              IQ-TREE 全输出（log, iqtree, treefile, etc）               results/trees/
                                          方便审计                                                    08_iqtree_concat.py

results/trees/.iqtree.done                IQ-TREE 完成哨兵                                           results/trees/
                                          标记树已构建                                               08_iqtree_concat.py


====================================================================================================================
【11：最终发布（Aphylo-ready 套件）】
--------------------------------------------------------------------------------------------------------------------
results/publish/aphylo_ready/strict_sco_msa/OG*.faa   严格发布用物种-only MSA                       results/publish/aphylo_ready/
                                                     完全匹配树叶节点                                 10_publish_aphylo_ready.py

results/publish/aphylo_ready/colmask/*                与 strict MSA 对应的掩码                        results/publish/aphylo_ready/
                                                     （若存在，则 1:1 匹配）                         10_publish_aphylo_ready.py

results/publish/aphylo_ready/species_tree.nwk         发布用物种树                                    results/publish/aphylo_ready/
                                                     aphylo / CAFE / codon downstream                 10_publish_aphylo_ready.py

results/publish/aphylo_ready/family.tsv               发布用家族矩阵（来自 OF）                       results/publish/aphylo_ready/
                                                     按叶顺序重排                                     10_publish_aphylo_ready.py

results/publish/aphylo_ready/pep2cds_resolved.tsv     严格 SCO 的蛋白→真实 CDS 映射                   results/publish/aphylo_ready/
                                                     完全匹配 fasta header                            10_publish_aphylo_ready.py

results/publish/aphylo_ready/all_pep2cds_resolved.tsv 全家族范围的蛋白→CDS 映射                      results/publish/aphylo_ready/
                                                     CAFE/PSG 联合用                                   10_publish_aphylo_ready.py

results/publish/aphylo_ready/pep2cds_unresolved.tsv   未匹配到 CDS 的蛋白列表                         results/publish/aphylo_ready/
                                                     QC 重要参考                                      10_publish_aphylo_ready.py

results/publish/aphylo_ready/ogs_selected.list        发布用 strict SCO 列表                           results/publish/aphylo_ready/
                                                     与原始列表同步                                    10_publish_aphylo_ready.py

results/publish/aphylo_ready/sco_filelist.txt         发布用 strict MSA 文件名列表                    results/publish/aphylo_ready/
                                                     （相对名，不含路径）                             10_publish_aphylo_ready.py

results/publish/aphylo_ready/ogs_selected.published.list   严格发布成功的 OG 列表                     results/publish/aphylo_ready/
                                                           排除物种集合不一致等情况                   10_publish_aphylo_ready.py

results/publish/aphylo_ready/sco_filelist.published.txt   与上方对应的 MSA 文件名                     results/publish/aphylo_ready/
                                                           最终上架集合                               10_publish_aphylo_ready.py

results/publish/aphylo_ready/og_species_protein.tsv   isoform 选择记录（折叠阶段产生）                results/publish/aphylo_ready/
                                                     发布时一并拷贝                                  10_publish_aphylo_ready.py

results/publish/aphylo_ready/matrix.tsv              覆盖度矩阵（report 结果）                       results/publish/aphylo_ready/
                                                     用于 QC 与 Methods                               10_publish_aphylo_ready.py

results/publish/aphylo_ready/excluded_reason.tsv     发布时被排除 OG 的原因                            results/publish/aphylo_ready/
                                                    （缺物种/空对齐/对应掩码缺失等）                 10_publish_aphylo_ready.py

results/publish/aphylo_ready/manifest.json           发布物元数据（机器可读）                          results/publish/aphylo_ready/
                                                     包含路径/统计/QC 信息                            10_publish_aphylo_ready.py

results/publish/aphylo_ready/QC_report.txt           人类可读 QC 报告                                  results/publish/aphylo_ready/
                                                     汇总 strict_sco、colmask、pep2cds 等              10_publish_aphylo_ready.py

results/publish/aphylo_ready/.done                   最终发布哨兵                                      results/publish/aphylo_ready/
                                                     标志 aphylo-ready 完整生成                        10_publish_aphylo_ready.py

====================================================================================================================


====================================================================================================
（END OF README.txt）
====================================================================================================