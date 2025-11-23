====================================================================================================
APHYLO PIPELINE —— 超详细说明文档README.txt
Version: stable | Format: pure text | Audience: researchers & pipeline users
Author: Li Ziqiang
====================================================================================================

本 README 说明当前 aphylo 项目中所有脚本的用途、输入依赖、配置字段与目录结构。

====================================================================================================
项目简介
====================================================================================================
aphylo假定已经有一个来自 phylo 流水线的“发布包接口”，包括：
  - 严格 SCO 的 AA MSA 与列掩码；
  - 统一的 pep→cds 映射表；
  - 物种树（用于 MCMCTree / codeml / CAFE）；
  - 基因家族矩阵 family.tsv 等。

在此基础上，aphylo 负责完成：
  1) 从发布包接口出发，构建 codon MSA 与 QC；
  2) codeml branch-site 正选择分析（PSG）及 PSG 富集输入；
  3) MCMCTree 两阶段运行 + 事后 ESS / finetune 建议；
  4) CAFE5 基因家族扩张/收缩分析与 per-species 富集输入；
  5) PSG × CAFE 联合整合；
  6) 最终报告目录搭建（统一拷贝关键 TSV）。


====================================================================================================
一、整体流程与推荐运行顺序
====================================================================================================

推荐严格顺序如下（括号内为脚本文件名）：

1) 接口体检与 codon 输入准备
   - 01_iface_check_publish.py
   - 02_build_codon_inputs.py
   - 03_pal2nal_like_and_colmask_sync.py
   - 04_qc_codon_alignments.py

2) PSG 模块
   - 05_define_foregrounds.py
   - 06_codeml_batch_branchsite.py
   - 07_codeml_aggregate.py

3) MCMCTree 模块
   - 08_clean_tree_and_noblize.py
   - 09_build_concat_clean_phy.py
   - 10_mcmctree_run_and_publish.py
   - mcmctree_postcheck.py（内部驱动 ess_report.R 和 finetune_suggest.py）

4) CAFE 模块
   - 11_cafe5_prepare_inputs.py
   - 12_cafe5_run_models.py
   - 13_cafe5_aggregate.py
   - 13b_cafe_build_enrich_inputs.py

5) 联合整合与报告
   - 14_joint_integration.py
   - 15_build_final_report.py


====================================================================================================
二、脚本逐个说明（功能 / 主要输入 / 特点 / 注意事项）
====================================================================================================

----------------------------------------
1) 01_iface_check_publish.py
----------------------------------------
功能：
  - 对 phylo 发布包接口进行完整体检，确认：
      * 物种树中的叶节点集合；
      * 严格 SCO AA-MSA 目录中的物种表头；
      * pep→cds 映射表中出现的物种集合；
    三者在物种集合与命名口径上完全一致。
  - 记录 publish 侧关键文件与 aphylo 所需接口的映射关系。

主要输入：
  - config.yaml：
      * publish_dir                       —— phylo 发布包根目录
      * inputs.sco_msa_dir / sco_msa_suffix
      * inputs.pep2cds_map
      * inputs.species_tree
      * inputs.family 或 inputs.family_tsv
      * paths.iface_dir
      * paths.logs_dir

特点：
  - 接口层的“健康检查”，不生成下游使用的新数据，只判断是否可以安全进入 aphylo。
  - 对不一致的地方给出明确条目（缺哪个物种、重复哪个 ID）。

注意事项：
  - 每次更新 phylo 流程输出或物种集合，都建议重新跑一次本脚本。
  - 如果 iface 检查不过，应先在 phylo/publish 侧修正，而不是在 aphylo 中硬适配。


----------------------------------------
2) 02_build_codon_inputs.py
----------------------------------------
功能：
  - 从严格 SCO AA-MSA 和 pep→cds 映射出发，为每个 OG 构建“已排好物种顺序的 CDS 输入”。
  - 统一清洗 ID（strip 前后缀、统一大小写、应用 alias_map），保证 AA 与 CDS 完全对应。

主要输入：
  - config.yaml：
      * publish_dir
      * inputs.sco_msa_dir / sco_msa_suffix      —— 严格 SCO AA-MSA
      * inputs.pep2cds_map                       —— pep→cds 映射表
      * inputs.cds_dir / cds_suffix              —— 按物种组织的 CDS FASTA 仓库
      * species.alias_map                        —— 物种别名映射
      * paths.bt_dir                             —— 作为本脚本的工作根目录
      * paths.logs_dir

特点：
  - 用正则从 MSA 文件名中提取 OGID（如 OG00001），杜绝双后缀误识别。
  - 所有路径写在 config.yaml，脚本内不硬编码 project 绝对路径。

注意事项：
  - 若某个 (OG, Species) 无法在 pep2cds_map 中找到对应 CDS，会直接报错并指明具体缺失条目。
  - 建议接口检查（01）通过之后再运行本脚本。


----------------------------------------
3) 03_pal2nal_like_and_colmask_sync.py
----------------------------------------
功能：
  - 结合 AA-MSA 与 02 生成的 ordered_cds，将氨基酸对齐回译为 codon MSA。
  - 同步列掩码：读取 AA 层面的 colmask，将其展开为核酸层面的 ntmask。

主要输入：
  - config.yaml：
      * publish_dir
      * inputs.sco_msa_dir / sco_msa_suffix      —— 严格 SCO AA-MSA
      * inputs.colmask_dir / colmask_suffix      —— AA 层面掩码（如有）
      * paths.bt_dir                             —— 02 步输出目录
      * paths.logs_dir
  - 02_build_codon_inputs.py 生成的 ordered CDS FASTA 与 order 表。

特点：
  - 通过统一的 OG 列表驱动回译，对每个 OG 只处理一次。
  - 若 AA 端存在列掩码，则严格按列同步到 nt 端，保证 mask 与 codon 对齐位置一一对应。

注意事项：
  - 若没有 colmask 目录，只生成 codon MSA，ntmask 相关部分自动跳过。
  - 建议抽查若干 OG，确认 codon MSA 与 AA-MSA 在物种顺序和长度上的一致性。


----------------------------------------
4) 04_qc_codon_alignments.py
----------------------------------------
功能：
  - 对 03 生成的 codon MSA 做质量控制，按照 config.codon 中的规则给每个 OG 打分与筛选。
  - 将“是否保留”与“丢弃原因”持久化，供后续 PSG 和 MCMCTree 使用同一 keep 列表。

主要输入：
  - config.yaml：
      * publish_dir
      * paths.codon_qc_dir
      * codon.min_taxa / min_codon_len
      * codon.max_gap_frac
      * codon.stopcodon_action
  - codon MSA 目录：
      * results/03_codon/codon_msa/OG*.codon.fna

特点：
  - 数量型阈值全部从 config 读取，便于不同项目之间统一调整。
  - 同时统计 taxa 数、有效长度、gap 比例、内部终止密码子数等多个指标，而不是单一标准。

注意事项：
  - stopcodon_action 建议先设为 “mask” 观察整体终止分布，再决定是否改为 “fail”。
  - PSG 和 MCMCTree 的信息量直接受 keep_og 数量影响，过严的 QC 会让后续分析变得过稀。


----------------------------------------
5) 05_define_foregrounds.py
----------------------------------------
功能：
  - 根据 config.codeml.foreground_sets 中的定义，在物种树上生成对应的前景集合。
  - 每个前景集合产生一个物种列表文件和一个打 #1 的树文件。

主要输入：
  - config.yaml：
      * publish_dir
      * inputs.species_tree
      * codeml.foreground_sets
      * codeml.add_phylip_header
      * species.alias_map
      * paths.codeml_dir

特点：
  - 当前模式为 terminals-only，前景只由叶节点组合而成，不自动构造内部节点前景。
  - 允许配置多个前景集合，便于针对不同生态或表型问题构建不同 PSG 对比。

注意事项：
  - foreground_sets 中的物种名要么直接与树叶名一致，要么在 alias_map 中有映射。
  - 若后续 codeml 报错“树上没有前景分支”，优先检查 .list 与 .nwk 中物种拼写。


----------------------------------------
6) 06_codeml_batch_branchsite.py
----------------------------------------
功能：
  - 在 QC 合格的 codon MSA 上，按前景集合批量运行 codeml branch-site 模型。
  - 对每个 OG × 前景，运行 ALT 和 NULL 两个模型，并在原始目录中保存 ctl / mlc 等。

主要输入：
  - config.yaml：
      * publish_dir
      * paths.codeml_raw_dir
      * paths.codeml_dir（sets/ 子目录下提供前景树）
      * codeml.threads
      * codeml.alt_template / null_template
      * binaries.codeml
  - 04 步输出的 keep_og.list
  - 03 步生成的 codon MSA

特点：
  - 所有运行线程数统一受 codeml.threads 控制，便于在本地/HPC 环境中限速。
  - 对单个 OG 的失败具有容错机制：记录失败但尽量不影响其它 OG。

注意事项：
  - 模板 ctl 文件路径和格式必须与脚本逻辑一致，尤其是 seqfile / treefile 的占位符命名。
  - 强烈建议先用少量 OG 做试运行，确认 mlc 中各字段能被 07 正确解析。


----------------------------------------
7) 07_codeml_aggregate.py
----------------------------------------
功能：
  - 聚合 06 的 branch-site ALT/NULL 结果，完成 PSG 判定与 FDR 校正。
  - 视配置生成 PSG 富集分析所需的 test/background 列表（CDS 层级）。

主要输入：
  - config.yaml：
      * publish_dir
      * paths.codeml_raw_dir
      * paths.codeml_agg_dir
      * inputs.pep2cds_map
      * inputs.all_members_tsv
      * psg.enable
      * psg.fdr_alpha
      * psg.cds_inputs_dir
  - 06 步 raw 目录结构下的 mlc 输出。

特点：
  - 采用混合 χ²₁ 分布和 BH-FDR，统计口径与主流 PSG 文献一致。
  - 使用 pep2cds_map 与 all_members_tsv 保证“OG → gene → cds” 映射过程透明。

注意事项：
  - 若之后要在 transcript 层级进行富集，应在 all_members_tsv 端先行汇总/映射，再统一到 cds 层级。
  - 如果某些前景集合理论上应有较多 PSG 却结果几乎为空，可回溯检查前景定义和 codon MSA 质量。


----------------------------------------
8) 08_clean_tree_and_noblize.py
----------------------------------------
功能：
  - 对 phylo 输出的物种树做定根和清洗，生成 MCMCTree 需要的两行树文件（有长度和无长度）。

主要输入：
  - config.yaml：
      * mcmctree.tree_in
      * mcmctree.out_tree
      * mcmctree.out_nobl
      * mcmctree.outgroup
      * mcmctree.outdir
      * paths.logs_dir
  - gotree 可执行程序（用于定根和简化树）。

特点：
  - 把分支支持值、注释等非必须信息全部清除，只保留分支长度与拓扑结构。
  - 可选在首行写入 “N 1”，让 MCMCTree 直接读取。

注意事项：
  - 外群物种列表必须是真实存在于树叶中的物种，否则 gotree 会失败。
  - 之后所有下游步骤都以 .nobl.trees 中的叶合集为准，不建议在此之后再改树。


----------------------------------------
9) 09_build_concat_clean_phy.py
----------------------------------------
功能：
  - 根据 QC 后的 codon MSA 和两行树叶序，生成 MCMCTree 需要的 concat.clean.phy。

主要输入：
  - config.yaml：
      * mcmctree.codon_dir
      * mcmctree.keep_list
      * mcmctree.codon_pos
      * mcmctree.nobl_trees
      * mcmctree.work_dir
      * paths.logs_dir
  - 来自 03 与 04 的 codon MSA 及 keep_og 列表。

特点：
  - 支持只使用第一、二密码子（codon_pos = "12"），减少同义位点影响。
  - 按树叶顺序排列物种行，避免 MCMCTree 再次 reorder。

注意事项：
  - 若合并时发现某个 OG 缺少特定物种，会在日志中给出明确提示，应结合 03 与 04 的日志排查。
  - concat.clean.phy 生成后建议用简单脚本查看前几行与树叶顺序是否匹配。


----------------------------------------
10) 10_mcmctree_run_and_publish.py
----------------------------------------
功能：
  - 在 mcmctree.work_dir 中运行 MCMCTree 两个阶段，并根据 out.txt / FigTree.tre 发布给 CAFE 使用的超时钟树。

主要输入：
  - config.yaml：
      * mcmctree.work_dir
      * mcmctree.seqfile
      * mcmctree.treefile
      * mcmctree.publish_ultrametric_tree
      * binaries.mcmctree
      * paths.logs_dir
  - work_dir 中事先放好的 mcmctree.ctl 与 mcmctree2.ctl 模板。

特点：
  - 通过 .stage1.done / .stage2.done / .stage3.done 记录执行进度，可以断点续跑。
  - Stage3 只负责从 FigTree.tre 提取并清洗 ultrametric_tree，便于 CAFE5 消费。

注意事项：
  - 运行前应确保 ctl 文件中的 seqfile / treefile 与 09 输出一致。
  - 若 out.txt 或 mcmc.txt 很短，多半是 finetune 或步长设置不合理，应结合 postcheck 再调整。


----------------------------------------
11) ess_report.R
----------------------------------------
功能：
  - 从 MCMCTree 的 mcmc.txt 中计算 ESS，并给出数值化建议。

主要输入：
  - WORK_DIR / QC_SUBDIR / MCMC_PATH（硬编码为 results/06_cafe/mcmctree 及其 qc_report 子目录）。

特点：
  - 只依赖 data.table 与 coda，方便在常见 R 环境中部署。
  - 支持使用 IPS 与 coda::effectiveSize 双重估计 ESS。

注意事项：
  - 安装 R 包时需要注意 R 版本与 Bioconductor 版本兼容。
  - KEEP_REGEX 如需调整，应直接在 ess_report.R 中修改正则表达式。


----------------------------------------
12) finetune_suggest.py
----------------------------------------
功能：
  - 读取 ESS 报告，根据最小 ESS 推断合理的 finetune 调整建议。

主要输入：
  - config.yaml：
      * mcmctree.work_dir
  - results/06_cafe/mcmctree/qc_report/ess.tsv
  - mcmctree.ctl（在 work_dir 中）

特点：
  - 默认只生成 finetune_suggestion.md 与 finetune_new_line.txt，不自动修改 ctl。
  - 能够根据 ESS 范围决定是整体放大还是缩小步长。

注意事项：
  - FORCE_CTL_NAME / AUTO_PATCH_CTL 等参数若被开启，会直接修改 ctl，使用前需确认是否接受自动改写。
  - 建议在接受新的 finetune 之前，先在短链测试中验证 ESS 改善情况。


----------------------------------------
13) mcmctree_postcheck.py
----------------------------------------
功能：
  - 统一执行 MCMCTree 的事后体检：
      * 从 out.txt 中提取节点年龄；
      * 调用 ess_report.R 计算 ESS；
      * 可选调用 finetune_suggest.py；
      * 生成汇总说明。

主要输入：
  - config.yaml（自动寻找 mcmctree.work_dir）
  - work_dir 中的 out.txt、mcmc.txt

特点：
  - 若没有显式提供 config.yaml 路径，会尝试在当前目录及父目录自动搜寻。
  - 各步骤出现错误时给出清晰的消息，而不是悄悄跳过。

注意事项：
  - 建议每次完整 MCMC 结束后都跑一次 postcheck，以确定 ESS 是否达标。


----------------------------------------
14) 11_cafe5_prepare_inputs.py
----------------------------------------
功能：
  - 将 phylo/aphylo 发布的 family 矩阵与 ultrametric_tree 规范化为 CAFE5 可直接使用的形式。

主要输入：
  - config.yaml：
      * publish_dir
      * inputs.family 或 inputs.family_tsv
      * inputs.ultrametric_tree
      * species.alias_map
      * paths.cafe_run_dir
      * paths.logs_dir

特点：
  - 使用 alias_map 校正物种名，并按 ultrametric_tree 的叶序重新排列列顺序。
  - 删除 Total 和非物种列，只保留真实物种计数列。

注意事项：
  - 若家族矩阵中存在树上没有的物种，会被删除；会影响可用物种的数量，需要结合 phylo 侧一起判断是否合理。
  - family 表的分隔符、表头命名需与脚本预期一致（例如第一列为家族 ID）。


----------------------------------------
15) 12_cafe5_run_models.py
----------------------------------------
功能：
  - 基于 family.tsv 和 ultrametric_tree，按 config.cafe5 中定义的多个模型运行 CAFE5。

主要输入：
  - config.yaml：
      * publish_dir
      * paths.cafe_run_dir
      * paths.logs_dir
      * cafe5.threads / models / k_cycles / p_alpha / max_autofix_rounds
      * cafe5.two_stage_large.*
      * cafe5.error_model.*
      * binaries.cafe5
  - 11 步准备的 family.tsv 与 cleaned_ultrametric_tree.nwk

特点：
  - 对 primary 集合支持多轮“极端 OG 自动剔除 + 重跑”，保证整体稳定。
  - 对 large 集合只进行单轮运行，节约资源。

注意事项：
  - 若 CAFE5 报错，优先查看模型子目录下的 run.log 或 stdout/stderr 文件。
  - 在调整模型参数（如 error_model）时，建议记录在 config.yaml 中以便追溯。


----------------------------------------
16) 13_cafe5_aggregate.py
----------------------------------------
功能：
  - 读取 12 步输出，完成家族 FDR 纠正与 λ / clade 统计整合。

主要输入：
  - config.yaml：
      * paths.cafe_run_dir
      * paths.cafe_agg_dir
      * report.fdr_alpha
      * cafe5.models
  - models/<model>/primary_global/ 与 large/ 内的 family/clade 结果文件。

特点：
  - FDR 阈值全局统一由 report.fdr_alpha 控制，避免不同脚本各自设置。
  - 兼容 family_results / clade_results 多种略有差异的列名。

注意事项：
  - 若 family 结果文件没有 pvalue 列，脚本会直接报错，需要检查 CAFE5 调用参数是否一致。
  - 对 high_fail_ogs.list 的处理会影响 no_highfail 版本的显著家族集合。


----------------------------------------
17) 13b_cafe_build_enrich_inputs.py
----------------------------------------
功能：
  - 将 CAFE 显著扩张/收缩家族映射到 per-species 基因集合，生成 GO/KEGG 富集分析输入。

主要输入：
  - config.yaml：
      * publish_dir
      * paths.cafe_agg_dir
      * cafe5.enrich_bridge.enable
      * cafe5.enrich_bridge.outputs_dir
      * cafe5.enrich_bridge.species_sets
      * report.fdr_alpha
      * inputs.all_members_tsv
  - 13 步的 cafe_significant_families_no_highfail.tsv
  - CAFE Base_clade_results, Base_change.tab, Base_branch_probabilities.tab

特点：
  - species_sets 配置可定义多种物种组合（单物种或多个物种合并）。
  - 使用 all_members_tsv 保证从 family 到 gene/cds 的映射与 PSG 一致。

注意事项：
  - 若 enable = false，本脚本不会生成任何富集输入，这在 config 中需明确选择。
  - species_sets 中物种名必须是“已经 alias 过的标准名”，否则会识别不到。


----------------------------------------
18) 14_joint_integration.py
----------------------------------------
功能：
  - 在家族/OG 层面把 PSG 与 CAFE 显著结果做交集与并集，形成一份联合视角的列表。

主要输入：
  - config.yaml：
      * paths.codeml_agg_dir
      * paths.cafe_agg_dir
      * paths.joint_dir
      * report.fdr_alpha
  - D_fdr_genes.tsv（来自 07）
  - cafe_significant_families.tsv（来自 13）

特点：
  - 只基于 significant 标记与 FDR 阈值做集合运算，不改变各自原始统计量。
  - 统计结果 integration_counts.tsv 方便快速了解 PSG/CAFE 与其交集的规模。

注意事项：
  - 若某一方没有任何显著条目（例如 PSG 全为 NS），integration_intersect.tsv 自然会为空，这通常是 upstream 的问题。


----------------------------------------
19) 15_build_final_report.py
----------------------------------------
功能：
  - 把 PSG、CAFE 与联合作用的关键表格集中复制到 report 目录，为后续制图与论文撰写提供统一入口。

主要输入：
  - config.yaml：
      * paths.codeml_agg_dir
      * paths.cafe_agg_dir
      * paths.joint_dir
      * paths.logs_dir
  - 07、13、13b、14 产生的各种 TSV 文件。

特点：
  - 纯“搬运工”，不修改任何表格内容。
  - 在 report/README.txt 中通过简单文字说明各表的来源与用途。

注意事项：
  - 若新增脚本产出新的关键表格，需要同步在 15 脚本中增加拷贝逻辑。
  - .report.done 可作为整套 aphylo 流程是否完整完成的总哨兵。


====================================================================================================
三、config.yaml 关键段与脚本之间的关系（简要）
====================================================================================================

1) publish_dir
   - 所有带 "<publish_dir>" 占位符的路径，都会在脚本中通过统一函数展开。
   - 与 phylo 的项目根目录一致时，接口最稳定。

2) inputs
   - sco_msa_dir / sco_msa_suffix       —— 01 / 02 / 03
   - colmask_dir / colmask_suffix       —— 03
   - pep2cds_map                        —— 01 / 02 / 07 / 13b
   - all_members_tsv                    —— 07 / 13b
   - species_tree                       —— 01 / 05 / 08
   - family / family_tsv                —— 01 / 11
   - ultrametric_tree                   —— 11（由 10 发布）
   - cds_dir / cds_suffix               —— 02

3) paths
   - iface_dir       —— 01
   - bt_dir          —— 02 / 03
   - codon_qc_dir    —— 04
   - codeml_dir      —— 05 / 06
   - codeml_raw_dir  —— 06
   - codeml_agg_dir  —— 07 / 14 / 15
   - cafe_run_dir    —— 11 / 12 / 10
   - cafe_agg_dir    —— 13 / 13b / 14 / 15
   - joint_dir       —— 14 / 15
   - logs_dir        —— 所有脚本
   - debug_dir       —— 10 可能使用

4) codon
   - 被 04 完整消费，用于 codon MSA QC。

5) codeml
   - 05 使用 foreground_sets；
   - 06 使用 alt_template / null_template / threads 等。

6) psg
   - 07 使用 enable / fdr_alpha / cds_inputs_dir。

7) mcmctree
   - 08 使用 tree_in / out_tree / out_nobl / outgroup / outdir；
   - 09 使用 codon_dir / keep_list / codon_pos / nobl_trees / work_dir；
   - 10 使用 work_dir / seqfile / treefile / publish_ultrametric_tree；
   - mcmctree_postcheck / ess_report / finetune_suggest 共用 work_dir 与 qc_report 子目录。

8) cafe5
   - 11 使用 family / ultrametric_tree；
   - 12 使用 models / threads / error_model / two_stage_large 等；
   - 13 使用 models；
   - 13b 使用 enrich_bridge 与 inputs.all_members_tsv。

9) report
   - fdr_alpha 影响 13 与 14；
   - selected_cafe_model 主要用于后续作图或报告侧的默认模型选择。


====================================================================================================
四、项目目录与结果目录树状结构（简略）
====================================================================================================

以下为一个典型 aphylo 工程的目录结构示意，风格与 phylo 项目类似，仅展示主要层级：

project/
├── config.yaml
├── README.txt
├── data/
│   ├── annotation/
│   ├── genomic/
│   ├── proteomes/
│   ├── cds/
│   └── maps/
│       └── pep2cds.tsv
└── results/
    ├── publish/
    │   └── aphylo_ready/
    │       ├── strict_sco_msa/
    │       ├── colmask/
    │       ├── ogs_selected.list
    │       ├── sco_filelist.txt
    │       ├── ogs_selected.published.list
    │       ├── sco_filelist.published.txt
    │       └── family.tsv
    ├── 01_iface/
    ├── 02_bt/
    ├── 03_codon/
    ├── 03_qc/
    ├── 04_codeml/
    ├── 05_cmlagg/
    ├── 06_cafe/
    │   ├── mcmctree/
    │   ├── input/
    │   └── models/
    ├── 07_cafeagg/
    │   └── enrich_inputs/
    └── 08_joint/
        └── report/

logs/


====================================================================================================
五、简要排查顺序（仅作提示）
====================================================================================================

1) 先看 01_iface：接口是否通过（物种树、MSA、pep2cds 是否完全一致）。
2) 再看 03_qc：codon MSA 保留的 OG 数量是否足以支撑 PSG 和 MCMCTree。
3) 用 mcmctree_postcheck 评估 ESS，必要时结合 finetune_suggest 调整链长与步长。
4) 确认 11 / 12 / 13 生成的 CAFE 结果正常，再进入 13b 做富集输入。
5) 最后检查 14 / 15 是否正常结束，report/ 下 TSV 是否齐全。



产物附录

====================================================================================================================
文件/目录                               描述                                                       路径
                                                                                                   来源脚本
====================================================================================================================

【01–04：Codon MSA 阶段】
--------------------------------------------------------------------------------------------------------------------
pep2cds_map.tsv                         PEP/CDS 对照表，用于所有 CDS 展开操作                      results/02_sync_pep_to_cds/
                                                                                                   02_sync_pep_to_cds.py

pep_clean/                               清洗后的 PEP FASTA                                       results/02_sync_pep_to_cds/pep_clean/
cds_clean/                               清洗后的 CDS FASTA                                       results/02_sync_pep_to_cds/cds_clean/
                                                                                                   02_sync_pep_to_cds.py

codon_msa/OG*.fa                         每个 OG 的 codon 对齐文件                                 results/03_codon/codon_msa/
                                                                                                   03_build_and_qc_codon_alignments.py

codon_msa_qc/*.tsv                       Codon MSA QC 报告（终止密码子/缺失）                      results/03_codon/qc/
                                                                                                   03_build_and_qc_codon_alignments.py

strict_sco_msa/                          严格 SCO codon MSA                                        results/03_codon/strict_sco_msa/
strict_sco_list.tsv                      严格 SCO 列表                                             results/03_codon/
OG.order.tsv                             物种顺序与 OG → 物种信息                                 results/03_codon/
                                                                                                   build_strict_sco_msa_dir.py

====================================================================================================================
【05：前景集合定义】
--------------------------------------------------------------------------------------------------------------------
FG.list                                  前景物种列表（terminals-only）                           results/04_codeml/sets/FG.list
FG.nwk                                   在树上对前景叶节点打 #1 的前景树                         results/04_codeml/sets/FG.nwk
.fgsets.done                             前景集合定义完成哨兵                                     results/04_codeml/sets/.fgsets.done
                                                                                                   05_define_foregrounds.py

====================================================================================================================
【06：codeml ALT/NULL 运行】
--------------------------------------------------------------------------------------------------------------------
raw/OG*/                                  每个 OG 的 codeml 运行目录（ctl/mlc/out）                results/04_codeml/raw/OG*/
                                                                                                   06_codeml_batch_branchsite.py

====================================================================================================================
【07：codeml 聚合（PSG 核心 + 富集输入）】
--------------------------------------------------------------------------------------------------------------------
D_fdr_genes.tsv                          PSG 基因显著性表（Q 校正后）                             results/05_cmlagg/D_fdr_genes.tsv
                                                                                                   07_codeml_aggregate.py

D_beb_sites.tsv                          BEB 正选择位点表                                         results/05_cmlagg/D_beb_sites.tsv
                                                                                                   07_codeml_aggregate.py

psg_cds_inputs/psg_<FG>/test.list        该前景的 PSG CDS 集合                                     results/05_cmlagg/psg_cds_inputs/
psg_cds_inputs/psg_<FG>/background.list  该前景背景 CDS 集合                                       results/05_cmlagg/psg_cds_inputs/
psg_cds_inputs/psg_<FG>/meta.tsv         PSG 富集输入元信息                                       results/05_cmlagg/psg_cds_inputs/
                                                                                                   07_codeml_aggregate.py

.cmlagg.done                              codeml 聚合完成哨兵                                      results/05_cmlagg/.cmlagg.done
                                                                                                   07_codeml_aggregate.py

====================================================================================================================
【08：树清理 & 去分支长度】
--------------------------------------------------------------------------------------------------------------------
species_calib.trees                      两行树（带分支长度）                                     results/06_cafe/mcmctree/
species_calib.nobl.trees                 两行树（无分支长度）                                     results/06_cafe/mcmctree/
                                                                                                   08_clean_tree_and_noblize.py

====================================================================================================================
【09：拼接序列 concat.clean.phy】
--------------------------------------------------------------------------------------------------------------------
concat.clean.phy                         严格 SCO codon MSA 拼接序列                              results/06_cafe/mcmctree/concat.clean.phy
                                                                                                   09_build_concat_clean_phy.py

====================================================================================================================
【10：MCMCTree 两阶段运行 + 发布树】
--------------------------------------------------------------------------------------------------------------------
mcmctree.ctl                              Stage1 控制文件                                          results/06_cafe/mcmctree/
mcmctree2.ctl                             Stage2 控制文件                                          results/06_cafe/mcmctree/

out.BV                                     Stage1 BV 输出                                           results/06_cafe/mcmctree/
in.BV                                      Stage2 输入 BV                                           results/06_cafe/mcmctree/

out.txt                                    节点年龄、HPD 等结果                                     results/06_cafe/mcmctree/
mcmc.txt                                   MCMC 轨迹文件（ESS 输入）                               results/06_cafe/mcmctree/
FigTree.tre                                含注释的超时钟树                                        results/06_cafe/mcmctree/

utree_for_cafe.nwk                         给 CAFE 使用的清洁超时钟树                              results/06_cafe/input/
                                                                                                   10_mcmctree_run_and_publish.py

.stage1.done /.stage2.done /.stage3.done   MCMCTree 运行哨兵                                       results/06_cafe/mcmctree/

====================================================================================================================
【MCMCTree QC：ESS + finetune】
--------------------------------------------------------------------------------------------------------------------
qc_report/node_ages.tsv                   从 out.txt 解析的节点年龄                                results/06_cafe/mcmctree/qc_report/
qc_report/ess.tsv                         ESS 数值（IPS + coda）                                   results/06_cafe/mcmctree/qc_report/
qc_report/ess_summary.md                  ESS 概览报告                                             results/06_cafe/mcmctree/qc_report/
qc_report/ess_recommendation.txt          推荐 MCMC 长度倍数                                       results/06_cafe/mcmctree/qc_report/
                                                                                                   ess_report.R

qc_report/finetune_suggestion.md          基于 ESS 的 finetune 建议                                results/06_cafe/mcmctree/qc_report/
qc_report/finetune_new_line.txt           finetune= 新行                                           results/06_cafe/mcmctree/qc_report/
                                                                                                   finetune_suggest.py + mcmctree_postcheck.py

====================================================================================================================
【11：CAFE5 输入准备】
--------------------------------------------------------------------------------------------------------------------
family.tsv                                按树叶顺序重排的家族矩阵                                 results/06_cafe/input/family.tsv
cleaned_ultrametric_tree.nwk              清洁后的超时钟树                                         results/06_cafe/input/
.prep.done                                CAFE 输入准备哨兵                                        results/06_cafe/
                                                                                                   11_cafe5_prepare_inputs.py

====================================================================================================================
【12：CAFE5 批量模型运行】
--------------------------------------------------------------------------------------------------------------------
models/<model>/primary_global/            primary 模型目录（Gamma/Base）                          results/06_cafe/models/<model>/
Gamma_family_results.txt                  family P（Gamma）                                        同上
Base_family_results.txt                   family P（Base）                                         同上
Gamma_results.txt / Base_results.txt      分支 λ                                                  同上
Base_clade_results.txt                    每个节点 Increase/Decrease                              同上
Base_change.tab                           Δ（child-parent）矩阵                                    同上
Base_branch_probabilities.tab             分支显著性概率                                           同上

flags/high_fail_ogs.list                  高失败家族                                               results/06_cafe/models/<model>/flags/
flags/extreme_ogs_round*.list             极端 OG 列表（多轮）                                    同上

.cafe.done                                CAFE 全流程完成哨兵                                     results/06_cafe/

                                                                                                   12_cafe5_run_models.py

====================================================================================================================
【13：CAFE 聚合汇总】
--------------------------------------------------------------------------------------------------------------------
cafe_significant_families.tsv             全局 FDR 后显著家族                                      results/07_cafeagg/
cafe_significant_families_no_highfail.tsv 剔除 high_fail 后显著家族                               results/07_cafeagg/

cafe_branch_summary.tsv                   分支 λ 摘要表                                            results/07_cafeagg/
cafe_family_universe.tsv                  家族宇宙：primary/large 覆盖情况                         results/07_cafeagg/
cafe_clade_summary.tsv                    每物种/内部节点的扩缩计数                               results/07_cafeagg/

inputs_used.tsv                           所有输入文件路径追踪                                    results/07_cafeagg/
                                                                                                   13_cafe5_aggregate.py

====================================================================================================================
【13b：CAFE per-species 富集输入】
--------------------------------------------------------------------------------------------------------------------
cafe_family_branch_status.tsv             family × species 的扩张/收缩状态                        results/07_cafeagg/

enrich_inputs/<tag>/background.list       背景基因                                                 results/07_cafeagg/enrich_inputs/<tag>/
enrich_inputs/<tag>/up.list               扩张基因                                                 同上
enrich_inputs/<tag>/down.list             收缩基因                                                 同上
enrich_inputs/<tag>/meta.tsv              元信息                                                   同上

                                                                                                   13b_cafe_build_enrich_inputs.py

====================================================================================================================
【14：PSG × CAFE 联合整合】
--------------------------------------------------------------------------------------------------------------------
integration_counts.tsv                    PSG vs CAFE 计数统计                                     results/08_joint/
integration_intersect.tsv                 交集（PSG ∩ CAFE）                                      results/08_joint/
integration_union.tsv                     并集（PSG ∪ CAFE）                                      results/08_joint/
.integration.done                         联合整合哨兵                                             results/08_joint/
                                                                                                   14_joint_integration.py

====================================================================================================================
【15：最终报告】
--------------------------------------------------------------------------------------------------------------------
report/tables/                            复制 PSG+CAFE+Integration 的全部关键表格                results/08_joint/report/tables/
report/figs/                              预留图输出目录                                           results/08_joint/report/figs/
README.txt                                最终报告说明                                             results/08_joint/report/README.txt
.report.done                              最终报告哨兵                                             results/08_joint/
                                                                                                   15_build_final_report.py

====================================================================================================================



====================================================================================================
END OF README.txt
====================================================================================================