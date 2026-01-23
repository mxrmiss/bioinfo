# bioinfo

一个面向科研日常的生信工具库与流程集合（以 **Snakemake + 脚本** 为主）。  
仓库按“子项目”组织：每个子目录通常就是一个相对独立的流程/工具箱，配套自己的 `config.yaml`、`scripts/`、以及更详细的说明文档。

---

## 这个仓库里有什么？

- **phylo/**  
  多物种系统发育基因组学全流程：从基因组与注释出发，提取 CDS/蛋白，跑 OrthoFinder，构建严格单拷贝（SCO）集合与拼接比对，并产出可复用的 *aphylo-ready* 发布包接口。

- **aphylo/**  
  基于 `phylo` 的发布包做更深层分析：构建/质控 codon MSA，进行 codeml（分支位点正选择等）与相关下游分析等。

- **rna-seq/**  
  转录组分析流程：围绕样本表/对比表组织计算，支持差异分析、富集分析、WGCNA（部分步骤输出“可复现原材料”，方便后续稳定作图）。

- **swissprot/**  
  Swiss-Prot（reviewed）同源注释流程：DIAMOND 比对 + 阈值过滤 + best hit 解析，生成结构化注释表与汇总报告。

- **snp_rbh/**  
  两物种 RBH（reciprocal best hits）同源对驱动的分歧景观：蛋白 RBH → MAFFT 对齐 → PAL2NAL 回译到 codon alignment → 统计 substitution/indel → 染色体与滑窗汇总并出图。

- **magic/**  
  “顶刊风格”可视化工具箱：标准化 R 作图脚本与模板，覆盖常见转录组分析图形需求。

- **tools/**  
  零散但常用的辅助脚本（质控、格式转换、快速检查、绘图小工具等）。

---

## 通用约定（大多数子项目遵循）

- **参数集中管理**：优先使用 `config.yaml` 统一配置，不依赖复杂命令行参数。
- **目录清晰可复现**：常见结构为 `data/ scripts/ results/ logs/ templates/`。
- **输出友好**：以 TSV/表格 + 图为主，方便直接进入论文与报告工作流。

---

## 从哪里开始？

建议直接进入你要用的子项目目录，先读它的说明文件：

- `phylo/README.phylo.txt`
- `aphylo/README.aphylo.txt`
- `rna-seq/README.rna-seq.txt`
- `swissprot/README/readme.md`
- `snp_rbh/README/readme.md`
- `magic/README_CN.txt`

---

## 备注

- 本仓库偏向“科研自用工具箱”，不包含大体积原始数据；请按各子项目说明准备输入数据。
- 若你用于论文/报告，建议在方法部分明确写出：子项目名称、关键软件（例如 Snakemake/DIAMOND/MAFFT/PAL2NAL 等）与核心参数配置（见对应 `config.yaml`）。
