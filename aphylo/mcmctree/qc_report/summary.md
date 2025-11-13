# MCMCTree QC Summary

- Generated at: 2025-11-10 13:02:32

- Work dir: `mcmctree`

## Artifacts

- ✅ `mcmctree/in.BV` (9.01 KB)
- ✅ `mcmctree/out.txt` (13.25 MB)
- ✅ `mcmctree/FigTree.tre` (1.05 KB)
- ✅ `mcmctree/mcmc.txt` (33.52 MB)
- ✅ `mcmctree/qc_report/node_ages.tsv`
- ✅ `mcmctree/qc_report/ess.tsv`, `mcmctree/qc_report/ess_summary.md`, `mcmctree/qc_report/ess_recommendation.txt`

## Node Ages & 95% HPD

- Parsed from `out.txt`: **14** nodes
- Approx. root by widest HPD (out.txt): mean=7.512, 95%HPD=(6.336, 9.935)
- Root~approx from `FigTree.tre`: mean≈7.768, 95%HPD=(6.129, 9.406)
- Node count check: ✅ parsed **14** ≈ expected **14** (taxa=15)

## MCMC Trace

- Lines (non-empty): 200002
- Header detected: Yes

## ESS Snapshot & Recommendation

- Min ESS: **318.75**; Median ESS: **1192.07**; Target: **250**
- Suggested `nsample` multiplier: **×1**
- Details: see `mcmctree/qc_report/ess_summary.md` and `mcmctree/qc_report/ess.tsv`

## Checklist

- [ ] 若更换对齐/模型/拓扑 → 重新生成 `in.BV` 再跑采样
- [ ] 如果 Min ESS < Target → 提高 `nsample`（按推荐倍数），必要时微调 `finetune`（目标接受率 0.2–0.4）
- [ ] 给 CAFE5：使用均值/中位数时间树（必要时转为仅含时间分支长度的 Newick）