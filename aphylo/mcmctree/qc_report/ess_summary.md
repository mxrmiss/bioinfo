# ESS Report (IPS + coda)

- Using config: `/home/mxrmiss/project/aphylo/config.yaml`
- File: `mcmctree/mcmc.txt`
- Parameters evaluated: **16**
- Target ESS: **250**

## Overview
- Min ESS (IPS): **318.75**
- Median ESS (IPS): **1150.85**
- 25% quantile (IPS): **610.10**
- Count < Target (IPS): **0**

## Worst parameters (Top-K by ESS, IPS)
| rank | param | ESS |
|---:|:------|---:|
| 1 | t_n28 | 318.75 |
| 2 | t_n24 | 482.72 |
| 3 | t_n27 | 527.09 |
| 4 | t_n25 | 528.88 |
| 5 | t_n23 | 637.17 |
| 6 | t_n26 | 707.65 |
| 7 | t_n22 | 952.15 |
| 8 | t_n16 | 1109.64 |
| 9 | t_n21 | 1192.07 |
| 10 | t_n20 | 1728.66 |
| 11 | t_n17 | 1866.16 |
| 12 | t_n19 | 1924.99 |
| 13 | t_n29 | 2307.26 |
| 14 | t_n18 | 3105.28 |
| 15 | sigma2 | 5915.44 |

## Method discrepancy (Top-K by |Δ%|)
| rank | param | ESS_ips | ESS_coda | Δ% |
|---:|:------|-------:|--------:|---:|
| 1 | t_n22 | 952.15 | 1339.88 | +40.7 |
| 2 | t_n21 | 1192.07 | 1609.49 | +35.0 |
| 3 | t_n27 | 527.09 | 675.71 | +28.2 |
| 4 | t_n23 | 637.17 | 804.47 | +26.3 |
| 5 | mu | 7936.25 | 9998.35 | +26.0 |
| 6 | sigma2 | 5915.44 | 7311.07 | +23.6 |
| 7 | t_n25 | 528.88 | 611.27 | +15.6 |
| 8 | t_n20 | 1728.66 | 1996.64 | +15.5 |
| 9 | t_n26 | 707.65 | 810.59 | +14.5 |
| 10 | t_n19 | 1924.99 | 2204.67 | +14.5 |
| 11 | t_n29 | 2307.26 | 2635.48 | +14.2 |
| 12 | t_n18 | 3105.28 | 3505.67 | +12.9 |
| 13 | t_n16 | 1109.64 | 1245.19 | +12.2 |
| 14 | t_n24 | 482.72 | 528.69 | +9.5 |
| 15 | t_n17 | 1866.16 | 2023.86 | +8.5 |

## Recommendation
- Suggested `nsample` multiplier: **×1**  (using IPS min ESS)
- Prefer two independent chains; IPS (Tracer风) 达标 + trace 平稳更审稿友好。

_Artifacts_: `ess.tsv`（含 ESS 与两口径细节）, `ess_recommendation.txt`。
