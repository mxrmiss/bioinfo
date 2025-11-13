# Finetune Recommendation (ESS-only)

- Generated at: 2025-11-10 13:02:32
- Work dir: `mcmctree`
- Ctl file: `mcmctree/mcmctree2.ctl`

## Current finetune

`finetune = 1: 8.00 12.00 1.00 5.00 12.00 5.00`

## Recommended finetune line

```
finetune = 1: 8.00 12.00 1.00 5.00 12.00 5.00
```

> 依据：all ESS ≥ 250, keep finetune；规则：对 ESS 最差 Top-3 参数按 mult = clamp((TARGET/ESS)^0.5, 0.8..1.25) 温和缩放到 [0.005, 2.0]。
