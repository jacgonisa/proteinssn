# Recombination vs homology — robust k-mer mosaics (internal span)

- family pairs tested (both families present): 91
- pairs with >=1 mosaic: 15
- mosaic elements: 619

## Spearman with homology (all pairs, zeros included; permutation p)

| metric | rho | p |
|---|---:|---:|
| raw | +0.35 | 0.0014 |
| loci | +0.35 | 0.0014 |
| per_opp | +0.35 | 0.0016 |
| obs_exp | +0.35 | 0.0016 |
| per_min | +0.35 | 0.0016 |
| per_opp_det | +0.34 | 0.0018 |

## Rate model: count ~ homology + log(detectability) + offset(log n_A n_B)

-   raw  NegBin: IRR per +10% homology = 54.75 (95% CI 34.65-86.50), p = 6e-66
-   raw Poisson: IRR per +10% homology = 4.08 (95% CI 3.67-4.55), p = 6.9e-144
-  loci  NegBin: IRR per +10% homology = 26.86 (95% CI 15.62-46.19), p = 1.3e-32
-  loci Poisson: IRR per +10% homology = 5.23 (95% CI 4.14-6.60), p = 5.6e-44
