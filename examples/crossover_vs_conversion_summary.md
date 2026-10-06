# Mosaic classes counted in independent variant sites (elements: 16454)

| min sites per run | class | observed | island-shuffle null | obs/null |
|---:|---|---:|---:|---:|
| 2 | crossover | 297 | 40.1 | 7.4 |
| 2 | conversion | 7 | 60.7 | 0.1 |
| 2 | complex | 0 | 28.7 | 0.0 |
| 3 | crossover | 180 | 12.7 | 14.2 |
| 3 | conversion | 2 | 15.0 | 0.1 |
| 3 | complex | 0 | 2.8 | 0.0 |
| 5 | crossover | 37 | 1.5 | 24.7 |
| 5 | conversion | 0 | 0.6 | 0.0 |
| 5 | complex | 0 | 0 | inf |

## Conversion-like calls supported by >=3 independent sites
- _Had-6b Chr4.4882538: ATHILA4c host, ATHILA9 tract ~1269 bp (13 independent sites)
- _Had-6b Chr4.4664029: ATHILA4c host, ATHILA1 tract ~806 bp (4 independent sites)
- (2-site only) 6a>2>6a, ~40 bp, Chr3 ~11.8 Mb in 4 accessions — one shared locus, weak.

## Detection power (conv_power.py; planted into clean host elements, >=3 sites)
| host>donor | crossover recall | conversion 100 | 250 | 500 | 1000 | 2000 bp |
|---|---:|---:|---:|---:|---:|---:|
| 6a>6b | 23% | 0% | 2% | 0% | 0% | 0% |
| 6>1 | 40% | 0% | 20% | 37% | 17% | 0% |
| 2>6b | 97% | 3% | 32% | 63% | 17% | 3% |
| 2>4c | 88% | 0% | 3% | 10% | 17% | 2% |
| 4>4c | 62% | 0% | 5% | 5% | 8% | 0% |
| 1>2 | 85% | 12% | 30% | 33% | 8% | 0% |
False-positive rate on clean host elements: 0% for all pairs.
