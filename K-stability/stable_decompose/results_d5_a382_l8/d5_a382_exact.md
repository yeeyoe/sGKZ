# Exact crease decomposition

- Crease edges: `256`
- Threshold: `5.136518752190988e-05`
- Approximate endpoints: `((4.878136671814176e-16, 0.0), (4.9999999999999964, 11.000000000000002))`
- Candidate boundary edge pairs: `[[0, 2], [0, 3], [4, 2], [4, 3]]`
- SVG diagram: `results_d5_a382_l8/d5_a382_crease.svg` (1 accepted crease(s))

## Selected algebraic branch

| Edges | `s` approx. | `t` approx. | `s` isolating interval | `t` isolating interval | Endpoint approximations | Endpoint coordinate intervals |
|---|---:|---:|---|---|---|---|
| `[0, 3]` | `0.449222552479224` | `0.00652993216254958` | `['712620455/1586341672', '891293813/1984080737']` | `['11602159/1776765625', '7755123/1187626886']` | (4.49222552479224, 0)<br>(4.89552108539921, 10.9412306105371) | (['3563102275/793170836', '8912938130/1984080737'], ['0', '0'])<br>(['2907026231/593813443', '8698193581/1776765625'], ['12994099639/1187626886', '19440002444/1776765625']) |

## Candidate systems

### Edges `[0, 2]`

- Continuity equation degrees: `16` and `16` (81 and 81 terms)
- Resultant degrees: `s=108`, `t=108`
- Common factor removed: `s*t + 4*s/7 - 211*t/70 - 4/7`
- Common factor is singular: `True`
- Common factor requires review: `False`
- Paired root boxes: `1`

- Full equations, resultants and root lists: [`d5_a382_exact_edges_0_2.txt`](d5_a382_exact_edges_0_2.txt)

| `s` approx. | `t` approx. | min poly degrees | `s` interval | `t` interval | Multiplicities |
|---:|---:|---|---|---|---|
| `0.886204256240215` | `0` | 28, 1 | `['5207855388/5876585845', '5300696785/5981348823']` | `['0', '0']` | `1, 10` |

### Edges `[0, 3]`

- Continuity equation degrees: `16` and `16` (81 and 81 terms)
- Resultant degrees: `s=108`, `t=108`
- Common factor removed: `s*t - 11*s/9 + 131*t/90 - 131/90`
- Common factor is singular: `True`
- Common factor requires review: `False`
- Paired root boxes: `1`

- Full equations, resultants and root lists: [`d5_a382_exact_edges_0_3.txt`](d5_a382_exact_edges_0_3.txt)

| `s` approx. | `t` approx. | min poly degrees | `s` interval | `t` interval | Multiplicities |
|---:|---:|---|---|---|---|
| `0.449222552479224` | `0.00652993216254958` | 28, 28 | `['712620455/1586341672', '891293813/1984080737']` | `['11602159/1776765625', '7755123/1187626886']` | `1, 1` |

### Edges `[2, 4]`

- Continuity equation degrees: `16` and `16` (81 and 81 terms)
- Resultant degrees: `s=108`, `t=108`
- Common factor removed: `s*t - 256*s/45 + 86*t/45 + 256/45`
- Common factor is singular: `True`
- Common factor requires review: `False`
- Paired root boxes: `1`

- Full equations, resultants and root lists: [`d5_a382_exact_edges_2_4.txt`](d5_a382_exact_edges_2_4.txt)

| `s` approx. | `t` approx. | min poly degrees | `s` interval | `t` interval | Multiplicities |
|---:|---:|---|---|---|---|
| `0.911354619264201` | `0.0256169596291567` | 28, 28 | `['2361281219/2590957646', '2523672729/2769144607']` | `['163312097/6375155341', '159410904/6222865879']` | `1, 1` |

### Edges `[3, 4]`

- Continuity equation degrees: `10` and `10` (29 and 31 terms)
- Resultant degrees: `s=43`, `t=43`
- Common factor removed: `s**4*t**4 - 4*s**3*t**4 + 6*s**2*t**4 - 4*s*t**4 + t**4`
- Common factor is singular: `True`
- Common factor requires review: `False`
- Paired root boxes: `0`

- Full equations, resultants and root lists: [`d5_a382_exact_edges_3_4.txt`](d5_a382_exact_edges_3_4.txt)

## Candidate screening

### Edges `[0, 2]`; s=0.886204256240215, t=0

- Accepted: `False`
- Rejection reason: `concavity_fail`
- Concavity: `fail`
- Left scan: `no_counterexample_found`
- Right scan: `unstable`

### Edges `[0, 3]`; s=0.449222552479224, t=0.00652993216254958

- Accepted: `True`
- Rejection reason: `None`
- Concavity: `pass`
- Left scan: `no_counterexample_found`
- Right scan: `no_counterexample_found`

### Edges `[2, 4]`; s=0.911354619264201, t=0.0256169596291567

- Accepted: `False`
- Rejection reason: `concavity_fail`
- Concavity: `fail`
- Left scan: `no_counterexample_found`
- Right scan: `unstable`

## Machine-readable output

The complete structured result is stored in the accompanying `.json` file.

