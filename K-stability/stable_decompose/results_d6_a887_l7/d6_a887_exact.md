# Exact crease decomposition

- Crease edges: `128`
- Threshold: `2.2930441839629417e-05`
- Approximate endpoints: `((-6.1574854654465e-15, 0.0), (25.000000000000007, 12.0))`
- Candidate boundary edge pairs: `[[0, 3], [0, 4], [5, 3], [5, 4]]`
- SVG diagram: `results_d6_a887_l7/d6_a887_crease.svg` (0 accepted crease(s))

## Selected algebraic branch

No branch passed the finite-element parameter window.

## Candidate systems

### Edges `[0, 3]`

- Continuity equation degrees: `18` and `18` (100 and 100 terms)
- Resultant degrees: `s=132`, `t=132`
- Common factor removed: `1`
- Common factor is singular: `False`
- Common factor requires review: `False`
- Paired root boxes: `1`

- Full equations, resultants and root lists: [`d6_a887_exact_edges_0_3.txt`](d6_a887_exact_edges_0_3.txt)

| `s` approx. | `t` approx. | min poly degrees | `s` interval | `t` interval | Multiplicities |
|---:|---:|---|---|---|---|
| `0.398947152472786` | `0.476112126688023` | 38, 38 | `['879807533/2205323506', '1054501515/2643211033']` | `['1927061855/4047495846', '180095001/378261739']` | `1, 1` |

### Edges `[0, 4]`

- Continuity equation degrees: `16` and `16` (81 and 81 terms)
- Resultant degrees: `s=108`, `t=108`
- Common factor removed: `s*t - 2*s + 9*t/4 - 9/4`
- Common factor is singular: `True`
- Common factor requires review: `False`
- Paired root boxes: `0`

- Full equations, resultants and root lists: [`d6_a887_exact_edges_0_4.txt`](d6_a887_exact_edges_0_4.txt)

### Edges `[3, 5]`

- Continuity equation degrees: `16` and `16` (81 and 81 terms)
- Resultant degrees: `s=108`, `t=108`
- Common factor removed: `s*t + 58*s/85 - 247*t/85 - 58/85`
- Common factor is singular: `True`
- Common factor requires review: `False`
- Paired root boxes: `1`

- Full equations, resultants and root lists: [`d6_a887_exact_edges_3_5.txt`](d6_a887_exact_edges_3_5.txt)

| `s` approx. | `t` approx. | min poly degrees | `s` interval | `t` interval | Multiplicities |
|---:|---:|---|---|---|---|
| `0.0945688907871198` | `0.0264517450229163` | 28, 28 | `['110521323/1168685834', '92280305/975799803']` | `['54045072/2043157151', '56906015/2151314212']` | `1, 1` |

### Edges `[4, 5]`

- Continuity equation degrees: `10` and `10` (29 and 31 terms)
- Resultant degrees: `s=43`, `t=43`
- Common factor removed: `s**4*t**4 - 4*s**3*t**4 + 6*s**2*t**4 - 4*s*t**4 + t**4`
- Common factor is singular: `True`
- Common factor requires review: `False`
- Paired root boxes: `0`

- Full equations, resultants and root lists: [`d6_a887_exact_edges_4_5.txt`](d6_a887_exact_edges_4_5.txt)

## Candidate screening

### Edges `[0, 3]`; s=0.398947152472786, t=0.476112126688023

- Accepted: `False`
- Rejection reason: `piece_unstable`
- Concavity: `pass`
- Left scan: `no_counterexample_found`
- Right scan: `unstable`

### Edges `[3, 5]`; s=0.0945688907871198, t=0.0264517450229163

- Accepted: `False`
- Rejection reason: `concavity_fail`
- Concavity: `fail`
- Left scan: `unstable`
- Right scan: `unstable`

## Machine-readable output

The complete structured result is stored in the accompanying `.json` file.

