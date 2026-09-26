# Quick Start

This guide starts with a fresh GitHub clone and shows how to run the command-line
programs and plotting scripts. All commands below are run from the repository
root, and generated files stay inside the repository.

## 1. Clone the repository

```bash
git clone https://github.com/yeeyoe/sGKZ.git
cd sGKZ
```

Prebuilt executables are expected at `build/shortest_gkz`,
`build/K-stability/k_stability`, `build/K-stability/k_stability_search`, and
`build/K-stability_highdim/k_stability_highdim`. If your clone does not contain
them, build the project using the main README. Binaries are platform-specific.

Create the output directories once:

```bash
mkdir -p results K-stability/K-results K-stability_highdim/nK_results
```

## 2. Run `shortest_gkz`

For a point set:

```bash
./build/shortest_gkz \
  --points examples/six_points.points \
  --output results/six_points.csv \
  --verbose 2> results/six_points.log
```

For a lattice polygon, choose a dilation level with `--k`:

```bash
./build/shortest_gkz \
  --polygon examples/unit_square.polygon --k 8 \
  --output results/square_k8.csv \
  --plot-prefix results/square_k8 \
  --verbose 2> results/square_k8.log
```

`points` and `twice_area` report the input size and hull area. `converged` and
`exact_certified` report numerical convergence and exact verification;
`iterations` and `active_size` describe the final active set. `gap` is the
Frank--Wolfe optimality gap, while `norm_squared` and `exact_norm_squared` are
the numerical and exact squared norms. The `output` line gives the CSV path.
Exit code `0` means success, `1` means an error, and `2` means the iteration
limit was reached before convergence.

## 3. Run `k_stability`

```bash
./build/K-stability/k_stability \
  --polygon examples/unit_square.polygon \
  --svg K-stability/K-results/unit_square.svg \
  > K-stability/K-results/unit_square.out
```

Evaluate one exact crease line with coefficients `a b c`:

```bash
./build/K-stability/k_stability \
  --polygon examples/unit_square.polygon \
  --check-line "1 0 -1/2" \
  > K-stability/K-results/unit_square_check.out
```

Add `--certify` for an exact rational certificate. Exit code `0` means a
certified witness, `2` means no counterexample (or failed certification), and
`1` indicates an error.

For the unit-square check, `ell_P(x,y) = 4` and `M_l=1/4` are expected.
`check_line_exact=true` confirms exact evaluation; a negative `M_l` is a
certified relative K-instability witness.

## 4. Run `k_stability_search`

Search reports are stored in `K-stability/K-results/`; the SQLite search state
is stored at `K-stability/k_stability_search.sqlite`:

```bash
./build/K-stability/k_stability_search \
  --d 6 --N 4 --M 4 --time-limit 300 \
  --database K-stability/k_stability_search.sqlite \
  --output-dir K-stability/K-results \
  > K-stability/K-results/k_stability_search.out
```

Run the same command again to resume the saved database. Add `--stop-on-first`
or `--verbose` as needed.

`generated`, `rejected`, `probes`, `confirms`, and `finals` are search counters.
`verified` and `unverified` count candidates after validation;
`best_twice_area` is the smallest exact area among verified unstable candidates,
and `best_verified_key` identifies that candidate. Exit code `0` means at least
one verified unstable candidate was found; `2` means none was found; `1` means an
error.

### Search keys and SQLite queries

Each candidate has a canonical key of the form
`dD|p=dx1:dy1;...;dxD:dyD|k=s1,...,sD`:

- `dD` is the number of polygon vertices.
- `p` lists the primitive edge directions in canonical cyclic order.
- `k` lists the positive integer edge lengths (steps) in that order.

The key printed by the search is the value to use in the database queries below.
The database contains a `candidates` table (polygon geometry and status) and an
`attempts` table (probe/confirm/final results and witnesses).

List the six smallest-volume (`twice_area`) candidates with `d=6` that are
certified unstable:

```bash
sqlite3 K-stability/k_stability_search.sqlite \
  "SELECT key, twice_area, singular_vertex_count
     FROM candidates
    WHERE d = 6 AND status = 'verified_unstable'
    ORDER BY CAST(twice_area AS INTEGER), key
    LIMIT 6;"
```

Inspect one candidate by replacing `CANDIDATE_KEY` with the exact key:

```bash
sqlite3 -header -column K-stability/k_stability_search.sqlite \
  "SELECT key, d, twice_area, status, directions, steps, vertices,
          normals, ell0, ell1, ell2, singular_vertex_count
     FROM candidates
    WHERE key = 'CANDIDATE_KEY';
```

Show every stored detector attempt and its witness for that candidate. Exact
columns are populated when rational certification succeeded; otherwise the
numeric witness columns may be populated:

```bash
sqlite3 -header -column K-stability/k_stability_search.sqlite \
  "SELECT profile, status, value, witness_ux, witness_uy, witness_t,
          exact_a, exact_b, exact_c, exact_value
     FROM attempts
    WHERE candidate_key = 'CANDIDATE_KEY'
    ORDER BY rowid;"
```

For a certified attempt, the exact witness is
`g(x,y) = max(a*x + b*y + c, 0)`, using `exact_a`, `exact_b`, and `exact_c`;
`exact_value` is its exact $M_\ell$ value. For an uncertified attempt, the
numeric convention is
`g(x,y) = max(witness_ux*x + witness_uy*y - witness_t, 0)`.

To render a verified candidate and its witness as an SVG, use the same key with
the Python helper:

```bash
python3 K-stability/plot_candidate.py \
  --database K-stability/k_stability_search.sqlite \
  --key 'CANDIDATE_KEY' \
  --output K-stability/K-results/CANDIDATE_KEY.svg
```

## 5. Run `k_stability_highdim`

Compute exact affine $ℓ_P$:

```bash
./build/K-stability_highdim/k_stability_highdim \
  --polytope K-stability_highdim/examples/unit_cube5.polytope --ell-only \
  > K-stability_highdim/nK_results/unit_cube5_ell.out
```

Check a PL function:

```bash
./build/K-stability_highdim/k_stability_highdim \
  --polytope K-stability_highdim/examples/unit_cube5.polytope \
  --check-pl K-stability_highdim/examples/unit_cube5_hinge.pl \
  > K-stability_highdim/nK_results/unit_cube5_check.out
```

Search one nonzero affine branch:

```bash
./build/K-stability_highdim/k_stability_highdim \
  --polytope K-stability_highdim/examples/unit_cube5.polytope \
  --pieces 1 --generations 80 --quadrature-samples 20000 \
  > K-stability_highdim/nK_results/unit_cube5_search.out
```

Exit codes are `0` for a certified unstable witness, `2` for no counterexample,
`3` for an unverified candidate, and `1` for an error.

`dimension`, `volume`, and `boundary_measure` describe the polytope. `ell_P`
lists its constant term and coordinate coefficients. In `--check-pl` mode,
`check_pl_M_l` is the exact invariant and `certified=true` means it is negative.
In `--pieces` mode, `sweep_min_M_l` is the best numerical value and
`certified_M_l` appears after exact certification.

## 6. Python scripts

Install NumPy and Matplotlib for PNG export:

```bash
python3 -m pip install numpy matplotlib
```

Plot CSV files produced with `--plot-prefix`:

```bash
python3 plot_results.py results/square_k8 --png
```

This reads the `_surface.csv`, `_triangles.csv`, and `_subdivision.csv` files
created by `--plot-prefix`, then writes HTML and (with `--png`) PNG plots beside
the prefix.

Plot the verbose iteration log:

```bash
python3 plot_iterations.py results/square_k8.log \
  --output-prefix results/square_k8_iterations --png
```

It reads `iteration=...` records from the verbose log and prints the generated
plot paths.

Render a verified search candidate from SQLite (replace `CANDIDATE_KEY` with the
exact `key=` printed by `k_stability_search`):

```bash
python3 K-stability/plot_candidate.py \
  --database K-stability/k_stability_search.sqlite \
  --key 'CANDIDATE_KEY' \
  --output K-stability/K-results/candidate.svg
```

Use `--save` instead of `--output` to save the SVG and vertex file under
`examples/`. Use `--help` with every Python script for its complete options.

Generated results belong in `results/`, `K-stability/K-results/`, and
`K-stability_highdim/nK_results/`; input examples are in the corresponding
`examples/` directories.
