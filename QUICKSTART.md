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

## 6. Python scripts

Install NumPy and Matplotlib for PNG export:

```bash
python3 -m pip install numpy matplotlib
```

Plot CSV files produced with `--plot-prefix`:

```bash
python3 plot_results.py results/square_k8 --png
```

Plot the verbose iteration log:

```bash
python3 plot_iterations.py results/square_k8.log \
  --output-prefix results/square_k8_iterations --png
```

Render a verified search candidate from SQLite (replace the key with the exact
`key=` printed by `k_stability_search`):

```bash
python3 K-stability/plot_candidate.py \
  --database K-stability/k_stability_search.sqlite \
  --key 'd6|p=...' \
  --output K-stability/K-results/candidate.svg
```

Use `--save` instead of `--output` to save the SVG and vertex file under
`examples/`. Use `--help` with every Python script for its complete options.

Generated results belong in `results/`, `K-stability/K-results/`, and
`K-stability_highdim/nK_results/`; input examples are in the corresponding
`examples/` directories.
