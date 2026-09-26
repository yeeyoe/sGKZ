# `unstable_Delzant`

Standalone smooth-polygon search CLI with a separate SQLite database. It reuses
`k_stability_core` only for exact moment calculations, relative DF evaluation,
numerical simple-convex scans, and rational witness certification.

Build from the repository root with CMake. The Release build files and
executable live in this directory:

```bash
cmake -S . -B K-stability/Delzant/build-release -DCMAKE_BUILD_TYPE=Release
cmake --build K-stability/Delzant/build-release --target unstable_Delzant -j2
K-stability/Delzant/build-release/K-stability/Delzant/unstable_Delzant \
  --d 10 --N 12 --M 12 --max-diameter 1000000 --time-limit 600
```

`--d`, `--time-limit`, and `--max-diameter` are required. The latter caps the
translation-invariant L-infinity diameter of a polygon after common-scale
normalization. Polygons exceeding it are rejected before exact moment or DF
computations. `--N` bounds fan-ray coordinates
(default 12), `--M` bounds the initially sampled positive edge lengths
(default 12), and `--seed` defaults to 1. The detector defaults to 360 direction
steps and 256 offset steps for its final tier. Every candidate first receives a
32x16 probe; candidates below a small negative normalized-DF cutoff use a
128x64 confirmation, and 5% of canonical candidate keys are deterministically
audited at final precision. Negative cutoffs avoid promoting floating-point
zero noise, which was common in the d=10 benchmark. Final stage thresholds and
profile are recorded in the database.
The database records detector profile, completed stage, normalized witness
values, and evaluation counts. The database is
`K-stability/Delzant/unstable_delzant.sqlite`; output files are placed in
`K-stability/Delzant/results/`. Override these paths with `--database` and
`--output-dir`.

The generator starts from the three-ray fan of `P^2` or a four-ray Hirzebruch
fan and inserts sums of adjacent rays. Each insertion preserves the smooth fan
condition. It then samples the first `d-2` positive edge lengths and solves the
two closure equations exactly for the final two lengths; candidates with a
nonpositive closing length are discarded. Candidates are normalized by their
common edge-length factor and canonicalized under cyclic rotation.

Candidate edge lengths are divided by their common gcd before sizing and
deduplication, so integral homotheties share the same canonical key. The first
candidate whose negative simple-convex witness is exactly certified
is saved as `unstable_dD.polygon`, `.svg`, and `.txt`, and the program exits 0.
When the time budget expires without a certificate, it exits 2; input or runtime
errors exit 1. Completed candidates and the random generator state are stored in
the SQLite database, so a later invocation resumes the stream. An existing
certified candidate for the same `d` is reported immediately.

The fan generator is a randomized search, not an exhaustive enumeration. A
timeout means only that this run did not find a certified unstable candidate; it
does not prove relative K-semistability.
