# sgkz

独立的动态维 shortest GKZ vector 求解器。支持 `--points FILE` 和
`--polytope FILE --k INTEGER` 两种模式，运行时从输入列数推断维数（`d>=2`）。

构建：

```bash
cmake -S sGKZ_ndim -B sGKZ_ndim/build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build sGKZ_ndim/build --target sgkz --parallel
ctest --test-dir sGKZ_ndim/build --output-on-failure
```

构建完成后，主程序位于 `sGKZ_ndim/build/sgkz`。

示例：

```bash
sGKZ_ndim/sgkz --points sGKZ_ndim/examples/simplex3.points
sGKZ_ndim/sgkz --polytope sGKZ_ndim/examples/unit_cube3.polytope --k 2 --ell
```

polytope 文件列出整数顶点；程序计算 `kP ∩ Z^d`。`--max-points` 默认值为
1000000，传入 0 可取消限制。结果 CSV 使用动态列，不生成绘图文件。
