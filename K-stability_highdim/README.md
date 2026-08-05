# K-stability_highdim

这是一个与旧 `K-stability/` 完全独立的动态维格 polytope 工具。它使用
`CGAL::Gmpq` 精确计算任意维数 $d\ge2$ 的 affine function $\ell_P$，并搜索
固定分支数的凸 PL 函数

$$
f(x)=\max\{0,L_1(x),\ldots,L_M(x)\}.
$$

其中 `1 <= M <= 8`，运行时推断维数，首版验收覆盖 2--5 维。

## 构建

独立构建需要 C++20、CGAL、Eigen3 和 Boost 1.89 或更新版本：

```bash
cmake -S K-stability_highdim -B build_highdim -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build_highdim --parallel
ctest --test-dir build_highdim --output-on-failure
```

根工程默认接入此子项目，也可以关闭：

```bash
cmake -S . -B build -DBUILD_K_STABILITY_HIGHDIM=OFF
```

## 输入

polytope 文件每行是一个整数点，维数由第一行列数推断；点不需要排序：

```text
0 0 0 0 0
1 0 0 0 0
...
```

支持空行、逗号和 `#` 注释。重复点和非满维输入会被拒绝。输入中的冗余非极点允许存在；它们只用于内部 Delaunay 剖分，输出中的 `hull_points` 是剔除冗余点后的凸包顶点数，`redundant_points` 会报告剔除数量。

PL 文件每行给一个非零 affine 分支的 $d+1$ 个系数，最后一个是常数项；系数可写成整数或 `p/q`：

```text
1 0 0 0 0 -1/2
0 1 0 0 0 -2/3
```

零函数分支隐含，不写入文件。

## 用法

只计算精确 $\ell_P$：

```bash
./build/K-stability_highdim/k_stability_highdim --polytope K-stability_highdim/examples/unit_cube5.polytope --ell-only
```

精确验算给定 PL 函数：

```bash
./build/K-stability_highdim/k_stability_highdim \
  --polytope K-stability_highdim/examples/unit_cube5.polytope \
  --check-pl K-stability_highdim/examples/unit_cube5_hinge.pl
```

自动搜索时必须指定非零分支数：

```bash
./build/K-stability_highdim/k_stability_highdim \
  --polytope K-stability_highdim/examples/unit_cube5.polytope --pieces 1
```

搜索使用固定 Halton 样本和 Boost differential evolution。输出会回显有效的
`population`、`generations`、`quadrature_samples`、`seed`、`threads` 和认证
分母上限。默认搜索参数可以用 `--population`、`--generations`、
`--quadrature-samples`、`--seed` 和 `--threads` 调整。数值候选的归一化值小于 `-1e-6` 后，程序自动按分母上限
`10, 100, ..., --certify-max-denom` 有理化全部分支，并用精确区域分解重算。

退出码：`0` 表示已经认证严格负 witness；`2` 表示未发现反例；`3` 表示数值上发现负候选但没有认证成功；`1` 表示输入或运行错误。

## 数学与限制

程序通过精确 Delaunay 单纯剖分计算内部二阶矩和边界格点测度。对每个 PL 分支，认证阶段在每个源单纯形内枚举其 dominance cell 的有理顶点，再用精确单纯剖分积分对应的 affine 分支。

高维中不存在二维 Donaldson“只需一个简单折痕即可完备判定”的同样事实。因此，认证成功的负 witness 是严格的相对 K-不稳定性证明；`no_counterexample_found` 只说明在指定分支数和有限启发式搜索中没有找到候选，不能作为高维半稳定性的证明。

## 已覆盖基准

- $d$ 维单位立方体：$\ell_P=2d$，$M_\ell(\max\{x_1-1/2,0\})=1/4$。
- 标准 $d$-单形：$\ell_P=d(d+1)$。
- 2、3、5 维精确矩、边界测度、PL 重复分支和解析认证回归。
