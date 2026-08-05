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

所有调用都必须提供 `--polytope FILE`，并且必须且只能选择以下三种模式之一：

- `--ell-only`：只计算并输出精确的 $\ell_P$，不执行 PL 搜索。
- `--check-pl FILE`：读取给定的 PL 函数并进行精确验算。`FILE` 中每行是一个非零 affine 分支；零函数分支隐含存在。
- `--pieces M`：搜索 `M` 个非零 affine 分支的凸 PL 函数，其中 `M` 必须满足 `1 <= M <= 8`。该模式还需要下面的搜索参数。

### 参数说明

| 参数 | 说明 |
| --- | --- |
| `--polytope FILE` | 必需。格 polytope 输入文件；维数由第一行的列数推断。 |
| `--ell-only` | 只计算 $\ell_P$。与 `--pieces`、`--check-pl` 互斥。 |
| `--check-pl FILE` | 精确验算给定 PL 函数。与 `--ell-only`、`--pieces` 互斥。 |
| `--pieces M` | 自动搜索 `M` 个非零分支，取值范围为 `1..8`。与 `--ell-only`、`--check-pl` 互斥。 |
| `--population N` | differential evolution 的种群大小。默认 `0`，此时自动取 `max(32, 8 * (M * d + M - 1))`；其中 `d` 是 polytope 维数。 |
| `--generations N` | differential evolution 的迭代代数，默认 `80`；必须为正数。 |
| `--quadrature-samples N` | 每个积分使用的 Halton 样本数，默认 `20000`；必须为正数。 |
| `--seed N` | 搜索随机种子，默认 `0`；用于控制确定性随机数生成器。 |
| `--threads N` | 搜索工作线程数，默认 `1`；必须为正数。 |
| `--certify-max-denom N` | 数值 witness 有理化时尝试的最大分母，默认 `1048576`；必须为正数。仅用于 `--pieces` 模式。 |
| `--verbose` | 将搜索摘要输出到标准错误。仅影响输出，不改变搜索参数。 |
| `--help` | 显示命令行帮助并退出。 |

`--population`、`--generations`、`--quadrature-samples`、`--seed`、`--threads`
和 `--certify-max-denom` 主要用于 `--pieces` 模式；在另外两种模式下不会改变
精确计算或验算结果。搜索得到的数值候选只有在归一化值小于 `-1e-6` 时才会进入
认证阶段，认证会按分母上限 `10, 100, ...` 逐步尝试，并使用精确区域分解重算。

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

## 输出说明

所有模式都会先输出 polytope 的基本信息。数值为整数或有理数时使用精确形式；
例如 `1/2` 表示有理数，而不是浮点近似。

| 输出字段 | 含义 |
| --- | --- |
| `dimension` | 输入 polytope 的维数 `d`。 |
| `input_points` | 输入文件中读取的点数，包含冗余点。 |
| `hull_points` | 去除冗余非极点后保留的凸包顶点数。 |
| `redundant_points` | 被识别为冗余非极点的点数。 |
| `simplices` | 内部 Delaunay 剖分中的 `d`-单纯形数量。 |
| `boundary_facets` | 边界 facet 数量。 |
| `volume` | polytope 的精确体积。 |
| `boundary_measure` | polytope 边界的精确格测度。 |
| `ell_P` | 精确 affine function $\ell_P(x)=c+a_1x_1+\cdots+a_dx_d$，输出顺序为常数项、各坐标系数。 |

### `--check-pl` 输出

| 输出字段 | 含义 |
| --- | --- |
| `check_pl_branches` | PL 文件中读取的非零 affine 分支数。 |
| `check_pl_M_l` | 给定 PL 函数的精确 $M_\ell$ 值。 |
| `certified` | 精确验算结果；`true` 表示 $M_\ell<0$。 |
| `relative_K_status` | `unstable` 表示已得到严格负的精确 witness；`no_counterexample_found` 表示该 PL 函数未给出反例。 |

### `--pieces` 输出

搜索模式使用固定 Halton 样本和 Boost differential evolution。除公共字段外，还会
输出：

| 输出字段 | 含义 |
| --- | --- |
| `search_pieces` | 搜索的非零分支数 `M`。 |
| `search_parameter_dimension` | 搜索参数维数，计算为 `M * d + M - 1`。 |
| `search_population` | 实际使用的种群大小；当 `--population 0` 时为自动计算值。 |
| `search_generations` | 实际使用的迭代代数。 |
| `search_quadrature_samples` | 每个积分使用的 Halton 样本数。 |
| `search_seed` | 实际使用的随机种子。 |
| `search_threads` | 实际使用的搜索线程数。 |
| `certify_max_denom` | 精确认证使用的最大分母上限。 |
| `search_evaluations` | 数值搜索评估目标函数的次数。 |
| `sweep_min_M_l` | 搜索阶段找到的最小数值 $M_\ell$。 |
| `sweep_min_M_l_normalized` | 按内部归一化尺度除后的数值，用于决定是否进入认证阶段。 |
| `witness_L1`, `witness_L2`, ... | 搜索得到的数值候选经有理化后的 affine 分支。 |
| `certified_M_l` | 认证阶段得到的精确 $M_\ell$ 值。 |
| `certified_L1`, `certified_L2`, ... | 通过精确认证的 affine 分支。 |

当 `sweep_min_M_l_normalized < -1e-6` 时，程序会按分母上限
`10, 100, ...` 逐步有理化候选，并用精确区域分解重新计算。否则直接输出
`relative_K_status=no_counterexample_found`，不会进行认证。

`relative_K_status` 和退出码对应如下：

| 状态 | 含义 | 退出码 |
| --- | --- | --- |
| `unstable` | 已认证严格负的 witness，证明相对 K-不稳定。 | `0` |
| `no_counterexample_found` | 在指定分支数和有限数值搜索中没有找到候选反例。 | `2` |
| `unverified_candidate` | 找到数值上为负的候选，但在给定分母上限内未能精确认证。 | `3` |
| 输入或运行错误 | 参数、文件或计算过程出错。 | `1` |

## 数学与限制

程序通过精确 Delaunay 单纯剖分计算内部二阶矩和边界格点测度。对每个 PL 分支，认证阶段在每个源单纯形内枚举其 dominance cell 的有理顶点，再用精确单纯剖分积分对应的 affine 分支。

高维中不存在二维 Donaldson“只需一个简单折痕即可完备判定”的同样事实。因此，认证成功的负 witness 是严格的相对 K-不稳定性证明；`no_counterexample_found` 只说明在指定分支数和有限启发式搜索中没有找到候选，不能作为高维半稳定性的证明。

## 已覆盖基准

- $d$ 维单位立方体：$\ell_P=2d$，$M_\ell(\max\{x_1-1/2,0\})=1/4$。
- 标准 $d$-单形：$\ell_P=d(d+1)$。
- 2、3、5 维精确矩、边界测度、PL 重复分支和解析认证回归。
