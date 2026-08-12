# Shortest GKZ 向量求解器

## 目录

- [依赖与编译](#依赖与编译)
- [快速调用](#快速调用)
- [输入文件](#输入文件)
  - [生成 Wang-Zhou 多边形](#生成-wang-zhou-多边形)
- [主算法](#主算法)
  - [第 1 步：读取和预处理点集](#第-1-步读取和预处理点集)
  - [第 2 步：regular triangulation oracle](#第-2-步regular-triangulation-oracle)
  - [第 3 步：初始化 active set](#第-3-步初始化-active-set)
  - [第 4 步：在 active hull 上完全校正](#第-4-步在-active-hull-上完全校正)
  - [第 5 步：Frank--Wolfe oracle 和停止判据](#第-5-步frank--wolfe-oracle-和停止判据)
  - [第 6 步：最终精确 QP 和全局认证](#第-6-步最终精确-qp-和全局认证)
  - [第 7 步：计算 $\ell_A$ 和写出结果](#第-7-步计算-ell_a-和写出结果)
  - [使用的外部算法汇总](#使用的外部算法汇总)
- [主程序参数](#主程序参数)
  - [输入模式参数](#输入模式参数)
  - [求解参数](#求解参数)
  - [文件输出参数](#文件输出参数)
- [终端输出解释](#终端输出解释)
- [CSV 输出解释](#csv-输出解释)
  - [`--output FILE`](#--output-file)
  - [`--plot-prefix PREFIX`](#--plot-prefix-prefix)
- [$\ell_A$ 与 relative Chow-semistability](#ell_a-与-relative-chow-semistability)
- [绘图脚本](#绘图脚本)
- [迭代历史折线图](#迭代历史折线图)
- [完整调用示例](#完整调用示例)
  - [例 1：一般六点点集，数值求解加精确认证](#例-1一般六点点集数值求解加精确认证)
  - [例 2：较大的 `examples/my`, $k=8$ 探索性计算](#例-2较大的-examplesmy-k8-探索性计算)
  - [例 3：允许较大的 exact active set](#例-3允许较大的-exact-active-set)
- [退出状态](#退出状态)
- [测试设计](#测试设计)
- [当前性能边界](#当前性能边界)
- [K-stability：$\ell_P$ 与相对 K-不稳定性检测](#k-stabilityell_p-与相对-k-不稳定性检测)

本目录包含一个二维求解器，用于计算 secondary polytope 的最小模长点

$$
\sigma_A=\operatorname*{argmin}_{g\in\Sigma(A)}\|g\|.
$$

这里 $A=(a_1,\ldots,a_N)\subset\mathbb Z^2$ 是标记整数点集，
$Q=\operatorname{conv}(A)$，$\Sigma(A)\subset\mathbb R^N$ 是 $A$ 的
secondary polytope。程序也支持论文中的特殊点集

$$
A_k=kP\cap\mathbb Z^2,
$$

其中输入文件给出格多边形 $P$，程序负责枚举 $A_k$。

对于三角剖分 $T$，本程序采用总质量归一化的 GKZ 向量。若
$\operatorname{area2}$ 表示二倍面积，则

$$
g_T(a_i)
=\frac{\sum_{\tau\in T,\,a_i\in\tau}\operatorname{area2}(\tau)}
{3\operatorname{area2}(Q)}.
$$

因此 $\sum_i g_T(a_i)=1$。程序同时保存分母归一化之前的整数面积分子，
以便最终进行精确有理数认证。

## 依赖与编译

macOS 上需要 Apple Command Line Tools、CMake、Ninja、CGAL 和 Eigen：

```bash
xcode-select --install
brew install cmake ninja cgal eigen
```

Release 构建及测试：

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

生成的主程序是：

```text
build/shortest_gkz
```

## 快速调用

一般点集模式：

```bash
./build/shortest_gkz \
  --points examples/six_points.points \
  --output build/six_points.csv \
  --verbose
```

$A_k=kP\cap\mathbb Z^2$ 模式：

```bash
./build/shortest_gkz \
  --polygon examples/unit_square.polygon \
  --k 8 \
  --output build/square_k8.csv \
  --plot-prefix build/square_k8
```

查看内置帮助：

```bash
./build/shortest_gkz --help
```

## 输入文件

必须且只能选择一种输入模式：

```text
shortest_gkz --points FILE [options]
shortest_gkz --polygon FILE --k INTEGER [options]
```

每个有效行包含两个整数坐标 `x y`。空行和 `#` 之后的内容被忽略。例如：

```text
# counterclockwise polygon vertices
0 0
4 1
3 2
0 3
-1 1
```

在 `--points` 模式中：

- 文件中的所有点组成标记点集 $A$；
- 点的输入顺序会保留，所有向量坐标和 CSV 行都按该顺序输出；
- 程序拒绝重复点和二维凸包面积为零的输入。

在 `--polygon FILE --k INTEGER` 模式中：

- 文件给出格多边形 $P$ 的顶点，顺时针或逆时针均可；
- 程序先计算顶点凸包，再将其放大为 $kP$；
- 当前实现扫描 $kP$ 的整数包围盒，并用凸多边形包含测试枚举
  $A_k=kP\cap\mathbb Z^2$；
- 内部点坐标是 $kP$ 中的整数坐标，绘图时再除以 $k$ 映回 $P$。

### 生成 Wang-Zhou 多边形

脚本 `generate_wang_zhou.py` 接受正整数 $a$，并默认写出主程序可直接读取的
`examples/Wang_Zhou_a{a}`。例如：

```bash
python3 generate_wang_zhou.py 5
./build/shortest_gkz --polygon examples/Wang_Zhou_a5 --k 1
```

不提供位置参数时，脚本会交互式询问 $a$：

```bash
python3 generate_wang_zhou.py
a = 5
```

每次运行都会覆盖同一个 $a$ 对应的 `examples/Wang_Zhou_a{a}`。参数 $a$ 必须是正整数；
这是因为主程序要求格点坐标，而 $a=0$ 会使给定的两个右侧顶点重合。
脚本按如下顺序写入九个顶点：

$$
(0,a+20),(1,a+20),(2,a+19),(3,a+17),(4,a+14),
(5,a+10),(7,a),(7,-a),(0,-a).
$$

## 主算法

主算法是 fully corrective Frank--Wolfe / active-set 方法。它不枚举
$\Sigma(A)$ 的所有顶点，也不枚举所有 regular triangulations。程序只在
需要时调用线性优化 oracle 产生新的 GKZ 顶点。

设当前 active set 为

$$
\mathcal V_t=\{v_1,\ldots,v_m\}\subset\operatorname{Vert}(\Sigma(A)).
$$

`active set` 指已经由 oracle 找到并且当前系数非零的 GKZ 向量集合，不是
$A$ 的子集。

### 第 1 步：读取和预处理点集

程序检查点的唯一性，使用单调链 convex hull 算法计算 $Q$，并以
128 位整数计算叉积和二倍面积。`--polygon` 模式还会枚举 $A_k$。

这部分是本项目自行实现的整数二维几何代码，没有调用外部凸包程序。

### 第 2 步：regular triangulation oracle

给定 height 向量 $h\in\mathbb R^A$，oracle 求

$$
v_h\in\operatorname*{argmin}_{v\in\Sigma(A)}\langle h,v\rangle.
$$

按照 secondary polytope 的基本性质，这等价于计算由 $h$ 诱导的 lower
regular triangulation，再返回其 GKZ 向量。

实现使用 `CGAL::Regular_triangulation_2`。传给 CGAL 的 weighted point 为

$$
\bigl(a_i,\,\|a_i\|^2-h_i\bigr),
$$

从而 CGAL 的 power lifting 对应 height $h_i$。数值主循环的核类型是
`CGAL::Exact_predicates_inexact_constructions_kernel`：谓词直接使用
`double` 过滤算术，每次调用比精确核快一个数量级以上。由于数值核在
secondary fan 墙面附近可能返回非最小化三角剖分（从而低估 gap），
程序在数值 gap 首次满足停止判据后切换到
`CGAL::Exact_predicates_exact_constructions_kernel` 继续迭代（精确
endgame）；只有精确 oracle 计算出的 gap 也满足判据时才报告收敛。
最终认证中的 rational height 会先减去 constant gauge、清除公分母
并约去公因子，再以 `CGAL::Gmpz` 整数 heights 传入精确核 oracle。

CGAL 返回有限三角形后，程序自行累加每个顶点相邻三角形的二倍面积，构造
GKZ 向量，并检查：

$$
\sum_{\tau\in T}\operatorname{area2}(\tau)
=\operatorname{area2}(Q),
\qquad
\sum_i g_T(a_i)=1.
$$

一次 oracle 调用只返回一个 minimizing refinement 和一个 GKZ 向量；程序不
在该调用中枚举所有具有相同线性目标值的 refining regular triangulations。
fully corrective 主循环会在后续调用中按需加入新的 GKZ 向量。

### 第 3 步：初始化 active set

程序先以零向量调用一次 oracle：

$$
v_1\in\operatorname*{argmin}_{v\in\Sigma(A)}\langle0,v\rangle.
$$

CGAL 在可能存在多个最优三角剖分时返回其中一个 refinement。初始化为

$$
\mathcal V_0=\{v_1\},\qquad x_0=v_1.
$$

### 第 4 步：在 active hull 上完全校正

对于固定的 $\mathcal V_t$，程序求其凸包上的最小模长点

$$
x_t=\sum_{i=1}^m\lambda_i v_i,
\qquad
\lambda_i\geq0,
\qquad
\sum_{i=1}^m\lambda_i=1,
$$

其中 $\lambda$ 解数值凸二次规划

$$
\min_{\lambda\geq0,\,\mathbf1^T\lambda=1}
\frac12\lambda^TG\lambda,
\qquad
G_{ij}=\langle v_i,v_j\rangle.
$$

这部分是本项目自行实现的 fully corrective active-set QP：

1. 对当前正系数 working set 建立 equality-constrained KKT 系统。
2. 使用 Eigen 的 `CompleteOrthogonalDecomposition` 求解。该秩揭示分解能
   处理相容的秩亏系统，因此数值阶段允许 active GKZ 向量仿射冗余。
3. 若求得的 affine minimizer 含负系数，从当前可行系数沿线段移动到
   simplex 边界，删除首先变为零的变量，再解较小的 KKT 系统。
4. 若 working set 外仍有变量违反 KKT 梯度条件，将其中梯度最小者加入
   working set，继续校正。
5. 完成校正后，系数不大于 `1e-15` 的 GKZ 向量从 active set 中删除。

加入一个新 oracle 顶点后，程序不是只沿新顶点做一次线搜索，而是重新优化
active set 中的所有系数；这就是 `fully corrective` 的含义。

### 第 5 步：Frank--Wolfe oracle 和停止判据

用当前候选点 $x_t$ 再调用 regular triangulation oracle：

$$
w_t\in\operatorname*{argmin}_{v\in\Sigma(A)}\langle x_t,v\rangle.
$$

定义 Frank--Wolfe gap

$$
\delta_t=\|x_t\|^2-\langle x_t,w_t\rangle.
$$

程序在

$$
\delta_t\leq
\texttt{absolute-tolerance}
+\texttt{tolerance}\cdot\|x_t\|^2
$$

时报告 `converged=true`，但有一个重要前提：数值主循环的 oracle 使用
`double` 谓词核，可能低估 gap，所以数值 gap 首次满足上式后，程序切换
到精确核 oracle 继续迭代（精确 endgame），只有精确 oracle 的 gap 也
满足上式才真正判收敛。判收敛前 gap 仍大时，将 $w_t$ 加入 active set，
扩充 Gram 矩阵，回到第 4 步。若精确认证或精确 oracle 返回与已有
active GKZ 向量精确相同的向量，但 gap 仍超过停止阈值，程序将其视为
数值 active QP 与 oracle 不一致并报错。重复判断直接比较整数面积分子，
不使用浮点容差；同一 $A$ 下所有 GKZ 向量具有共同分母
$3\operatorname{area2}(Q)$，所以这个判断等价于 GKZ 向量精确相等。

这个 gap 同时给出后验误差界

$$
\|x_t-\sigma_A\|\leq\sqrt{2\delta_t}.
$$

精确认证失败且 exact oracle 找到更小内积顶点时，该顶点会作为下一步
Frank--Wolfe 方向加入 active set（且不被系数阈值剪枝），主循环继续
迭代直到认证通过或达到迭代上限。

### 第 6 步：最终精确 QP 和全局认证

只有数值停止判据已经满足、没有使用 `--no-exact`，并且 active set 大小
没有超过 `--exact-max-active` 时，程序才执行这一阶段。

程序使用整数面积分子精确构造 Gram 矩阵，并调用 CGAL 的
`solve_quadratic_program` 求解

$$
\min_{\lambda_i\geq0,\,\sum_i\lambda_i=1}
\left\|\sum_i\lambda_i v_i\right\|^2.
$$

QP 的输入系数使用 `CGAL::Gmpz`，解系数是精确商数，再转换为
`CGAL::Gmpq`。该求解器允许 Gram 矩阵半正定，也允许 active GKZ 向量
仿射相关；它不要求最优凸组合系数唯一。程序随后检查

$$
\lambda_i\geq0,
\qquad
\sum_i\lambda_i=1,
$$

并精确构造

$$
\sigma=\sum_i\lambda_i v_i.
$$

最后把 $\sigma$ 的未归一化 rational height 清除公分母，得到等价的整数
heights，再调用一次 CGAL regular triangulation oracle，精确计算

$$
H_A(\sigma)=\min_{v\in\Sigma(A)}\langle\sigma,v\rangle.
$$

只有精确验证

$$
H_A(\sigma)=\|\sigma\|^2
$$

后，才输出 `exact_certified=true`。这个等式证明 $\sigma$ 是整个
$\Sigma(A)$ 上的最小模长点，而不只是 active hull 上的最小点。

数值主循环仍使用归一化后的 `double` GKZ 坐标。虽然数学上可以把整个
QP 乘以 $D^2$ 改写到整数面积分子空间，但未归一化分子的大小会随 $k^2$
增长，并会改变 Eigen KKT 残差和 correction tolerance 的数值尺度；当前
版本只把整数分子用于精确存储、重复判断和 exact QP。这样保留了数值 QP 的
稳定尺度，同时把必须精确的步骤延迟到认证阶段。

### 第 7 步：计算 $\ell_A$ 和写出结果

主算法结束后，程序以精确有理数解一个 $3\times3$ affine moment 系统，
得到

$$
\ell_A(x,y)=a+bx+cy
$$

及其限制 $\ell_A|_A$。然后进行 $\ell_A$ 与 shortest GKZ 的比较，写出
终端结果、普通 CSV 和可选绘图数据。

若要求绘图数据，程序在最终 height 上再次调用 regular triangulation
oracle，计算 lower convex envelope $\sigma_A^\vee$。有精确证书时，插值和
共面判断使用 `Gmpq`；只有写 CSV/SVG 时才转换为浮点数。程序自行把共面的
三角形合并为 $S(\sigma_A)$ 的子多面体。

### 使用的外部算法汇总

| 子任务 | 实现 |
| --- | --- |
| 二维 convex hull、格点枚举、面积与 GKZ 累加 | 本项目自行实现，整数运算 |
| regular triangulation / 线性优化 oracle（数值主循环） | CGAL `Regular_triangulation_2`，`Exact_predicates_inexact_constructions_kernel` |
| regular triangulation oracle（精确 endgame、认证、精确绘图路径） | CGAL `Regular_triangulation_2`，`Exact_predicates_exact_constructions_kernel` |
| 数值秩亏 KKT 线性系统 | Eigen `CompleteOrthogonalDecomposition` |
| fully corrective active-set QP 外层逻辑 | 本项目自行实现 |
| 最终精确半正定 QP | CGAL `solve_quadratic_program`；Gram 矩阵先用 128 位整数累加再转 `Gmpz` |
| 精确整数和有理数 | CGAL `Gmpz`、`Gmpq`，底层使用 GMP |
| lower envelope 插值、共面三角形合并 | 本项目自行实现 |
| 交互式 3D 页面 | Python 标准库生成 HTML，浏览器加载 Plotly |
| 可选静态 PNG | Matplotlib |

## 主程序参数

### 输入模式参数

| 参数 | 默认值 | 含义 |
| --- | ---: | --- |
| `--points FILE` | 无 | 从 `FILE` 读取一般标记点集 $A$。不能与 `--polygon` 同时使用。 |
| `--polygon FILE` | 无 | 从 `FILE` 读取格多边形 $P$。必须同时给出正整数 `--k`。 |
| `--k INTEGER` | 无 | 构造 $A_k=kP\cap\mathbb Z^2$ 的层数，只用于 `--polygon` 模式。 |

### 求解参数

| 参数 | 默认值 | 含义 |
| --- | ---: | --- |
| `--tolerance VALUE` | `1e-11` | Frank--Wolfe gap 的相对容差，必须非负。 |
| `--absolute-tolerance VALUE` | `1e-14` | Frank--Wolfe gap 的绝对容差，必须非负。 |
| `--prune-tolerance V` | `1e-15` | QP 后删除系数不超过 $V$ 的 active 顶点，必须非负。 |
| `--max-iterations N` | `500` | 最多进行 $N$ 次 active-set 扩张，必须为非负整数。 |
| `--exact-max-active N` | `128` | 最终 active set 超过 $N$ 时跳过精确 QP；`0` 表示不设上限。 |
| `--no-exact` | 关闭 | 完全跳过最终精确 QP 和精确 oracle 认证。 |
| `--projection` | 关闭 | 启用实验性的稳定暴露面投影；默认路径不受影响。 |
| `--projection-window N` | `32` | 自动进入投影模式时比较的 gap 历史窗口。 |
| `--projection-stall-ratio R` | `0.90` | 当前 gap 与窗口前 gap 的最小比值；越接近 `1` 越要求平台化。 |
| `--projection-relative-gap R` | `1e-2` | 只有当前 gap 不超过初始 gap 的该比例时才允许进入投影模式。 |
| `--projection-rank-stall-window N` | `64` | 后期普通顶点连续不增加 affine rank 的观测次数。 |
| `--projection-rank-tolerance R` | `1e-11` | 增量 affine QR 判断新方向的相对容差。 |
| `--verbose` | 关闭 | 每次 gap 检查向标准错误输出 `iteration`、`active`、`norm2` 和 `gap`。 |

`--max-iterations N` 的准确含义如下：初始化 active set 后，主循环检查
`iteration=0,1,\ldots,N`，但最多只进行 $N$ 次“加入新 GKZ 顶点并完全
校正”的更新。因此非收敛时最后输出 `iterations=N`。初始化、每次 gap
检查、精确认证和绘图数据还可能各自调用 oracle，所以该参数不是程序全部
regular triangulation 调用次数的上限。

`--exact-max-active` 只控制数值收敛后的精确认证，不限制数值 active set
增长，也不影响 `converged`。精确 QP 使用稠密 Gram 矩阵；变量数增大时，
有理数运算和整数位数增长可能非常昂贵，因此不建议在不了解资源开销时直接
使用 `--exact-max-active 0`。

### 实验性暴露面投影

`--projection` 在前期保持原 fully-corrective Frank--Wolfe 主循环。程序仍以
当前 QP 解 $x_t$ 调用普通 oracle 得到 $w_t$ 并计算唯一用于终止的 FW gap。
当 gap 已相对初始值足够小，且在 `--projection-window` 步内下降不明显时，
才新建稳定顶点集；它只接收此后普通 oracle 返回的 $w_t$，建立 affine hull
$L_t$，不包含触发前历史顶点，也不进行周期投影或额外 oracle 调用。

$$
p=\operatorname{Proj}_{L_t}(0).
$$

收集在 rank 达到二维 facet 上限 $|A|-4$，或连续
`--projection-rank-stall-window` 个后期顶点不增加 rank 时停止。仅此时计算
$p$，并以 $p$ 做一次精确 oracle 得到 $v_p$。同一轮的 $w_t$ 和去重后的
$v_p$ 一起加入 active set，只做一次 QP 校正，得到 final QP solution；随后用
普通 oracle 计算 final gap，算法停止。$p$ 从不作为 QP 顶点；程序不计算
`gap(p)`，也不检查 $p\in\Sigma(A)$。

因此 stable projection 仅是探索性高度向量，不是 shortest GKZ 或
$p\in\Sigma(A)$ 的证书。它的 lower envelope 给出

$$
\widehat\psi_k=\operatorname{area2}(P)k^3p^\vee-2k,
$$

用于观察 $-\Theta_P$ 的近似形状。若正常 FW 在收集完成前满足终止条件，或到达
最大迭代数仍未稳定，则不产生 stable projection，原终止/输出语义不变。

建议在探索 $\widehat\psi_k$ 时保留 `--no-exact`，让 stable rank 判据而非 exact
认证决定停止时刻：

```bash
./build/shortest_gkz --polygon examples/d5_a70 --k 8 \
  --max-iterations 10000 --no-exact --projection \
  --plot-prefix results/d5_a70 --verbose \
  2> results/d5_a70_projection.log

python3 plot_results.py results/d5_a70_stable
```

历史 `hexgon` 记录和当前 `d5_a70 --k 8` 基线都在约第 270 步首次满足默认的
32 步平台判据；这只用于给出默认值，不是对其他实例的收敛或几何保证。

### 文件输出参数

| 参数 | 默认值 | 含义 |
| --- | ---: | --- |
| `--output FILE` | 不写出 | 写出逐点 shortest GKZ 和 $\ell_A$ CSV。 |
| `--plot-prefix PREFIX` | 不写出 | 写出 surface、triangles、subdivision 和 $\ell_A$ 四个绘图 CSV。 |
| `--help` | - | 显示简要帮助并以状态码 `0` 退出。 |

`FILE` 和 `PREFIX` 都相对于启动程序时的当前目录解释。程序不会自动创建
父目录，应先执行例如 `mkdir -p build/results`。即使程序最终以状态码 `2`
退出，指定的 CSV 仍会保存最后一个数值近似。

## 终端输出解释

主程序把最终摘要写到标准输出。使用 `--verbose` 时，逐轮日志写到标准错误，
便于分别重定向。

| 字段 | 含义 |
| --- | --- |
| `points` | $|A|$，即 GKZ 向量的坐标数。 |
| `twice_area` | $\operatorname{area2}(Q)$。polygon 模式中 $Q=kP$，所以它等于 $k^2\operatorname{area2}(P)$。 |
| `level` | polygon 模式中的 $k$；points 模式不输出。 |
| `base_twice_area` | $\operatorname{area2}(P)$；只在 polygon 模式输出。 |
| `converged` | 精确 oracle 验证后的 gap 是否达到指定容差。它不等价于已有精确证书。 |
| `iterations` | 最后一次 gap 检查的编号，从 `0` 开始。 |
| `active_size` | 最终保留的正权重 GKZ 向量数。 |
| `norm_squared` | 数值候选点的 $\|x\|^2$。 |
| `gap` | 最后的 Frank--Wolfe gap $\delta(x)$。 |
| `l2_error_bound` | 后验界 $\sqrt{2\delta(x)}$。这是整个 GKZ 向量的 Euclidean norm 误差界。 |
| `exact_certified` | 是否完成精确 QP 并通过 exact oracle 全局等式。 |
| `exact_message` | 精确认证成功、失败或因 active set 过大而跳过的原因。使用 `--no-exact` 时不输出。 |
| `exact_norm_squared` | 精确的 $\|\sigma_A\|^2\in\mathbb Q$；仅在认证成功时输出。 |
| `projection_enabled` | 是否通过 `--projection` 请求 stable-projection 模式。 |
| `projection_start_iteration` | 开始收集稳定普通顶点的迭代号；`-1` 表示未触发。 |
| `stable_projection_available` | 是否完成 rank 判据并生成 stable projection。 |
| `stable_projection_norm_squared` | 可用时的 $\|p\|^2$。 |
| `stable_projection_rank` / `observations` | stable 顶点 affine hull 的数值 rank 与观测顶点数。 |
| `stable_projection_stop_reason` | `rank_limit`、`rank_stall` 或未完成原因。 |
| `final_qp_gap` / `final_qp_norm_squared` | stable 投影模式末轮 QP 解的普通 FW gap 与范数。 |
| `ell_A(x,y)` | $\ell_A$ 的精确有理数 affine 表达式。 |
| `ell_A_constant` | affine 斜率是否为零。 |
| `ell_A_sigma_comparison` | `exact` 表示在 $\mathbb Q$ 中比较，`numerical` 表示用数值容差比较。 |
| `ell_A_equals_shortest_GKZ` | $\ell_A|_A$ 是否与当前 shortest GKZ 结果相等。无证书时仅是数值判断。 |
| `relative_Chow_semistable` | 有精确证书时按 $\ell_A=\sigma_A$ 输出 `true/false`；无证书时输出 `unverified`。 |
| `ell_A_on_A` | 逐点列出 $\ell_A|_A$。大点集且指定文件输出时，终端只提示保存路径。 |
| `output` | `--output` 实际写出的文件路径。 |
| `plot_surface` 等 | `--plot-prefix` 实际写出的四个 CSV 路径。 |

判断“严格证明了 shortest GKZ”时，应同时看到：

```text
converged=true
exact_certified=true
exact_message=Exact convex QP and equality H(x) = ||x||^2 verified.
```

`converged=true, exact_certified=false` 只表示取得带 gap 误差界的数值近似。

## CSV 输出解释

### `--output FILE`

普通结果 CSV 的列为：

| 列 | 含义 |
| --- | --- |
| `x,y` | 点的内部整数坐标。polygon 模式中它表示 $kP$ 中的 $(x,y)$，对应 $P$ 中的 $(x/k,y/k)$。 |
| `sigma` | 数值算法得到的 shortest GKZ 坐标。即使有精确证书，这一列仍保留数值迭代结果。 |
| `sigma_exact` | 精确认证得到的有理数字符串；没有证书时为空。 |
| `ell_A` | $\ell_A|_A$ 的浮点值。 |
| `ell_A_exact` | $\ell_A|_A$ 的精确有理数字符串。 |

### `--plot-prefix PREFIX`

该参数生成：

| 文件 | 内容 |
| --- | --- |
| `PREFIX_surface.csv` | final QP solution 的列 `x,y,sigma,sigma_vee,psi`。坐标使用绘图坐标；polygon 模式下已经除以 $k$。 |
| `PREFIX_triangles.csv` | 列 `i,j,l`，是 lower regular triangulation 的三角形顶点行号。 |
| `PREFIX_subdivision.csv` | 列 `cell,vertex,x,y`，按子多面体及其边界顶点顺序保存 $S(\sigma_A)$。 |
| `PREFIX_ell.csv` | 列 `x,y,ell_A,ell_A_exact`，坐标与 surface CSV 一致。 |

stable projection 可用时，额外生成 `PREFIX_stable_surface.csv`、
`PREFIX_stable_triangles.csv` 与 `PREFIX_stable_subdivision.csv`。它们以 $p$ 为
height，`psi` 列为探索性的 $\widehat\psi_k$；可运行
`python3 plot_results.py PREFIX_stable` 绘图。

`surface.csv` 中：

- `sigma` 是点上的 height；有精确证书时，它由精确分数转换为 `double`；
- `sigma_vee` 是该 height 的 lower convex envelope 在点上的值；
- hidden point 可能满足 `sigma_vee < sigma`；
- 理论上已认证的 shortest GKZ 满足 $\sigma_A^\vee|_A=\sigma_A$，但 CSV
  最终写的是十进制浮点数，`sigma` 与 `sigma_vee` 仍可能有舍入量级差异；
- 一般 points 模式的 `psi` 为空；polygon 模式使用

$$
\psi_k=\operatorname{area2}(P)k^3\sigma_k^\vee-2k.
$$

这里二维情形中 $2V_P=\operatorname{area2}(P)$。

## $\ell_A$ 与 relative Chow-semistability

程序按照

$$
\frac1{|Q|}\int_Q h\,dx
=\sum_{u\in A}h(u)\ell_A(u)
$$

对所有 affine 函数 $h$ 解出 $\ell_A(x,y)=a+bx+cy$。polygon 模式的表达式
使用 $P$ 上的归一化坐标 $(x/k,y/k)$，这不会改变限制值
$\ell_A|_{A_k}$。

有精确 shortest GKZ 证书时，$\ell_A|_A$ 与 $\sigma_A$ 在 $\mathbb Q$ 中
逐坐标比较。没有证书时使用 $10^{-9}$ 的无穷范数容差作探索性比较，但
`relative_Chow_semistable` 保持为 `unverified`，不会把数值相等当作证明。

## 绘图脚本

先由 C++ 程序生成绘图 CSV：

```bash
./build/shortest_gkz \
  --polygon examples/unit_square.polygon \
  --k 8 \
  --plot-prefix build/square_k8
```

再运行：

```bash
python3 plot_results.py build/square_k8
```

points 模式默认输出：

- `build/square_k8_sigma_vee.html`：可旋转、缩放和平移的 3D 视图；
- `build/square_k8_subdivision.svg`：$S(\sigma_k)$ 的二维图。

对于 polygon 模式，默认输出：

- `build/square_k8_psi_k.html`：可旋转、缩放和平移的 $\psi_k$ 3D 视图；
- `build/square_k8_subdivision.svg`：$S(\sigma_k)$ 的二维图。

如果需要同时输出 `sigma_k^vee`，显式加入 `--sigma-vee`：

```bash
python3 plot_results.py build/square_k8 --sigma-vee
```

此时才会额外生成 `build/square_k8_sigma_vee.html`。

生成的每个 3D 视图都会在曲面上叠加 subdivision 每个子多面体边界的
lifted break lines：边界顶点的高度取自同一张 `surface.csv`，因此线段与
曲面使用同一组提升值。交互式 HTML 可以直接拖动旋转；PNG 是固定视角的
附加输出。

HTML 使用浏览器端 Plotly CDN，需要浏览器能够访问对应网络资源，但 Python
端不需要安装 Plotly。静态 PNG 需要 Matplotlib：

```bash
python3 -m venv .venv
.venv/bin/python -m pip install numpy matplotlib
.venv/bin/python plot_results.py build/square_k8 --png
```

绘图脚本参数：

| 参数 | 含义 |
| --- | --- |
| `PREFIX` | 必选位置参数，必须与 C++ 的 `--plot-prefix` 完全相同。 |
| `--sigma-vee` | 生成 $\sigma_A^\vee$ 3D 视图；polygon 模式默认关闭。 |
| `--no-sigma-vee` | 不生成 $\sigma_A^\vee$ 视图。 |
| `--no-psi` | 不生成 $\psi_k$ 视图。points 模式本来也不会生成。 |
| `--no-subdivision` | 不生成 subdivision 图。 |
| `--no-edges` | 不在 3D 曲面上叠加 lifted subdivision 边界线。 |
| `--z-scale FACTOR` | 仅改变 3D 图像的纵向视觉比例，将 z 方向显示拉伸 `FACTOR` 倍；默认 `1`。 |
| `--color` | 用函数值驱动 3D 曲面的 Viridis 颜色映射，并为不同 lifted 边界线使用不同颜色；默认关闭。关闭时曲面为浅蓝色、边界线为深红色。 |
| `--png` | 在 HTML/SVG 之外额外生成固定视角 PNG。 |
| `--show` | 生成 PNG 后打开 Matplotlib 窗口显示。 |

`--z-scale` 只改变 HTML/PNG 的 3D 场景纵横比，不会修改任何 CSV 中的
`sigma`、`sigma_vee` 或 `psi`，也不会修改 HTML 中曲面、颜色映射、hover
提示或 lifted break lines 的 z 坐标。若原始坐标范围为
$(\Delta x,\Delta y,\Delta z)$，绘图时只把场景的显示比例设置为
$(\Delta x,\Delta y,\texttt{FACTOR}\,\Delta z)$。因此这是视觉上的纵向
拉伸，而不是对函数值作变换。`FACTOR` 必须是有限的正数，小于 $1$ 会压缩
纵向显示，大于 $1$ 会拉伸纵向显示。

在 polygon 模式的默认设置下，唯一的 3D 曲面是 `psi_k`，所以
`--z-scale` 只针对 `psi_k`。若同时指定 `--sigma-vee`，同一个显示比例
也会分别作用于 `sigma_k^vee` 和 `psi_k` 两张曲面；它们的函数值仍不变。

例如将 z 方向的视觉比例拉伸 $10^4$ 倍并保留 lifted break lines：

```bash
python3 plot_results.py build/square_k8 --z-scale 10000
```

如只想看曲面而不显示提升后的 subdivision 边界：

```bash
python3 plot_results.py build/square_k8 --z-scale 10000 --no-edges
```

默认 3D 曲面使用固定的浅蓝色，不用颜色表示函数值；lifted break lines
使用固定的深红色。如需用函数值驱动颜色映射并使用多色 lifted 边界线：

```bash
python3 plot_results.py build/square_k8 --z-scale 10000 --color
```

只生成 subdivision：

```bash
.venv/bin/python plot_results.py build/square_k8 \
  --no-sigma-vee \
  --no-psi
```

## 迭代历史折线图

主程序使用 `--verbose` 时，每次 gap 检查都会把 `iteration`、active set
大小、`norm2` 和 `gap` 写入标准错误。先保存日志：

```bash
./build/shortest_gkz \
  --polygon examples/Wang_Zhou_a5 \
  --k 1 \
  --max-iterations 500 \
  --no-exact \
  --verbose \
  2> results/Wang_Zhou_a5.log
```

再生成交互式折线图：

```bash
python3 plot_iterations.py results/Wang_Zhou_a5.log
```

默认生成 `results/Wang_Zhou_a5_iterations.html`。图中包含两个共享 iteration
横轴的面板：上方显示 active set 大小，下方显示 Frank--Wolfe `gap` 的对数
坐标，并叠画实际的停机阈值曲线
`absolute_tolerance + tolerance * ||σ_t||^2`（蓝色实线）。HTML 使用 Plotly
CDN，不要求 Python 安装 Plotly。

启用 `--projection --verbose` 后，日志还会包含独立的 `projection event=...`
行；常规 `iteration=...` 行的格式不变。`plot_iterations.py` 会把 stable
projection 迭代标为绿色虚线，旧的无投影日志无需迁移。

迭代绘图参数：

| 参数 | 含义 |
| --- | --- |
| `LOG` | 必选位置参数，由 `shortest_gkz --verbose` 写出的标准错误日志。 |
| `--output-prefix PREFIX` | 指定输出路径但不写 `.html`/`.png` 后缀。 |
| `--tolerance VALUE` | 相对 gap 容差，用于画停机阈值线（默认：1e-11）。 |
| `--absolute-tolerance VALUE` | 绝对 gap 容差，用于画停机阈值线（默认：1e-14）。 |
| `--png` | 额外生成静态 PNG，需要 Matplotlib。 |
| `--show` | 生成 PNG 后打开 Matplotlib 窗口。 |

例如同时生成 HTML 和 PNG：

```bash
.venv/bin/python plot_iterations.py \
  results/Wang_Zhou_a5.log \
  --tolerance 1e-12 \
  --absolute-tolerance 1e-14 \
  --png
```

## 完整调用示例

### 例 1：一般六点点集，数值求解加精确认证

```bash
mkdir -p build/results
./build/shortest_gkz \
  --points examples/six_points.points \
  --tolerance 1e-12 \
  --absolute-tolerance 1e-14 \
  --max-iterations 200 \
  --output build/results/six_points.csv \
  --verbose \
  > build/results/six_points.out \
  2> build/results/six_points.iterations.log
status=$?
echo "exit_status=$status"
```

标准输出中的关键结果应包括：

```text
points=6
converged=true
active_size=3
exact_certified=true
exact_message=Exact convex QP and equality H(x) = ||x||^2 verified.
ell_A_sigma_comparison=exact
```

### 例 2：较大的 `examples/my`, $k=8$ 探索性计算

该例默认的 500 次更新不足，使用更大的迭代上限，并暂时关闭可能昂贵的
精确 QP：

```bash
mkdir -p build/results
./build/shortest_gkz \
  --polygon examples/my \
  --k 8 \
  --tolerance 1e-11 \
  --absolute-tolerance 1e-14 \
  --max-iterations 3000 \
  --no-exact \
  --output build/results/my_k8.csv \
  --plot-prefix build/results/my_k8
echo "exit_status=$?"
```

生成交互视图：

```bash
python3 plot_results.py build/results/my_k8
```

此时 `relative_Chow_semistable=unverified` 是预期行为，因为使用了
`--no-exact`。

### 例 3：允许较大的 exact active set

下面的参数允许至多 600 个 active GKZ 变量进入最终精确 QP：

```bash
./build/shortest_gkz \
  --polygon examples/my \
  --k 8 \
  --max-iterations 3000 \
  --exact-max-active 600
```

这不会使数值算法更快，也不保证精确阶段在合理时间内完成。该选项只是取消
默认的 128 变量保护阈值。

## 退出状态

程序退出码为：

- `0`：达到数值停止条件；是否有精确证书还要检查 `exact_certified`；
- `1`：参数、输入文件、输出文件或运行过程中发生错误；
- `2`：达到 `--max-iterations` 上限但未达到数值停止条件。

退出码必须在主程序结束后立即读取：

```bash
./build/shortest_gkz --points examples/six_points.points
status=$?
echo "$status"
```

若使用管道，例如 `| tee run.log`，zsh 中应读取 `pipestatus`，否则 `$?`
通常只是最后一个管道命令的状态：

```bash
./build/shortest_gkz --points examples/six_points.points | tee run.log
echo "shortest_gkz_status=$pipestatus[1]"
```

## 测试设计

自动测试分为五个相互独立的层次：

1. 几何不变量：检查格点枚举、多边形面积、GKZ 总质量和一阶矩。
2. 正方形穷尽检验：将 oracle 输出与四边形的两种三角剖分逐一比较，
   并要求求解器精确输出 $(1/4,1/4,1/4,1/4)$。
3. 仿射冗余测试：枚举严格凸六边形的全部 14 个三角剖分。对应 GKZ
   向量位于至多 3 维的 secondary polytope 中，故必然仿射冗余；测试要求
   精确 QP 仍能给出非负凸组合，并通过 exact oracle 的全局等式验证。
4. 论文回归测试：验证五点例子的 $\ell_A=(y+5)/27$，并将六点例子的
   shortest GKZ 每个坐标与 `main.tex` 中记录的精确分数比较。
5. 绘图数据测试：生成二维 $P_k$ 的 surface、triangles、subdivision 和
   $\ell_A$ CSV，同时检查归一化坐标、hidden point lower envelope 和
   $\psi_k$ 公式。若 CMake 找到 Python 解释器，CTest 还会运行
   `plot_results` 回归测试，检查高度尺度化、lifted edges 以及两个关闭开关。

完整验证命令：

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure

cmake -S . -B build-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build-debug --parallel
ctest --test-dir build-debug --output-on-failure

cmake -S . -B build-ubsan -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=undefined -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=undefined"
cmake --build build-ubsan --parallel
ctest --test-dir build-ubsan --output-on-failure
```

## 当前性能边界

当前版本以正确性和可验证性为优先：

- 每次 oracle 调用都从头计算 regular triangulation；
- 数值主循环使用 `double` 谓词核（EPICK），比精确核约快一个数量级；
  收敛判定和认证仍使用精确核（EPECK）；
- oracle 批量插入已使用 CGAL range insert（内部 spatial sort），并把
  点坐标与 $\|a_i\|^2$ 按内核类型缓存；
- GKZ 向量和 numerical Gram matrix 使用稠密存储；
- active-set QP 使用稠密 KKT 分解；
- 最终 exact QP 使用稠密整数 Gram matrix。Gram 条目先用 128 位整数
  累加再转换为 `Gmpz`（构造开销已可忽略），但 CGAL 的精确 QP 求解
  本身在 active set 大时可能非常慢；
- polygon 模式通过扫描整个整数包围盒枚举格点。

因此，大 $k$ 探索建议先使用 `--no-exact` 和适当增大的
`--max-iterations`，并依据 `gap` 与 `l2_error_bound` 判断数值精度。精确
认证应在确认 active set 规模可承受后开启。

## K-stability：$\ell_P$ 与相对 K-不稳定性检测

独立子目录 [`K-stability/`](K-stability/README.md) 实现
`paper/K-stability.tex` 的两个目标：精确计算连续 Donaldson $\ell_P$
（与本程序的离散 $\ell_A$ 是不同的数学对象），并用 Donaldson 简单凸函数
事实检测二维格点多边形的相对 K-不稳定性。生成独立二进制
`build/K-stability/k_stability`，不链接 `gkz_core`，只读 `--polygon`
文件；详见 [`K-stability/README.md`](K-stability/README.md)。
