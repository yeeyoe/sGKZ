# sGKZ 项目程序清单

本文按当前源码、CMake 构建目标和项目 README 汇总命令行程序、辅助脚本、测试程序及内部库。路径均相对于项目根目录。项目说明见 [README.md](README.md)、[QUICKSTART.md](QUICKSTART.md)、[K-stability/README.md](K-stability/README.md)、[K-stability/Delzant/README.md](K-stability/Delzant/README.md)、[K-stability_highdim/README.md](K-stability_highdim/README.md) 和 [sGKZ_ndim/README.md](sGKZ_ndim/README.md)。

## 二维 Shortest GKZ 求解器

| 程序 | 源码 / 构建路径 | 功能 |
|---|---|---|
| `shortest_gkz` | `src/main.cpp`；`build/shortest_gkz` | 二维 shortest GKZ 向量求解器。支持一般点集 `--points FILE` 和格点多边形 `--polygon FILE --k INTEGER` 两种输入，计算 secondary polytope 中的最小模长点。使用 fully-corrective Frank–Wolfe / active-set 算法，并可执行精确 QP 与全局认证；还可计算离散 affine 函数 `ell_A`、输出 CSV 和绘图数据。 |
| `plot_results.py` | `plot_results.py` | 读取 `shortest_gkz --plot-prefix` 输出的 CSV，绘制 lower envelope、`psi_k`、subdivision 和 lifted subdivision 边界。输出交互式 HTML/SVG，可选生成 PNG。 |
| `plot_iterations.py` | `plot_iterations.py` | 读取 `shortest_gkz --verbose` 日志，绘制 active-set 大小、Frank–Wolfe gap、停机阈值和 stable-projection 事件。输出交互式 HTML，可选生成 PNG。 |

二维求解器及离散 `ell_A` 的详细说明、参数和输出格式见 [README.md](README.md)。绘图脚本需要 Python；PNG 导出还需要 Matplotlib。

## 二维相对 K-stability

本模块位于 `K-stability/`，与 `gkz_core` 独立。它计算连续边界测度对应的 `ell_P`，并使用 Donaldson 简单凸函数检测相对 K-不稳定性。详见 [K-stability/README.md](K-stability/README.md)。

| 程序 | 源码 / 构建路径 | 功能 |
|---|---|---|
| `k_stability` | `K-stability/main.cpp`；`build/K-stability/k_stability` | 读取二维格点多边形，精确计算连续的 `ell_P`，扫描 `g=max{<x,u>-t,0}` 并评估相对 DF 不变量。支持 `--check-line` 检查指定折痕、`--certify` 精确认证 witness、`--svg` 输出多边形和折痕图。 |
| `k_stability_search` | `K-stability/search_main.cpp`；`build/K-stability/k_stability_search` | 面积优先生成并检测固定顶点数的严格凸整点多边形。按 probe、confirm、final 阶段验证并将候选、状态、witness 和 `Q_P(g)^2` 保存到 SQLite。支持从数据库恢复搜索。 |
| `k_stability_search --backfill-q` | 同上 | 数据库维护模式：为已有精确认证 witness 补算精确及浮点 `Q_P(g)^2` 字段。 |
| `plot_candidate.py` | `K-stability/plot_candidate.py` | 按 canonical key 从搜索数据库读取候选，绘制顶点、`ell_P=0`、witness 折痕和奇异顶点。可保存 SVG，也可一并导出多边形顶点文件。 |
| `plot_area_q.py` | `K-stability/plot_area_q.py` | 从搜索数据库绘制面积与候选级 `Q_P(g)^2` 的联合分布，并按 `Q_P(g)^2 / V_P` 排名前 `N` 输出候选。 |

`ell_P` 是连续多边形与边界测度的对象；根目录 GKZ 求解器里的 `ell_A` 是离散点集上的对象，二者不同。`no_counterexample_found` 只表示当前扫描未找到 witness，不能单独视为半稳定证明。

## Delzant / smooth polygon 搜索

| 程序 | 源码 / 构建路径 | 功能 |
|---|---|---|
| `unstable_Delzant` | `K-stability/Delzant/main.cpp`、`K-stability/Delzant/unstable_delzant.cpp`；Release 二进制通常为 `K-stability/Delzant/build-release/K-stability/Delzant/unstable_Delzant` | 生成满足 smooth fan 条件的 Delzant 格点多边形，并搜索相对 K-不稳定候选。候选、detector 阶段和随机状态保存到独立 SQLite 数据库；认证成功时输出多边形、SVG 和文本结果。搜索是随机的，不是穷举。 |

参数、数据库和输出目录说明见 [K-stability/Delzant/README.md](K-stability/Delzant/README.md)。

## 高维 K-stability

| 程序 | 源码 / 构建路径 | 功能 |
|---|---|---|
| `k_stability_highdim` | `K-stability_highdim/src/main.cpp`；`build/K-stability_highdim/k_stability_highdim` | 从 polytope 输入列数推断维度，精确计算任意维 `ell_P`。支持 `--ell-only`、`--check-pl FILE` 精确验算给定凸 PL 函数，以及 `--pieces M` 搜索 `max(0,L_1,...,L_M)` 并尝试有理精确认证。 |

搜索模式使用数值积分和 differential evolution；通过精确认证的负 witness 可证明相对 K-不稳定。有限搜索未找到反例不构成高维半稳定证明。参数和输出定义见 [K-stability_highdim/README.md](K-stability_highdim/README.md)。

## 动态维 Shortest GKZ 求解器

| 程序 | 源码 / 构建路径 | 功能 |
|---|---|---|
| `sgkz` | `sGKZ_ndim/src/main.cpp`；`sGKZ_ndim/build/sgkz` | 独立的动态维 shortest GKZ 求解器。支持 `--points FILE` 和 `--polytope FILE --k INTEGER`，从输入列数推断维数 `d>=2`，可枚举 `kP ∩ Z^d`、运行 shortest GKZ 求解、输出动态列 CSV，并通过 `--ell` 输出 affine `ell_A`。 |

与根目录二维程序不同，该子项目的 README 说明它不生成绘图文件；`--max-points` 默认限制格点数为 1,000,000，传 `0` 可取消限制。详见 [sGKZ_ndim/README.md](sGKZ_ndim/README.md)。

## 测试程序

以下是 CMake/CTest 使用的测试二进制，不属于日常数值计算 CLI：

| 测试程序 | 覆盖内容 |
|---|---|
| `gkz_tests` | 根目录二维几何、GKZ 向量、active-set QP、`ell_A`、精确认证及绘图 CSV 回归。 |
| `oracle_kernel_consistency_test` | 检查数值 regular-triangulation oracle 与精确 oracle 的支撑值一致性。 |
| `k_stability_tests` | 二维 `ell_P`、矩、边界测度、相对 DF、折痕扫描和认证。 |
| `k_stability_search_tests` | 候选生成、搜索状态、SQLite 持久化、detector profile 和 Q 值。 |
| `unstable_delzant_tests` | Delzant/smooth polygon 生成与搜索逻辑。 |
| `k_stability_highdim_tests` | 高维矩、边界测度、PL 函数计算和认证。 |
| `sgkz_tests` | 动态维 shortest GKZ 核心逻辑。 |

Python 回归测试位于 `tests/`：`plot_results_test.py`、`plot_iterations_test.py`、`plot_area_q_test.py`。测试命令通常使用 `ctest --test-dir build --output-on-failure`；动态维独立构建则使用 `ctest --test-dir sGKZ_ndim/build --output-on-failure`。

## 内部 CMake 库

这些是供可执行程序链接的实现库，不是独立命令行程序：

| 库 | 功能 |
|---|---|
| `gkz_core` | 根目录二维 Shortest GKZ 的几何、regular triangulation、QP 和精确认证。 |
| `k_stability_core` | 二维 `ell_P`、DF 计算、simple-convex 扫描与认证。 |
| `k_stability_search_core` | 二维面积优先候选生成、搜索检测流程与 SQLite 状态管理。 |
| `k_stability_highdim_core` | 高维 polytope 矩、Delaunay 剖分、PL 评估及认证。 |
| `sgkz_core` | `sGKZ_ndim` 子项目的动态维 shortest GKZ 核心。 |

## 其他脚本

`papers/mac.sh` 是论文构建脚本，依次运行 `pdflatex`、`biber` 和多次 `pdflatex`；它不参与 GKZ 或 K-stability 数值计算。
