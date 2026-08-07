# sGKZ 项目交接文档

更新日期：2026-08-04

## -1. 2026-08-04 K-stability 模块

新增自包含子目录 `K-stability/`（二进制 `k_stability`），实现
`paper/K-stability.tex` 的两个目标，详见 `K-stability/README.md`：

1. **精确 $\ell_P$**：3×3 矩系统（带符号扇形二阶矩 + 边界格点测度
   $\gcd$ 权重），`Gmpq` Cramer 精确求解。手算基准：单位正方形
   $\ell_P\equiv4$、直角三角形 $\equiv6$、梯形 $(0,0),(2,0),(1,1),(0,1)$
   为 $(54-24y)/13$，均已冻结为回归测试。注意 $\ell_P$（连续、边界测度）
   与主程序的 $\ell_A$（离散）是不同对象。
2. **Donaldson 简单凸函数测试**：$(\theta,t)$ 扫描 + golden 细化找
   $M_\ell(\max\{\langle x,u\rangle-t,0\})<0$ 的 witness；`--certify`
   用连分数有理逼近 + `Gmpq` 精确重算做认证；`--svg` 画多边形与折痕线。

实现要点与坑：

- 固定方向时 $M_\ell(t)$ 是**分段四次**（$\ell_P$ 非常数时），且折痕扫过
  与边平行的 breakpoint 时导数跳跃；不做多项式拟合，均匀采样 + breakpoint
  求值 + golden 细化，且 $t$ 括号必须夹在 $[w_{\min},w_{\max}]$
  （区间外 $g$ 仿射，对非真实 $\ell$ 无下界——曾因此漂移出界）。
- 输入解析合并共线连续顶点、拒绝非凸/回溯，顺时针自动反转
  （`examples/Wang_Zhou_a1` 是顺时针）。
- `WILL_FAIL` 会反转 `PASS_REGULAR_EXPRESSION` 的结果，两者不能组合；
  CLI 冒烟测试改用 `/bin/sh -c` 单行脚本同时查退出码与输出。

**经验结论**：仓库全部例子 + 五边形族 + 400 个随机格点凸多边形均未找到
反例（归一化最小值 $>-10^{-9}$）。程序的作用是精确 $\ell_P$ 与验证；搜索/认证
管线由错误 $\ell$ 的机械测试覆盖。Release 全部 12 个测试通过。

---

## 0. 2026-08-01 性能优化更新

本次完成了原第 13 节优先级 1–4 的优化，并引入一个此前未预料到的
数值核陷阱。要点如下，下文各节如与此处冲突以本节为准。

### 已实现的优化

1. **oracle 批量插入与坐标缓存**（原优先级 1、2，此前已部分完成）：
   `compute_oracle` 使用 CGAL range insert（内部 spatial sort），
   `RegularTriangulationOracle::Cache` 预存两种内核类型的
   $x,y,\|a_i\|^2$ 及 scaling 用的 `coordinate_scale`。
2. **双内核 oracle**（原优先级 3）：数值主循环使用
   `Exact_predicates_inexact_constructions_kernel`（EPICK，纯 `double`
   谓词）；精确 endgame、认证和精确绘图路径保持
   `Exact_predicates_exact_constructions_kernel`（EPECK）。
   Wang--Zhou $a=5$、$k=9$、200 次迭代从 11.6 秒降至 1.1 秒
   （约 10 倍；HANDOFF 原记录的 27.9 秒是批量插入之前的数据）。
3. **主循环不保存三角剖分面**：`minimize*` 新增 `keep_faces` 参数，
   迭代中不传 faces，只有绘图路径保留，active set 内存大幅下降。
4. **快照延后**：`SolverResult` 的 `active_vectors`/`sigma`/`coefficients`
   只在 return 前赋值，不再每轮深拷贝整个 active set。
5. **exact Gram 快路径**：Gram 条目先以 `__int128` 累加（带 $2^{120}$
   数量级预检，超出则回退 Gmpz 循环），`gmpz_from_wide` 改用
   `mpz_import`（去掉十进制字符串往返）。$m=100$、$n=13609$ 的 Gram
   构造从 12.1 秒降至 0.026 秒（约 460 倍）。大 active set 的认证瓶颈
   现在只剩 CGAL 精确 QP 求解本身（$m=241$ 时求解超过 3 分钟，
   固有开销，`--exact-max-active` 默认 128 的保护仍然必要）。

### EPICK 假收敛陷阱与精确 endgame

EPICK 在 secondary fan 墙面附近可能返回非最小化三角剖分，使计算
gap 偏小、过早判收敛。对 Wang--Zhou 等一般点集影响不大（六点例子
照常认证通过），但**强对称输入是极端情形**：unit square $k=8$ 的
$\sigma_A$ 是均匀向量，所有 GKZ 顶点对其内积相同，认证要求
$\sigma_A$ 精确落在 active hull 内；EPICK 收集的噪声顶点使认证失败。

当前的解决方案（`ShortestGkzSolver::solve`）：

- 数值 gap 首次满足停止判据后不直接收敛，切换到 EPECK oracle 继续
  迭代（**精确 endgame**），精确 gap 也满足判据才判 `converged=true`；
- 认证失败且 exact oracle 返回更小内积顶点（witness）时，把 witness
  加入 active set 继续迭代；
- **witness 加入时禁止 prune**：数值 QP 会给 witness 分配近零系数
  （当前候选点数值上已最优），若按 1e-15 阈值剪枝，精确信息永远
  无法积累，会空转到迭代上限（实测 $k=8$ 空转 222 轮仍失败）。

新增 `tests/oracle_kernel_consistency_test.cpp`：对随机 heights 比较
EPICK 与 EPECK oracle 的有理支撑值 $H(h)$ 精确相等（允许三角剖分
不同）。Release、Debug、UBSan 三套构建 8 个测试全部通过。

### 基准速查（Apple Silicon，Release）

| 输入 | 优化前 | 优化后 |
| --- | ---: | ---: |
| Wang--Zhou $a=5$, $k=9$, 200 轮, `--no-exact` | 11.6 s | 1.1 s |
| Wang--Zhou $a=5$, $k=16$, 5 轮, `--no-exact` | 31.2 s（旧 HANDOFF 记录） | 0.27 s |
| unit square $k=8$ 收敛+认证 | 0.15 s（43 轮） | 0.22 s（80 轮，含 endgame） |
| exact Gram $m=100$, $n=13609$ | 12.1 s | 0.026 s |

注意：本仓库的 git 元数据此前已损坏并被删除，目前不是 git 仓库；
`/tmp/sgkz_backup` 留有优化前源码备份（临时目录，重启即失）。

---

## 1. 项目位置与版本控制

当前项目根目录：

```text
/Users/yaoy/Documents/sGKZ
```

这是一个 Git 仓库：

- 当前分支：`main`
- 远程仓库：`https://github.com/yeeyoe/sGKZ.git`
- 跟踪分支：`origin/main`
- 本次交接前工作区干净
- 当前只有一个提交：`d658176 初始化仓库`

后续工作应以这个目录为准，不再使用旧的 `latex/code` 或 `latex/sGKZ`
目录。

## 2. 研究背景与记号偏好

论文的长期目标是研究极小元序列 $\eta_k^*$ 的 $C^0$ 预紧性，包括一致
有界性和等度连续性。论文重点考虑

$$
A_k=P_k=kP\cap\mathbb Z^2,
$$

以及 $k\to\infty$ 时 shortest GKZ vector $\sigma_k$ 的变化。

目前程序专注于二维点集 $A$ 的 secondary polytope $\Sigma(A)$，计算

$$
\sigma_A=\operatorname*{argmin}_{g\in\Sigma(A)}\|g\|.
$$

用户的记号偏好：

- 后续论文讨论使用 $\phi_k$，不再使用旧的 $\psi_k$ 记号；
- 但当前程序中 quantum minimizer 及绘图字段仍命名为 `psi` / `psi_k`；
- 用 $\vee_k(f)$ 表示相对于 $P_k$ 取 convex envelope；
- Markdown 数学公式必须使用 `$...$` 或 `$$...$$`；
- note 中不要使用环境里无法编译的平均积分命令；
- 生成普通 Markdown 文件，不使用 lark-markdown skill。

论文文件位于 `paper/`，主要包括：

- `paper/main.tex`
- `paper/easy question.tex`
- `paper/effective_algorithm_for_shortest_GKZ.md`

## 3. 程序结构

主要文件：

| 文件 | 作用 |
| --- | --- |
| `src/gkz.cpp` | 格点、多边形、regular triangulation oracle、数值 QP、精确认证、lower envelope 和绘图数据 |
| `src/main.cpp` | 命令行参数、运行摘要、退出状态 |
| `include/gkz/gkz.hpp` | 公共数据结构和接口 |
| `plot_results.py` | $\sigma_A^\vee$、$\psi_k$、subdivision 的 HTML/SVG/PNG 绘图 |
| `plot_iterations.py` | active size、log gap 与实际停机阈值曲线关于 iteration 的图 |
| `generate_wang_zhou.py` | 生成 Wang--Zhou 多边形输入 |
| `tests/` | C++ 与 Python 回归测试 |
| `README.md` | 中文使用文档和算法说明 |

依赖：C++20、CMake、Ninja、CGAL、Eigen、GMP；静态 Python 图片另需
Matplotlib。

## 4. GKZ 归一化

若 $Q=\operatorname{conv}(A)$，二维三角剖分 $T$ 的归一化 GKZ 向量为

$$
g_T(a_i)
=\frac{\sum_{\tau\in T,\,a_i\in\tau}\operatorname{area2}(\tau)}
{3\operatorname{area2}(Q)}.
$$

因此

$$
\sum_i g_T(a_i)=1.
$$

程序同时保存整数面积分子 `area_numerators`。数值阶段使用归一化后的
`double`，重复判断和精确认证使用整数面积分子。

## 5. 主算法

当前方法是 fully corrective conditional-gradient / active-set 算法。

### 5.1 Active set QP

设 active GKZ 向量为 $v_1,\ldots,v_m$。程序求解

$$
\min_{\lambda_i\geq0,\,\sum_i\lambda_i=1}
\frac12\lambda^TG\lambda,
\qquad
G_{ij}=\langle v_i,v_j\rangle,
$$

并令

$$
x_t=\sum_i\lambda_i v_i.
$$

数值 KKT 系统由 Eigen `CompleteOrthogonalDecomposition` 求解，可以处理
active GKZ 向量的仿射冗余。

### 5.2 Regular triangulation oracle 与 $H(x_t)$

程序把 $x_t=(h_i)$ 当作点集 $A$ 上的提升高度，调用 CGAL
`Regular_triangulation_2`。对 $a_i=(p_i,q_i)$ 插入 weighted point

$$
w_i=p_i^2+q_i^2-h_i.
$$

CGAL 的标准提升高度为 $p_i^2+q_i^2-w_i=h_i$，所以其 finite faces 给出
提升点下凸包诱导的 regular triangulation。

oracle 返回一个 GKZ 向量 $w_t$，程序计算

$$
H(x_t)=\langle x_t,w_t\rangle,
$$

以及 Frank--Wolfe gap

$$
\delta_t=\|x_t\|^2-H(x_t).
$$

停止条件为

$$
\delta_t\leq
\texttt{absolute-tolerance}
+\texttt{tolerance}\cdot\|x_t\|^2.
$$

若 gap 仍大，则把 $w_t$ 加入 active set 并完全校正所有系数。

### 5.3 CGAL 退化情形

若 lower hull 对应粗 subdivision，CGAL 返回一个 refinement triangulation。
已检查本机 CGAL 源码：退化 power-circle 判断使用全局 symbolic
perturbation，而不是分别对每个粗胞腔独立、任意地加入对角线。该 refinement
可视为由全局一致的无穷小扰动诱导的正则三角剖分。

仍建议后续对这个事实补充一个明确的程序回归测试或文档引用，尤其要验证：

- oracle 输出在复杂退化 subdivision 上确实对应 $\Sigma(A)$ 的顶点；
- symbolic perturbation 的结果与精确线性优化值一致；
- hidden weighted points 在所采用的 secondary-polytope 定义下如何解释。

### 5.4 精确认证

数值停止后，若 active set 不超过 `--exact-max-active`，程序使用
`CGAL::Gmpz/Gmpq`：

1. 用整数面积分子构造精确半正定 QP；
2. 精确求出凸组合系数与 $\sigma_A$；
3. 清除 $\sigma_A$ 的公分母，得到等价整数 heights；
4. 再调用 exact regular triangulation oracle；
5. 验证

$$
H(\sigma_A)=\|\sigma_A\|^2.
$$

只有出现

```text
converged=true
exact_certified=true
```

时，才可把结果视为 shortest GKZ 的严格证明。

退出状态：

- `0`：数值停止；是否精确认证还需查看 `exact_certified`；
- `2`：达到最大迭代数但尚未数值收敛，仍会写出当前近似；
- `1`：参数或运行时错误。

## 6. $\sigma^\vee$ 与 subdivision

主迭代中不显式计算整张 $x_t^\vee$。oracle 只返回最小化内积的 GKZ
向量，用于计算 $H(x_t)$。

停机并请求绘图数据后，程序再次以最终 $\sigma$ 为 height 调用 oracle，取得
lower regular triangulation。对每个三角形作重心坐标仿射插值，从而得到
$\sigma^\vee|_A$。

若有 exact certificate，使用 `Gmpq` 做精确插值；否则使用浮点插值。

当前 subdivision 恢复方式不是只比较相邻三角形，而是：

1. 计算每个三角形的提升平面；
2. 在全部三角形中全局寻找相同平面；
3. 合并同一平面的顶点；
4. 对这些顶点取二维凸包作为 cell 边界。

数值路径的共面容差固定为 `1e-8`；精确路径直接比较有理数平面。全局分组
在理论上依赖同一支撑平面的接触集为凸集。更稳健的实现可以改为沿三角形
邻接关系 flood-fill 共面面，再显式提取边界。

绘图中的 lifted break lines 不是所有 refinement triangle 的边，而是恢复出的
subdivision cell 边界。其 $z$ 坐标取同一张曲面的 $\sigma^\vee$ 或 $\psi_k$
值。

## 7. $\psi_k$ 计算与绘图

polygon 模式中程序计算

$$
\psi_k
=\operatorname{area2}(P)k^3\sigma_k^\vee-2k.
$$

默认三维图：

- 曲面为浅蓝色；
- subdivision 的 lifted break lines 为深红色；
- `--color` 才用函数值驱动曲面颜色；
- `--z-scale` 只改变三维场景的视觉纵向比例，不修改函数值；
- polygon 模式默认只生成 $\psi_k$ 三维图；
- `--sigma-vee` 可额外生成 $\sigma_k^\vee$ 三维图。

## 8. 尚未实现的绘图需求

用户要求增加一种新的二维 subdivision 图，用颜色显示 $\psi_k$ 大小，同时
保留清晰的 cell 边界。

推荐实现：

- 利用现有 refinement triangles 和顶点 $\psi_k$ 值，在每个 subdivision
  cell 内显示连续的分片仿射色场；
- 使用以 $0$ 为中心的发散色图，例如蓝色表示负值、浅色表示接近 $0$、
  红色表示正值；
- cell 边界使用深灰或黑色，并加一层细白色 halo，以免边界与红蓝色填充
  混淆；
- 只画真正的 subdivision cell 边界，不画 cell 内的 CGAL refinement
  对角线；
- 添加标注为 $\psi_k$ 的 colorbar；
- 建议独立输出为 `PREFIX_psi_subdivision.svg`，并设置独立开关。

不要简单地给每个 cell 只填一种颜色，除非明确规定使用 cell 重心处的值；
因为 $\psi_k$ 在每个 cell 上一般是仿射函数而不是常数。

## 9. 已完成的精度检查

此前 `surface.csv` 使用 17 位有效数字，而 `subdivision.csv` 只使用默认 6 位，
导致 Python 无法匹配同一个坐标。现已统一使用
`std::setprecision(17)`。

对 Wang--Zhou $k=9$ 的一次检查中：

- `surface.csv` 有 9001 个点；
- subdivision 的所有顶点都能在 surface 中找到；
- 所有数值均 finite；
- $\psi_k$ 公式重算的最大十进制差约为 $1.6\times10^{-15}$。

因此旧的

```text
ValueError: Subdivision vertex (...) is absent from surface
```

已经由输出精度统一解决。

需要严格区分：未收敛时 `sigma` 与 `sigma_vee` 的显著差异不是舍入误差。
例如 Wang--Zhou $k=9$、200 次迭代的一次运行结果为

```text
converged=false
active_size=186
gap=3.7063684984380815e-05
exact_certified=false
```

此时 `sigma` 只是 active hull 上的近似点，不保证已经凸，因此可能有
$\sigma^\vee<\sigma$。

## 10. 剩余数值风险

1. 数值 subdivision 使用固定 `1e-8` 共面容差；当平面系数很小时，这可能
   错误合并不同平面。
2. 数值 KKT 残差接受阈值固定为 `1e-8`，但 `correction_tolerance` 默认是
   `1e-14`，两者尺度不协调。
3. Apple Silicon 上 `long double` 与 `double` 都是 53 位二进制有效位，
   `long double` 累加没有额外精度。
4. ~~数值 oracle 对已经舍入为 binary64 的 height 做精确几何~~（2026-08-01
   起数值循环直接使用 `double` 谓词核，墙面附近可能翻转三角剖分；该风险
   已在 unit square $k=8$ 上实际发生，由精确 endgame + witness 反馈机制
   兜底，见第 0 节）。
5. 数值计算 $\psi_k$ 包含两个 $O(k)$ 项的相减；当前规模误差很小，但大 $k$
   时可能出现消去误差。exact-certified 路径使用有理数计算。
6. `plot_results.py` 使用 `round(x,12)` 匹配坐标；当前规模安全，极大 $k$
   时可能碰撞。长期应输出和使用格点索引。

## 11. 测试状态

此前 Release、Debug、UBSan 三套构建均通过全部 7 个自动测试：

```text
100% tests passed, 0 tests failed
```

重新验证命令：

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

## 12. 当前性能瓶颈

> 2026-08-01 更新：以下为双内核拆分**之前**的记录；当前数值见第 0 节
> 基准速查表。主瓶颈仍是每轮从零构造 CGAL
> `Regular_triangulation_2`，但数值循环已改用 `double` 谓词核。

主瓶颈是每轮从头构造 CGAL `Regular_triangulation_2`。当前运行基本单核，
外层迭代又存在

$$
x_t\longrightarrow w_t\longrightarrow
\text{active-set update}\longrightarrow x_{t+1}
$$

的数据依赖，不能直接并行连续迭代。

此前 Release 实测：

| 输入 | 格点数 | 迭代数 | 时间 |
| --- | ---: | ---: | ---: |
| Wang--Zhou $k=9$ | 9001 | 20 | 8.6 秒 |
| Wang--Zhou $k=9$ | 9001 | 60 | 15.6 秒 |
| Wang--Zhou $k=9$ | 9001 | 200 | 27.9 秒 |
| Wang--Zhou $k=9$，含绘图数据 | 9001 | 200 | 28.1 秒 |
| Wang--Zhou $k=16$ | 28209 | 5 | 31.2 秒 |

`user time` 与 `real time` 基本相等，说明只使用一个核心。$k=9$ 时绘图后
处理只增加约 0.2 秒，并非当前瓶颈。

## 13. 性能优化优先级

> 2026-08-01 更新：优先级 1–4 已完成，见第 0 节。以下为原始记录与
> 剩余方向。

建议先优化单个 oracle，再考虑并行化。

1. ~~**CGAL 批量插入。**~~（已完成）当前逐点调用 `triangulation.insert`，没有 location
   hint。CGAL 的 range insert 支持 `(Weighted_point, info)`，会 spatial
   sort 并在连续插入间传递 face hint。这是最直接的高收益修改。
2. ~~**缓存固定几何量。**~~（已完成）在 oracle 构造时预计算 $x,y,x^2+y^2$，避免每轮
   重复把整数通过字符串转换为 CGAL exact number。
3. ~~**拆分 numerical/exact kernel。**~~（已完成，但务必阅读第 0 节的
   EPICK 假收敛陷阱与 endgame 设计）数值主循环可研究使用
   `Exact_predicates_inexact_constructions_kernel`，精确认证继续使用当前
   exact kernel。必须增加退化数据和随机数据的一致性测试。
4. ~~**迭代 oracle 不保存 faces。**~~（已完成）active set 只需要 GKZ 数值和整数面积
   分子；完整 triangulation 只在最终绘图时需要。
5. **多 refinement 并行。** 对多个全局微扰 height 并行调用 oracle，一轮
   加入多个不同 GKZ 顶点。若大量重复则会浪费时间；MacBook Air 建议先用
   3--4 个 worker，而不是同时占满 10 核。
6. active set 很大后，再并行 Gram 新列和 exact Gram 构造。

当前 16GB 内存对 $k\leq16$、active set 为几百的单进程数值计算通常足够。
真正可能造成内存和时间爆炸的是大 active set 上的精确稠密 QP **求解**
（不是 Gram 构造，后者已优化到可忽略），而不是当前数值 oracle。

## 14. 推荐的下一轮工作顺序

> 2026-08-01 更新：原第 1–4 项已完成，见第 0 节。剩余顺序如下。

1. 实现二维 $\psi_k$ 彩色 subdivision 图及测试（第 8 节）。
2. 改进数值共面容差、KKT residual 和坐标索引（第 10 节）。
3. 若继续追求 oracle 速度：评估 flip-based 增量 regular
   triangulation（从上一轮三角剖分出发做加权边翻转），比从零重建
   更接近「变化量」复杂度；CGAL 不支持 weight 更新，需自行实现。
4. 最后再评估多 refinement 并行算法。

任何算法优化都应继续保留 exact certificate 作为最终正确性标准。
