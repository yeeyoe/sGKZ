# K-stability：精确计算 $\ell_P$ 并检测相对 K-不稳定性

本目录实现 `paper/K-stability.tex` 的两个目标，与主程序 `shortest_gkz`
完全隔离：不链接 `gkz_core`，只依赖 CGAL（`Gmpq` 精确有理数）。

设 $P\subset\mathbb R^2$ 是格点多边形，边界上带格点测度
$\mathrm{d}\sigma$（在原语法向量 $v_i$ 的边 $F_i$ 上
$\mathrm{d}\sigma|_{F_i}=\mathrm{d}s/|v_i|$；等价地，边 $p\to q$ 的
$\mathrm{d}\sigma$ 长度是格点长度 $\gcd(|q_x-p_x|,|q_y-p_y|)$），即顶点差向量的分量的 $\gcd$。

## 目标 1：精确计算 $\ell_P$

$\ell_P$ 是唯一满足下式的仿射函数：对所有仿射函数 $h$，

$$
\int_{\partial P}h\,\mathrm{d}\sigma=\int_P h\,\ell_P\,\mathrm{d}x.
$$

取 $h=1,x,y$ 得 $3\times3$ 线性系统

$$
\begin{pmatrix}
V & \int_P x & \int_P y\\
\int_P x & \int_P x^2 & \int_P xy\\
\int_P y & \int_P xy & \int_P y^2
\end{pmatrix}
\begin{pmatrix}a\\b\\c\end{pmatrix}
=
\begin{pmatrix}
|\partial P|_{\mathrm{d}\sigma}\\
\int_{\partial P}x\,\mathrm{d}\sigma\\
\int_{\partial P}y\,\mathrm{d}\sigma
\end{pmatrix},
\qquad \ell_P=a+bx+cy.
$$

格点多边形的所有矩都是有理数：面积分用从原点出发的带符号三角形扇公式，
边界积分对仿射 $h$ 有精确公式
$\int_{p q}h\,\mathrm{d}\sigma=\gcd(|\Delta x|,|\Delta y|)\,\frac{h(p)+h(q)}2$。
程序用 `CGAL::Gmpq` 以 Cramer 法则精确求解，输出分数形式的 $\ell_P$。

手算基准（回归测试冻结）：

| 多边形 | $\ell_P$ |
| --- | --- |
| 单位正方形 | $4$ |
| 三角形 $(0,0),(1,0),(0,1)$ | $6$ |
| 梯形 $(0,0),(2,0),(1,1),(0,1)$ | $\dfrac{54-24y}{13}$ |

## 目标 2：Donaldson 简单凸函数测试

Donaldson 事实（$n=2$）：$P$ 相对 K-不稳定当且仅当存在简单凸函数

$$
g=\max\{ax+by+c,\,0\}
$$

使相对 DF 不变量

$$
M_\ell(g)=\int_{\partial P}g\,\mathrm{d}\sigma-\int_P g\,\ell_P\,\mathrm{d}x<0.
$$

$M_\ell$ 只依赖折痕线 $\{ax+by+c=0\}$ 及其法向一侧。程序以单位法向
$u=(\cos\theta,\sin\theta)$ 和偏移 $t$ 参数化
$g=\max\{\langle x,u\rangle-t,0\}$：

1. **方向扫描**：$\theta\in[0,2\pi)$ 均匀网格（默认 720），另加所有边法向
   与顶点对法向作为种子方向。
2. **偏移采样**：固定 $u$ 时，对每个 $t$ 用半平面裁剪
   $P^+=P\cap\{s\ge0\}$（Sutherland–Hodgman）在 double 下计算
   $M_\ell=K-J$，其中 $K$ 逐条原始边对子段积分，$J$ 由 $P^+$ 的二阶矩
   组合。$t$ 在投影范围 $[w_{\min},w_{\max}]$ 上均匀采样（默认 512）并
   显式评估所有顶点投影 breakpoint；区间外 $M\equiv0$。
3. **细化**：最优小区间内 golden-section 细化，再对前 5 个候选做
   3 轮交替 $(\theta,t)$ 细化（$t$ 括号始终夹在 $[w_{\min},w_{\max}]$）。
4. **判定**：报告归一化值
   $M_\ell/(|\partial P|_{\mathrm{d}\sigma}\cdot\sup_P|\ell_P|\cdot\mathrm{diam} P)$；
   小于 $-10^{-6}$ 判 `unstable`，否则 `no_counterexample_found`。
5. **可选认证**（`--certify`）：对数值 witness 的 $(\cos\theta^*,
   \sin\theta^*,-t^*)$ 做连分数有理逼近（分母上限 $10,10^2,\dots$ 递增），
   逐档用 `Gmpq` 精确重算 $M_\ell$；有理折痕与格点多边形的交点参数是
   有理数，因此精确评估没有舍入。首个负值即认证 witness。

严格性层级：找到 witness（尤其 `certified=true`）是相对 K-不稳定的严格
证明；`no_counterexample_found` 只是数值证据，不是半稳定性的证明。

## 用法

```bash
k_stability --polygon FILE [options]
```

| 参数 | 默认值 | 含义 |
| --- | ---: | --- |
| `--polygon FILE` | 必选 | 格点多边形顶点文件（每行 `x y`，`#` 注释，逗号视为空白；顺时针输入自动反转）。 |
| `--theta-steps N` | `720` | 折痕方向数。 |
| `--t-steps N` | `512` | 每方向偏移采样数。 |
| `--no-refine` | 关 | 跳过局部细化。 |
| `--certify` | 关 | 用有理数精确认证 witness。 |
| `--certify-max-denom N` | `1048576` | 认证时有理逼近的分母上限。 |
| `--svg FILE` | 不写出 | 写出多边形 + witness 折痕线（含 $P\cap\{s\ge0\}$ 阴影）的 SVG。 |
| `--check-line "a b c"` | 无 | 直接精确/数值评估 $M_\ell(\max\{ax+by+c,0\})$；`p/q` 或整数走精确路径。 |
| `--verbose` | 关 | 逐方向输出最小值到标准错误。 |

例子：

```bash
./build/K-stability/k_stability --polygon examples/unit_square.polygon
./build/K-stability/k_stability --polygon examples/my --svg /tmp/my.svg
./build/K-stability/k_stability --polygon examples/unit_square.polygon \
  --check-line "1 0 -1/2"        # 精确输出 M_l=1/4
```

## 输出字段

标准输出为逐行 `key=value`：

| 字段 | 含义 |
| --- | --- |
| `vertices` / `twice_area` | 顶点数 / 二倍面积。 |
| `boundary_length_dsigma` | $|\partial P|_{\mathrm{d}\sigma}$，精确有理数。 |
| `ell_P(x,y)` | $\ell_P$ 的精确分数系数。 |
| `ell_P_constant` | $\ell_P$ 是否为常数。 |
| `search_evaluations` | $M_\ell$ 求值次数。 |
| `sweep_min_M_l` / `sweep_min_M_l_normalized` | 搜索得到的（归一化）最小值。 |
| `witness_theta` / `witness_t` / `witness_g` | 最优候选的参数与函数形式。 |
| `relative_K_status` | `unstable` 或 `no_counterexample_found`。 |
| `certified` 等 | `--certify` 时的认证 witness 与精确 `certified_M_l<0`。 |

退出码：`0` 找到（且若要求则已认证）witness；`2` 未发现反例或认证失败；
`1` 参数/输入/运行错误。

## 面积优先搜索

`k_stability_search` 在固定顶点数 `d` 下随机生成严格凸整点多边形。前
`d-1` 条边的方向和步长按极角递增生成，闭合边自动计算且不受 `N,M`
限制；整体步长公因子会被约去。候选按精确 `twice_area` 与 Donaldson
probe 分数进入两个 frontier，依次使用 `probe`、`confirm`、`final` profile，
最终只有 `df_simple_exact(...) < 0` 才计为 `verified_unstable`。

搜索状态、候选的方向/步长/顶点/facet normals/`ell_P`、各 detector profile
和认证 witness 均保存在 SQLite 中；同一 profile 的 `unverified` 候选在重启
后会跳过。搜索不内置文献多边形或种子。

```bash
./build/K-stability/k_stability_search --d 6 --N 4 --M 4 \
  --time-limit 3600 --database search.sqlite
```

### 输入参数

必需参数：

| 参数 | 含义 |
| --- | --- |
| `--d D` | 固定顶点数，`D >= 3`。实际搜索通常从 `D >= 6` 开始。 |
| `--N N` | 初始受约束方向的坐标上界，`|p_x|,|p_y| <= N`。 |
| `--M M` | 初始前 `d-1` 条边的步长上界，`1 <= k_i <= M`。闭合边不受此限制。 |
| `--time-limit SEC` | 全局搜索时间预算，单位为秒，必须为正数。 |

可选参数：

| 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `--database FILE` | `k_stability_search.sqlite` | SQLite 状态及候选记录文件。 |
| `--output-dir DIR` | `.` | 结果报告目录；写入 `k_stability_search_result.txt`。 |
| `--shell-seconds SEC` | `60` | 每个 `(N,M)` 自动扩展 shell 的时间片。 |
| `--beam-width K` | `48` | 每批随机/beam 候选生成数量。 |
| `--seed S` | `1` | 随机搜索种子；相同数据库会优先恢复已保存的随机状态。 |
| `--stop-on-first` | 关闭 | 第一个精确认证候选出现后立即结束。默认持续搜索到时间预算耗尽。 |
| `--certify-max-denom Q` | `1048576` | 数值 witness 有理化时的最大分母。该值改变会产生新的 detector profile。 |
| `--verbose` | 关闭 | 输出额外进度信息。 |

自动 shell 按 `(N,M)`, `(N,2M)`, `(2N,2M)`, `(2N,4M), ...` 扩展，直到达到
`--time-limit`。重新运行同一个数据库会恢复 shell、随机状态和已有候选；同一
detector profile 下的 `unverified` 候选不会重复检测。

### 输出说明

程序向标准输出写逐行 `key=value`：

| 字段 | 含义 |
| --- | --- |
| `database` | 实际使用的 SQLite 路径。 |
| `report_file` | 实际写出的结果报告路径。 |
| `generated` / `rejected` | 本次运行接受/拒绝的候选生成数量。 |
| `probes` / `confirms` / `finals` | 本次运行各检测 profile 的调用数量。 |
| `verified` / `unverified` | 数据库中当前已认证/未认证候选数量。 |
| `skipped` | 因已有相同 detector profile 记录而跳过的检测数量。 |
| `have_verified` | 是否已经找到精确认证的不稳定候选。 |
| `best_twice_area` | 当前数据库中已认证候选的最小二倍面积；实际面积为该值的一半。 |
| `first_verified_key` | 首个认证候选的 canonical key。 |
| `best_verified_key` | 当前最小面积认证候选的 canonical key。 |

只有 `verified_unstable` 候选才会记录数值 witness、精确有理 witness 及精确
`M_l<0`；`unverified` 不能解释为半稳定或稳定。

### 输出文件位置

SQLite 持久化输出在 `--database FILE` 指定的文件中，文字结果报告在
`--output-dir DIR/k_stability_search_result.txt`：

- 相对路径相对于启动命令时的当前工作目录解析；
- 文件不存在时自动创建父目录和数据库；
- `candidates` 表保存方向序列、步长序列、顶点、facet normals、`ell_P`、面积和状态；
- `attempts` 表保存每个候选/profile 的检测结果，只有认证记录含 witness 字段；
- `state` 表保存 shell、随机游标等断点状态。
- 找到认证候选时，报告文件保存其完整几何信息、面积、边界测度、`ell_P`、
  profile 和精确 witness；没有找到时报告内容为 `没找到`。

SQLite 运行期间可能同时出现同名的 `-wal` 和 `-shm` 临时文件；它们与主
数据库放在同一目录。除上述 SQLite 和结果报告外，程序不会生成多边形文本
或 SVG 文件。

返回码 `0` 表示数据库中已有精确认证候选，`2` 表示时间预算内尚未找到，
`1` 表示参数、输入或运行错误。

## 测试

```bash
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

`k_stability_tests` 覆盖：三个手算 $\ell_P$ 基准、扇形矩公式、边界测度、
$\ell_P$ 定义性质（任意仿射 $h$ 精确 $M_\ell(h)=0$）、正方形手算回归
$M_\ell(\max\{x-\frac12,0\})=\frac14$、double/exact 一致性、退化折痕、
解析器校验、连分数逼近、半稳定集成扫描，以及（用错误 $\ell$ 的）搜索 +
认证 + SVG 机械测试。
