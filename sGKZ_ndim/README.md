# `sgkz`：动态维 Shortest GKZ 求解器

本目录包含一个独立的动态维 Shortest GKZ vector 求解器。程序从输入文件每一行的列数自动推断维数 `d`，支持任意 `d >= 2`。它使用 regular triangulation、active-set / fully-corrective Frank--Wolfe 方法寻找 secondary polytope 中欧氏模长最小的 GKZ 向量，并在条件允许时使用精确有理数二次规划进行全局认证。

## 目录结构

* `src/main.cpp`：命令行程序和参数解析。
* `src/gkz_ndim.cpp`：点配置、格点枚举、regular triangulation、Shortest GKZ 和精确认证实现。
* `include/gkz_ndim/gkz_ndim.hpp`：公开 C++ 接口和求解器选项。
* `tests/gkz_ndim_tests.cpp`：核心逻辑测试。
* `examples/`：点集、单纯形和多胞体示例输入。

## 依赖

需要：

* CMake 3.20 或更高版本；
* 支持 C++20 的编译器；
* CGAL（包含 `CGAL::Gmpq`）；
* Eigen3；
* Ninja（可选，下面的命令使用 Ninja）。

在 macOS 上可以用 Homebrew 安装常见依赖：

```bash
brew install cmake ninja cgal eigen
```

## 构建

在仓库根目录 `/Users/mac/Documents/sGKZ` 执行：

```bash
cmake -S sGKZ_ndim -B sGKZ_ndim/build \
  -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build sGKZ_ndim/build --target sgkz --parallel
```

可执行文件为 `sGKZ_ndim/build/sgkz`。同时构建并运行测试：

```bash
cmake --build sGKZ_ndim/build --parallel
ctest --test-dir sGKZ_ndim/build --output-on-failure
```

如果系统没有 Ninja，可省略 `-G Ninja`，或指定其他 CMake generator。

## 输入文件

输入文件是空白分隔的整数矩阵，每一行是一个点。空行会被忽略；所有非空行必须有相同列数，且点必须互不相同、张成满维凸包。

例如，`examples/simplex3.points` 表示三维标准单纯形的四个顶点：

```text
0 0 0
1 0 0
0 1 0
0 0 1
```

维数由第一行的列数推断，不需要额外的 `--dimension` 参数。 `--points` 模式直接把文件中的所有点作为点配置。

`--polytope` 模式的文件同样每行一个整数顶点，但程序会先构造 `kP`，枚举其中的整数点 `kP ∩ Z^d`，再把这些格点作为点配置。例如，`unit_cube3.polytope` 是单位立方体的 8 个顶点；使用 `--k 2` 时，求解对象是 `2P` 中的全部整数格点。

## 基本用法

```bash
# 直接读取点配置
sGKZ_ndim/build/sgkz \
  --points sGKZ_ndim/examples/simplex3.points

# 从多胞体顶点枚举 kP 中的格点
sGKZ_ndim/build/sgkz \
  --polytope sGKZ_ndim/examples/unit_cube3.polytope \
  --k 2

# 计算 affine ell_A，并把逐点结果写成 CSV
sGKZ_ndim/build/sgkz \
  --points sGKZ_ndim/examples/simplex3.points \
  --ell \
  --output sGKZ_ndim/build/simplex3.csv
```

查看完整帮助：

```bash
sGKZ_ndim/build/sgkz --help
```

`--points` 与 `--polytope` 必须且只能指定一个；`--polytope` 必须同时指定正整数 `--k`。

## 命令行参数

### 输入和输出

| 参数 | 默认值 | 说明 |
|---|---:|---|
| `--points FILE` | 无 | 读取点配置文件。 |
| `--polytope FILE` | 无 | 读取整数多胞体顶点，并枚举 `kP ∩ Z^d`。 |
| `--k INTEGER` | `0` | 多胞体缩放倍数；仅和 `--polytope` 一起使用，必须大于 0。 |
| `--output FILE` | 不写文件 | 写出逐点动态列 CSV。父目录必须已经存在。 |
| `--ell` | 关闭 | 计算 affine 函数 `ell_A`，同时加入标准输出和 CSV。 |
| `--max-points N` | `1000000` | `--polytope` 模式最多允许枚举的格点数。`0` 表示不限制。 |
| `--verbose` | 关闭 | 将每轮 active-set 的诊断信息写到 stderr。 |
| `--help` | - | 打印帮助并退出。 |

### 数值求解选项

| 参数 | 默认值 | 说明 |
|---|---:|---|
| `--tolerance VALUE` | `1e-11` | 相对 Frank--Wolfe gap 停止阈值。 |
| `--absolute-tolerance VALUE` | `1e-14` | 绝对 gap 停止阈值；实际阈值为 `max(absolute, tolerance * max(1, norm_squared))`。 |
| `--correction-tolerance VALUE` | `1e-14` | active-set 二次规划修正和进入判定的容差。 |
| `--prune-tolerance VALUE` | `1e-15` | 删除过小 active 系数的阈值。 |
| `--max-iterations N` | `500` | active-set 扩展和精确认证最多使用的迭代次数。 |
| `--max-correction-steps N` | `10000` | 每次 active-set QP 最多修正步数。 |
| `--exact-max-active N` | `128` | active 向量数量不超过该值时尝试精确认证；`0` 表示不限制。 |
| `--no-exact` | 关闭 | 跳过精确 QP 和全局认证，仅返回数值结果。 |

所有数值阈值和步数参数都不能为负数。

## 标准输出

程序每行输出一个 `key=value` 字段，例如：

```text
dimension=3
points=4
volume=1/6
converged=true
iterations=0
active_size=1
norm_squared=...
gap=0
l2_error_bound=0
exact_certified=true
exact_norm_squared=...
```

字段含义：

* `dimension`：点的维数 `d`。
* `points`：点配置中的点数；`--polytope` 模式下是枚举得到的格点数。
* `volume`：点配置凸包的几何体积，即行列式体积除以 `d!`，以精确有理数打印。
* `converged`：浮点 active-set 求解是否达到停止阈值。
* `iterations`：数值求解执行的外层迭代次数。
* `active_size`：最终 active GKZ 向量数量。
* `norm_squared`：当前 Shortest GKZ 候选的平方范数。
* `gap`：Frank--Wolfe 对偶 gap；越接近 0，候选越接近最优。
* `l2_error_bound`：由 gap 给出的 `L2` 误差上界 `sqrt(2 * gap)`。
* `exact_certified`：是否通过精确有理数 QP 和 exact oracle 验证了全局最优性。
* `exact_norm_squared`：仅在 `exact_certified=true` 时输出的精确平方范数。
* `ell_A`：仅在使用 `--ell` 时输出，内容为 `d+1` 个 affine 系数，顺序为常数项、`x1`、...、`xd`。
* `output`：仅在使用 `--output` 时输出实际写入的 CSV 路径。

数值求解未收敛且没有精确认证时，程序退出码为 `2`；参数或输入错误退出码为 `1`；成功求解或成功认证退出码为 `0`。

## CSV 输出

使用 `--output result.csv` 后，文件第一行包含动态列：

```text
x1,x2,...,xd,sigma,sigma_exact
```

使用 `--ell` 时还会追加：

```text
ell_A,ell_A_exact
```

每一行对应一个输入点：

* `x1` 到 `xd`：点坐标；
* `sigma`：数值算法得到的 Shortest GKZ 坐标；
* `sigma_exact`：精确认证得到的有理数坐标。未认证时该列为空；
* `ell_A`：affine `ell_A` 的浮点值；
* `ell_A_exact`：affine `ell_A` 的精确有理数值。

CSV 中使用小数点表示浮点列，精确列使用整数或 `numerator/denominator` 表示。

## 典型运行与日志

下面的命令把标准输出、CSV 和 stderr 日志分别保存：

```bash
mkdir -p sGKZ_ndim/build/results
sGKZ_ndim/build/sgkz \
  --polytope sGKZ_ndim/examples/unit_cube3.polytope \
  --k 2 --ell --verbose \
  --output sGKZ_ndim/build/results/cube3_k2.csv \
  > sGKZ_ndim/build/results/cube3_k2.out \
  2> sGKZ_ndim/build/results/cube3_k2.log
```

`--verbose` 日志中的 `iteration`、`active` 和 `gap` 描述数值 active-set 过程，不会写入 CSV。若格点数超过 `--max-points`，程序会报错；高维或较大的 `k` 可能导致枚举和 exact certification 显著变慢，此时可提高 `--max-points`、使用 `--no-exact`，或设置合适的 `--exact-max-active`。

## C++ 接口

除命令行程序外，也可以链接 `sgkz_core`：

```cpp
#include <gkz_ndim/gkz_ndim.hpp>

auto configuration = sgkz::PointConfiguration::from_points_file("input.points");
sgkz::ShortestGkzSolver solver;
sgkz::SolverResult result = solver.solve(configuration);
```

`PointConfiguration::from_polytope_file(path, k, max_points)` 对应命令行的 `--polytope` 模式；`SolverOptions` 可用于在 C++ 中设置与命令行相同的容差、迭代上限和 exact certification 选项。

## 常见问题

* **`points are not full dimensional`**：输入点没有张成其列数对应的维数，或输入点数不足。
* **`input rows have different dimensions`**：输入文件各行列数不一致。
* **`duplicate point`**：输入中存在重复点。
* **`lattice point count exceeds --max-points`**：多胞体模式枚举的格点数超过上限；确认资源足够后提高上限，或使用 `0` 取消限制。
* **`exact certification skipped by active-set limit`**：active 向量数量超过 `--exact-max-active`；使用 `--exact-max-active 0` 可取消该限制，但计算可能更慢。

