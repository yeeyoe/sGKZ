# 计算 shortest GKZ vector 的有效算法

本文讨论如何在不枚举 secondary polytope 的全部顶点、也不穷尽所有 triangulations 的前提下，计算

$$
\sigma_k
=
\operatorname*{argmin}_{g\in\Sigma(P_k)}\|g\|.
$$

推荐方案是：

> 使用 fully corrective Frank--Wolfe/Wolfe minimum-norm-point algorithm，并以 regular triangulation 作为 linear minimization oracle。

这种方法的基本版本每次迭代只计算一个由当前 height vector 诱导的 regular triangulation；也可以在同一个 exposed subdivision 中生成一小批 regular refinements。两种版本都具有严格的最优性判据和误差证书。

## 1. 记号

令

$$
P_k=P\cap k^{-1}M,
\qquad
N_k=|P_k|,
\qquad
\Sigma_k:=\Sigma(P_k).
$$

对于 $f:P_k\to\mathbb R$，本文用

$$
\vee_k(f)
$$

表示关于 $P_k$ 所作的凸包络。

若 $T$ 是 $(P,P_k)$ 的一个 triangulation，则相应的归一化 GKZ vector 为

$$
g_T(u)
=
\frac{1}{(n+1)V_P}
\sum_{\substack{\Delta\in T\\u\in\operatorname{vert}(\Delta)}}
|\Delta|,
\qquad
u\in P_k.
$$

secondary polytope 定义为

$$
\Sigma_k
=
\operatorname{conv}
\left\{
g_T:
T\text{ 是 }(P,P_k)\text{ 的 triangulation}
\right\}.
$$

欧氏范数采用

$$
\|g\|^2
=
\sum_{u\in P_k}g(u)^2.
$$

## 2. Linear minimization oracle

定义

$$
H_k(f)
:=
\frac{1}{V_P}
\int_P\vee_k(f)\,dx.
$$

secondary polytope 的 support-function 公式是

$$
\boxed{
H_k(f)
=
\min_{g\in\Sigma_k}\langle f,g\rangle.
}
\tag{2.1}
$$

因此，给定 $f\in\mathbb R^{P_k}$，可以在不知道 $\Sigma_k$ 全部顶点的情况下求解右侧的 linear minimization problem：

1. 考虑 lifted points

   $$
   (u,f(u))\in M_{\mathbb R}\times\mathbb R,
   \qquad u\in P_k.
   $$

2. 计算这些点的 lower convex hull。
3. lower faces 给出 regular subdivision $S(f)$。
4. 对 lower faces 作任意 regular triangulation $T_f\prec S(f)$。
5. 计算 $g_{T_f}$。

由 support-function 公式，

$$
g_{T_f}
\in
\operatorname*{argmin}_{g\in\Sigma_k}
\langle f,g\rangle.
\tag{2.2}
$$

若 lower hull 有非单纯形面，可以使用 symbolic perturbation 选出一个 refining triangulation。扰动只用于选择 refinement；objective value 仍使用原来的 $f$。

一次在 $\Sigma_k$ 上的线性优化只需要计算一次 regular triangulation，而不需要枚举其他 triangulations。

## 3. Shortest vector 的变分判据

$\sigma_k$ 是 $\Sigma_k$ 的最小范数点，当且仅当

$$
\langle\sigma_k,g-\sigma_k\rangle
\geq0,
\qquad
\forall g\in\Sigma_k.
\tag{3.1}
$$

等价地，

$$
H_k(\sigma_k)=\|\sigma_k\|^2.
\tag{3.2}
$$

对任意当前可行点 $x\in\Sigma_k$，令 oracle 返回

$$
s(x)
\in
\operatorname*{argmin}_{g\in\Sigma_k}
\langle x,g\rangle.
$$

定义 Frank--Wolfe gap

$$
\delta(x)
:=
\|x\|^2-\langle x,s(x)\rangle
=
\|x\|^2-H_k(x).
\tag{3.3}
$$

因为 $x\in\Sigma_k$，所以 $\delta(x)\geq0$。由 (3.1)，

$$
\boxed{
\delta(x)=0
\quad\Longleftrightarrow\quad
x=\sigma_k.
}
\tag{3.4}
$$

因此，regular-triangulation oracle 产生一个由单个 triangulation 给出的 vertex GKZ vector $g_T$；把它与当前可行点比较，便得到严格的最优性证书。需要强调：oracle 的输出和整个 fully corrective algorithm 的输出不是同一个对象。

## 4. Fully corrective algorithm

维护一个由已经发现的 GKZ vectors 组成的 active set：

$$
\mathcal A_t
=
\{g_1,\ldots,g_m\}.
$$

这里每个

$$
g_i=g_{T_i}
$$

都是某个 triangulation $T_i$ 对应的 vertex GKZ vector。但是算法的当前候选解以及最终输出是

$$
x_t
=
\sum_{i=1}^m\lambda_i g_{T_i},
\qquad
\lambda_i\geq0,
\qquad
\sum_{i=1}^m\lambda_i=1.
\tag{4.0}
$$

因此 $x_t$ 是 $\Sigma_k$ 中的一般点，通常不是顶点，也通常不对应任何单个 triangulation。

这正是 fully corrective step 的作用：oracle 每轮只提供一个顶点，但 active-set quadratic program 在已经发现的顶点的凸包中寻找最小范数点。若 $\sigma_k$ 位于一条边、一个高维 face 的相对内部，算法会用多个 triangulation GKZ vectors 的非平凡凸组合逼近并最终表示 $\sigma_k$。

由 Carathéodory theorem，

$$
\dim\Sigma_k=N_k-n-1
$$

意味着 $\sigma_k$ 总可以表示成至多

$$
N_k-n
$$

个 vertex GKZ vectors 的凸组合；但它不必等于其中任何一个。

### 4.1 初始化

任选一个 pulling triangulation $T_0$，令

$$
\mathcal A_0=\{g_{T_0}\}.
$$

### 4.2 Active-set quadratic program

在当前 active set 的凸包中计算最小范数点：

$$
\min_{\substack{
\lambda_i\geq0\\
\sum_{i=1}^m\lambda_i=1
}}
\frac12
\left\|
\sum_{i=1}^m\lambda_ig_i
\right\|^2.
\tag{4.1}
$$

记其解为

$$
x_t
=
\sum_{i=1}^m\lambda_i g_i.
\tag{4.2}
$$

这里 $x_t$ 才是 fully corrective algorithm 的候选输出，而不是下一步 oracle 返回的 $s_t$。把所有满足 $\lambda_i=0$ 的 active vectors 删除。

### 4.3 Oracle step

以 $x_t:P_k\to\mathbb R$ 为 lifting height，计算 lower hull 和一个 refining regular triangulation $T_t$。令

$$
s_t=g_{T_t}.
$$

计算

$$
\delta_t
=
\|x_t\|^2-\langle x_t,s_t\rangle.
\tag{4.3}
$$

- 若 $\delta_t\leq\varepsilon$，停止并输出凸组合 $x_t$ 以及系数 $(\lambda_i)$；
- 否则，把 $s_t$ 加入 active set，重新求解 (4.1)。

伪代码如下：

    T0 = pulling_triangulation(P_k)
    active = [GKZ(T0)]

    repeat:
        solve λ = argmin 1/2 ||Σ_i λ_i active[i]||²
        subject to λ_i ≥ 0 and Σ_i λ_i = 1

        x = Σ_i λ_i active[i]
        remove active[i] with λ_i = 0

        T = regular_triangulation(P_k, heights=x)
        s = GKZ(T)

        gap = ||x||² - <x,s>

        if gap ≤ tolerance:
            return x and the convex coefficients λ

        active.append(s)

这可以理解为 Wolfe minimum-norm-point algorithm、fully corrective Frank--Wolfe、column generation，或作用于 dual functional 的 bundle method。

### 4.4 非顶点 shortest GKZ vector 的情形

假设 exact shortest vector 满足

$$
\sigma_k
=
\sum_{i=1}^r\lambda_i^*g_{T_i},
\qquad
\lambda_i^*>0,
\qquad
r\geq2.
\tag{4.4}
$$

fully corrective algorithm 的目标正是发现足够多的 $g_{T_i}$，然后由 quadratic program 恢复系数 $\lambda_i^*$。当 active set 已经包含一个能够表示 $\sigma_k$ 的子集时，(4.1) 的解就是 $\sigma_k$；随后 oracle 验证

$$
\|\sigma_k\|^2
=
\min_{g\in\Sigma_k}\langle\sigma_k,g\rangle,
$$

算法终止。

所以算法不假设 $\sigma_k$ 是 secondary polytope 的顶点。相反，它专门允许并处理 $\sigma_k$ 位于 positive-dimensional face 内部的情况。

### 4.5 Batched face oracle：一次加入多个 refining triangulations

令当前 fully corrective candidate 为 $x_t$，由其 lower hull 得到 regular subdivision

$$
S_t:=S(x_t).
$$

由 subdivision-face correspondence，

$$
\begin{aligned}
F_t
&:=
\operatorname*{argmin}_{g\in\Sigma_k}
\langle x_t,g\rangle\\
&=
F(S_t).
\end{aligned}
\tag{4.5}
$$

若

$$
T_t^{(1)},\ldots,T_t^{(r)}
\prec
S_t
$$

是 $S_t$ 的多个 regular refining triangulations，则相应的 GKZ vectors 全部满足

$$
g_{T_t^{(j)}}\in F_t
$$

以及

$$
\langle x_t,g_{T_t^{(j)}}\rangle
=
H_k(x_t),
\qquad
1\leq j\leq r.
\tag{4.6}
$$

因此它们都是当前 linear minimization problem 的最优解，可以在一次 oracle step 中同时加入 active set：

$$
\mathcal A_{t+1}
=
\mathcal A_t
\cup
\left\{
g_{T_t^{(1)}},\ldots,g_{T_t^{(r)}}
\right\}.
\tag{4.7}
$$

然后在扩大的凸包上重新求解 fully corrective quadratic program。

这里“可以产生多个”并非无条件成立。若 $S_t$ 本身已经是一个 regular triangulation，或者它虽有非单纯形 cells、但只有一个 regular triangulating refinement，则这一轮只能得到一个不同的 regular GKZ vertex。只有当 exposed face $F_t$ 至少含有两个不同 vertices 时，才存在至少两个不同的 regular refining triangulations 可供批量加入。

记

$$
C_t=\operatorname{conv}(\mathcal A_t),
\qquad
\widehat C_t
=
\operatorname{conv}
\left(
\mathcal A_t
\cup
\{g_{T_t^{(1)}},\ldots,g_{T_t^{(r)}}\}
\right).
$$

因为 $C_t\subseteq\widehat C_t\subseteq\Sigma_k$，

$$
\min_{g\in\Sigma_k}\|g\|
\leq
\min_{g\in\widehat C_t}\|g\|
\leq
\min_{g\in C_t}\|g\|.
\tag{4.8}
$$

更精确地，对 batch 中任意一个 vertex 令

$$
C_t^{(j)}
=
\operatorname{conv}
\left(C_t\cup\{g_{T_t^{(j)}}\}\right).
$$

由于 $C_t^{(j)}\subseteq\widehat C_t$，有

$$
\min_{g\in\widehat C_t}\|g\|
\leq
\min_{g\in C_t^{(j)}}\|g\|,
\qquad 1\leq j\leq r.
$$

所以加入多个 refinements 不会使当前 primal approximation 变差，而且本轮 full correction 的结果不会劣于只加入 batch 中任意一个 vertex 的结果。它可能借助不同 vertices 的横向分量相消而获得更小的范数，从而减少后续 lower-hull oracle 的调用轮数；这里“可能”而非“必然”，因为新增 vertices 也可能在下一次 quadratic program 中得到零系数。

但是，“扩大 active hull”不等于“立刻包住 $\sigma_k$”。即使把 $F_t$ 的全部 vertices 都加入，所得最大 active hull 也只是

$$
\operatorname{conv}(C_t\cup F_t).
$$

因此，这一轮通过枚举当前 exposed face 能够包住 $\sigma_k$ 的精确条件是

$$
\sigma_k
\in
\operatorname{conv}(C_t\cup F_t).
$$

条件 $\sigma_k\in F_t$ 是一个更强但方便的充分条件，并非必要条件，因为旧 active hull $C_t$ 也可以参与表示 $\sigma_k$。对实际选取的有限 batch，精确条件相应地是 $\sigma_k\in\widehat C_t$。当 $x_t$ 已经接近 $\sigma_k$，或者已经识别出 optimal face 时，这种 batched strategy 尤其有效。

#### 产生多个 regular refinements 的方法

1. **多次 symbolic perturbation。**  
   取若干 generic vectors $r_j$，使用 infinitesimal heights

   $$
   x_t+\varepsilon r_j,
   \qquad
   0<\varepsilon\ll1.
   $$

   选择 $\varepsilon$ 足够小，使所得 regular triangulation 仍细化 $S_t$。不同的 $r_j$ 往往给出不同的 regular refinements。

2. **在 exposed face 内作 zero-cost flips。**  
   从一个 regular refinement 出发，对满足当前 circuit equality 的 circuits 作 bistellar flips，但只保留仍然细化 $S_t$ 且仍为 regular 的结果。因为这些合格的 flips 不改变

   $$
   \langle x_t,g_T\rangle,
   $$

   它们沿着 secondary polytope 的 exposed face $F_t$ 移动。

3. **利用对称群。**  
   若 $x_t$ 是 $G$-invariant，而 $T$ 是一个 minimizing triangulation，则

   $$
   \gamma T,
   \qquad
   \gamma\in G,
   $$

   也都是 minimizing triangulations。可以同时加入其 orbit 中的 GKZ vertices，或直接加入 orbit average。

4. **对 maximal cells 作局部 refinements。**  
   可以在 $S_t$ 的 maximal cells 内产生多个 local triangulations，再将它们 gluing 成 global refinements。必须保证共享 faces 上的 triangulations 相容；任意独立选择 local triangulations 未必能够 gluing。

任意 triangulation $T\prec S_t$ 都给出一个满足 (4.6) 的 feasible GKZ vector，但它未必是 $\Sigma_k$ 的 vertex，因为非正则 triangulation 的 GKZ vector 可能位于 face 内部。若目标是快速扩大 active hull，优先生成 regular refinements。

#### Batch selection

$F_t$ 可能包含指数多个 regular triangulations，所以不能枚举整个 exposed face。推荐每轮只加入一个小 batch，例如 $2$ 到 $10$ 个，并优先选择彼此差异较大的 GKZ vectors。

可以按以下准则贪心选择：

- 与现有 active hull 的距离较大；
- 新方向 $g_T-x_t$ 与已选择方向近似线性无关；
- 在 symmetry orbits 中没有重复；
- 能改善 active Gram matrix 的条件数。

加入 batch 后，fully corrective QP 会自动把无用 vertices 的系数压到零，再通过 active-set pruning 删除它们。

因此推荐的实际版本是：

> batched fully corrective Frank--Wolfe：每次 lower-hull 计算得到一个 exposed subdivision，在其中采样少量多样化的 regular refinements，将多个 GKZ vertices 同时加入 active set，再作 full correction。

## 5. Primal-dual 解释和误差证书

原问题是

$$
p^*
:=
\min_{g\in\Sigma_k}
\frac12\|g\|^2
=
\frac12\|\sigma_k\|^2.
\tag{5.1}
$$

定义 dual functional

$$
D_k(f)
:=
H_k(f)-\frac12\|f\|^2.
\tag{5.2}
$$

对任意 $f$ 和任意 $g\in\Sigma_k$，

$$
H_k(f)
\leq
\langle f,g\rangle
\leq
\frac12\|f\|^2+\frac12\|g\|^2.
$$

因此 $D_k(f)\leq p^*$。在 $f=\sigma_k$ 处，由 (3.2)，

$$
D_k(\sigma_k)
=
\frac12\|\sigma_k\|^2
=
p^*.
$$

所以 dual problem 是

$$
\max_{f\in\mathbb R^{P_k}}D_k(f),
$$

而其唯一 maximizer 是 $\sigma_k$。

对于算法产生的可行点 $x\in\Sigma_k$，

$$
\begin{aligned}
\frac12\|x\|^2-D_k(x)
&=
\|x\|^2-H_k(x)\\
&=
\delta(x).
\end{aligned}
$$

因此

$$
0
\leq
\frac12\|x\|^2-\frac12\|\sigma_k\|^2
\leq
\delta(x).
\tag{5.3}
$$

另一方面，由 $\sigma_k$ 的投影变分不等式，

$$
\begin{aligned}
\frac12\|x\|^2-\frac12\|\sigma_k\|^2
&=
\langle\sigma_k,x-\sigma_k\rangle
+\frac12\|x-\sigma_k\|^2\\
&\geq
\frac12\|x-\sigma_k\|^2.
\end{aligned}
$$

结合 (5.3)，得到

$$
\boxed{
\|x-\sigma_k\|
\leq
\sqrt{2\delta(x)}.
}
\tag{5.4}
$$

所以算法具有可计算、可验证的 $\ell^2$ 误差上界。

## 6. 利用 $P_k$ 的格点结构

### 6.1 使用整数点集 $kP\cap M$

实际计算时，可以使用

$$
A_k:=kP\cap M
$$

代替分数坐标 $P_k$。simplex volumes 可以通过整数行列式计算，从而：

- 使用 exact orientation predicates；
- 避免 lower-hull 判定中的浮点误差；
- 精确计算 rational GKZ vectors；
- 在算法末期恢复 exact rational answer。

### 6.2 动态 regular triangulation

相邻迭代的 height vectors 往往比较接近。可以保留上一轮 triangulation，并在 heights 改变后通过 bistellar flips 恢复 regularity，而不是每次重新计算整个 lower convex hull。

在二维，主要操作包括：

- quadrilateral diagonal flips；
- 点的插入和删除；
- circuit inequalities 的局部检查。

对稠密格点集，这通常是最重要的工程加速。

### 6.3 沿整除关系 warm start

如果已经求出 level $k$ 的 active representation，而目标 level 是 $kl$，则

$$
J_{k,l}\big(\Sigma(P_k)\big)
$$

是 $\Sigma(P_{kl})$ 的一个 face。

把 coarse triangulations 看成忽略新增点的 fine triangulations，就得到 level $kl$ 的可行 initial active set。

因此适合按整除链计算，例如

$$
1,\ 2,\ 4,\ 8,\ldots.
$$

每一级都可以复用前一级的 active triangulations、convex coefficients 和 triangulation data structure。

### 6.4 使用对称群

设有限群 $G$ 保持 $P$。由于 $\sigma_k$ 唯一，

$$
\sigma_k\in(\mathbb R^{P_k})^G.
$$

若 invariant height vector 的 oracle 返回 $g_T$，则

$$
\overline g_T
:=
\frac{1}{|G|}
\sum_{\gamma\in G}\gamma g_T
$$

仍属于 $\Sigma_k$，仍是 linear minimizer，而且是 $G$-invariant。

可以只保留 $P_k$ 的 $G$-orbits 上的变量。此时范数为

$$
\|g\|^2
=
\sum_{\mathcal O}
|\mathcal O|\,g(\mathcal O)^2.
$$

### 6.5 去掉 height vector 的仿射部分

对任意仿射函数 $h$，

$$
S(f+h|_{P_k})=S(f).
$$

而且所有 $g\in\Sigma_k$ 都有相同的总质量和一阶矩，所以

$$
\operatorname*{argmin}_{g\in\Sigma_k}
\langle f+h,g\rangle
=
\operatorname*{argmin}_{g\in\Sigma_k}
\langle f,g\rangle.
$$

因此，在调用 lower-hull oracle 前，可以：

1. 从 $x$ 中减去其最佳仿射拟合；
2. 再乘以一个正数作 rescaling。

这不会改变 regular subdivision，但可以显著改善数值条件。对很大的 $k$，这是重要的稳定化步骤。

### 6.6 Active face identification

如果后期 oracle 反复返回细化同一个 regular subdivision $S$ 的 triangulations，可以推测已经识别出包含 $\sigma_k$ 的 face

$$
F(S)\subset\Sigma_k.
$$

此后可以：

- 只在 $S$ 的 maximal cells 内作 flips；
- 只处理共享 faces 上的兼容条件；
- 避免产生不细化 $S$ 的全局 triangulations。

这可能把后期计算分解为若干较小的局部 secondary-polytope 问题。

## 7. 复杂度

设

$$
N_k=|P_k|\sim V_Pk^n.
$$

穷尽所有 triangulations 的数量通常呈指数增长。

本文算法每轮只需要：

1. 一次 $N_k$ 个 lifted points 在 $\mathbb R^{n+1}$ 中的 lower-hull/regular-triangulation 计算；
2. 一次 active-set size 为 $m$ 的小型 quadratic program；
3. 一次 GKZ vector 累加。

在 $n=2$ 时，oracle 是三维 lower convex hull。一般实现的复杂度约为

$$
O(N_k\log N_k)
$$

加上输出 triangulation 的复杂度；典型 triangulation 具有 $O(N_k)$ 个 triangles。

在 $n=3$ 时，需要计算四维 lower hull，最坏复杂度更高，但仍避免了枚举全部 triangulations。

该方案不声称 exact computation 在最坏情况下具有 strongly polynomial complexity。实际实现应结合：

- fully corrective updates；
- away steps；
- active-set pruning；
- warm starts；
- dynamic regular triangulations。

## 8. Exact rational certificate

若最终需要严格的 rational $\sigma_k$，可以采用两阶段策略。

第一阶段用浮点算法找出较小的 active triangulation set

$$
\{T_1,\ldots,T_m\}.
$$

第二阶段：

1. 精确计算 rational GKZ vectors $g_{T_i}$；
2. 用 exact rational arithmetic 求解

   $$
   \min_{\substack{
   \lambda_i\geq0\\
   \sum_i\lambda_i=1
   }}
   \frac12
   \left\|
   \sum_i\lambda_i g_{T_i}
   \right\|^2;
   $$

3. 得到 candidate

   $$
   x=\sum_i\lambda_i g_{T_i};
   $$

4. 对 height vector $x$ 调用一次 exact regular-triangulation oracle；
5. 精确计算

   $$
   H_k(x)
   =
   \min_{g\in\Sigma_k}\langle x,g\rangle;
   $$

6. 验证

   $$
   H_k(x)=\|x\|^2.
   $$

由 (3.4)，一旦最后的等式成立，就可以严格断言

$$
x=\sigma_k.
$$

这个 certificate 不需要知道 $\Sigma(P_k)$ 的其他顶点。

## 9. 推荐实现顺序

第一版程序建议按以下顺序实现：

1. 枚举 $A_k=kP\cap M$；
2. 实现 pulling triangulation；
3. 实现由 triangulation 计算 normalized GKZ vector；
4. 实现 lifted lower-hull oracle；
5. 实现 active-set quadratic program；
6. 实现 gap $\delta(x)$ 和终止判据；
7. 加入 active-set pruning；
8. 加入对称群约化；
9. 加入 level $k\to kl$ 的 warm start；
10. 最后实现 dynamic flips 和 exact rational verification。

第一版不必立即实现 dynamic triangulation。即使每轮重新计算 lower hull，oracle 与 fully corrective QP 也已经从根本上消除了遍历整个 secondary polytope 的瓶颈。

## 10. 结论

shortest GKZ vector 的计算可以改写为一个带 regular-triangulation oracle 的 minimum-norm problem：

$$
\boxed{
\sigma_k
=
\operatorname*{argmin}_{g\in\Sigma(P_k)}
\frac12\|g\|^2.
}
$$

每次迭代只访问一个由 lower hull 产生的 GKZ vector。最优性由

$$
\boxed{
\delta(x)
=
\|x\|^2
-
\min_{g\in\Sigma(P_k)}
\langle x,g\rangle
}
$$

控制，而且

$$
\boxed{
\|x-\sigma_k\|
\leq
\sqrt{2\delta(x)}.
}
$$

推荐的总体方案是：

> fully corrective Frank--Wolfe/Wolfe minimum-norm-point algorithm  
> $+$ regular-triangulation oracle  
> $+$ lattice exact predicates  
> $+$ divisibility warm start  
> $+$ symmetry reduction  
> $+$ dynamic bistellar flips。

它不能消除 secondary polytope 在最坏情形下的组合复杂性，但避免了预先枚举全部 triangulations，并且能够输出严格可验证的近似误差或 exact rational certificate。
