# Optimal test function 的 piecewise-affine 问题：研究思路

## 1. 问题和目前真正已知的内容

令

$$
u:=-\Theta _P.
$$

则 $u$ 是 $P$ 上连续的凸函数，并且是

$$
\Psi(u)=\frac12\int_Pu^2\,dx+\int_{\partial P}u\,d\sigma
$$

的唯一极小元。原问题等价于：$u$ 是否为有限个仿射函数的最大值？下文把这种性质称为 polyhedral 或 piecewise affine；这里先不要求仿射函数的系数有理。

已有结果可以整理为：

1. 对任意凸函数 $f$，
   $$
   \int_{\partial P}f\,d\sigma-\int_Pf\Theta_P\,dx\geq0.
   $$
2. 记
   $$
   \lambda_u:=\sigma_{\partial P}+u\,dx=\sigma_{\partial P}-\Theta_P\,dx.
   $$
   则 $\int_Pf\,d\lambda_u\geq0$ 对所有凸 $f$ 成立；对所有仿射 $f$ 等号成立。
3. 还有
   $$
   \int_Pu\,d\lambda_u=\int_{\partial P}u\,d\sigma+\int_Pu^2\,dx=0.
   $$
4. $u$ 满足 Alexandrov 意义下的齐次 Monge--Ampère 方程：对紧集 $K\Subset\operatorname{int}P$，
   $$
   |\partial u(K)|=0.
   $$
5. 若假设 $\Theta_P$ 已经是 piecewise affine，则其最大仿射区域给出 Székelyhidi 的半稳定分解；这只是条件性结论，不能反过来当作一般证明。

Székelyhidi 的文章明确把 piecewise-linearity 说成未知、甚至可能不成立的问题；Sheng--Yao 的文章也把它列为中心开放问题。因此，下面的目标是把开放问题拆成若干可验证的引理，而不是把 HMA 直接误读成 polyhedral 性。

## 2. 一个关键约化：$u$ 是边界迹的凸包络

设

$$
g:=u|_{\partial P}.
$$

定义边界数据的凸包络

$$
E_g(x):=\sup\{\ell(x):\ell\text{ 仿射且 }\ell\leq g\text{ on }\partial P\}.
$$

在现有的 HMA 结论上，可以证明

$$
u=E_g.
$$

证明思路如下。取 $x\in\operatorname{int}P$ 及 $u$ 在 $x$ 处的支撑仿射函数 $\ell$，令

$$
C_\ell:=\{y\in P:u(y)=\ell(y)\}.
$$

Li--Zhou 使用的极点引理表明：若 $u$ 是 HMA 广义解，则 $C_\ell$ 的每个极点都在 $\partial P$。有限维 Carathéodory 定理于是给出 $x\in\operatorname{conv}(C_\ell\cap\partial P)$。因此 $\ell$ 是一个受边界迹 $g$ 约束的仿射函数，得到 $u\leq E_g$。反过来，$u$ 本身是一个边界值为 $g$ 的凸函数，而 $E_g$ 是所有这样的凸函数的上确界，所以 $u\geq E_g$。两边相等。

这个约化给出一个很有用的等价性：

$$
\boxed{
u\text{ piecewise affine}
\Longleftrightarrow
g=u|_{\partial P}\text{ piecewise affine}.
}
$$

其中“$\Leftarrow$”来自下面的标准事实：若 $g$ 在每个 facet 的有限多面体剖分上仿射，则条件 $\ell\leq g$ 只需在有限个边界顶点检查；可行的仿射函数集合是有限个线性不等式给出的 polyhedron，其支撑函数 $E_g$ 是有限个仿射函数的最大值。

所以真正的核心问题不是先研究内部的 HMA，而是：

> **Facet-trace 问题：** 能否证明 $u|_F$ 在每个 facet $F$ 上是 piecewise affine？

这也解释了为什么“零 Monge--Ampère 测度”本身不够：例如在矩形上 $u(x,y)=x^4$ 是凸函数，且 $\det D^2u=0$，但显然不是 piecewise affine。额外的边界变分条件必须被实质性地使用。

## 3. 变分条件的一个精确证书

下面的命题适合同时做证明和找反例。

**KKT/凸序证书。** 若凸连续函数 $v$ 满足

$$
\lambda_v:=\sigma_{\partial P}+v\,dx,
$$

并且

$$
\int_P f\,d\lambda_v\geq0\quad(\forall f\text{ convex}),
\qquad
\int_Pv\,d\lambda_v=0,
$$

那么 $v$ 就是 $\Psi$ 的唯一极小元。事实上，对任意凸 $w$，

$$
\Psi(w)-\Psi(v)
=\frac12\|w-v\|_{L^2(P)}^2+\int_Pw\,d\lambda_v-\int_Pv\,d\lambda_v\geq0.
$$

反过来，极小元必满足这些条件。这把反例问题变成了一个明确的构造问题：找一个非 polyhedral 的凸函数 $v$，使得 $\lambda_v$ 在凸函数序下为正，并且对 $v$ 达到等号。

## 4. 正向路线 A：先证明 facet trace lemma

这是最直接、也最值得先攻的路线。

### A1. 局部化到一个 facet

固定 facet $F$，取坐标使得 $F=\{x_n=0\}$，$P$ 在 $x_n\geq0$ 一侧。对 $F$ 上的切向凸函数 $h$，构造只在 $F$ 附近改变边界值的凸包络扰动 $u_{\varepsilon,h}$。将

$$
\left.\frac{d}{d\varepsilon}\right|_{0+}\Psi(u_{\varepsilon,h})\geq0
$$

写出来，目标是得到一个 $F$ 上的低维变分不等式。这里需要处理两个技术点：

- 凸包络对边界数据的方向导数；
- 体积分项 $\int_Pu\,\dot u\,dx$ 在法向投影下产生的响应测度。

如果得到的 facet 响应测度是分片多项式或分片常数，那么低维问题有希望强迫 $D^2(g|_F)$ 只集中在有限个位置。

### A2. 二维情形的目标引理

在 $n=2$ 时，facet 是区间。最理想的形式是证明

$$
D^2(g|_F)=\sum_{j=1}^N a_j\,\delta_{s_j},
\qquad a_j\geq0,
$$

即边界迹的二阶分布导数是有限原子测度；这就等价于 $g|_F$ 分段仿射。再结合第 2 节的凸包络约化，立即得到二维的 polyhedral 性。

可以先只证明以下弱版本：若 $g|_F$ 在某个开区间上严格非仿射，则能构造一个凸扰动，严格降低 $\Psi$。这会直接排除连续曲率区间。

### A3. 从边界到内部的有限性

若每个 facet 上的迹都已分片仿射，有限个 facet 的有限剖分合起来就是有限边界剖分；凸包络 $E_g$ 因而是有限个仿射函数的最大值。这里不再需要额外的内部正则性定理。

## 5. 正向路线 B：利用 martingale decomposition 强迫有限个 cell

由 Cartier--Fell--Meyer 型测度主化定理，$\lambda_u$ 可以写成

$$
\lambda_u=\int_P(T_x-\delta_x)\,d\nu(x),
$$

其中 $T_x$ 是重心为 $x$ 的概率测度。对 $u$ 使用 Jensen 不等式：

$$
\int u\,d\lambda_u
=\int_P\left(\int u\,dT_x-u(x)\right)d\nu(x)=0.
$$

因此对 $\nu$-几乎处处的 $x$，$u$ 必须在 $\operatorname{conv}(\operatorname{supp}T_x)$ 上仿射。也就是说，martingale transport 的支撑正好落在 $u$ 的 affine cells 中。

一个可能的有限性定理是：

> 若 $P$ 的边界测度只由有限个平面 facet 组成，且 $u$ 连续、梯度有界、$\lambda_u$ 满足上面的等号条件，则存在一个有限的 martingale decomposition，其支撑凸包给出有限个多面体 cell。

这条表述目前不能直接引用，正是需要证明的新引理。主要障碍有两个：

1. Cartier--Fell--Meyer 分解一般不唯一，也不保证有限支撑；
2. 边界测度是有限 facet 测度，并不自动排除连续的一参数 transport foliation。

可操作的办法是选取一个极值 martingale transport，证明其每个 transport component 至少包含一个外边界 facet（这与 Sheng--Yao 的 no-island 观察一致），然后用 facet 对的有限组合排除无限分支。

## 6. 正向路线 C：二维 developable 凸曲面的分类

在二维，$\operatorname{MA}(u)=0$ 意味着图像是 developable convex surface。非仿射部分应由一族线段（rulings）组成。可按下列方式具体化：

1. 取一族 affine contact segments $[p(s),q(s)]$，其端点位于 $\partial P$；
2. 由于 $P$ 是多边形，在一个子区间内 $p(s)$、$q(s)$ 分别落在固定的两个边上；
3. 用
   $$
   X(s,t)=(1-t)p(s)+tq(s)
   $$
   参数化 ruling foliation；
4. 把 $\lambda_u$ 的 martingale 等号和边界测度的常密度代入，得到关于端点参数和 affine slope 的一维方程。

如果能证明该方程迫使 $p(s)$、$q(s)$ 或 slope 在每个组合区间上为常数/仿射函数，就能排除弯曲 ruling，只剩有限条 crease，进而证明 $u$ 分片仿射。

这一方法的优点是几何对象有限：只需逐一处理固定的边对。建议先做三角形和四边形，再尝试一般多边形。

## 7. 正向路线 D：一维切片和 ridge ansatz

为了快速判断反例是否可能，先研究

$$
u(x)=\phi(\ell(x)),
$$

其中 $\ell$ 是线性函数、$\phi$ 是一维凸函数。令

$$
A(t):=\mathcal H^{n-1}(P\cap\{\ell=t\}),
$$

并令 $B(t)$ 为边界测度沿 $\ell$ 的推前密度，则泛函的主体变成

$$
J(\phi)=\frac12\int_I A(t)\phi(t)^2\,dt+\int_I B(t)\phi(t)\,dt
$$

外加端点项。

在二维多边形中，组合区间上 $A(t)$ 是仿射函数，而 $B(t)$ 是常数。若进一步要求 $\lambda_u$ 沿每条切片分别成为一维 martingale measure，则质量条件给出大致形式

$$
\phi(t)=-\frac{B(t)}{A(t)}.
$$

当 $B$ 为正常数、$A$ 为正的仿射函数时，右侧是凹函数而不是凸函数。因此最简单的 ridge ansatz 通常只会产生仿射或分段仿射极小元。这是一个有用的负面证据：寻找反例时不能只试 $u=\phi(\ell)$，需要允许 transport fibers 之间相互耦合。

在一维变量表示中，还可以把

$$
\phi(t)=a+bt+\int_I(t-s)_+\,d\mu(s),
\qquad \mu\geq0,
$$

代入 $J$。KKT 条件转化为一个关于 $s$ 的积分势函数：在 $\operatorname{supp}\mu$ 上等于零，在其余位置具有固定符号。由于 $A,B$ 分片多项式，若能证明该势函数在每个组合区间不可能有连续零区间，就可得到 $\mu$ 有限原子，从而得到分段仿射性。

## 8. 反例路线：直接构造凸序证书

如果猜想不成立，最有效的反例不是只给出一个 HMA 函数，而是直接构造第 3 节的证书：

1. 选取非 polyhedral 的凸函数 $v$，例如一个 developable envelope 或一族弯曲 rulings 的凸包络；
2. 计算
   $$
   \lambda_v=\sigma_{\partial P}+v\,dx;
   $$
3. 构造显式 martingale decomposition，证明
   $$
   \int f\,d\lambda_v\geq0
   $$
   对所有凸 $f$ 成立；
4. 检查
   $$
   \int v\,d\lambda_v=0.
   $$

对 $v=\phi(\ell)$，可以尝试沿切片做 martingale disintegration；上节说明二维多边形的最简单切片模型很可能被凸性排除。更有希望的是：

- 让 transport segment 的方向随位置变化；
- 在三维中使用二维 transport sheets；
- 选择有大量平行 facet 的 prism 或截角多面体，使不同切片之间的质量和重心条件能够耦合。

反例搜索的关键是同时验证所有凸函数，而不是只测试仿射函数或 hinge 函数。有限网格上的凸二次规划可以先发现候选，再用 martingale decomposition 做解析证书。

## 9. 有限维数值实验方案

数值实验的目的不是替代证明，而是快速判断二维是否存在弯曲 ruling，以及寻找最小的候选多面体。

### 9.1 离散变量和约束

取 $P$ 的三角剖分，节点为 $x_i$。给每个节点设置函数值 $u_i$ 和离散次梯度 $p_i$，使用支撑平面约束

$$
u_j\geq u_i+p_i\cdot(x_j-x_i)
$$

（至少对相邻节点，验证时再对全部节点检查）。这样得到离散凸函数。

### 9.2 离散泛函

用精确的单元质量矩阵和边界质量矩阵，写成

$$
\Psi_h(u)=\frac12\mathbf u^{\mathsf T}M\mathbf u+\mathbf b^{\mathsf T}\mathbf u.
$$

这是带线性不等式约束的严格凸二次规划，故离散极小元唯一。

### 9.3 应观察的量

- 相邻三角形梯度跳跃的支持是否在网格加密后稳定为有限条线；
- 活跃 crease 的数量是否稳定，还是随网格线性增长；
- 次梯度像是否收敛到有限个二维区域，还是收敛到一条曲线；
- 将离散 $\Theta_h$ 代入对偶问题后，最大凸函数违约量是否趋于零。

若 crease 数量稳定且位置趋于固定，可以猜测一个有限 subdivision，再用精确积分验证。若 crease 数量持续增长、梯度像趋于曲线，则应转向反例路线。

## 10. 有限 subdivision 的解析验证算法

给定一个候选多面体剖分 $P=\bigcup_iQ_i$，可以完全有限地检查它是否产生 $\Theta_P$：

1. 在每个 $Q_i$ 上设未知 affine density $\theta_i$；
2. 用局部 affine moment 条件确定 $\theta_i$：
   $$
   \int_{Q_i}a\,\theta_i\,dx
   =\int_{\partial Q_i\cap\partial P}a\,d\sigma
   \quad(\forall a\text{ affine});
   $$
3. 检查相邻 cell 上的 $\theta_i$ 是否拼成全局凹函数；
4. 检查每个 $(Q_i,d\sigma_i)$ 的半稳定不等式
   $$
   \int_{\partial Q_i\cap\partial P}f\,d\sigma
   \geq\int_{Q_i}f\theta_i\,dx
   \quad(\forall f\text{ convex});
   $$
5. 检查
   $$
   \int_{\partial P}\Theta\,d\sigma=\int_P\Theta^2\,dx.
   $$

若这些条件成立，把 cell 不等式相加即可得到全局 KKT 证书；由唯一性，候选函数就是 $\Theta_P$。在低维中，第 4 步可通过离散 convex-order LP，再对候选结构做精确证明。

这条算法路线还有一个实际用途：枚举小尺度 lattice dilation 的 regular subdivisions，寻找稳定的 active subdivision。它和 secondary polytope 的组合结构相容，但要注意：固定 dilation 的有限性并不自动推出原问题的有限性，需要证明 active subdivision 在尺度上稳定。

## 11. 建议的优先级

### 第一阶段：马上可做

1. 把第 2 节的“边界凸包络引理”写成正式命题，并补齐对 HMA contact set 的引用和 Carathéodory 论证。
2. 对三角形、矩形、一般四边形做精确二维数值实验。
3. 证明一维 ridge ansatz 的 KKT 条件，并确认其极小元确实分段仿射。

### 第二阶段：最有希望的正向结果

1. 证明二维 facet-trace lemma；哪怕先只处理三角形。
2. 对固定边对的 ruling family 推导端点参数的 ODE。
3. 将二维结果推广到“边界迹分片仿射 $\Rightarrow$ 内部凸包络分片仿射”的一般维数结论。

### 第三阶段：决定猜想真假的分叉

- 如果所有二维实验都显示有限 active creases，并且 ruling ODE 排除连续族，应集中证明二维定理。
- 如果出现 crease 数量随网格增长的稳定现象，应立刻寻找 martingale 证书，而不是继续假设 polyhedral 性。
- 若二维成立而三维出现连续 transport sheet，可能的最终结论是“二维成立、高维一般不成立”。

## 12. 需要避免的几个逻辑陷阱

1. **HMA 不推出 piecewise affine。** $u(x,y)=x^4$ 已足以否定这个推理。
2. **PL 逼近不推出极限 PL。** 凸函数的 PL 逼近可以收敛到任意连续凸函数。
3. **只检查仿射测试函数不够。** 仿射函数只给出 $n+1$ 个矩条件。
4. **条件性的半稳定分解不能反用。** Sheng--Yao/Székelyhidi 的 cell 分解是在假设 $\Theta_P$ PL 后得到的。
5. **实系数 PL 和有理 PL 要区分。** 前者回答当前分析问题；后者才直接对应有限 exponent 的 toric test configuration。

## 13. 相关文献

- Székelyhidi, *Optimal test-configurations for toric varieties*, J. Differential Geom. 80 (2008), DOI: [10.4310/jdg/1226090485](https://doi.org/10.4310/jdg/1226090485)。文章明确指出最优 destabilizer 是否 piecewise linear 未知，并在条件性假设下建立 semistable decomposition。
- Donaldson, *Scalar Curvature and Stability of Toric Varieties*, J. Differential Geom. 62 (2002), DOI: [10.4310/jdg/1090950195](https://doi.org/10.4310/jdg/1090950195)。
- Li--Zhou, *K-stability and polystable degenerations of polarized spherical varieties*, arXiv: [2111.04269](https://arxiv.org/abs/2111.04269)。其中的 contact-set / HMA 方法是第 2 节边界凸包络约化的来源。
- Sheng--Yao, *Properties of optimal test-function for toric varieties*，本地参考文件：`piecewise_affiness.refs/Sheng-Yao_Properties of optimal test-function for toric varieties.tex`。其中包含连续性、有界次梯度、HMA 以及条件性半稳定分解。
- Gale--Klee--Rockafellar, *Convex functions on convex polytopes*, Proc. Amer. Math. Soc. 19 (1968), DOI: [10.2307/2035330](https://doi.org/10.2307/2035330)。可用于边界分片仿射数据的凸包络有限性。

本笔记没有使用 `.refs` 中与当前问题无直接关系的大量 bib 条目。目录中出现的 Yao 关于 secondary polytopes 的 2025--2026 条目值得继续核对，但在未取得正文前，不把其可能结论当作已知定理。
