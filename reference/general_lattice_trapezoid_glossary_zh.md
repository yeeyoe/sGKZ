# 一般格点多面体上的 lattice trapezoid estimate：术语与分解图

本文说明 `main.tex` 中 Lemma `general lattice trapezoid estimate` 的几何分解。设 $P\subset\mathbb R^n$ 是满维格点多面体，$P_k=P\cap k^{-1}\mathbb Z^n$，$g\geq0$ 是凸函数。这里的目标不是重新证明该引理，而是说明证明中每个区域和每类格子的定义，以及它们分别产生哪一项误差。

![一般格点多面体网格分解示意](general_lattice_trapezoid_mesh.svg)

图是二维横截面：二维时，$\partial^{>1}P$ 就是顶点。高维时，它表示余维至少二的面之并；图中的深色角区应理解为这些面的一个固定小邻域。

## 基本边界数据

### $\partial^{>1}P$：余维至少二的骨架

定义为

$$
\partial^{>1}P=\bigcup\{G\subset P:\operatorname{codim}G\geq2\}.
$$

它是两个或更多 facet 同时相交的位置。远离它的边界点只属于一个 facet，因而局部只有一个法向方向，能够作一维 trapezoid 配对。靠近它时，多个 facet 同时出现，局部网格不再可由单一法向的半盒子描述。

它在证明中的作用是划出唯一不能通过普通 facet 计算完全控制的区域；所有这部分贡献保留为角点缺陷，而不把它误认为小的连续模误差。

### Primitive facet defining function $l_F$

对每个 facet $F$，选取 primitive integral affine defining function $l_F$，使得

$$
P\subset\{l_F\geq0\},\qquad F=P\cap\{l_F=0\},
$$

且 $l_F$ 的线性部分在对偶格中是 primitive 的。primitive 意味着其在格子上的取值步长恰为 $1$；因此 $l_F(x)$ 自然标记从 $F$ 向内的整数层。

它的作用是把 regular facet 附近的法向方向离散化为 $0,1/k,2/k,\ldots$。这正是将高维局部和式化成一维 trapezoid 配对的坐标。

### Integral inward normal $v_F$

由 primitive 性，可取 $v_F\in\mathbb Z^n$ 使

$$
l_F(v_F)=1.
$$

这里取向使 $v_F$ 指向 $P$ 的内部。沿 $v_F$ 的平移把 $l_F$ 的层数增加一，因此对 $y\in F$，$y+t v_F$ 给出 facet collar 中的法向参数化。

它的作用是保证这个参数化与全局格子相容：在尺度 $1/k$ 上，相邻的法向层正好相距 $v_F/k$，不需要 Delzant 条件。

### Integral-affine chart：整仿射坐标图

在只接触 facet $F$ 的区域，选择 $\ker l_F\cap\mathbb Z^n$ 的整基，并以 $l_F$ 作为第一坐标。得到的坐标变换属于整仿射变换，故它把 $k^{-1}\mathbb Z^n$ 仍送到标准的 $k^{-1}$ 网格。

它的作用是把 collar 中的普通格子写成标准盒子。盒子的凸插值积分给出顶点权重；在 facet 上，外侧半盒子不存在，故该顶点总权重是 $1/2$，这就是边界主项系数的来源。

## 几何区域

以下先取小数 $\eta,\delta>0$。它们只依赖 $P$ 和预先指定的 $K_F\Subset\operatorname{relint}F$，并先固定，再令 $k\to\infty$。

### Truncated facet region $F_F^\eta$

若 facets 由 $l_G$ 标记，定义

$$
F_F^\eta=\{y\in F:l_G(y)\geq3\eta\ \text{for every }G\ne F\}.
$$

它是从 $F$ 的边缘切掉固定宽度后的中间部分。通过将 $\eta$ 取得足够小，可令给定的紧集 $K_F$ 包含于 $F_F^\eta$。

它的作用是确保该区域只看见一个 facet。于是建立 collar 时不会触及相邻 facet，也不会把角点现象混进 $1/2$ 的主项。

### Facet collar $\mathcal C_F$

定义

$$
\mathcal C_F=\{y+t v_F:y\in F_F^\eta,\ 0\leq t\leq\delta\}.
$$

这是从截断 facet 沿整内法向推进固定厚度得到的带状区域；图中浅蓝色条带即为两个二维 collar 的示意。

它的作用是承载常规 facet 的主贡献。collar 内除接口附近的少数格子外，都可由完整标准盒子组成，因此能精确记录内点权重 $1$ 和 facet 顶点权重 $1/2$。

### Corner neighborhood $U_\eta(\partial^{>1}P)$

取 $\partial^{>1}P$ 的一个固定小开邻域，例如其距离小于 $\eta$ 的点集，再按需要略作多面体化。图中深色区域表示它。

它的作用是隔离多 facet 相交的区域。一般格点多面体上不试图在这里恢复单一 facet 的 $1/2$ 权重；尺度 $1/k$ 的相关格点被归入 $Z_k$，留下显式缺陷项。

### Interior core $\mathcal C_0$

定义为

$$
\mathcal C_0=\overline{P\setminus\left(\bigcup_F\mathcal C_F\cup U_\eta(\partial^{>1}P)\right)}.
$$

它是远离原始边界、只剩有限个人工截断面的紧内部区域；图中淡灰色中心区域表示它。

它的作用是让剩余部分能够在固定的全局整坐标中用标准网格分割。由于它离 $\partial P$ 有正距离，所有由它产生的接口误差都可在一个固定紧集 $Q\Subset\operatorname{int}P$ 上估计。

## 人工切割及其离散化

### Artificial interface：人工接口

collar 的内端和 core 的外边界是人为选定的固定有理多面体超曲面。它们不是 $P$ 的原始边界；图中的虚线就是一个接口。

它的作用是把不同整仿射坐标图下的网格分块拼接起来。由于接口位置通常不经过 $k^{-1}\mathbb Z^n$，不能直接把它当作网格面。

### Mesh interface $H_k$

设 $H$ 是一个 artificial interface。在选定图中，将 $H$ 沿一个坐标方向移动到最近的格点超平面，得到 $H_k$。移动距离为 $O(1/k)$。

它的作用是使 collar 侧与 core 侧各自都可由真正的格子盒组成；但 $H$ 与 $H_k$ 之间会留下薄层，必须单独处理。

### Slab $\operatorname{Slab}_k(H)$

定义为 $H$ 与 $H_k$ 所夹的闭区域：

$$
\operatorname{Slab}_k(H)=\overline{\text{the region between }H\text{ and }H_k}.
$$

它的厚度为 $O(1/k)$，而切向面积有界，故其中只含 $O(k^{n-1})$ 个尺度 $1/k$ 的单元。

它的作用是产生连续模项。将 slab 三角剖分后，与相邻标准盒子所用的参考顶点权重相比，两套权重的总质量相同，因此常数抵消；所有顶点距离为 $O(1/k)$，每个单元的差由 $\omega_{g,Q}(C/k)$ 控制，合计为

$$
Ck^{n-1}\omega_{g,Q}(C/k).
$$

### One-interface cell：单接口单元

这是只接触一个人工接口的 slab 单元，或与之相邻、需要在两套局部网格间转换的单元。它的所有顶点仍落在固定内部紧集内。

它的作用是贡献上面的连续模项，而不是 $C^0$ 项：因为两种插值系数的总质量相同，可以先减去一个固定顶点处的函数值，再使用连续模。

### Junction cell：交汇单元

这是同时接触两个或更多人工接口的单元，例如二维图中虚线端点附近的少量单元。交汇集合的固定几何维数至多为 $n-2$，故其 $1/k$ 邻域内只有 $O(k^{n-2})$ 个这样的单元。

它的作用是贡献

$$
Ck^{n-2}\sup_Q g.
$$

这里不使用权重总质量的精细消去：接口交汇处的剖分有限但复杂，保留粗的 $C^0$ 控制正是该项的来源。

## 角点格子与缺陷

### Corner lattice set $Z_k$

取足够大的常数 $c_P$，定义

$$
Z_k=\{u\in P_k:\operatorname{dist}(u,\partial^{>1}P)\leq c_P/k\}.
$$

常数应大到使任何被设为一旁的 corner cell 的所有格点顶点都属于 $Z_k$。由于余维至少二，$|Z_k|=O(k^{n-2})$。

它的作用是精确记录一般多面体中尚未处理的格点值，而不以 $\sup_Qg$ 或内部连续模替代它们。注意 $Z_k$ 位于原始边界附近，通常不包含于 $Q$。

### Corner defect $D_k(g)$

在 $B_k$ 的归一化中，定义

$$
D_k(g)=k^{1-n}\sum_{u\in Z_k}g(u).
$$

它的作用是说明一般格点多面体估计的准确边界：Cauchy--Schwarz 与 $|Z_k|=O(k^{n-2})$ 给出 $D_k(g)\leq C\|g\|_{k,2}$，足以得到 coercivity；但有界离散 $L^2$ 并不推出 $D_k(g_k)\to0$，所以不能由此单独推出无缺陷的边界 liminf。

## $Q$ 与 $Q'$ 从哪里来

把所有 artificial interfaces、slabs 和 junction cells 的固定几何支撑合成一个闭集 $S\Subset\operatorname{int}P$。这是可能的，因为 collars 在距原始边界正距离的地方截断，而角点邻域已被移除。取 $\varepsilon>0$ 足够小，使

$$
Q=\operatorname{conv}(S)+\varepsilon\overline B,
\qquad
Q'=\operatorname{conv}(S)+2\varepsilon\overline B
$$

仍满足 $Q\Subset Q'\Subset\operatorname{int}P$。其中 $\overline B$ 是单位闭球。对充分大的 $k$，所有普通 interface cell 的顶点都在 $Q$；$Q'$ 只是为了在内点 $C^0$ 与 Lipschitz 估计中给 $Q$ 留出固定距离。

## 这些对象如何合成估计

对普通盒子进行凸插值积分时，内点格点的总系数为 $1$，regular facet 上格点的总系数为 $1/2$。后者只收到内侧半数盒子的贡献，给出

$$
\frac12\sum_F\sum_{u\in K_F\cap P_k}g(u).
$$

随后保留三种彼此不同的负误差：

$$
Ck^{n-2}\sup_Qg,
\qquad
Ck^{n-1}\omega_{g,Q}(C/k),
\qquad
C\sum_{u\in Z_k}g(u).
$$

它们的来源依次是 junction cells、one-interface cells 与 corner cells。因而引理中的下界是

$$
\begin{aligned}
\sum_{u\in P_k}g(u)-k^n\int_Pg\,\mathrm dx
\geq{}&
\frac12\sum_F\sum_{u\in K_F\cap P_k}g(u)
-Ck^{n-2}\sup_Qg\\
&-Ck^{n-1}\omega_{g,Q}(C/k)
-C\sum_{u\in Z_k}g(u).
\end{aligned}
$$

乘以 $2k^{1-n}$ 后，这三项仍应分开看：它们分别变为 $Ck^{-1}\sup_Qg$、$C\omega_{g,Q}(C/k)$ 和 $C D_k(g)$。只有在证明 coercivity 时，才分别用内点 $C^0$ 估计、内点 Lipschitz 估计和离散 $L^2$ 控制它们。

