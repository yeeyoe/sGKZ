# 凸函数求积紧性引理

记

$$
\lVert f\rVert _{k,2}^{2}=\frac{1}{k^{n}}
\sum_{u\in P_{k}}f(u)^{2},\qquad
B_{k}(f)=\frac{2}{k^{n-1}}\left(\sum_{u\in P_{k}}f(u)
-k^{n}\int_{P}f\,\mathrm{d}x\right).
$$

于是 $\Psi_k(f)=\frac12\lVert f\rVert_{k,2}^2+B_k(f)$。下面给出凸函数求积引理的紧性形式。注意：不要求控制函数在余维至少为二的面上的取值。

## 引理（凸函数求积紧性）

设 $P$ 为 Delzant 多面体。存在只依赖于 $P$ 的常数 $C_P$，使得对每个充分大的 $k$ 和每个 $f\in\mathcal C(P)$，

$$
B_k(f)\geq-C_P\lVert f\rVert_{k,2},\qquad
\lVert f\rVert_{L^2(P)}\leq C_P\lVert f\rVert_{k,2}.

$$

下文只用到 $k\to\infty$，所以有限多个较小的 $k$ 无关紧要。

再设 $f_k\in\mathcal C(P)$ 且 $\sup_k\lVert f_k\rVert_{k,2}<\infty$。取一列子列后，存在凸函数 $f\in L^2(P)$，使得

$$
f_k\rightharpoonup f\quad\text{在 }L^2(P)\text{ 中弱收敛},
\qquad
f_k\to f\quad\text{在 }\operatorname{int}P\text{ 上局部一致收敛}。
$$

把 $f$ 取为 $P$ 上的下半连续延拓，并允许边界积分为 $+\infty$。对该子列，

$$
\liminf_{k\to\infty}\lVert f_k\rVert_{k,2}^{2}
\geq\lVert f\rVert_{L^2(P)}^{2},\qquad
\liminf_{k\to\infty}B_k(f_k)
\geq\int_{\partial P}f\,\mathrm{d}\sigma.

$$

若进一步满足

$$
\lim_{k\to\infty}\lVert f_k\rVert_{k,2}^{2}
=\lVert f\rVert_{L^2(P)}^{2},

$$

则 $f_k\to f$ 在 $L^2(P)$ 中强收敛。

## 证明

证明分四步。写 $P$ 的 facets 为 $F_1,\ldots,F_N$，取 primitive integral affine 函数 $l_r$，使得

$$
P=\{l_r\geq0,\ 1\leq r\leq N\},\qquad F_r=\{l_r=0\}.
$$

对 $\delta>0$，令

$$
P^\delta=\{x\in P:l_r(x)\geq\delta\text{ 对所有 }r\}.
$$

改变欧氏范数只会改变常数，因此以下所有常数只依赖于 $P$。Delzant 条件只在后面的格点梯形估计中使用。

### 第一步：两个强制性估计

固定一个顶点为格点的 $P$ 的三角剖分 $\mathcal T$。其标准 $k$ 重 edgewise subdivision 给出 $\mathcal T_k$，并满足：每个顶点在 $P_k$ 中；每个单形直径至多为 $C/k$；经 $k$ 倍缩放后只出现有限种仿射形状；每个顶点至多属于 $C$ 个单形。

若 $\Delta\in\mathcal T_k$ 的顶点为 $u_0,\ldots,u_n$，且 $x=\sum_i\lambda_i u_i$，凸性给出

$$
h(x)\leq\sum_{i=0}^n\lambda_i h(u_i).

$$

因此

$$
h_+(x)^2\leq\left(\sum_i\lambda_i|h(u_i)|\right)^2
\leq\sum_i\lambda_i h(u_i)^2.
$$

在 $\Delta$ 上积分并对所有单形求和，得到

$$
\int_P h_+^2\,\mathrm{d}x
\leq\frac{C}{k^n}\sum_{u\in P_k}h(u)^2.

$$

若 $g\geq0$ 为凸函数，$E\subset P$ 可测，令 $E_k$ 为所有与 $E$ 相交的 $\mathcal T_k$ 单形之并。$E_k$ 的每个顶点距 $E$ 至多为 $C/k$，故同样有

$$
\int_Eg^2\,\mathrm{d}x
\leq\frac{C}{k^n}
\sum_{\substack{u\in P_k\\\operatorname{dist}(u,E)\leq C/k}}g(u)^2.

$$

下面把一般凸函数化为非负情形。取 $p_0\in\operatorname{int}P$ 与 $r>0$，使 $B(p_0,5r)\Subset P$。对充分大的 $k$，取 $p_k\in P_k$，满足 $|p_k-p_0|\leq C/k$，并令

$$
Z_k=\{z\in\mathbb Z^n:|z|\leq rk\}.
$$

于是 $|Z_k|\geq ck^n$，且对 $z\in Z_k$，点 $p_k\pm z/k$ 和 $p_k+2z/k$ 都在 $P_k$ 中。凸性给出

$$
h(p_k)\leq\frac{h(p_k+z/k)+h(p_k-z/k)}2,
\qquad
h(p_k)\geq2h(p_k+z/k)-h(p_k+2z/k).
$$

对 $Z_k$ 取平均。对 $j\in\{-1,1,2\}$，由 Cauchy--Schwarz 不等式和映射 $z\mapsto p_k+jz/k$ 的单射性，

$$
\frac1{|Z_k|}\sum_{z\in Z_k}|h(p_k+jz/k)|
\leq |Z_k|^{-1/2}\left(\sum_{u\in P_k}h(u)^2\right)^{1/2}
\leq C\lVert h\rVert_{k,2}.
$$

故

$$
|h(p_k)|\leq C\lVert h\rVert_{k,2}.

$$

取 $q_k\in\partial h(p_k)$。若 $q_k\ne0$，令 $e_k=q_k/|q_k|$，并令

$$
S_k=P_k\cap B(p_k+3re_k,r).
$$

当 $k$ 充分大时，$B(p_k+3re_k,r)\subset B(p_0,5r)\subset P$，所以 $S_k$ 就是该球中的所有 $k^{-1}\mathbb Z^n$ 格点。固定半径的球至少包含 $ck^n$ 个这样的格点，且该下界与 $e_k$ 无关。对 $u\in S_k$，次梯度不等式和 (7) 给出

$$
h(u)\geq h(p_k)+\langle q_k,u-p_k\rangle
\geq-C\lVert h\rVert_{k,2}+2r|q_k|.
$$

若 $|q_k|\leq C_0r^{-1}\lVert h\rVert_{k,2}$，所需估计已经成立；否则 $h(u)\geq r|q_k|$，于是

$$
k^n\lVert h\rVert_{k,2}^2=\sum_{u\in P_k}h(u)^2
\geq |S_k|r^2|q_k|^2\geq ck^nr^2|q_k|^2.
$$

两种情形下均有 $|q_k|\leq C\lVert h\rVert_{k,2}$。令

$$
\ell_k(x)=h(p_k)+\langle q_k,x-p_k\rangle,
\qquad g_k=h-\ell_k\geq0.
$$

则

$$
\lVert\ell_k\rVert_{C^0(P)}+\lVert\ell_k\rVert_{k,2}+\lVert g_k\rVert_{k,2}
\leq C\lVert h\rVert_{k,2}.

$$

这是因为 (7)、$q_k$ 的估计和 $P$ 的有界性控制 $\lVert\ell_k\rVert_{C^0(P)}$；$|P_k|\leq Ck^n$ 控制 $\lVert\ell_k\rVert_{k,2}$；最后对 $g_k=h-\ell_k$ 使用离散三角不等式。

由于 $g_k\geq0$，可将 (5) 应用于 $g_k$。结合 (8)，

$$
\begin{aligned}
\lVert h\rVert_{L^2(P)}
&\leq\lVert g_k\rVert_{L^2(P)}+\lVert\ell_k\rVert_{L^2(P)}\\
&\leq C\bigl(\lVert g_k\rVert_{k,2}+\lVert\ell_k\rVert_{C^0(P)}\bigr)
\leq C\lVert h\rVert_{k,2}.
\end{aligned}
$$

即

$$
\lVert h\rVert_{L^2(P)}\leq C\lVert h\rVert_{k,2}.

$$

还要用到内部估计：若 $K\Subset K'\Subset\operatorname{int}P$，且 $h$ 为凸函数，则

$$
\sup_K|h|+\operatorname{Lip}_K(h)
\leq C_{K,K'}\lVert h\rVert_{L^2(K')}.

$$

为说明细节，取 $r>0$，使以 $K$ 中任一点为中心、半径 $4r$ 的球都包含于 $K'$。Jensen 不等式给出

$$
h(x)\leq |B(x,r)|^{-1}\int_{B(x,r)}h
\leq C_r\lVert h\rVert_{L^2(K')}.
$$

这是上界。若 $h(x)=-M$，并记 $K$ 的 $2r$ 邻域上的上述上界为 $A$，则对 $z\in B(x,r)$，

$$
h\left(\frac{x+z}{2}\right)\leq\frac{h(x)+h(z)}2\leq\frac{-M+A}{2}.
$$

若 $M>2A$，当 $z$ 变化时，在 $B(x,r/2)$ 上有 $|h|\geq M/4$，故 $M\leq C_r\lVert h\rVert_{L^2(K')}$；若 $M\leq2A$，结论显然。因此 $|h|$ 在 $K$ 的 $2r$ 邻域上一致有界。最后，对 $q\in\partial h(x)$，在 $x\pm re_i$ 处使用次梯度不等式：

$$
h(x)\pm r\langle q,e_i\rangle\leq h(x\pm re_i).
$$

$C^0$ 界控制 $q$ 的所有分量，再由线段上的次梯度不等式得到 (10) 的 Lipschitz 估计。

下面证明格点梯形估计。若给定每个 facet 的紧多面体 $K_F\Subset\operatorname{relint}F$，则存在
$Q\Subset Q'\Subset\operatorname{int}P$，使每个非负凸函数 $g$ 满足

$$
\sum_{u\in P_k}g(u)-k^n\int_Pg\,\mathrm{d}x
\geq\frac12\sum_F\sum_{u\in K_F\cap P_k}g(u)
-Ck^{n-2}\bigl(\sup_Qg+\operatorname{Lip}_Q(g)\bigr).

$$

下面的二维示意图展示 collar/core 分解、人工分隔面、one-layer slab、普通 interface cells 和 junction cells：

![二维 Delzant mesh 示意图](/Users/yaoy/Documents/sGKZ/paper/delzant_mesh_claim.svg)

当 $n=1$ 时，最后的误差项不存在。

右端的结构来自以下局部计数：内部格点邻接 $2^n$ 个 boxes，而 facet 相对内部格点只邻接 $2^{n-1}$ 个向内 boxes。每个 box 对每个顶点贡献 $2^{-n}$，所以总权重分别为 $1$ 和 $1/2$。坐标变换界面有 $O(k^{n-1})$ 个单元；普通界面单元中参考权重质量相等，常数项抵消，每个单元留下 $O(k^{-1}\operatorname{Lip}_Q g)$，总计 $O(k^{n-2}\operatorname{Lip}_Q g)$。两个界面的交集余维至少二，只产生 $O(k^{n-2})$ 个 junction cells，每个可用 $O(\sup_Qg)$ 控制。

#### Delzant 网格构造

先证明存在权重 $a_k(u)\geq0$，使

$$
a_k(u)\leq1,qquad a_k(u)\leq\frac12\quad(u\in K_F\cap P_k),

$$

并且

$$
k^n\int_Pg\,\mathrm{d}x
\leq\sum_{u\in P_k}a_k(u)g(u)+R_k(g),
\qquad
|R_k(g)|\leq Ck^{n-2}\bigl(\sup_Qg+\operatorname{Lip}_Q(g)\bigr).

$$

当 $n=1$ 时取 $R_k=0$。

令 $G$ 为余维 $r$ 的面，$F_1,\ldots,F_r$ 为包含它的 facets。Delzant 条件保证 primitive affine 函数 $l_{F_1},\ldots,l_{F_r}$ 可延拓为局部 integral affine 坐标

$$
(x_1,\ldots,x_n)=(l_{F_1},\ldots,l_{F_r},y_{r+1},\ldots,y_n).
$$

坐标变换的线性部分属于 $\mathrm{GL}(n,\mathbb Z)$，所以保持 $k^{-1}\mathbb Z^n$ 和 Lebesgue 测度；局部地 $P$ 由 $x_1,\ldots,x_r\geq0$ 给出。

按余维递减构造一个 $\partial P$ 的固定 collar：顶点附近取坐标正交锥；边的未覆盖部分取边坐标图与向内法向小区间的乘积；依次处理更高维的面。处理面 $G$ 时只剩下 $\operatorname{relint}G$ 的紧子集，有限次细分后可用一个坐标图覆盖。重叠处保留先前的法向坐标，只细分 $G$ 的切向格，因此得到 half-open integral-affine box 分解。余维 $r$ 面上的点恰邻接 $2^{n-r}$ 个向内 boxes。取 collar 足够宽，使所有 $K_F$ 落在 facet boxes 内并避开人为切向分隔面。

沿内部固定的 rational polyhedral 超曲面截断 collar，核心部分使用一个固定 integral affine 坐标图。经过任意小的 rational 位移，可使各人为切割处于一般位置：任意两个不同切割要么不交，要么交集余维至少二；相对开集上重合的片段合并为一个。取 $Q\Subset\operatorname{int}P$ 包含这些内部超曲面的固定邻域，再取 $Q'$ 满足 $Q\Subset Q'\Subset\operatorname{int}P$。

在层级 $k$ 中用 $x_i=j/k$ 细分各 box。每个人为分隔面替换为两侧最近的网格面之并，并以一层 slab 填充其间区域。多个阶梯状替换相交处，用与交集相交的网格单元之并扩大；这仍在余维二交集的 $C/k$ 邻域内，且具有格点顶点。沿切向进一步细分，并对每个 slab 多面体作 pulling lattice triangulation，所得满维格单形称为 slab cells。

每个 slab cell 的直径至多为 $C/k$。若其顶点为 $v_0,\ldots,v_n$，则

$$
|C|=\frac1{n!}|\det(v_1-v_0,\ldots,v_n-v_0)|.
$$

由于 $k(v_i-v_0)$ 为整数，行列式是非零整数；因此

$$
\frac1{n!}k^{-n}\leq|C|\leq Ck^{-n}.

$$

令 $\mathcal I_k$ 为接触恰好一个分隔面的单元，令 $\mathcal J_k$ 为其余单元。固定超曲面的 $C/k$ 邻域体积为 $O(k^{-1})$，由 (14) 得

$$
|\mathcal I_k|\leq Ck^{n-1}.

$$

$\mathcal J_k$ 中的单元位于两个不同分隔面交集的 $C/k$ 邻域内，该交集余维至少二，故

$$
|\mathcal J_k|\leq Ck^{n-2}.

$$

充分大 $k$ 时，这些单元及参考系数所涉及的格点都在 $Q$ 内。

对普通 box

$$
R=\prod_{i=1}^n[a_i,a_i+k^{-1}],
$$

令 $v_\varepsilon=(a_i+\varepsilon_i/k)_i$、$t_i=k(x_i-a_i)$，并设

$$
\lambda_\varepsilon(x)=\prod_{i=1}^nt_i^{\varepsilon_i}(1-t_i)^{1-\varepsilon_i}.
$$

则 $x=\sum_\varepsilon\lambda_\varepsilon v_\varepsilon$，$\sum_\varepsilon\lambda_\varepsilon=1$，且 $\int_R\lambda_\varepsilon=2^{-n}k^{-n}$。凸性给出

$$
k^n\int_Rg\,\mathrm{d}x
\leq2^{-n}\sum_{\varepsilon\in\{0,1\}^n}g(v_\varepsilon).

$$

对普通 boxes 求和：内部点的权重为 $2^n2^{-n}=1$；余维 $r$ 面相对内部点的权重为

$$
2^{n-r}2^{-n}=2^{-r}.

$$

特别地，facet 权重为 $1/2$，且不会超过 $1$；$K_F$ 避开人为分隔面，所以其中格点充分大时恰有权重 $1/2$。

若 $C\in\mathcal I_k\cup\mathcal J_k$ 的顶点为 $v_0,\ldots,v_n$，对单形重心坐标积分，得到非负系数 $\beta_{C,u}$，支撑在 $C$ 的顶点上，并满足

$$
k^n\int_Cg\,\mathrm{d}x\leq\sum_{u\in\operatorname{vert}(C)}\beta_{C,u}g(u),
\qquad
\sum_u\beta_{C,u}=k^n|C|\leq C.

$$

下面定义参考系数。把各坐标图的 box mesh 向分隔面另一侧延伸一层，用两个延伸 mesh 之间的中间 polyhedral 超曲面切开 slab。两侧分别归相应坐标图所有，共同面各分一半；多个分隔面相交时，固定坐标图顺序并按字典序分配 half-open sectors。再细分，使每个被分配片段 $E$ 包含于其所有者坐标图的一个 box $R$ 中。

若 $u$ 是 $R$ 的顶点，令 $\lambda_{R,u}$ 为 (17) 中的 multilinear 坐标，并设

$$
q_{E,u}=k^n\int_E\lambda_{R,u}(x)\,\mathrm{d}x;
$$

对其他格点令 $q_{E,u}=0$，定义

$$
\bar\beta_{C,u}=\sum_{E\text{ 是 }C\text{ 的被分配片段}}q_{E,u}.

$$

由 $\lambda_{R,u}\geq0$ 和 $\sum_{u\in\operatorname{vert}(R)}\lambda_{R,u}=1$，有

$$
\sum_u\bar\beta_{C,u}=k^n|C|=\sum_u\beta_{C,u},
\qquad
\sum_u\bigl(|\beta_{C,u}|+|\bar\beta_{C,u}|\bigr)\leq C.

$$

第一式通过先对 $u$、再对分配片段 $E$ 求和得到；第二式由第一式、非负性及 (14) 得到。由于 box 只延伸一层，其顶点距 $C$ 至多为 $C/k$；扩大 $Q$ 后，所有支撑 $\bar\beta_{C,\bullet}$ 的格点在充分大 $k$ 时均位于 $Q$。

把 $(R,u)$ 称为一个 box--corner incidence。普通 box 使用整个 incidence，贡献 $2^{-n}$；被分配片段 $E\subset R$ 使用的系数满足

$$
0\leq q_{E,u}\leq k^n\int_R\lambda_{R,u}\,\mathrm{d}x=2^{-n}.
$$

half-open ownership 规则也施加于 incidences：普通 box 未使用的 incidence 恰分配给一个相邻 sector，且不重复分配；junction 处由固定字典序确定。因此每个格点至多使用 $2^n$ 个 incidence，每个容量为 $2^{-n}$，所以总权重至多为 $1$。余维 $r$ 面上只有 $2^{n-r}$ 个向内 incidences，权重至多为 $2^{-r}$。这证明了 (12)。

固定 $C$ 的顶点 $u_C$。由 (21) 常数项相消；所有出现的格点距 $u_C$ 至多为 $C/k$，故

$$
\begin{aligned}
\left|\sum_u(\beta_{C,u}-\bar\beta_{C,u})g(u)\right|
&=\left|\sum_u(\beta_{C,u}-\bar\beta_{C,u})\bigl(g(u)-g(u_C)\bigr)\right|\\
&\leq\frac{C}{k}\operatorname{Lip}_Q(g).
\end{aligned}

$$

由 (15)，普通界面的总误差为 $Ck^{n-2}\operatorname{Lip}_Q(g)$。对 junction cell，不使用常数相消，只用 (21) 和 $Q$ 上的上界，单元误差至多为 $C\sup_Qg$；由 (16)，总误差至多为 $Ck^{n-2}\sup_Qg$。

令 $a_k(u)$ 为普通 box 系数与参考系数 $\bar\beta_{C,u}$ 的总和（注意不是用 $\beta_{C,u}$）。把每个 $\beta$ 写成 $\bar\beta+(\beta-\bar\beta)$，对所有普通 boxes 和 slabs 的求积不等式求和，再用 (15)、(16)、(22)，得到 (13)。

从 $\sum_{u\in P_k}g(u)$ 中减去 (13)。由 $g\geq0$ 及 (12)，

$$
\begin{aligned}
\sum_{u\in P_k}g(u)-k^n\int_Pg\,\mathrm{d}x
&\geq\sum_{u\in P_k}(1-a_k(u))g(u)-R_k(g)\\
&\geq\frac12\sum_F\sum_{u\in K_F\cap P_k}g(u)
-Ck^{n-2}\bigl(\sup_Qg+\operatorname{Lip}_Q(g)\bigr).
\end{aligned}
$$

这就是 (11)。当 $n=1$ 时没有人为界面，故 $R_k=0$。Delzant 条件只用于保证面坐标 unimodular，从而坐标 boxes 的顶点在 $P_k$ 中、体积为 $k^{-n}$，并且余维 $r$ 面有 $2^{n-r}$ 个向内 incidences。

回到 (8) 中的分解 $h=g_k+\ell_k$。在 (11) 中取 $K_F=\varnothing$ 并乘以 $2/k^{n-1}$。由 (5)、(8)、(10)，

$$
\sup_Qg_k+\operatorname{Lip}_Q(g_k)
\leq C\lVert g_k\rVert_{L^2(Q')}
\leq C\lVert g_k\rVert_{k,2}
\leq C\lVert h\rVert_{k,2}.
$$

所以

$$
B_k(g_k)\geq-\frac{C}{k}\lVert h\rVert_{k,2}
\geq-C\lVert h\rVert_{k,2}.

$$

对仿射函数 $\ell_k$，加权 Ehrhart 多项式的前两项给出

$$
B_k(\ell_k)=\int_{\partial P}\ell_k\,\mathrm{d}\sigma
+\mathrm O(k^{-1})\lVert\ell_k\rVert_{C^0(P)}.
$$

而 $|\int_{\partial P}\ell_k\,\mathrm{d}\sigma|\leq C\lVert\ell_k\rVert_{C^0(P)}$，故 $B_k(\ell_k)\geq-C\lVert h\rVert_{k,2}$。由 $B_k$ 的线性性和 (23)，$B_k(h)\geq-C\lVert h\rVert_{k,2}$。结合 (9) 即得 (1)。

### 第二步：边界项的下极限

由 (9)，$\{f_k\}$ 在 $L^2(P)$ 中有界，因此一列子列弱收敛于 $f\in L^2(P)$。对非负 $\chi\in C_c^\infty(\operatorname{int}P)$ 和任意向量 $\xi$，凸性给出

$$
\int_Pf_k\,\partial_{\xi\xi}\chi\,\mathrm{d}x\geq0.
$$

弱收敛取极限后，$f$ 的分布 Hessian 非负；分布意义下凸性的标准刻画给出 $f$ 的凸代表元，仍记为 $f$。

对 $K\Subset K'\Subset\operatorname{int}P$，由 (10) 得到一致的 $C^0$ 和 Lipschitz 界。取穷尽 $K_1\Subset K_2\Subset\cdots\Subset\operatorname{int}P$，反复应用 Arzelà--Ascoli 定理并取对角子列，得到在内部局部一致收敛的子列。局部一致收敛蕴含分布收敛，与 $L^2$ 弱极限比较可知局部一致极限就是 $f$。

固定每个 facet 的 $K_F\Subset\operatorname{relint}F$。将 (11) 应用于非负凸函数列 $g_k$，得到

$$
B_k(g_k)\geq\sum_F\frac1{k^{n-1}}
\sum_{y\in K_F\cap P_k}g_k(y)
-\frac{C}{k}\bigl(\sup_Qg_k+\operatorname{Lip}_Q(g_k)\bigr).

$$

系数为 $1$ 是因为 (11) 的 facet 权重 $1/2$ 被 $B_k$ 中的因子 $2$ 抵消。

取 $p\in\operatorname{int}P$，选择 $q_k\in\partial f_k(p)$。若 $B(p,2r)\Subset P$，则

$$
|\langle q_k,e_i\rangle|
\leq r^{-1}\bigl(|f_k(p+re_i)|+|f_k(p-re_i)|+2|f_k(p)|\bigr),
$$

所以 $q_k$ 有界。再取子列，支撑仿射函数

$$
\ell_k(x)=f_k(p)+\langle q_k,x-p\rangle
$$

一致收敛于仿射函数 $\ell$；于是 $g_k=f_k-\ell_k\geq0$ 在内部局部一致收敛于 $g=f-\ell$。

对 facet $F$，primitive 性给出整向量 $v_F$，满足 $l_F(v_F)=1$。由于 $K_F$ 与其他 facets 分离，可取 $\delta>0$，使

$$
\{y+tv_F:y\in K_F,\ 0\leq t\leq2\delta\}
\Subset P\setminus\bigcup_{F'\neq F}F'.
$$

令 $m_k=\lfloor k\delta\rfloor$。当 $k$ 充分大时，$y+jv_F/k\in P_k$（$y\in K_F\cap P_k$，$0\leq j\leq2m_k$）。对这样的 $y$，沿线段的凸性给出

$$
g_k(y)\geq2g_k(y+m_kv_F/k)-g_k(y+2m_kv_F/k).

$$

右端位于内部紧集上。局部一致收敛、facet 上的 Riemann 和定理以及 (24) 因而给出

$$
\liminf_{k\to\infty}B_k(g_k)
\geq\sum_F\int_{K_F}
\bigl(2g(y+\delta v_F)-g(y+2\delta v_F)\bigr)\,\mathrm{d}\sigma_F.
$$

对 $\phi(t)=g(y+tv_F)$，割线斜率单调性给出

$$
2\bigl(\phi(t)-\phi(t/2)\bigr)\leq\phi(2t)-\phi(t),
$$

即 $2\phi(t/2)-\phi(t)\geq2\phi(t)-\phi(2t)$。因此当二进制 $\delta\downarrow0$ 时，$2g(y+\delta v_F)-g(y+2\delta v_F)$ 单调递增至下半连续边界值 $g(y)$。在固定 $K_F$ 上，这些函数的负部由内部紧集上的 $g(y+2\delta v_F)$ 控制；加上一个常数后可以使用单调收敛定理。先令二进制 $\delta\downarrow0$，再取 $K_F\nearrow\operatorname{relint}F$，再次使用单调收敛定理，得到

$$
\liminf_{k\to\infty}B_k(g_k)\geq\sum_F\int_Fg\,\mathrm{d}\sigma_F.
$$

仿射 Ehrhart 公式和 $\ell_k\to\ell$ 一致收敛给出 $B_k(\ell_k)\to\int_{\partial P}\ell\,\mathrm{d}\sigma$。由于 $f_k=g_k+\ell_k$，

$$
\liminf_{k\to\infty}B_k(f_k)
\geq\int_{\partial P}f\,\mathrm{d}\sigma.

$$

若右端为 $+\infty$，对 $K_F$ 取递增穷尽即可得到左端也是 $+\infty$。

### 第三步：二次项的下极限

固定 $\delta>0$，使 $P^\delta$ 有非空内部。局部一致收敛和通常的 Riemann 和定理给出

$$
\lim_{k\to\infty}\frac1{k^n}\sum_{u\in P_k\cap P^\delta}f_k(u)^2
=\int_{P^\delta}f^2\,\mathrm{d}x.
$$

省略的求和项均非负，所以

$$
\liminf_{k\to\infty}\lVert f_k\rVert_{k,2}^2
\geq\int_{P^\delta}f^2\,\mathrm{d}x.
$$

令 $\delta\downarrow0$ 并用单调收敛定理，得到

$$
\liminf_{k\to\infty}\lVert f_k\rVert_{k,2}^2
\geq\lVert f\rVert_{L^2(P)}^2.
$$

结合 (26)，即得 (2)。

### 第四步：没有 $L^2$ 边界缺陷

现在假设 (3)。对固定 $\delta>0$，内部 Riemann 和收敛给出

$$
\begin{aligned}
\limsup_{k\to\infty}\frac1{k^n}
\sum_{u\in P_k\setminus P^\delta}f_k(u)^2
&\leq\lVert f\rVert_{L^2(P)}^2-\int_{P^\delta}f^2\,\mathrm{d}x\\
&=\int_{P\setminus P^\delta}f^2\,\mathrm{d}x.
\end{aligned}

$$

对每个 $k$，将 (8) 应用于 $f_k$，写成

$$
f_k=g_k+\ell_k,\qquad g_k\geq0,
\qquad\lVert\ell_k\rVert_{C^0(P)}\leq C\lVert f_k\rVert_{k,2}\leq C.
$$

令

$$
E=P\setminus P^{2\delta},\qquad
U_{k,\delta}=\{u\in P_k:\operatorname{dist}(u,E)\leq C/k\}.
$$

当 $k$ 充分大时，$U_{k,\delta}\subset P_k\setminus P^\delta$。边界带 $E$ 的体积至多为 $C\delta$，并且

$$
|U_{k,\delta}|\leq C(k^n\delta+k^{n-1}).
$$

将局部采样估计 (6) 应用于 $g_k$，再用 $g_k^2\leq2f_k^2+2\ell_k^2$，得到

$$
\begin{aligned}
\int_Ef_k^2\,\mathrm{d}x
&\leq2\int_Eg_k^2\,\mathrm{d}x+2|E|\lVert\ell_k\rVert_{C^0(P)}^2\\
&\leq\frac{C}{k^n}\sum_{u\in U_{k,\delta}}f_k(u)^2
+C\bigl(k^{-n}|U_{k,\delta}|+|E|\bigr)\lVert\ell_k\rVert_{C^0(P)}^2.
\end{aligned}
$$

对固定 $\delta$，当 $k$ 充分大时，最后一个系数至多为 $C\delta$。结合 (27)，

$$
\limsup_{k\to\infty}\int_{P\setminus P^{2\delta}}f_k^2\,\mathrm{d}x
\leq C\int_{P\setminus P^\delta}f^2\,\mathrm{d}x+C\delta.
$$

在 $P^{2\delta}$ 上收敛一致，故在那里 $L^2$ 强收敛。因此

$$
\begin{aligned}
\limsup_{k\to\infty}\lVert f_k-f\rVert_{L^2(P)}^2
&\leq2\limsup_{k\to\infty}\int_{P\setminus P^{2\delta}}f_k^2\,\mathrm{d}x
+2\int_{P\setminus P^{2\delta}}f^2\,\mathrm{d}x\\
&\leq2C\int_{P\setminus P^\delta}f^2\,\mathrm{d}x+2C\delta
+2\int_{P\setminus P^{2\delta}}f^2\,\mathrm{d}x.
\end{aligned}
$$

由于 $f\in L^2(P)$，令 $\delta\downarrow0$，右端趋于零。因此 $f_k\to f$ 在 $L^2(P)$ 中强收敛，引理得证。$\square$
