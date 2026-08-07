# K-stability of Toric Varieties：文献综述与 relatively K-unstable polytope 例子

> 整理日期：2026-08-05
> 来源：arXiv、zbMATH、J-STAGE、Project Euclid、期刊官网、会议/讨论班公告等公开网络资料。
> 用途：为本项目（`K-stability/` 中 Wang–Zhou 多边形相对 K-不稳定性的精确计算，见 `generate_wang_zhou.py` 与 `paper/K-stability.tex`）提供文献背景和可直接引用的例子清单。

---

## 1. 背景与定义速览

设 `(X_Δ, L_Δ)` 是带 Delzant 多胞体 `P ⊂ ℝⁿ` 的极化 toric 流形（或 toric Fano 流形，取 `L = -K_X`）。K-stability 在 toric 情形有纯组合的表述：

- **Donaldson 的 polytope K-稳定性**：toric 检验构形对应 `P` 上的分片线性（PL）凸函数 `u`；Donaldson 定义了 `L`-泛函（约化的 Futaki 不变量）`L(u)`，并提出 polytope 稳定性与 cscK 度量存在性的对应（Donaldson 猜想）。相关参考文献：[Donaldson 2002](#donaldson-2002)、[Donaldson 2008](#donaldson-2008)。
- **相对 K-稳定性（toric sense）**：令 `L(u)` 为相对（约化）Futaki 不变量，则
  - `(X_Δ, L_Δ)` 称为**相对 K-半稳定**（toric sense），若对所有 PL 凸函数 `u` 都有 `L(u) ≥ 0`；
  - 称为**相对 K-稳定**（toric sense），若此外 `L(u) = 0` 当且仅当 `u` 是仿射线性函数；
  - `P`（或 `X_Δ`）称为**相对 K-不稳定**，若存在 PL 凸函数 `u` 使 `L(u) < 0`（在适当归一化下）。

  这是 [Yotsutani–Zhou 2019](#yotsutani-zhou-2019) 等文献中采用的表述。

- **Donaldson 的简单凸函数事实（n = 2）**：在二维情形，相对 K-不稳定当且仅当存在“简单凸函数”

  ```
  g = max{ax + by + c, 0}
  ```

  使得相对 DF 不变量 `< 0`。本项目 `K-stability/k_stability.cpp` 正是对这一判据做精确（`Gmpq`）计算。

- **极值度量与相对 K-稳定性**：[Stoppa–Székelyhidi 2011](#stoppa-szekelyhidi-2011) 证明：若极化流形在某 Kähler 类中容许极值度量（extremal Kähler metric），则它相对其自同构群的极大环面是 K-多稳定的（K-polystable）。因此**相对 K-不稳定 ⇒ 不存在极值度量**，这是 [Hwang–Sato–Yotsutani 2024](#hwang-sato-yotsutani-2024) 反例的逻辑主线。

- **Abreu 方程**：toric 流形上的极值度量方程可化为 Delzant 多胞体上的四阶实 PDE（Abreu 方程），见 [Abreu 1998](#abreu-1998) 与 [Guillemin 1994](#guillemin-1994)。

- **相对 Ding 稳定性与 Mabuchi 常数**：[Yao 2022](#yao-2022) 在 toric Fano 情形引入相对 Ding 稳定性；Mabuchi 常数可通过 moment polytope 的组合数据计算（Yao 的公式），[Nitta–Saito–Yotsutani 2023](#nitta-saito-yotsutani-2023) 用它完成了 ≤ 4 维 toric Fano 流形的 uniform relative Ding stability 分类。

- **Folklore 猜想**：长期以来人们猜测“每个 toric Fano 流形在 `c₁` 中都容许极值度量”（Mabuchi 提出的问题）。[Hwang–Sato–Yotsutani 2024](#hwang-sato-yotsutani-2024) 在维数 10 构造了反例，回答了 Mabuchi 的问题。

---

## 2. 奠基文献

### 2.1 <a id="abreu-1998"></a>Abreu（1998）
**M. Abreu, *Kähler geometry of toric varieties and extremal metrics*, Internat. J. Math. 9 (1998), no. 6, 641–651.**

建立了 toric Kähler 度量的势函数理论，把 Calabi 极值度量方程化为多胞体上的 Abreu 方程。是 toric 极值度量问题的出发点。

### 2.2 <a id="guillemin-1994"></a>Guillemin（1994）
**V. Guillemin, *Kähler structures on toric varieties*, J. Differential Geom. 40 (1994), no. 2, 285–309.**

给出了 toric 流形上标准 Kähler 势（Guillemin 势）的显式公式，Abreu 方程中边界条件的标准来源。

### 2.3 <a id="donaldson-2002"></a>Donaldson（2002）
**S. K. Donaldson, *Scalar curvature and stability of toric varieties*, J. Differential Geom. 62 (2002), no. 2, 289–349.** [zbMATH: Zbl 1074.53059](https://zbmath.org/1074.53059)

奠基性论文。定义极化簇的 K-稳定性并猜测与 cscK 度量存在性等价；在 toric 情形把 K-稳定性完全转化为 Delzant 多胞体上 PL 凸函数的 `L`-泛函非负性，即 **polytope K-stability**。本文中“简单凸函数”破坏不稳定的思想也源于此。

### 2.4 <a id="donaldson-2008"></a>Donaldson（2008）
**S. K. Donaldson, *Extremal metrics on toric surfaces: a continuity method*, J. Differential Geom. 79 (2008), no. 3, 389–432.** [zbMATH: Zbl 1151.53030](https://zbmath.org/1151.53030)

用连续性方法结合“blow-up”论证研究 Abreu 方程的可解性，在某种解析“M-条件”下给出 toric 曲面极值度量的存在性。是 [Wang–Zhou 2011](#wang-zhou-2011) 之前最重要的存在性工作之一。

### 2.5 Ross–Thomas（2007）
**J. Ross, R. Thomas, *A study of the Hilbert–Mumford criterion for the stability of projective varieties*, J. Algebraic Geom. 16 (2007), no. 2, 201–255.**

系统发展检验构形与 Chow 权（Chow weight）理论；toric 情形的相对 Chow 稳定性约化直接用到其方法（见 [Yotsutani–Zhou 2019](#yotsutani-zhou-2019)）。

### 2.6 <a id="szekelyhidi-thesis"></a>Székelyhidi（2006，博士论文）
**G. Székelyhidi, *Extremal metrics and K-stability*, Ph.D. thesis, Imperial College London (2006).**

引入相对 K-稳定性的现代表述（对自同构群的极大环面），为极值度量与稳定性对应奠定框架。

### 2.7 <a id="stoppa-szekelyhidi-2011"></a>Stoppa–Székelyhidi（2011）
**J. Stoppa, G. Székelyhidi, *Relative K-stability of extremal metrics*, J. Eur. Math. Soc. (JEMS) 13 (2011), no. 4, 899–909.** [arXiv:0912.4095](https://arxiv.org/abs/0912.4095)

证明：若极化流形容许极值度量，则它相对于极大自同构环面是 K-多稳定的。这是“相对 K-不稳定 ⇒ 无极值度量”这一常用推断的依据。

---

## 3. 相对 K-稳定性与 toric 极值度量

### 3.1 <a id="zhou-zhu-2006"></a>Zhou–Zhu（2006）
**Bin Zhou, Xiaohua Zhu, *Relative K-stability and modified K-energy on toric manifolds*, arXiv:math/0603237.** [arXiv](https://arxiv.org/abs/math/0603237)

研究 toric 流形上与 Calabi 极值度量相伴的相对 K-稳定性与修正 K-能量。给出纯多胞体意义的充分条件，使两者同时成立；特别适用于 Futaki 不变量为零的 toric Fano 流形，并在 toric Fano 曲面（toric del Pezzo）上验证。

### 3.2 <a id="wang-zhou-2011"></a>Wang–Zhou（2011）★ 用户指定论文
**Xu-Jia Wang, Bin Zhou, *On the existence and nonexistence of extremal metrics on toric Kähler surfaces*, Adv. Math. 226 (2011), no. 5, 4429–4455.** [zbMATH: Zbl 1220.53088](https://zbmath.org/1220.53088) · [DOI: 10.1016/j.aim.2010.12.008](https://doi.org/10.1016/j.aim.2010.12.008) · [ANU 开放存档](https://openresearch-repository.anu.edu.au/items/f9111a67-3436-4a0f-82d7-b1de95e72da1/)

**主要内容：**
- 证明**每个 toric Kähler 曲面都存在某个 Kähler 类，其中含 Calabi 极值度量**（用 Arezzo–Pacard–Singer 的 blow-up 方法等）。
- 同时给出**不存在性**的例子：找到一个具有 **9 个 `T_ℂ²` 不动点** 的 toric Kähler 曲面，其某 Kähler 类是 **K-不稳定的**，因而该类中**没有极值度量**。
- 证明 toric 曲面的 K-稳定性可由**简单分片线性函数**（simple PL functions，即形如 `max{ax+by+c, 0}`）刻画。

> 本项目 `generate_wang_zhou.py` 中生成的 9 顶点 Wang–Zhou 多边形
> `(0,a+20), (1,a+20), (2,a+19), (3,a+17), (4,a+14), (5,a+10), (7,a), (7,-a), (0,-a)`
> 正是该论文所研究的（相对）K-不稳定多胞体族。

### 3.3 <a id="wang-zhou-2014"></a>Wang–Zhou（2014）
**Xu-Jia Wang, Bin Zhou, *K-stability and canonical metrics on toric manifolds*, Bull. Inst. Math. Acad. Sinica (N.S.) 9 (2014), no. 1, 85–110.** [论文页](https://www.math.sinica.edu.tw/bulletin/journals/704) · [zbMATH: Zbl 1290.32027](https://zbmath.org/1290.32027)

讨论 toric Kähler 流形上的 K-稳定性，并给出一个 **8 个 `T_ℂ²` 不动点** 的 toric Kähler 曲面不稳定例子（另一类不稳定多胞体例子）。

### 3.4 <a id="chen-li-sheng-2014"></a>Chen–Li–Sheng（2014）
**Bingyan Chen, An-Min Li, Li Sheng, *Uniform K-stability for extremal metrics on toric varieties*, J. Differential Equations 257 (2014), no. 5, 1487–1500.** [arXiv:1109.5228](https://arxiv.org/abs/1109.5228) · [zbMATH: Zbl 1300.53067](https://zbmath.org/1300.53067)

证明 toric 情形 Abreu 方程解的存在性蕴含 **uniform K-稳定性**（以及相对强 K-稳定性），并证明 Mabuchi 泛函的 properness 是极值度量存在的必要条件。

### 3.5 Chen–Li–Sheng（2018）
**Bingyan Chen, An-Min Li, Li Sheng, *Extremal metrics on toric surfaces*, Adv. Math. 340 (2018), 363–405.** [arXiv:1008.2607](https://arxiv.org/abs/1008.2607) · [zbMATH: Zbl 1411.53054](https://zbmath.org/1411.53054)

研究 Delzant 多胞体上带指定标量曲率的 Abreu 方程；作为应用，零 Futaki 不变量的 K-稳定 toric 曲面容许 cscK 度量（Donaldson 猜想在曲面情形的进展）。

### 3.6 <a id="szekelyhidi-2008"></a>Székelyhidi（2008）
**Gábor Székelyhidi, *Optimal test-configurations for toric varieties*, J. Differential Geom. 80 (2008), no. 3, 501–523.** [arXiv:0709.2687](https://arxiv.org/abs/0709.2687)

对 K-不稳定的 toric 簇证明存在**最优去稳定凸函数**（optimal destabilising convex function）；若该函数为 PL，则给出类似 Harder–Narasimhan 滤过的“半稳定分解”。还证明 Calabi 泛函下确界等于所有去稳定检验构形的归一化 Futaki 不变量的上确界。

### 3.7 Hisamoto（2016）
**Tomoyuki Hisamoto, *Stability and coercivity for toric polarizations*, arXiv:1610.07998.** [arXiv](https://arxiv.org/abs/1610.07998)

证明 toric 极化流形 **uniform K-稳定（toric sense）⇔ K-能量泛函模极大环面作用 coercive**。把分析条件与组合稳定性完全对应起来。

### 3.8 Apostolov–Calderbank–Gauduchon–Tønnesen-Friedman（Ambitoric geometry II）
**V. Apostolov, D. M. J. Calderbank, P. Gauduchon, C. W. Tønnesen-Friedman, *Ambitoric geometry II: Extremal toric surfaces and Einstein 4-orbifolds*, arXiv:1302.6979.** [arXiv](https://arxiv.org/abs/1302.6979)

对第二 Betti 数 `b₂(M) = 2` 的 toric 4-orbifold 完全解决极值度量存在性问题：存在极值度量 **iff** 有理 Delzant 多胞体在 toric 意义下是**相对 K-多稳定**的（Székelyhidi 意义）。是多胞体相对稳定性的一个“等价性”实例。

### 3.9 Delcroix–Jubert（2023）
**Thibaut Delcroix, Simon Jubert, *An effective weighted K-stability condition for polytopes and semisimple principal toric fibrations*, Annales Henri Lebesgue 6 (2023), 117–?; arXiv:2202.02996.** [arXiv](https://arxiv.org/abs/2202.02996) · [AHL](https://ahl.centre-mersenne.org/item/AHL_2023__6__117_0/) · [zbMATH: Zbl 1531.14063](https://zbmath.org/7697374)

给出带标记多胞体（labelled polytopes）上可有效检验的 **weighted uniform K-stability** 充分条件，并由此得到半单主 toric 纤维化上极值度量的低维新例子。属于“多胞体稳定性判据”方向的现代发展。

### 3.10 Fujita（2020）
**Kento Fujita, *Notes on K-semistability of toric polarized varieties*, Springer INdAM Series（2020）.** [zbMATH: Zbl 1439.14131](https://zbmath.org/1439.14131)

对 toric 极化簇从环面不变素除子系统构造“basic blowup type”检验构形，给出 toric K-半稳定性的组合判据表述。

---

## 4. 相对 algebro-geometric stabilities（★ 用户指定论文）

### 4.1 <a id="yotsutani-zhou-2019"></a>Yotsutani–Zhou（2019 / arXiv 2023 修订版）
**Naoto Yotsutani, Bin Zhou, *Relative algebro-geometric stabilities of toric manifolds*, Tohoku Math. J. (2) 71 (2019), no. 4, 495–524.** [arXiv:1602.08201](https://arxiv.org/abs/1602.08201)（v3, 2023） · [J-STAGE](https://www.jstage.jst.go.jp/article/tmj/71/4/71_1576724790/_article/-char/en) · [DOI: 10.2748/tmj/1576724790](https://doi.org/10.2748/tmj/1576724790) · [zbMATH: Zbl 1451.53099](https://zbmath.org/1451.53099)

**主要内容：**
- 给出 toric Fano 流形在 toric 意义下**相对 K-稳定与不稳定的判据**（基于 moment polytope 上的约化 Futaki / `L`-泛函与归一化势函数）。
- 用 Hilbert–Mumford 判据两种方式约化相对 Chow 稳定性：
  1. 极大环面作用与权多胞体（Ono 的策略，贴合 Székelyhidi 的相对 GIT 稳定性）；
  2. `ℂ*`-作用与 toric 退化相伴的 Chow 权（Donaldson、Ross–Thomas）。
- 应用：**部分确定全部 toric Fano 三维流形的相对 K-稳定性**（论文表 6），并给出**相对 K-稳定（toric sense）但渐近相对 Chow 不稳定**的反例（基于 [Nill–Paffenholz](#nill-paffenholz-2011) 的例子）。

> ⚠️ **重要勘误提示**：arXiv v3（2023-05-16）附录指出，**已发表版本**（Tohoku Math. J. 2019，对应正文第 1–5 节）存在错误，其中表 6 关于相对 K-稳定性的结果不可靠（“inconclusive results”）。引用其三维分类结论时务必使用 arXiv v3 或后续文献（如 [Nitta–Saito 2023](#nitta-saito-2023)）的修订意见。

### 4.2 <a id="nitta-saito-2023"></a>Nitta–Saito（2023）
**Yasufumi Nitta, Shunsuke Saito, *A note on the Yotsutani–Zhou condition for relative K-instability*, Kodai Math. J. 46 (2023), no. 2, 219–227.** [J-STAGE](https://www.jstage.jst.go.jp/article/kodaimath/46/2/46_219/_article/-char/en) · [Project Euclid](https://projecteuclid.org/journals/kodai-mathematical-journal/volume-46/issue-2/A-note-on-the-Yotsutani-Zhou-condition-for-relative-K/10.2996/kmj46205.short) · [zbMATH: Zbl 1530.14090](https://zbmath.org/1530.14090)

Yotsutani–Zhou 曾给出 toric Fano 流形**相对 K-不稳定**的充分条件；本文给出应用该条件时的**简单障碍**（obstruction）。即：并非所有看起来“形如 YZ 判据要求”的多胞体都能用 YZ 条件判定，需先检查该障碍。这是使用 YZ 判据时必须注意的细节。

### 4.3 Ono（2010）
**Hajime Ono, *Algebro-geometric semistability of polarized toric manifolds*, arXiv:1009.0087.** [arXiv](https://arxiv.org/abs/1009.0087)

给出积分 Delzant 多胞体对应的 toric 流形在极大环面作用下 **Chow 半稳定**的充要条件；不借助 Riemann–Roch 与检验构形即证明：渐近（相对）Chow 半稳定 ⇒（相对）K-半稳定（toric 退化意义下，Ross–Thomas 的结论）。

### 4.4 <a id="yotsutani-2023"></a>Yotsutani（2023）
**Naoto Yotsutani, *Asymptotic Chow semistability implies Ding polystability for Gorenstein toric Fano varieties*, Mathematics 11 (2023), no. 19, 4114.** [arXiv:1711.10113](https://arxiv.org/abs/1711.10113) · [MDPI](https://www.mdpi.com/2227-7390/11/19/4114)

- 定理 1.3：Gorenstein toric Fano 簇若渐近 Chow 半稳定，则关于 toric 检验构形 **Ding 多稳定**（推广了已知的定理 1.2 至 Gorenstein 奇异情形）。
- 命题 1.4：利用 Mabuchi 常数的**可加性**（基于 [Ono–Sano–Yotsutani 2023](#ono-sano-yotsutani-2023)），系统构造**无穷多个**例子说明**相对 K-稳定性与相对 Ding 稳定性不同**。
- 用 [Yotsutani–Zhou](#yotsutani-zhou-2019) 的组合判据验证 Gorenstein toric del Pezzo 曲面的相对 Chow 稳定性：16 个同构类中（据 arXiv 早期版本命题 5.4）5 个渐近相对 Chow 多稳定、4 个渐近相对 Chow 不稳定。

### 4.5 <a id="ono-sano-yotsutani-2023"></a>Ono–Sano–Yotsutani（2023）
**Hajime Ono, Yuji Sano, Naoto Yotsutani, *Bott manifolds with vanishing Futaki invariants for all Kähler classes*, arXiv:2305.05924.** [arXiv](https://arxiv.org/abs/2305.05924)

证明：对所有 Kähler 类 Futaki 不变量都为零的 Bott 流形恰为射影直线的乘积；其公式支撑了 Mabuchi 常数的可加性与上述无穷多相对 K-/Ding 稳定性分离例子。

### 4.6 Lee–Li–Sturm–Wang（2019）
**King Leung Lee, Zhiyuan Li, Jacob Sturm, Xiaowei Wang, *Asymptotic Chow stability of toric del Pezzo surfaces*, Math. Res. Lett. 26 (2019), no. 6, 1769–?; arXiv:1711.10099.** [arXiv](https://arxiv.org/abs/1711.10099) · [zbMATH: Zbl 1467.14110](https://zbmath.org/1467.14110)

研究 Odaka–Spotti–Sun 的 Kähler–Einstein Fano 簇模空间中出现的 toric del Pezzo 曲面的有效 Chow 稳定性（toric 多胞体的 GIT 计算）。

### 4.7 <a id="nill-paffenholz-2011"></a>Nill–Paffenholz（2011）
**Benjamin Nill, Andreas Paffenholz, *Examples of non-symmetric Kähler–Einstein toric Fano manifolds*, Beitr. Algebra Geom. 52 (2011), no. 2, 297–304.** [arXiv:0905.2054](https://arxiv.org/abs/0905.2054)

构造 7 维与 8 维的**非对称** toric Fano 流形，仍容许 Kähler–Einstein 度量（回答 Batyrev–Selivanova 的问题）。该例子被 [Yotsutani–Zhou](#yotsutani-zhou-2019) 用作“相对 K-稳定但渐近相对 Chow 不稳定”的反例来源。

---

## 5. 相对 Ding 稳定性

### 5.1 <a id="yao-2022"></a>Yao（2022）
**Yi Yao, *Mabuchi solitons and relative Ding stability of toric Fano varieties*, Int. Math. Res. Not. IMRN 2022, no. 24, 19790–?; arXiv:1701.04016.** [arXiv](https://arxiv.org/abs/1701.04016) · [IMRN](https://academic.oup.com/imrn/article-abstract/2022/24/19790/6374223) · [zbMATH: Zbl 1510.53084](https://zbmath.org/1510.53084)

作为 Futaki 不变量非零 Fano 流形上 Kähler–Einstein 度量的推广（Mabuchi 孤子），引入 toric Fano 情形的**相对 Ding 稳定性**。在不稳定情形，确定**最大去稳定子**——moment polytope 上的一个简单凸函数，并建立连接 Ding 能量与 Ding 不变量（Berman–Ding 不变量）的 **Moment-Weight 等式**。NSY 计算 Mabuchi 常数即基于本文。

### 5.2 <a id="nitta-saito-yotsutani-2023"></a>Nitta–Saito–Yotsutani（2023）
**Yasufumi Nitta, Shunsuke Saito, Naoto Yotsutani, *Relative Ding and K-stability of toric Fano manifolds in low dimensions*, European J. Math. 9 (2023), no. 2, art. 44; arXiv:1712.01131.** [arXiv](https://arxiv.org/abs/1712.01131) · [zbMATH: Zbl 1517.32072](https://zbmath.org/1517.32072)

- 主要定理（定理 1.1）：完整列出所有 **≤ 4 维 toric Fano 流形**的 **uniform relative Ding stability**（稳定与不稳定），关键量是可由 moment polytope 组合数据计算的 **Mabuchi 常数**（Yao 的公式）。表 1–3 给出全部列表。
- 推论 1.6 与 1.9：利用 Bott tower 结构，澄清**相对 K-稳定性与相对 Ding 稳定性之间的差别**（存在相对 K-多稳定但相对 Ding 不稳定的 toric Fano 流形）。
- 具体数据点（据 arXiv 表格）：例如某编号 38、记作 `G_6` 的 4 维 toric Fano 例子被判为 relative Ding unstable，对应 Mabuchi 常数 `6431616/4388521`（具体编号以原文表为准）。

---

## 6. 近期进展

### 6.1 <a id="hwang-sato-yotsutani-2024"></a>Hwang–Sato–Yotsutani（2024）★ 最重要的 relatively K-unstable 例子
**DongSeon Hwang, Hiroshi Sato, Naoto Yotsutani, *A toric Fano manifold that does not admit an extremal Kähler metric*, arXiv:2411.17574.** [arXiv](https://arxiv.org/abs/2411.17574) · [Semantic Scholar](https://www.semanticscholar.org/paper/Toric-Fano-manifolds-that-do-not-admit-extremal-Hwang-Sato/10871e7ddee55dc09a45ede38af9eea273263d93) · [讨论班公告（BIMSA/PKU）](https://www.bimsa.cn/talk/30422.html)

- 主要结果：存在 **10 维 toric Fano 流形，其在 toric 意义下相对 K-不稳定**；由 [Stoppa–Székelyhidi](#stoppa-szekelyhidi-2011)，它在 `c₁` 中**不含极值度量**，回答 Mabuchi 的问题。
- 方法：对该流形的 moment polytope **显式构造去稳定凸函数**，从而验证相对 K-不稳定性（使用相对 K-不稳定判据）。
- 高维推广：讨论表明通过适当取乘积可构造 **所有更高维数（n ≥ 11）** 的类似例子（见 2025 KMS 年会摘要）。

### 6.2 Lee–Yotsutani（2024）
**King Leung Lee, Naoto Yotsutani, *Chow stability of λ-stable toric varieties*, arXiv:2405.06883.** [arXiv](https://arxiv.org/abs/2405.06883)

对极化 toric 簇定义 **λ-稳定性**（uniform K-稳定性的自然推广，源于 [Yotsutani–Zhou](#yotsutani-zhou-2019) 命题 5.1.2 的思想）；作为应用，证明 Futaki–Ono 不变量为零的 K-半稳定光滑 toric 极化簇渐近 Chow 多稳定。

### 6.3 综述文章
- **Bin Zhou, *Extremal metrics on toric manifolds——existence and K-stability*（极值度量与 K-稳定性）, Scientia Sinica Mathematica 44 (2014), no. 1.** [sciengine: 10.1360/012013-144](https://www.sciengine.com/SSM/doi/10.1360/012013-144)。中文综述：toric 曲面每个都含极值度量的 Kähler 类（Arezzo–Pacard–Singer 方法），同时给出大量不稳定 Kähler 类反例。
- **An-Min Li, Li Sheng, *Extremal Kähler metrics of toric manifolds*, Chinese Ann. Math. Ser. B 44 (2023).** 近年 toric 极值度量与 K-稳定性进展的综述。
- **Harold Blum, Yuchen Liu, Chenyang Xu, Ziquan Zhuang, *K-stability of Fano varieties: an algebro-geometric approach*, EMS Surv. Math. Sci. 9 (2022), 411–464; arXiv:1911.02461.** 一般 Fano 情形 K-稳定性的现代综述，toric 是其中最重要的具体类之一。（一般背景参考）

---

## 7. Relatively K-unstable polytope 例子清单（核心）

### 7.1 判定框架（toric sense）
对带 Delzant 多胞体 `P` 的（Fano）toric 流形，相对 K-不稳定性归结为：**存在 PL 凸函数 `u` 使约化 Futaki / 相对 `L`-泛函 `L(u) < 0`**。在 `n = 2` 时只需检验简单凸函数 `g = max{ax+by+c, 0}`（Donaldson）。

下面按文献整理已知的（相对）不稳定多胞体例子与相关反例。

### 7.2 例 1：Wang–Zhou 9 顶点多边形（曲面情形，经典例子）
来源：[Wang–Zhou 2011](#wang-zhou-2011)。

存在具有 **9 个 `T_ℂ²` 不动点** 的 toric Kähler 曲面，其某个 Kähler 类 K-不稳定，类内无极值度量。其 moment polytope 即本仓库 `generate_wang_zhou.py` 生成的多边形族：

```
(0, a+20), (1, a+20), (2, a+19), (3, a+17), (4, a+14), (5, a+10), (7, a), (7, -a), (0, -a)
```

这是**本项目数值/精确计算验证的对象**（见 `K-stability/README.md`：计算 `ℓ_P` 并对简单凸函数检验相对 DF 不变量）。Wang–Zhou 同时证明 K-稳定性可由简单 PL 函数刻画，故该族是检验计算程序的天然 benchmark。

另见 [Wang–Zhou 2014](#wang-zhou-2014) 的 **8 不动点** 不稳定例子（另一族）。

### 7.3 例 2：Yotsutani–Zhou 判据与 toric Fano 三维流形（表 6）
来源：[Yotsutani–Zhou 2019](#yotsutani-zhou-2019)。

YZ 给出 toric Fano 流形相对 K-稳定/不稳定的判据，并据以（部分）确定 toric Fano 三维流形的相对 K-稳定性（表 6）。

⚠️ 注意：已发表版本表 6 的若干结论已被作者在 arXiv v3 附录中声明不可靠（第 1–5 节有错误）；且 [Nitta–Saito 2023](#nitta-saito-2023) 指出 YZ 不稳定判据存在适用障碍。因此**引用三维分类表时应以 arXiv v3 及后续文献为准**。

### 7.4 例 3：相对 K-稳定但渐近相对 Chow 不稳定（反向反例）
来源：[Yotsutani–Zhou 2019](#yotsutani-zhou-2019) 第 5.2 节，基于 [Nill–Paffenholz](#nill-paffenholz-2011)。

[Nill–Paffenholz](#nill-paffenholz-2011) 的 7、8 维非对称 toric Fano 流形（容许 Kähler–Einstein 度量）被 YZ 用作例子：其 polytope **相对 K-稳定（toric sense）但渐近相对 Chow 不稳定**。这提醒我们：相对 K-稳定性并不自动蕴含渐近相对 Chow 稳定性，两条稳定线之间有真实差别。

### 7.5 例 4：10 维 toric Fano 流形（相对 K-不稳定，无极值度量）★
来源：[Hwang–Sato–Yotsutani 2024](#hwang-sato-yotsutani-2024)。

目前最直接、最醒目的 **relatively K-unstable polytope** 例子：10 维 toric Fano 流形 `X`，其 moment polytope 上存在显式构造的**去稳定凸函数**，使 `X` 在 toric 意义下相对 K-不稳定；由 [Stoppa–Székelyhidi](#stoppa-szekelyhidi-2011)，`X` 在 `c₁` 中不容许极值度量（回答 Mabuchi 的问题）。通过取乘积可得所有 `n ≥ 11` 的类似例子。

### 7.6 例 5：低维相对 Ding 不稳定例子（≤ 4 维，含 Bott tower）
来源：[Nitta–Saito–Yotsutani 2023](#nitta-saito-yotsutani-2023)、[Yotsutani 2023](#yotsutani-2023)。

- NSY 表 1–3 列出全部 ≤ 4 维 toric Fano 的 uniform relative Ding stability 与 Mabuchi 常数；其中有相对 Ding 不稳定者（如 4 维例子 `G_6`，Mabuchi 常数 `6431616/4388521`）。
- 推论 1.9：相对 Ding 不稳定的 toric Fano 流形可由 Bott tower 结构给出；18 类低阶（stage ≤ 4）Fano Bott 流形中恰有一个 twist-one 者相对 Ding 不稳定。
- Yotsutani 2023 用 Mabuchi 常数的可加性构造**无穷多**相对 K-多稳定但相对 Ding 不稳定的 toric Fano 流形（命题 1.4/1.5）。

### 7.7 例 6：相对 Ding 不稳定情形下的“最大去稳定子”
来源：[Yao 2022](#yao-2022)。

对相对 Ding 不稳定的 toric Fano 簇，Yao 确定其**最大去稳定子**——moment polytope 上的一个简单凸函数，并建立 Moment-Weight 等式（Calabi 型能量下确界 = Berman–Ding 不变量）。这说明“不稳定 polytope 的显式去稳定函数”可以非常具体，甚至可取为简单凸函数。

### 7.8 例 7：K-不稳定 toric 簇的最优去稳定函数与 HN 分解
来源：[Székelyhidi 2008](#szekelyhidi-2008)。

一般性定理：对任何 K-不稳定的 toric 簇，存在**最优去稳定凸函数**；若它是 PL 的，则诱导出半稳定块的 Harder–Narasimhan 型分解。这给出了“不稳定多胞体 → 去稳定凸函数 → 分解”的统一图景，也暗示不稳定多胞体通常由简单的 PL 函数即可检测。

---

## 8. 汇总表

| 文献 | 类型 | 是否给出 relatively K-unstable / 不稳定多胞体例子 | 关键点 |
| --- | --- | --- | --- |
| Donaldson 2002 | 奠基 | 理论框架 | polytope K-稳定性、`L`-泛函 |
| Wang–Zhou 2011 | 曲面存在性/不存在性 | ✅ 9 不动点不稳定多胞体 | 简单 PL 函数刻画；本项目对象 |
| Wang–Zhou 2014 | 曲面 K-稳定性 | ✅ 8 不动点不稳定多胞体 | toric Kähler 流形 K-稳定性 |
| Zhou–Zhu 2006 | 相对 K-稳定性 | 充分条件 | 修正 K-能量 properness |
| Chen–Li–Sheng 2014/2018 | 存在性与稳定性 | 必要条件方向 | Abreu 方程解 ⇒ uniform K-稳定 |
| Székelyhidi 2008 | 理论 | ✅ 一般不稳定情形 | 最优去稳定函数、HN 分解 |
| Stoppa–Székelyhidi 2011 | 理论 | 推断工具 | 极值度量 ⇒ 相对 K-多稳定 |
| Yotsutani–Zhou 2019 | 相对代数几何稳定性 | ✅ 判据 + Fano 3 维表（注意 errata） | 相对 K-/Chow 稳定性约化 |
| Nitta–Saito 2023 | 判据注记 | ⚠️ 使用障碍 | YZ 不稳定判据的适用性 |
| Ono 2010 | Chow 半稳定 | 理论 | Chow ⇒ K-半稳定（toric 退化） |
| Yotsutani 2023 | Gorenstein Fano | ✅ del Pezzo 相对 Chow 不稳定类 | 相对 K-与相对 Ding 分离的无穷族 |
| Ono–Sano–Yotsutani 2023 | Bott 流形 | 辅助结果 | Mabuchi 常数可加性 |
| Nill–Paffenholz 2011 | 反例来源 | ⚠️ 反向例子 | KE 存在但相对 Chow 不稳定 |
| Yao 2022 | 相对 Ding | ✅ 不稳定情形最大去稳定子 | Moment-Weight 等式 |
| Nitta–Saito–Yotsutani 2023 | 低维分类 | ✅ ≤ 4 维 Ding 不稳定表 | Mabuchi 常数分类 |
| Hwang–Sato–Yotsutani 2024 | 反例 | ✅ 10 维（及 n ≥ 11） | 显式去稳定凸函数；回答 Mabuchi 问题 |
| Lee–Yotsutani 2024 | λ-稳定性 | 理论 | λ-稳定 ⇒ 渐近 Chow 多稳定 |
| Delcroix–Jubert 2023 | weighted K-稳定性 | 有效判据 | labelled polytopes |
| Apostolov 等（Ambitoric II） | 4-orbifold 等价性 | ✅ 相对 K-多稳定 ⇔ 极值度量 | `b₂ = 2` toric 情形完全解决 |

---

## 9. 快速索引（按 arXiv 编号）

| arXiv | 标题 | 作者 |
| --- | --- | --- |
| [math/0603237](https://arxiv.org/abs/math/0603237) | Relative K-stability and modified K-energy on toric manifolds | Zhou, Zhu |
| [0709.2687](https://arxiv.org/abs/0709.2687) | Optimal test-configurations for toric varieties | Székelyhidi |
| [0912.4095](https://arxiv.org/abs/0912.4095) | Relative K-stability of extremal metrics | Stoppa, Székelyhidi |
| [0905.2054](https://arxiv.org/abs/0905.2054) | Examples of non-symmetric Kähler–Einstein toric Fano manifolds | Nill, Paffenholz |
| [1008.2607](https://arxiv.org/abs/1008.2607) | Extremal metrics on toric surfaces | Chen, Li, Sheng |
| [1009.0087](https://arxiv.org/abs/1009.0087) | Algebro-geometric semistability of polarized toric manifolds | Ono |
| [1109.5228](https://arxiv.org/abs/1109.5228) | Uniform K-stability for extremal metrics on toric varieties | Chen, Li, Sheng |
| [1302.6979](https://arxiv.org/abs/1302.6979) | Ambitoric geometry II: Extremal toric surfaces and Einstein 4-orbifolds | Apostolov, Calderbank, Gauduchon, Tønnesen-Friedman |
| [1602.08201](https://arxiv.org/abs/1602.08201) | Relative algebro-geometric stabilities of toric manifolds | Yotsutani, Zhou |
| [1610.07998](https://arxiv.org/abs/1610.07998) | Stability and coercivity for toric polarizations | Hisamoto |
| [1701.04016](https://arxiv.org/abs/1701.04016) | Mabuchi solitons and relative Ding stability of toric Fano varieties | Yao |
| [1711.10099](https://arxiv.org/abs/1711.10099) | Asymptotic Chow stability of toric del Pezzo surfaces | Lee, Li, Sturm, Wang |
| [1711.10113](https://arxiv.org/abs/1711.10113) | Asymptotic Chow semistability implies Ding polystability for Gorenstein toric Fano varieties | Yotsutani |
| [1712.01131](https://arxiv.org/abs/1712.01131) | Relative Ding and K-stability of toric Fano manifolds in low dimensions | Nitta, Saito, Yotsutani |
| [2202.02996](https://arxiv.org/abs/2202.02996) | An effective weighted K-stability condition for polytopes and semisimple principal toric fibrations | Delcroix, Jubert |
| [2305.05924](https://arxiv.org/abs/2305.05924) | Bott manifolds with vanishing Futaki invariants for all Kähler classes | Ono, Sano, Yotsutani |
| [2405.06883](https://arxiv.org/abs/2405.06883) | Chow stability of λ-stable toric varieties | Lee, Yotsutani |
| [2411.17574](https://arxiv.org/abs/2411.17574) | A toric Fano manifold that does not admit an extremal Kähler metric | Hwang, Sato, Yotsutani |

---

## 10. 结论与对本项目的建议

1. **最直接的 relatively K-unstable polytope 例子**按文献时间线是：Wang–Zhou（9/8 不动点曲面，2011/2014）→ Yotsutani–Zhou 判据（Fano 3 维，2019，需注意 errata）→ Hwang–Sato–Yotsutani（10 维 Fano，2024）。
2. 本仓库计算的 Wang–Zhou 9 顶点多边形族正是上述第一个例子的具体实现；在论文引用时建议同时引用 Wang–Zhou 2011（例子的来源）与 Wang–Zhou 2014（简单 PL 函数刻画），并可用 Székelyhidi 2008 / Yao 2022 的“去稳定函数”语言补充理论背景。
3. 若论文涉及 YZ 判据或表 6，务必引用 arXiv:1602.08201 v3 并说明已发表版本存在错误（作者附录声明）；同时可引用 Nitta–Saito 2023 说明判据适用性的限制。
4. 若需要“相对 K-稳定 ⇔ 极值度量”的正面等价性声明，可引用 Apostolov 等（Ambitoric II，`b₂ = 2` toric 4-orbifold）或 Chen–Li–Sheng 2014/2018（必要条件/曲面充分条件）。
