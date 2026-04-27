# PRD: Implement and optimize Kazhdan-Lusztig core in pyw

## Goal

**Embed** the functionality of `refs/Kazhdan-Lusztig/Algebra.py` into the `pyw` project. **Optimize** performance by caching and other methods.

Reuse the existing `pyw` infrastructure whenever possible.

Resort to `sagemath` framework when no tool is available in `pyw`


## CRITICAL

Read `refs/Kazhdan-Lusztig/Algebra.py` line by line CAREFULLY. Extract the structure of that code. Extrac the workflow of that code.
- map out which line of code can be implemented using existing tools in `pyw`
- map out which line of code can be implemented using `sagemath`
- map out which line of code needs to be implemented from scratch (which is unlikely since that code was implemented in `sagemath` in the first place)

A schematic workflow is listed later in this document.



## API

- `kazhdan_lusztig.py` 下 `KazhdanLusztigPolynomials` 负责 KL 核心计算和缓存机制
  - `.P(x, y)`: 计算 KL 多项式
  - `.Q(x, y)`: 计算 KL 逆多项式 of $\widehat W_{\Lambda}$, compute from coxeter3 by default
  - `.Qtilde(x_coset, y_coset)`: 计算 Qtildes of $\widehat W_{\Lambda}/\widehat W_\Lambda^0$
- `character.py` 下 `KazhdanLusztigCharacter` 负责组装特征标
  - `KazhdanLusztigCharacter.character(lambda)`: 组装 $L_\lambda$ 的特征标
  $$
  \operatorname{ch}(L_\lambda)

  = \sum_{\substack{[\widehat w'] \in \widehat W_\Lambda/ \widehat W_\Lambda^0\\ [\widehat w] \le [\widehat w']}}

  \widetilde Q_{[\widehat w], [\widehat w']}(1)

  \frac{e^{\widehat w' \cdot \widehat \Lambda} }{
  \sum_{\widehat w \in \widehat W} (-1)^{\widehat w} e^{\widehat w(\widehat \rho) - \widehat \rho} 
  } \\
  = \frac{1}{\sum_{\widehat w \in \widehat W} (-1)^{\widehat w} e^{\widehat w(\widehat \rho) - \widehat \rho}} \sum_{\substack{[\widehat w'] \in \widehat W_\Lambda/ \widehat W_\Lambda^0\\ [\widehat w] \le [\widehat w']}}

  \widetilde Q_{[\widehat w], [\widehat w']}(1)
  e^{\widehat w' \cdot \widehat \Lambda}
  $$


## Current Progress (2026-04-20)

- 已逐段拆解 `refs/Kazhdan-Lusztig/Algebra.py`，并把其职责分成：初始化/权转换/affine 枚举/缓存/Bruhat 与 coset 辅助/`Q`/`Qtilde`/`\Lambda` 与稳定子群/character 组装。
- 已对照 `pyw` 现有基础设施，确认可复用的核心在 `pyw/core/bruhat.py`、`pyw/core/kazhdan_lusztig.py`、`pyw/core/weyl_group.py`、`pyw/core/affine_lie_algebra.py`、`pyw/core/affine_weight.py`。
- 已完成一轮关键语义校正：
  - `KazhdanLusztigPolynomials.Q(x, y)` 现在表示 **Weyl 群元素级** inverse KL；
  - `KazhdanLusztigPolynomials.Q_tilde(x_coset, y_coset)` 现在表示 **陪集级 / parabolic** inverse KL；
  - `parabolic_Q_tilde(...)` 保留为明确的 coset-level 入口；
  - `affine_bounded_Q(...)` 成为 canonical 名称，`affine_bounded_Q_tilde(...)` 仅保留为 legacy alias。
- 已同步修正文档字符串、测试与缓存语义：磁盘缓存 `value_kind` 改为 `Q_at_one`，并兼容旧的 `Q_tilde_at_one` 载荷。
- 已完成定向验证：
  - `sage -python -m pytest pyw/tests/test_kazhdan_lusztig.py -q --no-cov`
  - 结果：**36 passed**
- 已新增 `pyw/core/affine_kl.py`，将 `order` 截断、dominant `\Lambda` 搜索、bounded stabilizer、quotient representative 收口为 workflow/context 层。
- 已在 `pyw/core/kazhdan_lusztig.py` 补充 bounded affine coset-level `Qtilde` helper，用于 `\widehat W_\Lambda / \widehat W_\Lambda^0` 的有界工作流。
- 已在 `pyw/core/character.py` 新增第一版 `KazhdanLusztigCharacter`，当前通过 `VermaCharacter` 进行 **grade-only** character 组装。
- 已新增测试：
  - `pyw/tests/test_affine_kl_context.py`
  - `pyw/tests/test_affine_kl_character.py`
- 已完成当前实现阶段的定向验证：
  - `sage -python -m pytest pyw/tests/test_kazhdan_lusztig.py pyw/tests/test_affine_kl_context.py pyw/tests/test_affine_kl_character.py -q --no-cov`
  - 结果：**44 passed**
- 已对 `demos/Algebra.py:GetLambda` 做逐函数上游对照，确认其核心搜索语义为：
  - 使用 `weight_lattice(extended=True).weyl_group(prefix='w')` 生成 `self.WW`
  - 直接扫描 `self.WW[:50]`
  - 用 `w.action(llambda + rho) - rho` 构造候选
  - 用 `RemoveDelta(weight).coefficients()` 的 `>= -1` 条件接受候选
- 已在 `pyw` 中开始按 legacy 语义回对 `GetLambda`：
  - 修正了 `D_4^{(1)}` 的 `marks` 计算
  - 修正了 legacy extended weight-lattice 组装中的 `FiniteFamily` 取值逻辑
  - 对 `-2\hat\omega_0` 的 D4 vacuum 例子，`old Algebra` 与 `pyw` 现已在以下中间量上逐项一致：
    - `llambda`
    - `rho`
    - `llambda + rho`
    - 前若干个 affine Weyl prefix 的 `w.action(llambda + rho) - rho`
- 当前仍在继续收尾 `GetLambda` 的三组 D4 基准权对照：
  - $-2\hat\omega_0$
  - $-2\hat\omega_1$
  - $-\hat\omega_2$
- 已补全 `AffineLieAlgebra` / `finite_lie_algebra` 的一组常用标量接口：
  - `.dim`
  - `.num_roots`
  - `.num_positive_roots`
  - `.num_negative_roots`
- 已将几处不必要的内部转发压平，避免 affine 情况下再绕回 `finite_lie_algebra` 一层获取：
  - `alpha_0()`
  - `rho()`
  - `positive_roots()`
  - `negative_roots()`
  - `theta()` 的最高根缓存路径
- 已完成这一轮定向验证：
  - `sage -python -m pytest pyw/tests/test_affine_lie_algebra.py -q --no-cov`
  - 结果：**45 passed**

## Next Concrete Step

- 在当前已经澄清的 `Q` / `Qtilde` 语义基础上，继续实现 PRD 中尚未完成的部分：
  - 完成 `GetLambda` 的 legacy 对齐，确保 D4 benchmark 权在 `pyw` 中得到与 `Algebra.py` 相同的 `\widehat\Lambda`；
  - 将当前以 integral/full-affine 情形为主的 workflow 推广到一般 `\widehat W_\Lambda`；
  - 从 `\widehat \Delta^{\mathrm{re}}_+(\widehat\Lambda)` 中真正恢复 simple roots 与 subgroup 生成元；
  - 细化 `\widehat W_\Lambda^0` 与 coset 最大代表元的有界构造；
  - 视需要从当前 grade-only character 扩展到更丰富的 character 表示。


## Computation workflow

- 用户提供 `alg`, `lambda_hat`, `order`
- 构造 affine Weyl group (as semi-direct product $\widehat W = Q^\vee \rtimes W$) 的系列 elements 
  $$
  t_{\alpha^\vee} = s_\alpha s_{ \alpha + \delta}, \qquad \alpha^\vee \in Q^\vee
  $$
- 寻找 $\widehat \Lambda = \widehat w_\lambda \cdot \widehat{\lambda}$，且 $\widehat \Lambda + \hat \rho$ 是 dominant 的。其中 $\widehat w_\lambda$ 是 $\widehat W$ 中的元素。 确定元素 $\widehat w_\lambda$ 使得 $\widehat w \cdot \widehat \Lambda = \widehat \lambda$ (注意是作用在 $\widehat \Lambda$ 上)。
- 确定 $\widehat \Delta^\text{re}_+(\widehat \Lambda) = \{\alpha \in \widehat \Delta^\text{re}_+ \ | \ (\widehat \Lambda, \alpha^\vee) \in \mathbb{Z} \}$，提取里头的 simple roots 以及对应的 reflection 生成的 subgroup $\widehat W_\Lambda$。

  When $\widehat \Lambda$ is integral (but not necessarily dominant), $\widehat W_\Lambda = \widehat W$. This applies to $\widehat{\mathfrak{so}}(8)_{-2}$, $(\mathfrak{e}_6)_{-3}$, for example.
- 计算 $\widehat W_\Lambda^0$ which is the subgroup keeping $\widehat \Lambda$ fixed
- 获取 inverse KL 多项式 `Q` for $\widehat W_{\Lambda}$
- 获取 inverse KL 多项式 `Qtilde` for $\widehat W_\Lambda/\widehat W_\Lambda^0$
  $$
  \widetilde Q_{[x], [y]} = \sum_{z \in [y]} Q_{\bar x, z}(-1)^{\ell (\bar x)} (-1)^{\ell(z)} \ .
  $$
  Here $\bar x$ means the maximal (not minimal) representative of the coset element $[x]$.
- 组装特征标: 列举商群元素 $[\widehat w'] \in \widehat W_\Lambda/ \widehat W_\Lambda^0$ 满足 Bruhat order $[\widehat w] \le [\widehat w']$，计算 $\widetilde Q_{[\widehat w], [\widehat w']}(1)$，以及 $e^{\widehat w' \cdot \widehat \Lambda}$，最后除以 Weyl denominator $\sum_{\widehat w \in \widehat W} (-1)^{\widehat w} e^{\widehat w(\widehat \rho) - \widehat \rho}$。
