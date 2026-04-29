# 背景

从任意 affine weight $\widehat \lambda = (\lambda; k, n)$ 出发，我们需要计算哪些 translations 让 $\widehat \lambda$ 的 $\delta$ 分量系数 (也称为 $\widehat \lambda$ 的 level/grade，或者 $n$-value) 变化量 $\Delta n$ 满足不等式
$$
0 \le - \Delta n \le \text{upper bound}
$$
`characters.py` 里面已经有 `_translations_by_n_shift_impl`，但设计不合理：它接受两个 upper bound 地位的参数：`order`, `max_neg_shift`。这不是正确的设计。

# 目标
- 修正 `_translations_by_n_shift_impl` 的设计，使得它接受一个 upper bound 参数。
- 修正引用 `_translations_by_n_shift_impl` 的业务逻辑。比如，KazhdanLusztigCharacter 计算中，$-\Delta n$ 应当满足
  $$
  0 \le -\Delta n \le \text{order} + n(\widehat \Lambda + \widehat \rho) - n(\widehat \lambda)
  $$
  其中 $\widehat \lambda$ 是表示的最高权，`order` 是特征标计算准确到 $q^\text{order}$ 的幂次，$\widehat \Lambda = w.\widehat \lambda$ 使得 $\widehat \Lambda + \widehat \rho$ 是 dominant 的。
  
  因此应该使用正确的 upper bound ``order + (Lambda_hat + rho_hat).grade - lambda_hat.grade`` 来调用 `_translations_by_n_shift_impl`。
- IntegrableModuleCharacter 也要正确调整 upper bound。