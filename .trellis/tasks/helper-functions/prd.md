# PRD: add helper functions to facilitate computations

## weights and roots
- [x] 所有的 (affine) weight, (affine) root 都配备一些 translator，写成 simple root basis, fundamental weight basis, etc，写成相应的基底的线性组合
- [x] AffineLieAlgebra 的 finite_weyl_group 添加列举所有元素的方法 `.list(prefix="s")` (output list of prefix products `[1, s1, s2, s1*s2, ...]`) 
- [x] finite_weyl_group 通过 `associated_reflection(root)` 得到的输出能否改为标准 sage 的 `weight_lattice.weyl_group(prefix="s").from_reduced_word` 的输出？换言之，没有必要单独设立 `FiniteWeylGroupElement`。直接用回 `sage` 的外尔群元表示即可。
- [ ] symbol 与 weights, roots, coroots, coweights 等的乘积


## weyl group
- [x] `AffineLieAlgebra`  的 affine root 的 `.associated_reflection()` 现在输出一个数组，但我想你输出一个 `sage` 的外尔群元表示 (比如输出为 simple reflection 的积)，这样就可以直接用来做计算了。
- [x] `AffineWeylGroupSemidirectElement` 添加一个方法 `to_simple_reflection_basis_sage()`，输出对应的 `sage` 的外尔群元表示 (比如输出为 simple reflection 的积)
- [x] `AffineWeylGroupSemidirectElement` 的 `reduced_word()` 应该把 translation 部分也转化为 `s_{i = 0, ..., r}` 的乘积

