# background

给定一个 affine weight $\hat \lambda$ (level $k$) 以及一个 finite coroot $\alpha^\vee$，可以计算它的 translation $t_{\alpha^\vee} \widehat \lambda$，
$$
\begin{align}
t_{\alpha^\vee} \widehat\lambda
= & \ (\lambda + k \alpha^\vee; k ; n + \frac{|\lambda|^2 - |\lambda + k \alpha^\vee|^2}{2k}) \\

= & \ (\lambda + k \alpha^\vee; k ; n - \frac{2}{|\alpha|^2}[(\lambda, \alpha) + k] ) \\
= & \ \Big(\lambda + k \alpha^\vee; k ; n - \big[(\lambda, \alpha^\vee) + \frac{1}{2}k(\alpha^\vee, \alpha^\vee) \big] \Big)
\end{align}
$$
特别地，
$$
\Delta n = - \big[(\lambda, \alpha^\vee) + \frac{1}{2}k(\alpha^\vee, \alpha^\vee) \big]
$$

# task: list translations by nshift

我希望列举所有的 $\alpha^\vee$ 使得 $0 \le - \Delta n \le \text{order}$. 请你
- 将这个不等式写成关于 $\alpha^\vee$ 分量或者 Dynkin labels 的不等式
- 实现一个函数，严格、高效地列举所有满足这个不等式的 $\alpha^\vee$
  
  可以参考 `_translations_by_n_shift` 这个旧函数，但这个函数是通过有限盒子穷举的笨方法实现的，有可能漏解，效率也很差。当你找到更好的方法之后，这个旧函数就可以被废弃了。

  