# Kazhdan-Lusztig

# Kazhdan-Lusztig Readme

　　We use $\widehat{\mathfrak{g}}$ to denote the affine Lie algebra of the finite Lie algebra $\widehat{\mathfrak{g}}$, and $r \coloneqq \operatorname{rank}\mathfrak{g}$, $\operatorname{rank}(\widehat{\mathfrak{g}}) = r + 1$.

- Download and install `SageMath`

  Use mainly the jupyter notebook environment for computation
- Navigate to folder `~/.sage/local/lib/python3.10/site-packages/coxeter3_sage`​, inside the file `coxeter3.py`, add the following codes following the comments

  > 文件也有可能在 `~/Library/SageMath-10-4/lib/python3.12/site-packages/coxeter3_sage/coxeter3.py`
  >
  > 可以在 Jupyter notebook 中通过 `from coxeter3_sage import Coxeter3`​ 然后 `?Coxeter3` 查看文件所在位置
  >
  > ![image](assets/image-20240812012209-8mqlem2.png)
  >

  ```python
  # at the beginning of the file, for using SR function at the end
  from sage.all import *


  # inside the Coxeter3 class definition, add the following method
  class Coxeter3:
  	... 

  	def invpol(self, x, y):
  	    # returns invpol result from coxeter3
  	        x = self.__ensure_element(x)
  	        y = self.__ensure_element(y)
  	        if x.bruhat_le(y):
  	            wordx = self.__convert_element_to_coxeter_input(x)
  	            wordy = self.__convert_element_to_coxeter_input(y)
  	            self._process.sendline("invpol")
  	            self._process.expect("first : ")
  	            self._process.sendline(wordx)
  	            self._process.expect("second : ")
  	            self._process.sendline(wordy)
  	            self._process.sendline("")
  	            self._process.expect("coxeter : ")
  	            before = self._process.before
  	            r = before.splitlines()[-1]
  	            return SR(str(r)[2:-1])
  	        else:
  	            return 0
  ```

- ​`Algebra.py`​ provides a class `Alg` with useful methods and properties.

  - It can only handle integral weights, unfortunately
  - The key method is `Alg.Kazhdan_Lusztig_numerator(weight, order)`​ it produces information on the numerator of the Kazhdan-Lusztig formula, up to the given `order`
  - ​`Alg.Kazhdan_Lusztig_denominator` gives the denominator
  - ​`Alg.save_numerator`​ and `Alg.save_denominator`​ save the generated numerator and denominator to file in the `numerators`​ and `denominators`​ subfolder, which will then processed by Mathematica file `Kazhdan-Lusztig.nb`
  - ​`.dat` files are the stored results of inverse KL polynomials, drastically accelerate some computation

　　‍

# Affine Lie algebra conventions

> # Finite weights and roots
>
> - The fundamental weights of $\mathfrak{g}$ are denoted as $\omega_{i = 1, ..., r}$
> - The simple roots are $\alpha_{i = 1, ..., r}$. They form a integral basis for the set of roots $\Delta$.
> - 对每个根都可以计算其 co-version
>
>   $$
>   \alpha^\vee = \frac{2\alpha}{|\alpha|^2}
>   $$
>
>   其中
>
>   $$
>   (\alpha_i^\vee, \omega_j) = \delta_{ij}
>   $$
>
>   > 当 Lie algebra 是 simply-laced 的，全部 $|\alpha|^2 = 2$，因此
>   >
>   > $$
>   > \alpha^\vee = \alpha, \qquad
>   > \alpha_i^\vee = \alpha_i
>   > $$
>   >
> - 每一个 **fundamental** weight 都可以计算它的 co-version，定义为
>
>   $$
>   (\omega_i^\vee, \alpha_j) = \delta_{ij}
>   $$
>
>   > 注意 weight 的 co-version **不是通过**
>   >
>   > $$
>   > \omega^\vee = \frac{2\omega}{|\omega|^2}
>   > $$
>   >
>
>   > 当 Lie algebra 是 simply laced 的，由于 $\alpha_i = \alpha_i^\vee$，因此
>   >
>   > $$
>   > \omega_i = \omega_i^\vee
>   > $$
>   >
> - The finite Weyl vector is $\rho = \sum_{i = 1}^r \omega_i$
> - The heighest root $\theta$, and it's related o the dual Coxeter number $h^\vee$,
>
>   $$
>   \theta = \sum_{i = 1}^r a_i^\vee \alpha_i^\vee = \sum_{i = 1}^r a_i \alpha_i \ , \qquad
>   h^\vee = 1 + \sum_{i = 1}^r a^\vee_i \ .
>   $$
> - A finite weight can be denoted using Dynkin labels
>
>   $$
>   \lambda = \sum_{i = 1}^r \lambda_i \omega_i, \qquad
>   \sim 
>   \qquad
>   [\lambda_1, ..., \lambda_r] \ .
>   $$
>
> # Affine weights and roots
>
> - An affine weight can be denoted as (the <span data-type="text" style="background-color: var(--b3-font-background8);">extended weight space</span> notation, in an $r + 2$ dimensional vector space)
>
>   $$
>   \widehat{\lambda}
>
>   = (\lambda; k; n) \ .
>   $$
>
>   where $k$ is the level (eigenvalue of the central element $K$, and $n$ is the eigenvalue of $- L_0$).
>
>   > For example, $H^I_n \sim (0;0;n)$, $E^\alpha_n \sim (\alpha; 0; n)$.
>   >
>
>   <span data-type="text" style="background-color: var(--b3-font-background8);">Inner product</span> of affine weights is given by
>
>   $$
>   (\widehat \lambda, \widehat \mu) = (\lambda, \mu) + k_\lambda n_\mu + k_\mu n_\lambda \ .
>   $$
>
>   An affine coweight with respect to a affine is defined by
>
>   $$
>   \widehat \lambda^\vee \coloneqq \frac{2 \widehat\lambda}{(\widehat \lambda, \widehat \lambda)} \ .
>   $$
> - The affine fundamental weights are denoted $\widehat{\omega}_{i = 0, 1, ..., r}$, and explicity
>
>   $$
>   \widehat{\omega}_0 = ([0]; 1; 0), \qquad 
>   \widehat{\omega}_i = (\omega_i, a_i^\vee, 0) \ .
>   $$
> - The affine simple roots are denoted $\widehat{\alpha}_{i = 0, 1, ..., r}$,
>
>   $$
>   \widehat{\alpha}_0 = (- \theta; 0; 1), \qquad
>   \widehat{\alpha}_{i = 1, ..., r} = (\alpha_i; 0; 0) \ .
>   $$
>
>   > Sometimes we just denote an affine root $(\alpha; 0; 0)$ as $\alpha$.
>   >
> - Define the ***basic imaginary root*** $\delta = (0; 0; 1)$. It's related to the simple roots by
>
>   $$
>   \delta = \sum_{i = 0}^r a_i \widehat\alpha_i \ .
>   $$
>
>   Affine roots of the form $n \delta$ are called ***imaginary roots***.
> - An affine weight can also be written in terms of ***affine Dynkin labels***,
>
>   $$
>   [\lambda; λ_0 + \sum_{i = 1}^r \lambda_i a^\vee_i; n = 0] \sim [\lambda_0, \lambda_1, \cdots,\lambda_r] \ .
>   $$
>
>   > The affine Dynkin label notation always implies $n = 0$. One can include $n \ne 0$ by adding $\delta$:
>   >
>   > $$
>   > [\lambda_0, \lambda_1, \cdots, \lambda_r] + n \delta \ .
>   > $$
>   >
>
>   > Ignoring $n$, we often denote $\widehat \lambda$ by just $\lambda$ and expliclitly states its $k$ value.
>   >
>   > $$
>   > (\lambda, k) \sim [\lambda_0 = k - (\lambda, \theta), \lambda_1, ..., \lambda_r] \ .
>   > $$
>   >
>   > See Di Francesco's eq (14.57)
>   >
> - The ***affine roots*** of $\widehat{\mathfrak{g}}$ are given by $(\alpha; 0; n) = \alpha + n \delta$, for all $\alpha \in \Delta$, $n \in \mathbb{Z}$.
>
>   Among them, the ***affine positive roots*** are
>
>   $$
>   \widehat \Delta_+ = \{\alpha + n \delta \ | \ n \in \mathbb{Z}_{> 0}, \alpha \in \Delta \text{ or } \alpha = 0\} \cup \{\alpha \in \Delta_+\} \ .
>   $$
> - ***Imaginary roots*** defined to be $n \delta$, while other affine roots are all ***real roots***.
>
>   > All real roots $\widehat \alpha$ have multiplicity $\text{mult}(\widehat\alpha) = 1$, but $\text{mult}(n \delta) = r$.
>   >
>
>   Therefore, within the affine positive roots, we define the set of ***affine positive real roots***
>
>   $$
>   \widehat \Delta^\text{re}_+ = \{\alpha + n \delta \ | \ n \in \mathbb{Z}_{> 0}, \alpha \in \Delta\} \cup \{\alpha \in \Delta_+\} \ .
>   $$
> - The affine Weyl vector $\widehat{\rho}$
>
>   $$
>   \widehat{\rho} = \sum_{i = 0}^{r} \widehat{\rho}_i
>   $$
>
> # Affine Weyl group
>
> For any affine root $\widehat{\alpha} = (\alpha; 0; m)$, the associated affine Weyl reflection is defined as
>
> $$
> s_{\widehat{\alpha}} \lambda = \lambda - (\widehat{\lambda}, \widehat{\alpha}^\vee) \widehat{\alpha} \ .
> $$
>
> - When $m = 0$ (so $\widehat \alpha$ is a **finite** root)
>
>   $$
>   s_{\widehat \alpha} \widehat \lambda = (s_\alpha \lambda; k; n) \ , \qquad
>
>   \widehat \lambda = (\lambda; k; n) \ .
>   $$
> - In general, for $\widehat \lambda = (\lambda; k; n)$,
>
>   $$
>   \begin{align}
>   s_{\widehat \alpha} \widehat \lambda
>
>   = & \ (s_\alpha \lambda - km\alpha^\vee; k; n - [(\lambda, \alpha) + km] \frac{2m}{|\alpha|^2}) \\
>
>   = & \ (s_\alpha(\lambda + k m \alpha^\vee); k ; n + \frac{|\lambda|^2 - |\lambda + k m \alpha^\vee|^2}{2k})  
>   \end{align}
>   $$
>
>   > See Di Francesco's eq (14.64)
>   >
> - The $\delta$ component of $\widehat \lambda,$ namely, the $n$-value of $\widehat \lambda$ is shifted by the ammount
>
>   $$
>   \Delta n= - \Big[(\lambda, \alpha) + km \Big] \frac{2m}{|\alpha|^2}
>
>   = \frac{|\lambda|^2 - |\lambda + k m \alpha^\vee|^2}{2 k}
>    \ .
>   $$
> - All $s_{\widehat \alpha}$ generate the affine Weyl group $\widehat W$.
> - All elemenets in $\widehat W$ can be generated by $s_{\widehat \alpha_{i = 0, 1, 2, ..., r}}$.
>
>   - The 0-th affine simple reflection is
>
>     $$
>     s_{\widehat \alpha_0} \widehat \lambda = \widehat \lambda - [(\lambda, - \theta^\vee) + k] \widehat \alpha_0 \\
>     = (s_{-\theta}\lambda + k \theta; k ; n + (\lambda, \theta) - k) \ .
>     $$
>   - The $i = 1, ..., r$-th affine simple reflections are
>
>     $$
>     s_{\widehat \alpha_i} \widehat \lambda = \widehat \lambda - (\lambda, \alpha_i^\vee)\widehat \alpha_i
>
>     = (s_\alpha \lambda; k; n) \ .
>     $$
>   - The $i = 1, ..., r$ simple reflections generate the finite Weyl group $W \le \widehat W$. Finite Weyl group $W$ does not change the $n$-value of an affine weight $(\lambda; k; n)$.
>   - Each element $\widehat w$ in $\widehat W$ can be written as some product of the simple affine Weyl reflections. The length of the **shortest** such product of $\ell$ simple reflections defines the ***length function*** $\ell(\widehat w)$.
>
>     > The shortest such product expression of $\widehat w$ is called a ***reduced expression***<span data-type="text" style="background-color: var(--b3-font-background8);"> </span>for $\widehat w$.
>     >
> - Define the dot-action of any element $\widehat w \in \widehat W$
>
>   $$
>   \widehat w \cdot \widehat \lambda \coloneqq \widehat w (\widehat \lambda + \widehat \rho) - \widehat \rho \ .
>   $$
>
> ### Affine Weyl group as semidirect product
>
> The **Affine Weyl group** $\widehat W$ can be understood as a semi-product group $\widehat W = Q^\vee \rtimes W$, where $W$ is the **finite Weyl group** of $\mathfrak{g}$, and $Q^\vee$ is the **finite coroot lattice**. The $Q^\vee$ elements correspond to Weyl translations,
>
> $$
> t_{\alpha^\vee} = s_\alpha s_{ \alpha + \delta}
> $$
>
>> Here $\alpha = (\alpha; 0; 0)$，$\alpha^\vee$ 属于 finite coroot lattice
>>
>
> By direct computation, one has
>
> $$
> \begin{align}
> t_{\alpha^\vee} \widehat\lambda
> = & \ (\lambda + k \alpha^\vee; k ; n + \frac{|\lambda|^2 - |\lambda + k \alpha^\vee|^2}{2k}) \\
>
> = & \ (\lambda + k \alpha^\vee; k ; n - \frac{2}{|\alpha|^2}[(\lambda, \alpha) + k] ) \\
> = & \ \Big(\lambda + k \alpha^\vee; k ; n - \big[(\lambda, \alpha^\vee) + \frac{1}{2}k(\alpha^\vee, \alpha^\vee) \big] \Big)
> \end{align}
> $$
>
>> See Di Franceso's eq (14.69)
>>
>
>> 也有的文献将 $t$ 的下标拓展到 finite coweight lattice $P^\vee$
>>
>> $$
>> t_\beta \widehat\lambda = (\lambda + k \beta; k; n - [(\lambda, \beta) + \frac{1}{2}k(\beta, \beta)]) = \widehat\lambda + k \beta - [(\lambda, \beta) + \frac{1}{2}k(\beta, \beta)]\delta, \qquad
>> \beta \in P^\vee
>> $$
>>
>
> Note that $\delta$ in $s_{\alpha + \delta}$ provides the $m = 1$ effect, participating in the shift
>
> $$
> \Delta n = - [(\lambda, \alpha) + k] \frac{2}{|\alpha|^2} \ .
> $$
>
>> Note also that $\delta$ can be obtained by
>>
>> $$
>> \delta = \sum_{i = 0}^r a_i \widehat \alpha_i \ .
>> $$
>>
>
>> For an arbitrary finite weight $\lambda$, and a **finite positive** root $\alpha$,
>>
>> $$
>> (\lambda, \alpha) \frac{2}{|\alpha|^2} = (\lambda, \alpha^\vee) = (\lambda, \sum_{i = 1}^r m_i\alpha_i^\vee)
>>
>> = \sum_{i = 1}^r \lambda_i m_i\ .
>> $$
>>
>> The coefficients $m_i\ge 0$, and cannot be all zero.
>>
>> - When $\lambda_i \ge 0$, namely when $\lambda$ is dominant, $(\lambda, \alpha) \ge 0$.
>> - When the affine weight $\widehat \lambda$ is integral dominant, then $\lambda_0 \ge 0$, $\lambda_i \ge 0$, and
>>
>>   $$
>>   k + (\lambda, \alpha) = \lambda_0 + (\lambda, \theta + \alpha) \ .
>>   $$
>>
>>   Since $\theta$ is highest, $\theta + \alpha$ must be non-negative linear combination of the simple roots $\alpha_i$, and therefore $(\lambda, \theta + \alpha) \ge 0$. Hence, $\Delta n \le 0$.
>>
>
> More generally,
>
> $$
> t_{\alpha^\vee}^m \widehat \lambda
> = (\lambda + m k \alpha^\vee; k ; n + \frac{|\lambda|^2 - |\lambda + k m \alpha^\vee|^2}{2k}) \ .
> $$
>
> - Weyl translations mutually commute,
>
>   $$
>   t_{\alpha^\vee} t_{\beta^\vee} = t_{\beta^\vee} t_{\alpha^\vee} = t_{\alpha^\vee + \beta^\vee} \ .
>   $$
> - All Weyl translations can be generated by the simple translations associated with the finite simple corrots
>
>   $$
>   t_{\alpha^\vee_1}, ..., t_{\alpha^\vee_r} \ .
>   $$
> - Any affine Weyl element $\widehat w \in \widehat W$ can be written as
>
>   $$
>   \widehat w = w \prod_{i = 1}^r t_{\alpha_i^\vee}^{m_i} = w t_{\alpha^\vee}, \qquad
>   \alpha^\vee = \sum_{i = 1}^r m_i \alpha^\vee_i \qquad
>   w \in W \ .
>   $$
> - Note that the finite Weyl element $w$ does not change the $n$-value of an affine weight $(\lambda; k; n)$: change of $n$-value comes solely from the simple translation,
>
>   $$
>   \Delta n = \frac{|\lambda|^2 - |\lambda + k m \alpha^\vee|^2}{2k}
>
>   = - \Bigg[
>   \sum_{i = 1}^r m_i \lambda_i
>
>   + \frac{k}{2} \sum_{i, j = 1}^r m_i m_j (\alpha_i^\vee, \alpha_j^\vee)
>
>   \Bigg] \ .
>   $$
>
>   > Consider translation $t_{\alpha_i^\vee}$ and $t_{\alpha_i^\vee}^{-1}$, namely $m_i = \pm 1$, $m_{j \ne i} = 0$. Then
>   >
>   > $$
>   > \Delta n = \frac{|\lambda|^2 - |\lambda + k m \alpha^\vee|^2}{2k}
>   >
>   > = - \Bigg[
>   > \pm\lambda_i
>   >
>   > + \frac{2k}{|\alpha_i|^2}
>   >
>   > \Bigg] \ .
>   > $$
>   >
>   > When the Lie algebra is **ADE-type**, $|\alpha_i|^2 = 2$ and the above gives
>   >
>   > $$
>   > \Delta n = - (\lambda_i + k) = -(\pm \lambda_i + \lambda_0 + (\lambda, \theta)) \ .
>   > $$
>   >
>   > Apparently, for **affine integral dominant** weight $\widehat \lambda$,
>   >
>   > $$
>   > \pm \lambda_i + \lambda_0 + (\lambda, \theta) \ge 0 \Rightarrow \Delta n \le 0 \ .
>   > $$
>   >
>   > Hence, any **translations**, and in fact the **entire affine Weyl group** can only **decrease** the $n$-value of an affine integral dominant $\widehat \lambda$.
>   >

# Kazhdan-Lusztig formula

　　Consider an affine Lie algebra $\widehat{\mathfrak{g}}_k$ at level $k$, with highest weight state $\widehat \lambda = (\lambda; k; 0)$. In terms of Dynkin labels,

$$
\widehat\lambda = [k - (\lambda, \theta)] \widehat\omega_0 + \sum_{i = 1}^r \lambda_i \widehat\omega_i \ , \qquad
\lambda \sim [\lambda_1, ..., \lambda_r] \ .
$$

　　Take $\widehat \lambda$ as the highest weight, one can construct two modules $M(\lambda)$ and $L(\lambda)$

- The level-$k$ Verma module $M_{\lambda}$. Its character is given by

  $$
  \operatorname{ch}M(\lambda) 
  = e^{\widehat \lambda} \prod_{\widehat\alpha > 0} (1 - e^{- \widehat\alpha})^{- \text{mult}(\widehat\alpha)} \ .
  $$

  Here $\text{mult}$ is the multiplicity of an affine root.

  By Weyl-Kac, or the Macdonald-Weyl denominator formula, we have

  $$
  \prod_{\widehat\alpha > 0} (1 - e^{- \widehat\alpha})^{- \text{mult}(\widehat\alpha)} \ .


  = \frac{1}{\sum_{\widehat w \in \widehat W} (-1)^{\widehat w} e^{\widehat w(\widehat \rho) - \widehat \rho} }\ .
  $$

  > See Di Francesco eq (14.150)
  >

  > The $e^{-\delta}$ factors in the character correspond to the parameter $q$. Module $M(\lambda)$ is obtained by acting $J^a_{n < 0}$ on the highest weight state, and therefore the $n$ value decreases. Consequently, the character $\operatorname{ch}M(\lambda)$ is a increasing power series in $q$.
  >
- The irreducible module $L_{\lambda}$ by removing nulls from $M_{ \lambda}$ (if there is any). As long as $k + h^\vee > 0$, one can apply the Kazhdan-Lusztig formula,

  $$
  \operatorname{ch}L(\widehat \lambda) = \sum_{\widehat \mu \le \widehat \lambda} m_{\widehat \lambda, \widehat \mu} \operatorname{ch}M(\widehat \mu) \ .
  $$

  where the weights $\widehat \mu$ labels the so-called primitive null vectors that appear in the reducible $M(\widehat \lambda)$, and $m_{\widehat \lambda, \widehat \mu}$ could be positive or negative, reflecting the intricate subtraction of the nulls.

  All the weights $\widehat \mu$ that can appear in the above sum must satisfy the equality

  $$
  |\widehat \mu + \widehat \rho| ^2 = | \widehat \lambda + \widehat \rho|^2 \ .
  $$

  When $k + h^\vee > 0$, the above equality is fairly restrictive: all the solutions to the above equation lies on the **affine dot-Weyl orbit**, namely, for any solution $\widehat \mu$ there must exists a $\widehat w$ such that

  $$
  \widehat \mu = \widehat w \cdot \widehat \lambda \ .
  $$

  Moreover, among all the solutions, there is a **unique** $\widehat \Lambda$ such that $\widehat \Lambda + \widehat \rho$ is dominant (having non-negative Dynkin labels).

  > Note that it is AFTER adding $\widehat \rho$, that $\widehat \Lambda + \widehat \rho$ is dominant.
  >

  Since all the solutions lie on the same orbit, one can denote all the **potential** solutions $\widehat \mu$ by $\widehat w \cdot \widehat \Lambda$, and exhaust all $\widehat w \in \widehat W$.

  Fix a special $\widehat w$ such that $\widehat w \cdot \widehat \Lambda = \widehat \lambda,$ and therefore

  $$
  \operatorname{ch}(L_\lambda)

  = \sum_{\widehat w' \in \widehat W} m_{\widehat w, \widehat w'}

  \operatorname{ch}M(\widehat \mu = \widehat w'\cdot \widehat \Lambda) \ .
  $$

  **Kazhdan-Lusztig formula** precisely tells us what the $m_{\widehat w, \widehat w'}$ are.

  - Let $\widehat \Delta^\text{re}_+(\Lambda) = \{\alpha \in \Delta^\text{re}_+ \ | \ (\Lambda, \alpha^\vee) \in \mathbb{Z} \}$ be a subset of positive real roots, and consider the subgroupo $\widehat W_\Lambda$ generated by the corresponding affine Weyl reflections.

    > When $\Lambda$ is integral (but not necessarily dominant), $\widehat W_\Lambda = \widehat W$. This applies to $\widehat{\mathfrak{so}}(8)_{-2}$, $(\mathfrak{e}_6)_{-3}$, for example.
    >
  - In $\widehat W_\Lambda$, there could be a subgroup $\widehat W^0_\lambda$ that leaves $\Lambda$ itself **invariant**: consider the coset $\widehat W_\Lambda / \widehat W^0_\Lambda$. The level-$k$ character $\operatorname{ch}(L_\lambda)$ receives only contributions from $[\widehat w'] \in \widehat W_\Lambda / \widehat W^0_\Lambda$. More precisely, one should pick the **shortest ** representative for $[\widehat w']$ based on the length function $\ell(\widehat w')$.
  - In $\widehat W$, and therefore in $\widehat W_\Lambda$, there is a Bruhat (partial)-order $\le$.

    It induces an order in the coset $\widehat W_\Lambda/\widehat W_\Lambda^0$ by comparing the Bruhat order of the **shortest/minimal** representatives.

    > Two elements $\widehat w, \widehat w' \in \widehat W$, if the reduced expression for $\widehat w$ can be obtained by **dropping** some simple reflections from that of $\widehat w'$, then we say in **Bruhat order** $\widehat w \le \widehat w'$.
    >
    > Obviously $\widehat w \le \widehat w' \Rightarrow \ell(\widehat w) \le\ell(\widehat w')$.
    >
  - With the inverse **Kazhdan-Lusztig inverse polynomial** $\tilde Q$,

    $$
    m_{\widehat w, \widehat w'} = \tilde Q_{\widehat w, \widehat w'}(1) \ ,
    $$
  - To summarize,

    $$
    \operatorname{ch}(L_\lambda)

    = \sum_{\substack{[\widehat w'] \in \widehat W_\Lambda/ \widehat W_\Lambda^0\\ [\widehat w] \le [\widehat w']


    }} \widetilde Q_{[\widehat w], [\widehat w']}(1) \operatorname{ch} M(\widehat \mu = \widehat w'\cdot \widehat \Lambda) .
    $$

    > The **leading **​**$q$**​ **-power** of $M(\widehat \mu = \widehat w' \cdot \widehat \Lambda)$ is given by $- n(\widehat w' \cdot \widehat \Lambda)$.
    >
    > As $\widehat w'$ changes, since $\widehat w' \cdot \widehat \Lambda = \widehat w'(\widehat \Lambda + \rho) - \rho$, the $n$-value of $\widehat w ' \cdot \widehat \Lambda$ will only **decrease. ** Therefore $\operatorname{ch}(L_\lambda)$ has the lowest $q$-power given by $n(\widehat \Lambda= \operatorname{id} \cdot \widehat \Lambda)$, and followed by **higher powers** of $q$: these higher powers of $q$ are due to the **translations**.
    >
    > Hence we shall exploit the semidirect product of $\widehat W$ to implement the sum over cosets: we shall construct the infinite group $\widehat W$ by first work out translations that shift $n(\widehat \Lambda + \widehat \rho)$ by certain amount, and then take the product with the finite Weyl group.
    >

    > However, the finite Weyl group $W$ is huge for $E_6$, so constructing $\widehat W$ is a computational heavy task.
    >

    > The inverse KZ polynomial for the **coset **​**$\widehat W_\Lambda/\widehat W^0_\Lambda$** is derived from the KZ polynomial of $\widehat W_\Lambda$,
    >
    > $$
    > \widetilde Q_{[x], [y]} = \sum_{z \in [y]} Q_{\bar x, z}(-1)^{\ell (\bar x)} (-1)^{\ell(z)} \ .
    > $$
    >
    > See the eq (C.28) in paper [Cordova, Gaiotto, Shao] Infrared Computations of Defect Schur Indices. Here $\bar x$ means the **maximal (not minimal)**  representative of the coset element $[x]$.
    >

    Writing out the $\operatorname{ch}M$, we have

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

    ‍

　　‍

　　‍

　　panyw5@mail.sysu.edu.cn
