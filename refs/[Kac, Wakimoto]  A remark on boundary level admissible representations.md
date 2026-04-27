# A remark on boundary level admissible representations

Victor G. Kac* and Minoru Wakimoto†

Recently a remarkable map between 4-dimensional superconformal field theories and vertex algebras has been constructed [BLLPRV15]. This has lead to new insights in the theory of characters of vertex algebras. In particular it was observed that in some cases these characters decompose in nice products [XYY16], [Y16].

The purpose of this note is to explain the latter phenomena. Namely, we point out that it is immediate by our character formula [KW88], [KW89] that in the case of a boundary level the characters of admissible representations of affine Kac-Moody algebras and the corresponding  $W$ -algebras decompose in products in terms of the Jacobi form  $\vartheta_{11}(\tau, z)$ .

We would like to thank Wenbin Yan for drawing our attention to this question.

Let  $\mathfrak{g}$  be a simple finite-dimensional Lie algebra over  $\mathbb{C}$ , let  $\mathfrak{h}$  be a Cartan subalgebra of  $\mathfrak{g}$ , and let  $\Delta \subset \mathfrak{h}^*$  be the set of roots. Let  $Q = \mathbb{Z}\Delta$  be the root lattice and let  $Q^{*} = \{h\in \mathfrak{h}\mid \alpha (h)\in \mathbb{Z}$  for all  $\alpha \in \Delta \}$  be the dual lattice. Let  $\Delta_{+}\subset \Delta$  be a subset of positive roots, let  $\{\alpha_{1},\ldots ,\alpha_{\ell}\}$  be the set of simple roots and let  $\rho$  be half of the sum of positive roots. Let  $W$  be the Weyl group. Let  $(\cdot |.)$  be the invariant symmetric bilinear form on  $\mathfrak{g}$ , normalized by the condition  $(\alpha |\alpha) = 2$  for a long root  $\alpha$ , and let  $h^\vee$  be the dual Coxeter number ( $= \frac{1}{2}$  eigenvalue of the Casimir operator on  $\mathfrak{g}$ ). We shall identify  $\mathfrak{h}$  with  $\mathfrak{h}^*$  using the form  $(\cdot |.)$ .

Let  $\widehat{\mathfrak{g}} = \mathfrak{g}[t,t^{-1}] + \mathbb{C}K + \mathbb{C}d$  be the associated to  $\mathfrak{g}$  affine Kac-Moody algebra (see [K90] for details), let  $\widehat{\mathfrak{h}} = \mathfrak{h} + \mathbb{C}K + \mathbb{C}d$  be its Cartan subalgebra. We extend the symmetric bilinear form  $(\cdot |.)$  from  $\mathfrak{h}$  to  $\widehat{\mathfrak{h}}$  by letting  $(\mathfrak{h}|\mathbb{C}K + \mathbb{C}d) = 0$ ,  $(K|K) = 0$ ,  $(d|d) = 0$ ,  $(d|K) = 1$ , and we identify  $\widehat{\mathfrak{h}}^*$  with  $\widehat{\mathfrak{h}}$  using this form. Then  $d$  is identified with the  $0^{th}$  fundamental weight  $\Lambda_0 \in \widehat{\mathfrak{h}}^*$ , such that  $\Lambda_0|_{\mathfrak{g}[t,t^{-1}] + \mathbb{C}d} = 0$ ,  $\Lambda_0(K) = 1$ , and  $K$  is identified with the imaginary root  $\delta \in \widehat{\mathfrak{h}}^*$ . Then the set of real roots of  $\widehat{\mathfrak{g}}$  is  $\hat{\Delta}^{\mathrm{re}} = \{\alpha + n\delta | \alpha \in \Delta, n \in \mathbb{Z}\}$  and the subset of positive real roots is  $\hat{\Delta}_{+}^{\mathrm{re}} = \Delta_{+} \cup \{\alpha + n\delta | \alpha \in \Delta, n \in \mathbb{Z}_{\geq 1}\}$ . Let  $\hat{\rho} = h^{\vee}\Lambda_0 + \rho$ . Let

$$
\hat {\Pi} _ {u} = \{u \delta - \theta , \alpha_ {1}, \ldots , \alpha_ {\ell} \},
$$

where  $\theta \in \Delta_{+}$  is the highest root, so that  $\hat{\Pi}_1$  is the set of simple roots of  $\widehat{\mathfrak{g}}$ . For  $\alpha \in \hat{\Delta}^{\mathrm{re}}$  one lets  $\alpha^{\vee} = 2\alpha / (\alpha | \alpha)$ . Finally, for  $\beta \in Q^{*}$  define the translation  $t_{\beta} \in \operatorname{End} \widehat{\mathfrak{h}}^{*}$  by

$$
t _ {\beta} (\lambda) = \lambda + \lambda (K) \beta - ((\lambda | \beta) + \frac {1}{2} \lambda (K) | \beta | ^ {2}) \delta .
$$

Given  $\Lambda \in \widehat{\mathfrak{h}}^*$  let  $\hat{\Delta}^{\Lambda} = \{\alpha \in \hat{\Delta}^{\mathrm{re}} | (\Lambda |\alpha^{\vee}) \in \mathbb{Z}\}$ . Then  $\Lambda$  is called an admissible weight if the following two properties hold

(i)  $(\Lambda + \widehat{\rho} | \alpha^{\vee}) \notin \mathbb{Z}_{\leq 0}$  for all  $\alpha \in \hat{\Delta}_{+}$ ,  
(ii)  $\mathbb{Q}\hat{\Delta}^{\Lambda} = \mathbb{Q}\hat{\Delta}$

If instead of (ii) a stronger condition holds:

(ii)'  $\varphi (\hat{\Delta}^{\Lambda}) = \hat{\Delta}$  for a linear isomorphism  $\varphi :\widehat{\mathfrak{h}}^{*}\to \widehat{\mathfrak{h}}^{*}$

then  $\Lambda$  is called a principal admissible weight. In [KW89] the classification and character formulas for admissible weights is reduced to that for principal admissible weights. The latter are described by the following proposition.

Proposition 1. [KW89] Let  $\Lambda$  be a principal admissible weight and let  $k = \Lambda(K)$  be its level. Then

(a)  $k$  is a rational number with denominator  $u\in \mathbb{Z}_{\geq 1}$ , such that

$$
k + h ^ {\vee} \geq \frac {h ^ {\vee}}{u} a n d \operatorname * {g c d} (u, h ^ {\vee}) = \operatorname * {g c d} (u, r ^ {\vee}) = 1, \tag {1}
$$

where  $r^{\vee} = 1$  for  $\mathfrak{g}$  of type  $A$ - $D$ - $E$ ,  $= 2$  for  $\mathfrak{g}$  of type  $B$ ,  $C$ ,  $F$ , and  $= 3$  for  $\mathfrak{g} = G_2$ .

(b) All principal admissible weights are of the form

$$
\Lambda = (t _ {\beta} y). (\Lambda^ {0} - (u - 1) (k + h ^ {\vee}) \Lambda_ {0}), \tag {2}
$$

where  $\beta \in Q^{*},y\in W$  are such that  $(t_{\beta}y)\hat{\Pi}_u\subset \hat{\Delta}_{+},\Lambda^0$  is an integrable weight of level  $u(k + h^{\vee}) - h^{\vee}$ , and dot denotes the shifted action:  $w.\Lambda = w(\Lambda +\widehat{\rho}) - \widehat{\rho}$ .

(c) For  $\mathfrak{g} = s\ell_N$  all admissible weights are principal admissible.

Recall that the normalized character of an irreducible highest weight  $\widehat{\mathfrak{g}}$ -module  $L(\Lambda)$  of level  $k \neq -h^{\vee}$  is defined by

$$
\mathrm {c h} _ {\Lambda} (\tau , z, t) = q ^ {m _ {\Lambda}} \mathrm {t r} _ {L (\Lambda)} e ^ {2 \pi i h}
$$

where

$$
h = - \tau d + z + t K, z \in \mathfrak {h}, \tau , t \in \mathbb {C}, \operatorname {I m} \tau > 0, q = e ^ {2 \pi i \tau}, \tag {3}
$$

and  $m_{\Lambda} = \frac{|\Lambda + \widehat{\rho}|^2}{2(k + h^{\vee})} - \frac{\dim \mathfrak{g}}{24}$  (the normalization factor  $q^{m_{\Lambda}}$  "improves" the modular invariance of the character).

In [KW89] the characters of the  $\widehat{\mathfrak{g}}$ -modules  $L(\Lambda)$  for arbitrary admissible  $\Lambda$  were computed, see Theorem 3.1, or formula (3.3) there for another version in case of a principal admissible  $\Lambda$ . In order to write down the latter formula, recall the normalized affine denominator for  $\widehat{\mathfrak{g}}$ :

$$
\hat {R} (h) = q ^ {\frac {\dim \mathfrak {g}}{2 4}} e ^ {\widehat {\rho} (h)} \prod_ {n = 1} ^ {\infty} (1 - q ^ {n}) ^ {\ell} \prod_ {\alpha \in \Delta_ {+}} (1 - e ^ {\alpha (z)} q ^ {n}) (1 - e ^ {- \alpha (z)} q ^ {n - 1}).
$$

In coordinates (3) this becomes:

$$
\hat {R} (\tau , z, t) = (- i) ^ {| \Delta_ {+} |} e ^ {2 \pi i h ^ {\vee} t} \eta (\tau) ^ {\frac {1}{2} (3 \ell - \dim \mathfrak {g})} \prod_ {\alpha \in \Delta_ {+}} \vartheta_ {1 1} (\tau , \alpha (z)), \tag {4}
$$

where

$$
\vartheta_ {1 1} (\tau , z) = - i q ^ {\frac {1}{1 2}} e ^ {- \pi i z} \eta (\tau) \prod_ {n = 1} ^ {\infty} (1 - e ^ {- 2 \pi i z} q ^ {n}) (1 - e ^ {2 \pi i z} q ^ {n - 1})
$$

is one of the standard Jacobi forms  $\vartheta_{ab}$ ,  $a, b = 0$  or 1 (see e.g., Appendix to [KW14]), and  $\eta(\tau)$  is the Dedekind eta function.

For a principal admissible  $\Lambda$ , given by (2), formula (3.3) from [KW89] becomes in coordinates (3):

$$
\left. \left(\hat {R} \operatorname {c h} _ {\Lambda}\right) (\tau , z, t) = \left(\hat {R} \operatorname {c h} _ {\Lambda^ {0}}\right) \left(u \tau , y ^ {- 1} (z + \tau \beta), \frac {1}{u} (t + (z | \beta) + \frac {\tau | \beta | ^ {2}}{2})\right). \right. \tag {5}
$$

It follows from (5) that if  $\Lambda^0 = 0$  in (2) (so that  $\mathrm{ch}_{\Lambda^0} = 1$ ), which is equivalent to

$$
k + h ^ {\vee} = \frac {h ^ {\vee}}{u} \text {a n d} \operatorname {g c d} (u, h ^ {\vee}) = \operatorname {g c d} (u, r ^ {\vee}) = 1, \tag {6}
$$

the (normalized) character  $\mathrm{ch}_{\Lambda}$  turns into a product. The level  $k$ , defined by (6), is naturally called the boundary principal admissible level in [KRW03], see formula (3.5) there. We obtain from Proposition 1, (4) and (5)

Proposition 2. (a) All boundary principal admissible weights are of level  $k$ , given by (6), and are of the form

$$
\Lambda = (t _ {\beta} y). (k \Lambda_ {0}), \tag {7}
$$

where  $\beta \in Q^{*},y\in W$  are such that  $(t_{\beta}y)\hat{\Pi}_u\subset \hat{\Delta}_{+}$ . In particular,  $k\Lambda_0$  is a principal admissible weight of level (6).

(b) If  $\Lambda$  is of the form (7), then

$$
\mathrm {c h} _ {\Lambda} (\tau , z, t) = e ^ {2 \pi i (k t + \frac {h ^ {\vee}}{u} (z | \beta))} q ^ {\frac {h ^ {\vee}}{2 u} | \beta | ^ {2}} \left(\frac {\eta (u \tau)}{\eta (\tau)}\right) ^ {\frac {1}{2} (3 \ell - \dim \mathfrak {g})} \prod_ {\alpha \in \Delta_ {+}} \frac {\vartheta_ {1 1} (u \tau , y (\alpha) (z + \tau \beta))}{\vartheta_ {1 1} (\tau , \alpha (z))}.
$$

Remark 1. For the vacuum module  $L(k\Lambda_0)$  of the boundary principal admissible level  $k$  the character formula from Proposition 2(b) becomes

$$
\mathrm {c h} _ {k \Lambda_ {0}} (\tau , z, t) = e ^ {2 \pi i k t} \left(\frac {\eta (u \tau)}{\eta (\tau)}\right) ^ {\frac {1}{2} (3 \ell - \dim \mathfrak {g})} \prod_ {\alpha \in \Delta_ {+}} \frac {\vartheta_ {1 1} (u \tau , \alpha (z))}{\vartheta_ {1 1} (\tau , \alpha (z))}.
$$

Example 1. Let  $\mathfrak{g} = s\ell_2$ , so that  $h^\vee = 2$ . Then the boundary levels are  $k = \frac{2}{u} - 2$ , where  $u$  is a positive odd integer, and all admissible weights are

$$
\Lambda_ {k, j} := t _ {- \frac {j}{2} \alpha_ {1}}. (k \Lambda_ {0}) = (k + \frac {2 j}{u}) \Lambda_ {0} - \frac {2 j}{u} \Lambda_ {1}, j = 0, 1, \ldots , u - 1,
$$

and the character formula from Proposition 2(b) becomes:

$$
\mathrm {c h} _ {\Lambda_ {u, j}} = e ^ {2 \pi i (k t - \frac {j}{u} z)} q ^ {\frac {j ^ {2}}{2 u}} \frac {\vartheta_ {1 1} (u \tau , z - j \tau)}{\vartheta_ {1 1} (\tau , z)}. \tag {8}
$$

For  $u = 3$  and 5 some of these formulas were conjectured in [Y16].

Example 2. Let  $\mathfrak{g} = s\ell_N$ , so that  $h^\vee = N$ , let  $N > 1$  be odd, and let  $u = 2$ . Then the boundary admissible level is  $k = -\frac{N}{2}$ , and the boundary admissible weights of the form  $t_\beta.(k\Lambda_0)$  are:

$$
\Lambda_ {N, p} = - \frac {N}{2} \Lambda_ {p}, p = 0, 1, \dots ,, N - 1,
$$

where  $\Lambda_p$  are the fundamental weights of  $\widehat{\mathfrak{g}}$ . Letting  $z = \sum_{i=1}^{N-1} z_i \bar{\Lambda}_i$ , where  $\bar{\Lambda}_i$  are the fundamental weights of  $\mathfrak{g}$ , the character formula from Proposition 2 (b) becomes:

$$
\begin{array}{l} \mathrm {c h} _ {\Lambda_ {N, p}} (\tau , z, t) = i ^ {p (N - p)} e ^ {- \pi i N t} \left(\frac {\eta (2 \tau)}{\eta (\tau)}\right) ^ {- \frac {(N - 1) (N - 2)}{2}} \\ \times \frac{\prod\limits_{\substack{1\leq i\leq j <   p\\ \text{Or} p <   i\leq j <   N}}\vartheta_{11}(2\tau,z_{i} + \ldots +z_{j})\prod\limits_{1\leq i\leq p\leq j <   N}\vartheta_{01}(2\tau,z_{i} + \ldots +z_{j})}{\prod\limits_{1\leq i\leq j <   N}\vartheta_{11}(\tau,z_{i} + \ldots +z_{j})}, \\ \end{array}
$$

where

$$
\vartheta_ {0 1} (\tau , z) = \prod_ {n = 1} ^ {\infty} (1 - q ^ {n}) (1 - e ^ {2 \pi i z} q ^ {n - \frac {1}{2}}) (1 - e ^ {- 2 \pi i z} q ^ {n - \frac {1}{2}}).
$$

This follows from Proposition 2(b) by applying to  $\vartheta_{11}$  an elliptic transformation (see e.g. [KW14], Appendix). In particular

$$
\mathrm {c h} _ {- \frac {N}{2} \Lambda_ {0}} = e ^ {- \pi i N t} \left(\frac {\eta (2 \tau)}{\eta (\tau)}\right) ^ {- \frac {(N - 1) (N - 2)}{2}} \prod_ {1 \leq i \leq j <   N} \frac {\vartheta_ {1 1} (2 \tau , z _ {i} + \ldots + z _ {j})}{\vartheta_ {1 1} (\tau , z _ {i} + \ldots + z _ {j})}.
$$

The latter formula was conjectured in [XYY16].

Remark 2. For principal admissible weights  $\Lambda = (t_{\beta}y).(k\Lambda_0)$  and  $(t_{\beta^{\prime}}y^{\prime}).(k\Lambda_{0})$  of boundary level  $k = \frac{h^{\vee}}{u} - h^{\vee}$  the  $S$ -transformation matrix  $(a(\Lambda, \Lambda'))$ , given by [KW89], Theorem 3.6, simplifies to

$$
a (\Lambda , \Lambda^ {\prime}) = | Q / u h ^ {\vee} Q ^ {*} | ^ {- \frac {1}{2}} \varepsilon (y y ^ {\prime}) \prod_ {\alpha \in \Delta_ {+}} 2 \sin \frac {\pi i u (\rho | \alpha)}{h ^ {\vee}} e ^ {- 2 \pi i \left((\rho | \beta + \beta^ {\prime}) + \frac {h ^ {\vee} (\beta | \beta^ {\prime})}{u}\right)}.
$$

Remark 3. If  $\mathfrak{g} = s\ell_2$  and  $k$  is as in Example 1, then

$$
a (\Lambda_ {k, j}, \Lambda_ {k, j ^ {\prime}}) = (- 1) ^ {j + j ^ {\prime}} e ^ {- \frac {2 \pi i j j ^ {\prime}}{u}} \frac {1}{\sqrt {u}} \sin \frac {u \pi}{2}.
$$

One can compute fusion coefficients by Verlinde's formula:

$$
N _ {\Lambda_ {k, j _ {1}}, \Lambda_ {k, j _ {2}}, \Lambda_ {k, j _ {3}}} = (- 1) ^ {j _ {1} + j _ {2} + j _ {3}} \mathrm {i f} j _ {1} + j _ {2} + j _ {3} \in u \mathbb {Z}, \mathrm {a n d} = 0 \mathrm {o t h e r w i s e}.
$$

Example 3. Let  $\mathfrak{g} = sl_3$ , so that  $h^\vee = 3$ , and let  $u$  be a positive integer, coprime to 3. Then all (principal) admissible weights have level  $k = \frac{3}{u} - 3$  and are of the form (7), where

$$
\beta = - (- 1) ^ {p} (k _ {1} \bar {\Lambda} _ {1} + k _ {2} \bar {\Lambda} _ {2}), y = r _ {\theta} ^ {p}, p = 0 \mathrm {o r} 1, k _ {i} \in \mathbb {Z}, k _ {i} \geq \delta_ {p, 1}, k _ {1} + k _ {2} \leq u - \delta_ {p, 0}.
$$

Denote this weight by  $\Lambda_{u;k_1,k_2}^{(p)} = (t_\beta y).(k\Lambda_0)$ . Using Remark 2, one computes the fusion coefficients by Verlinde's formula:

$$
N _ {\Lambda_ {u; k _ {1}, k _ {2}} ^ {(p)} \Lambda_ {u; k _ {1} ^ {\prime}, k _ {2} ^ {\prime}} ^ {(p ^ {\prime})} \Lambda_ {u; k _ {1} ^ {\prime \prime}, k _ {2} ^ {\prime \prime}}} = (- 1) ^ {p + p ^ {\prime} + p ^ {\prime \prime}} \mathrm {i f} (- 1) ^ {p} k _ {i} + (- 1) ^ {p ^ {\prime}} k _ {i} ^ {\prime} + (- 1) ^ {p ^ {\prime \prime}} k _ {i} ^ {\prime \prime} \in u \mathbb {Z} \mathrm {f o r} i = 1, 2,
$$

and  $= 0$  otherwise.

Remark 4. If  $\Lambda$  is an arbitrary admissible weight, then  $\hat{\Delta}^{\Lambda}$  decomposes in a disjoint union of several affine root systems. Then  $\Lambda$  has boundary level if restrictions of it to each of them has boundary level, and formula (3.4) from [KW89] shows that  $\mathrm{ch}_{\Lambda}$  decomposes in a product of the corresponding boundary level characters. Note also that all the above holds also for twisted affine Kac-Moody algebras [KW89].

Remark 5. The product character formula for boundary level affine Kac-Moody superalgebras holds as well, see [GK15], formula (2).

Recall that to any  $sl_2$ -triple  $\{f,x,e\}$  in  $\mathfrak{g}$ , where  $[x,f] = -f$ ,  $[x,e] = e$ , one associates a  $W$ -algebra  $W^{k}(g,f)$ , obtained from the vacuum  $\widehat{\mathfrak{g}}$ -module of level  $k$  by quantum Hamiltonian reduction, so that any  $\widehat{\mathfrak{g}}$ -module  $L(\Lambda)$  of level  $k$  produces either an irreducible  $W^{k}(g,f)$ -module  $H(\Lambda)$  or zero. The characters of  $L(\Lambda)$  and  $H(\Lambda)$  are related by the following simple formula ([KRW03] or [KW14]):

$$
\left( \begin{array}{c} W \\ R \operatorname {c h} _ {H (\Lambda)} \end{array} \right) (\tau , z) = \left( \begin{array}{c} \hat {R} \operatorname {c h} _ {\Lambda} \end{array} \right) (\tau , - \tau x + z, \frac {\tau}{2} (x | x)). \tag {9}
$$

Here  $z\in \mathfrak{h}^f$  , the centralizer of  $f$  in  $\mathfrak{h}$  , and

$$
R (\tau , z) = \eta (\tau) ^ {\frac {3}{2} l - \frac {1}{2} \dim \left(\mathfrak {g} _ {0} + \mathfrak {g} _ {1 / 2}\right)} \prod_ {\alpha \in \Delta_ {+} ^ {0}} \vartheta_ {1 1} (\tau , \alpha (z)) \left(\prod_ {\alpha \in \Delta_ {1 / 2}} \vartheta_ {0 1} (\tau , \alpha (z))\right) ^ {1 / 2}, \tag {10}
$$

where  $\mathfrak{g} = \oplus_{j}\mathfrak{g}_{j}$  is the eigenspace decomposition for ad  $x$ ,  $\Delta_j \subset \Delta$  is the set of roots of root spaces in  $\mathfrak{g}_j$  and  $\Delta_+^0 = \Delta_+ \cap \Delta_0$  (we assume that  $\Delta_j \subset \Delta_+$  for  $j > 0$ ). If  $k$  is a boundary level (6), we obtain from Proposition 2(b) and formulas (9), (10) the following character formula for  $H(\Lambda)$  if  $\Lambda$  is a principal admissible weight (7) ( $z \in \mathfrak{h}^f$ ):

$$
\begin{array}{l} \mathrm {c h} _ {H (\Lambda)} (\tau , z) = (- i) ^ {| \Delta + |} q ^ {\frac {h ^ {\vee}}{2 u} | \beta - x | ^ {2}} e ^ {\frac {2 \pi i h ^ {\vee}}{u} (\beta | z)} \\ \times \frac {\eta (u \tau) ^ {\frac {3}{2} \ell - \frac {1}{2} \dim \mathfrak {g}}}{\eta (\tau) ^ {\frac {3}{2} \ell - \frac {1}{2} \dim \left(\mathfrak {g} _ {0} + \mathfrak {g} _ {1 / 2}\right)}} \frac {\prod_ {\alpha \in \Delta_ {+}} \vartheta_ {1 1} (u \tau , y (\alpha) (z + \tau \beta - \tau x))}{\prod_ {\alpha \in \Delta_ {+} ^ {0}} \vartheta_ {1 1} (\tau , \alpha (z)) \left(\prod_ {\alpha \in \Delta_ {1 / 2}} \vartheta_ {0 1} (\tau , \alpha (z))\right) ^ {1 / 2}}. \tag {11} \\ \end{array}
$$

Remark 6. A formula, similar to Proposition 2(b) and to formula (11), holds if  $\mathfrak{g}$  is a basic Lie superalgebra; one has to replace the character by the supercharacter, dim by  $sdim$ , and the factor  $\vartheta_{ab}$ , corresponding to a root  $\alpha$ , by its inverse if this root is odd. Also, the character is obtained from the supercharacter by replacing  $\vartheta_{ab}$  by  $\vartheta_{a,b+1 \bmod 2}$  if the root  $\alpha$  is odd.

Remark 7. An example of (11) is the minimal series representations of the Virasoro algebra with central charge  $c = 1 - \frac{3(u - 2)^2}{u}$ , obtained by the quantum Hamiltonian reduction from the boundary admissible  $\hat{sl}_2$ -modules from Example 1. For  $j = u - 1$  one gets 0, for  $u = 3$  and  $j = 0, 1$  one gets the trivial representation, but for all other  $j$  and  $u \geq 5$  the characters are the product sides of the Gordon generalizations of the Rogers-Ramanujan identities (the latter correspond to  $u = 5$ ). Another example is the minimal series representations of the  $N = 2$  superconformal algebras, see [KRW03], Section 7.

# References

[BLLPRV15] C. Beem, M. Lemos, P. Liendo, W. Peelaers, L. Rastelli, B.C. van Rees, Infinite chiral symmetry in four dimensions. Comm. Math. Phys. 336 (2015), no. 3, 13591433.  
[GK15] M. Gorelik, V. G. Kac, Characters of (relatively) integrable modules over affine Lie superalgebras, Jpn. J. Math. 10 (2015), No. 2, 135-235.  
[K90] V. G. Kac, Infinite-dimensional Lie algebras, Third edition, Cambridge University press, 1990.  
[KRW03] V. G. Kac, S.-S. Roan, M. Wakimoto, Quantum reduction of affine superalgebras, Comm. Math. Phys. 241 (2003), 307-342.  
[KW88] V. G. Kac, M. Wakimoto, Modular invariant representations of infinite-dimensional Lie algebras and superalgebras, Proc. Nat. Acad. Sci. USA 85 (1988), pp 4956-4960.  
[KW89] V. G. Kac, M. Wakimoto, Classification of modular invariant representations of affine algebras, Adv. Ser. Math. Phys. 7, World Sci. 1989, pp 138-177.  
[KW14] V. G. Kac, M.Wakimoto, Representations of affine superalgebras and mock theta functions, Transf. Groups 19 (2014), 387-455.  
[XYY16] D. Xie, W. Yan, S.-T. Yau, Chiral algebra of Argyres-Douglas theory from M5 brane, arXiv:1604.02155  
[Y16] W. Yan, Observations on characters of some Kac-Moody algebras, 2016 preprint.