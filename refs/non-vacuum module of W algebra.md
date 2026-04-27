# 1 Non-vacuum module of W algebra at boundary admissible level

# 1.1 Argyres–Douglas theories and associated VOAs

In general, an AD theory is characterized by its Higgs field

$$
\Phi (z) = \left(\frac {T _ {k}}{z ^ {2 + \frac {k}{b}}} + \sum_ {- b \leq l <   k} \frac {T _ {l}}{z ^ {2 + \frac {l}{b}}} + \dots\right) d z. \tag {1}
$$

These AD theories are expected to correspond to certain W algebras at boundary admissible levels. To make this correspondence precise, we introduce the rational number

$$
\nu = 1 + \frac {k}{b} = \frac {k + b}{b} = \frac {u}{m}, u \equiv \frac {k + b}{\operatorname * {g c d} (k + b , b)}, m \equiv \frac {b}{\operatorname * {g c d} (k + b , b)} \qquad \qquad (2)
$$

For a simple Lie algebra $\mathfrak { g }$ , when $m$ takes the values listed in Table 1 [], the corresponding AD theory is associated with a W-algebra at boundary admissible level.

Table 1: Regular number $m$ for simple Lie algebra   

<table><tr><td>g</td><td>m</td></tr><tr><td>An</td><td>n+1</td></tr><tr><td>Dn</td><td>m even, 2n-2m odd
m even, 2n/m even</td></tr><tr><td>E6</td><td>12, 9, 6, 3</td></tr><tr><td>E7</td><td>18, 14, 6, 2</td></tr><tr><td>E8</td><td>2, 3, 5, 6, 10, 15, 30
4, 8, 12, 24
20</td></tr></table>

In this case, the boundary admissible level is given by

$$
k _ {F} = - h ^ {\vee} + \frac {h ^ {\vee}}{k + b}, \quad (h ^ {\vee}, k + b) = 1. \tag {3}
$$

In addition to the irregular data, one may also introduce a regular puncture labeled by a nilpotent orbit $f$ . The resulting AD theory is therefore specified by

$$
(J ^ {b} [ k ], f). \tag {4}
$$

We first consider the AD theory with only an irregular puncture $J ^ { b } [ k ]$ . This theory corresponds to the affine Kac–Moody (AKM) algebra $L _ { - k _ { F } } ( { \mathfrak { g } } )$ at boundary admissible level. The boundary admissible weights are given by []

$$
\Lambda = \left(t _ {\beta} y\right). \left(k _ {F} \Lambda_ {0}\right), \tag {5}
$$

where $\beta \in Q ^ { * }$ and $y \in W$ satisfy $( t _ { \beta } y ) \hat { \Pi } _ { u } \subset \hat { \Delta } _ { + }$ for

$$
\hat {\Pi} _ {u} = \{u \delta - \theta , \hat {\alpha} _ {1}, \dots , \hat {\alpha} _ {r} \}. \tag {6}
$$

The character of the vacuum module is []

$$
\left(\frac {\eta (u \tau)}{\eta (\tau)}\right) ^ {\frac {1}{2} (3 r - \dim \mathfrak {g})} \prod_ {\alpha \in \Delta_ {+}} \frac {\vartheta_ {1} (\alpha (z) , u \tau)}{\vartheta_ {1} (\alpha (z) , \tau)} \tag {7}
$$

while the character of a non-vacuum module is []

$$
e ^ {2 \pi i \left(\frac {h ^ {\vee}}{u} (z | \beta)\right)} q ^ {\frac {h ^ {\vee}}{2 u} | \beta | ^ {2}} \left(\frac {\eta (u \tau)}{\eta (\tau)}\right) ^ {\frac {1}{2} (3 r - \dim \mathfrak {g})} \prod_ {\alpha \in \Delta_ {+}} \frac {\vartheta_ {1} (y (\alpha) (z + \tau \beta) , u \tau)}{\vartheta_ {1} (\alpha (z) , \tau)}. \tag {8}
$$

The conformal dimension of $L ( \Lambda )$ is

$$
h _ {\Lambda} = \frac {(\Lambda , \Lambda + 2 \hat {\rho})}{2 (\kappa + h ^ {\vee})}. \tag {9}
$$

Closing the regular puncture corresponds to performing the quantum Drinfeld–Sokolov reduction, which produces a $W$ -algebra $W ^ { k _ { F } } ( { \mathfrak { g } } , f ) \ [ ]$ . This construction requires choosing an ${ \mathfrak { s l } } _ { 2 }$ -triple $( x , e , f )$ in $\mathfrak { g }$ , guaranteed by the Jacobson–Morozov theorem [], satisfying

$$
[ x, e ] = e, [ x, f ] = - f, [ e, f ] = 2 x \tag {10}
$$

The Lie algebra then decomposes into eigenspaces of ad $x$

$$
\mathfrak {g} = \oplus \mathfrak {g} _ {j} \quad \mathfrak {g} _ {j} = \{[ x, g _ {j} ] = j g _ {j} \}. \tag {11}
$$

Under the quantum Drinfeld–Sokolov reduction, an admissible AKM module is mapped either to zero or to an admissible W algebra module [1, 2]

$$
\Psi^ {\pm}: V _ {k _ {F}} (\mathfrak {g}) - \operatorname {m o d} \rightarrow W _ {k _ {F}} (\mathfrak {g}, f) - \operatorname {m o d}. \tag {12}
$$

In fact, there are two versions of the quantum reduction for AKM modules, denoted by $\Psi ^ { \pm }$ The $\Psi ^ { + }$ reduction yields the standard VOA; we denote the corresponding $\mathcal { W }$ -algebra by chW .

The vacuum character chW of the W algebra is given by

$$
(- i) ^ {| \Delta + |} q ^ {\frac {h ^ {\vee}}{2 u} | x | ^ {2}} \frac {\eta (u \tau) ^ {\frac {3}{2} r} - \frac {1}{2} \dim \mathfrak {g}}{\eta (\tau) ^ {\frac {3}{2} l - \frac {1}{2} \dim \left(\mathfrak {g} _ {0} + \mathfrak {g} _ {\frac {1}{2}}\right)}} \frac {\prod_ {\alpha \in \Delta_ {+}} - \vartheta_ {1} (\alpha (z - \tau x) , u \tau)}{\prod_ {\alpha \in \Delta_ {+} ^ {0}} - \vartheta_ {1} (\alpha (z) , \tau) \left(\prod_ {\alpha \in \Delta_ {\frac {1}{2}}} \vartheta_ {4} (\alpha (z) , \tau)\right) ^ {\frac {1}{2}}} \tag {13}
$$

where $\Delta _ { + } ^ { 0 } = \Delta _ { + } \cap \Delta _ { 0 }$ . The character of a non-vacuum W algebra module takes the form

$$
\begin{array}{l} \mathrm {c h} _ {W} = (- i) ^ {| \Delta_ {+} |} q ^ {\frac {h ^ {\vee}}{2 u} | \beta - x | ^ {2}} e ^ {\frac {2 \pi i h ^ {\vee}}{u} (\beta | z)} \\ \times \frac {\eta (u \tau) ^ {\frac {3}{2} r - \frac {1}{2} \dim (\mathfrak {g} _ {0} + \mathfrak {g} _ {\frac {1}{2}})}}{\eta (\tau) ^ {\frac {3}{2} r - \frac {1}{2} \dim (\mathfrak {g} _ {0} + \mathfrak {g} _ {\frac {1}{2}})}} \frac {\prod_ {\alpha \in \Delta_ {+}} - \vartheta_ {1} (y (\alpha) (z + \tau \beta - \tau x) , u \tau)}{\prod_ {\alpha \in \Delta_ {+} ^ {0}} - \vartheta_ {1} (\alpha (z) , \tau) \left(\prod_ {\alpha \in \Delta_ {\frac {1}{2}}} \vartheta_ {4} (\alpha (z) , \tau)\right) ^ {\frac {1}{2}}}. \tag {14} \\ \end{array}
$$

Finally, the characters of the AKM algebra and the associated $W$ -algebra are related by

$$
(R _ {W} \mathrm {c h} _ {\Psi^ {+} (\Lambda)}) (\tau , z) = (\hat {R} \mathrm {c h} _ {L (\Lambda)}) (\tau , - \tau x + z, \frac {\tau}{2} (x | x)), \tag {15}
$$

where $H ( \Lambda )$ and $L ( \Lambda )$ denote W algebra and AKM modules, respectively. The prefactors are

$$
R _ {W} (\tau , z) = \eta (\tau) ^ {\frac {3}{2} l - \frac {1}{2} \dim \left(\mathfrak {g} _ {0} + \mathfrak {g} _ {1 / 2}\right)} \prod_ {\alpha \in \Delta_ {+} ^ {0}} - \vartheta_ {1} (\tau , \alpha (z)) \left(\prod_ {\alpha \in \Delta_ {1 / 2}} \vartheta_ {2} (\tau , \alpha (z))\right) ^ {1 / 2}, \tag {16}
$$

$$
\hat {R} (\tau , z) = (- i) ^ {| \Delta_ {+} |} \eta (\tau) ^ {\frac {1}{2} (3 \ell - \dim \mathfrak {g})} \prod_ {\alpha \in \Delta_ {+}} - \vartheta_ {1} (\tau , \alpha (z)).
$$

There are situations in which two distinct AKM modules reduce to the same W algebra module when they are related by the action of a certain Weyl group [].

# 1.2 Modularity of W algebra

In this subsection, we review the modular properties of the $\mathcal { W }$ algebra at boundary admissible level. The general form of the character of a module associated with a W algebra or an AKM algebra is

$$
\mathcal {I} = \operatorname {T r} _ {\Lambda} q ^ {L _ {0} + c _ {4 d} / 2} \mathbf {y} ^ {\alpha}. \tag {17}
$$

Here the lowest eigenvalue of $L _ { 0 }$ defines the conformal dimension $h _ { \Lambda }$ , and $c _ { 4 d }$ denotes the central charge of the corresponding four-dimensional theory. For a W algebra at boundary admissible level, the associated four-dimensional central charge is given by

$$
\frac {- 1}{1 2} \left(\dim \mathfrak {g} _ {0} - \frac {1}{2} \dim \mathfrak {g} _ {\frac {1}{2}} - \frac {1 2}{\kappa + h ^ {\vee}} | \rho - (k + h ^ {\vee}) x | ^ {2}\right), \quad \kappa = - h ^ {\vee} + \frac {h ^ {\vee}}{h ^ {\vee} + k}. \tag {18}
$$

Let us first consider the AKM case. For a boundary admissible weight $\Lambda$ , the conformal dimension is

$$
h _ {\Lambda} = \frac {(\Lambda , \Lambda + 2 \hat {\rho})}{2 (\kappa + h ^ {\vee})}. \tag {19}
$$

Consequently, the leading power in the $q$ -expansion of the Schur index is

$$
h _ {\Lambda} + \frac {c _ {4 d}}{2}. \tag {20}
$$

The modular T and $\mathbb { S }$ matrices for two VOA weights $\Lambda = ( t _ { \beta } y ) \cdot ( \kappa \Lambda _ { 0 } )$ and $\Lambda ^ { \prime } = ( t _ { \beta ^ { \prime } } y ^ { \prime } ) \cdot ( \kappa \Lambda _ { 0 } )$ are given by

$$
\mathbb {T} _ {\Lambda , \Lambda^ {\prime}} = e ^ {2 \pi i \left(h _ {\Lambda} - \frac {c}{2 4}\right)} \delta_ {\Lambda , \Lambda^ {\prime}},
$$

$$
\mathbb {S} _ {\Lambda , \Lambda^ {\prime}} = \left| \frac {P ^ {\vee}}{u h ^ {\vee} Q ^ {\vee}} \right| ^ {- \frac {1}{2}} \epsilon \left(y y ^ {\prime}\right) \prod_ {\alpha \in \Delta_ {+}} \left(2 \sin \frac {\pi u (\rho , \alpha)}{h ^ {\vee}}\right) e ^ {- 2 \pi i \left(\left(\rho , \beta + \beta^ {\prime}\right) + \frac {h ^ {\vee} \left(\beta , \beta^ {\prime}\right)}{u}\right)}. \tag {21}
$$

The modular behavior of the character of the $\Psi _ { f } ^ { - }$ -reduction was studied in [3–5]. It is known to work in the case of a good even grading, that is, when all eigenvalues of $\operatorname { a d } x$ are even. In this setting, one can conjugate the associated ${ \mathfrak { s l } } _ { 2 }$ -triple to another triple such that the resulting element $f ^ { \prime }$ becomes a regular nilpotent element in a Levi subalgebra ${ \mathfrak { l } } \subset { \mathfrak { g } }$ .

In this paper, we restrict our attention to the $\mathcal { W }$ -algebras obtained via the $\Psi _ { f } ^ { + }$ reduction. To discuss the modularity of $\mathcal { W }$ -algebra characters, we first introduce the notation of [6]. Define

$$
R _ {W} ^ {-} (\tau , z) \equiv R _ {W} (\tau , z + x),
$$

$$
R _ {W} ^ {*} (\tau , z) \equiv R _ {W} (\tau , z + \tau x + x). \tag {22}
$$

The modular transformation properties of $R _ { W }$ , $R _ { W } ^ { - }$ , and $R _ { W } ^ { * }$ are as follows:

Observe that

$$
(\hat {R} \mathrm {c h} _ {\Lambda}) (\tau , z, t) = (i) ^ {| \Delta_ {+} |} e ^ {2 \pi i (h ^ {\vee} + k) t} e ^ {2 \pi i \frac {h ^ {\vee}}{u} (z | \beta)} q ^ {\frac {h ^ {\vee}}{2 u} | \beta | ^ {2}} \eta (u \tau) ^ {\frac {1}{2} (3 \ell - \dim \mathfrak {g})} \prod_ {\alpha \in \Delta_ {+}} \vartheta_ {1} (u \tau , y (\alpha) (z + \tau \beta)), \tag {23}
$$

The modular property of $\hat { R } \mathrm { c h } _ { \Lambda }$ is

Motivated by the modular properties of RˆchΛ together with those of $R _ { W }$ , $R _ { W } ^ { - }$ , and $R _ { W } ^ { * }$ , we expect that a natural basis for the $\mathcal { W }$ -algebra characters is given by

$$
\mathrm {c h} _ {W}, \frac {\hat {R} \mathrm {c h} _ {L (\Lambda)} (\tau , - \tau x + x + z , \frac {\tau}{2} | x | ^ {2})}{R _ {W} ^ {-}}, \frac {\hat {R} \mathrm {c h} _ {L (\Lambda)} (\tau , x + z , 0)}{R _ {W} ^ {*}}. (2 4)
$$

# 1.3 Example: sl3 at level $\textstyle k = - 3 + { \frac { 3 } { 4 } }$

We now illustrate the example of ${ \mathfrak { g } } = { \mathfrak { s l } } _ { 3 }$ at level $\begin{array} { r } { k = - 3 + \frac { 3 } { 4 } } \end{array}$ . In this case, the vacuum character of the affine Kac–Moody algebra $L _ { - 3 + \frac { 3 } { 4 } } ( \mathfrak { s l } _ { 3 } )$ is given by

$$
\left(\frac {\eta (u \tau)}{\eta (\tau)}\right) ^ {- 1} \prod_ {\alpha \in \Delta_ {+}} \frac {\vartheta_ {1} (\alpha (z) , u \tau)}{\vartheta_ {1} (\alpha (z) , \tau)} \tag {25}
$$

The total number of simple modules is $u ^ { r } = 1 6$ .

The corresponding boundary admissible weights are summarized in Table 2 [].

Table 2: Admissible weight for $L _ { - 3 + \frac { 3 } { 4 } } ( \mathfrak { s l } _ { 3 } )$   

<table><tr><td>[tβy]</td><td>Λ</td><td>[tβy]</td><td>Λ</td></tr><tr><td>1</td><td>-9/4 Λ0</td><td>t-ω2</td><td>-3/2 Λ0 - 3/4 Λ2</td></tr><tr><td>t-2ω2</td><td>-3/4 Λ0 - 3/2 Λ2</td><td>t-3ω2</td><td>-9/4 Λ2</td></tr><tr><td>t-ω1</td><td>-3/2 Λ0 - 3/4 Λ1</td><td>t-ω1-ω2</td><td>-3/4 Λ0 - 3/4 Λ1 - 3/4 Λ2</td></tr><tr><td>t-ω1-2ω2</td><td>-3/4 Λ1 - 3/2 Λ2</td><td>t-2ω1</td><td>-3/4 Λ0 - 3/2 Λ1</td></tr><tr><td>t-2ω1-ω2</td><td>-3/2 Λ1 - 3/4 Λ2</td><td>t-3ω1</td><td>-9/4 Λ1</td></tr><tr><td>tω1+ω2sθ</td><td>1/4 Λ0 - 5/4 Λ1 - 5/4 Λ2</td><td>tω1+2ω2sθ</td><td>-1/2 Λ0 - 5/4 Λ1 - 1/2 Λ2</td></tr><tr><td>tω1+3ω2sθ</td><td>-5/4 Λ0 - 5/4 Λ1 + 1/4 Λ2</td><td>t2ω1+ω2sθ</td><td>-1/2 Λ0 - 1/2 Λ1 - 5/4 Λ2</td></tr><tr><td>t2ω1+2ω2sθ</td><td>-5/4 Λ0 - 1/2 Λ1 - 1/2 Λ2</td><td>t3ω1+ω2sθ</td><td>-5/4 Λ0 + 1/4 Λ1 - 5/4 Λ2</td></tr></table>

Using these data, the lowest powers appearing in the $q$ -expansion are

$$
h _ {\Lambda} + \frac {2}{2} = \left\{1, \frac {1}{4}, 0, \frac {1}{4}, \frac {1}{4}, - \frac {1}{4}, - \frac {1}{4}, 0, - \frac {1}{4}, \frac {1}{4}, - \frac {1}{4}, - \frac {1}{4}, \frac {1}{4}, - \frac {1}{4}, 0, \frac {1}{4} \right\}, \tag {26}
$$

which precisely matches the lowest dimensions appearing in the expansion of the corresponding Schur index. For this case, we can explicitly determine the modular $T$ and $S$ matrices, which

are given by

$$
T = \left( \begin{array}{c c c c c c c c c c c c c c c} 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & i & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & i & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & i & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & - i & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & - i & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & - i & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & i & 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & - i & 0 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & - i \end{array} \right) \tag {27}
$$

and

$$
S = \left( \begin{array}{c c c c c c c c c c c c c c c} \frac {1}{4} & \frac {1}{4} & \frac {1}{4} & \frac {1}{4} & \frac {1}{4} & \frac {1}{4} & \frac {1}{4} & \frac {1}{4} & \frac {1}{4} & \frac {1}{4} & - \frac {1}{4} & - \frac {1}{4} & - \frac {1}{4} & - \frac {1}{4} & - \frac {1}{4} & - \frac {1}{4} \\ \frac {1}{4} & - \frac {1}{4} & \frac {1}{4} & - \frac {1}{4} & - \frac {i}{4} & \frac {i}{4} & - \frac {i}{4} & - \frac {1}{4} & \frac {1}{4} & \frac {i}{4} & \frac {i}{4} & - \frac {i}{4} & \frac {i}{4} & - \frac {1}{4} & \frac {1}{4} & - \frac {i}{4} \\ \frac {1}{4} & \frac {1}{4} & \frac {1}{4} & \frac {1}{4} & - \frac {1}{4} & - \frac {1}{4} & - \frac {1}{4} & \frac {1}{4} & \frac {1}{4} & - \frac {1}{4} & \frac {1}{4} & \frac {1}{4} & - \frac {1}{4} & - \frac {1}{4} & \frac {1}{4} \\ \frac {1}{4} & - \frac {1}{4} & \frac {1}{4} & - \frac {1}{4} & \frac {i}{4} & - \frac {i}{4} & \frac {i}{4} & - \frac {1}{4} & \frac {1}{4} & - \frac {i}{4} & - \frac {i}{4} & \frac {i}{4} & - \frac {i}{4} & - \frac {1}{4} \\ \frac {1}{4} & - \frac {i}{4} & - \frac {1}{4} & \frac {i}{4} & - \frac {1}{4} & \frac {i}{4} & \frac {1}{4} & \frac {1}{4} & - \frac {i}{4} & - \frac {1}{4} & \frac {i}{4} & - \frac {1}{4} & - \frac {i}{4} & - \frac {i}{4} \\ \frac {1}{4} & i / 2 = 0. 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 0. 0 p t \\ \frac {1}{4} & - i / 2 = 0. 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8. \\ \frac {1}{4} & - i / 2 = i / (2) \\ \frac {1}{4} & i / (2) \\ - i / (2) & i / (2) \\ - i / (2) & i / (2) \\ - i / (2) & i / (2) \\ - i / (2) & i / (2) \\ - i / (2) & i / (2) \\ - i / (2) & i / (2) \\ - i / (2) & i / (2) \\ - i / (2) & i / (2) \\ - [ i, j ] = [ i, j ] = [ i, j ] = [ i, j ] = [ i, j ] = [ i, j ] = [ i, j ] = [ i, j ] = [ i, j ] = [ i, j ] = [ i, j ] = [ i, j ] = [ i, j ] = [ i, j ] = [ i, j ] = [ i, j ] = [ i, j ] = [ i, k ] = [ i, k ] = [ i, k ] = [ i, k ] = [ i, k ] = [ i, k ] = [ i, k ] = [ i, k ] = [ i, k ] = [ i, k ] = [ i, k ] = [ i, k ] = [ i, k ] = [ i, k ] = [ i, k ] = [ i, k ] = [ i, k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = [ k ] = \( {\left[ {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left| {\left|     .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .      } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   } +   }\right) , \\ S = (28)
$$

which is exactly the one given by the formula (21), as expected.

We next introduce a regular puncture corresponding to the nilpotent orbit $f = \lfloor 2 , 1 \rfloor$ . After performing the quantum Drinfeld–Sokolov reduction, the admissible weights of the resulting $W$ -algebra are listed in Table 3.

Among these modules, we identify the cases that are isomorphic.

From the $W$ -algebra character formula (14), one can readily determine the non-vacuum modules of the W -algebra,

$$
y (\alpha) (z) \in \mathbb {Z} \text {a n d} y (\alpha) (\beta - x) = \mathbb {Z} u. \tag {29}
$$

Table 3: Admissible weights for $W ^ { - 3 + \frac { 3 } { 4 } } ( \mathfrak { s l } _ { 3 } , [ 2 , 1 ] )$   

<table><tr><td>[tβy]</td><td>Λ</td><td>[tβy]</td><td>Λ</td></tr><tr><td>1</td><td>-9/4 Λ0</td><td>t-ω2</td><td>-3/2 Λ0 - 3/4 Λ2</td></tr><tr><td>t-2ω2</td><td>-3/4 Λ0 - 3/2 Λ2</td><td>t-ω1</td><td>-3/2 Λ0 - 3/4 Λ1</td></tr><tr><td>t-ω1-ω2</td><td>-3/4 Λ0 - 3/4 Λ1 - 3/4 Λ2</td><td>t-2ω1</td><td>-3/4 Λ0 - 3/2 Λ1</td></tr><tr><td>tω1+ω2sθ</td><td>1/4 Λ0 - 5/4 Λ1 - 5/4 Λ2</td><td>tω1+2ω2sθ</td><td>-1/2 Λ0 - 5/4 Λ1 - 1/2 Λ2</td></tr><tr><td>tω1+3ω2sθ</td><td>-5/4 Λ0 - 5/4 Λ1 + 1/4 Λ2</td><td>t2ω1+ω2sθ</td><td>-1/2 Λ0 - 1/2 Λ1 - 5/4 Λ2</td></tr><tr><td>t2ω1+2ω2sθ</td><td>-5/4 Λ0 - 1/2 Λ1 - 1/2 Λ2</td><td>t3ω1+ω2sθ</td><td>-5/4 Λ0 + 1/4 Λ1 - 5/4 Λ2</td></tr></table>

<table><tr><td>tβ</td><td>1</td><td>t-ω2</td><td>t-2ω2</td><td>t-ω1</td><td>t-ω1-ω2</td><td>t-2ω1</td></tr><tr><td>tβsθ</td><td>tω1+ω2sθ</td><td>t2ω1+ω2sθ</td><td>t3ω1+ω2sθ</td><td>tω1+2ω2sθ</td><td>t2ω1+2ω2sθ</td><td>tω1+3ω2sθ</td></tr></table>

Here the indepent modules give the six characters

$$
\begin{array}{l} \mathrm {c h} _ {- \frac {9}{4} \Lambda_ {0}} = - \frac {i q ^ {3 / 1 6} \vartheta_ {1} (- \tau , q ^ {4}) \vartheta_ {1} \left(- \frac {\tau}{2} - 3 z _ {1} , q ^ {4}\right) \vartheta_ {1} \left(3 z _ {1} - \frac {\tau}{2} , q ^ {4}\right)}{\eta \eta (q) \eta (q ^ {4}) \vartheta_ {4} (3 z _ {1} , q)}, \\ \mathrm {c h} _ {- \frac {3}{2} \Lambda_ {0} - \frac {3}{4} \Lambda_ {2}} = - \frac {i q ^ {1 3 / 1 6} e ^ {\frac {3}{2} i \pi z _ {1}} \vartheta_ {1} (- 2 \tau , q ^ {4}) \vartheta_ {1} \left(- \frac {3 \tau}{2} - 3 z _ {1} , q ^ {4}\right) \vartheta_ {1} \left(3 z _ {1} - \frac {\tau}{2} , q ^ {4}\right)}{\eta \eta (q) \eta (q ^ {4}) \vartheta_ {4} (3 z _ {1} , q)}, \\ \operatorname {c h} _ {- \frac {3}{4} \Lambda_ {0} - \frac {3}{2} \Lambda_ {2}} = - \frac {i q ^ {3 1 / 1 6} e ^ {3 i \pi z _ {1}} \vartheta_ {1} (- 3 \tau , q ^ {4}) \vartheta_ {1} \left(- \frac {5 \tau}{2} - 3 z _ {1} , q ^ {4}\right) \vartheta_ {1} \left(3 z _ {1} - \frac {\tau}{2} , q ^ {4}\right)}{\eta \eta (q) \eta \left(q ^ {4}\right) \vartheta_ {4} \left(3 z _ {1} , q\right)}, \tag {30} \\ \operatorname {c h} _ {- \frac {3}{2} \Lambda_ {0} - \frac {3}{4} \Lambda_ {1}} = - \frac {i q ^ {1 3 / 1 6} e ^ {- \frac {3}{2} i \pi z _ {1}} \vartheta_ {1} (- 2 \tau , q ^ {4}) \vartheta_ {1} \left(- \frac {\tau}{2} - 3 z _ {1} , q ^ {4}\right) \vartheta_ {1} \left(3 z _ {1} - \frac {3 \tau}{2} , q ^ {4}\right)}{\eta \eta (q) \eta (q ^ {4}) \vartheta_ {4} (3 z _ {1} , q)}, \\ \operatorname {c h} _ {- \frac {3}{4} \Lambda_ {0} - \frac {3}{4} \Lambda_ {1} - \frac {3}{4} \Lambda_ {2}} = - \frac {i q ^ {2 7 / 1 6} \vartheta_ {1} (- 3 \tau , q ^ {4}) \vartheta_ {1} \left(- \frac {3 \tau}{2} - 3 z _ {1} , q ^ {4}\right) \vartheta_ {1} \left(3 z _ {1} - \frac {3 \tau}{2} , q ^ {4}\right)}{\eta   \eta (q)   \eta (q ^ {4})   \vartheta_ {4} (3 z _ {1} , q)}, \\ \mathrm {c h} _ {- \frac {3}{4} \Lambda_ {0} - \frac {3}{2} \Lambda_ {1}} = - \frac {i q ^ {3 1 / 1 6} e ^ {- 3 i \pi z _ {1}} \vartheta_ {1} (- 3 \tau , q ^ {4}) \vartheta_ {1} \left(- \frac {\tau}{2} - 3 z _ {1} , q ^ {4}\right) \vartheta_ {1} \left(3 z _ {1} - \frac {5 \tau}{2} , q ^ {4}\right)}{\eta \eta (q) \eta (q ^ {4}) \vartheta_ {4} (3 z _ {1} , q)}. \\ \end{array}
$$

However, these six characters do not form the complete basis of the $S L ( 2 , \mathbb { Z } )$ orbit of the Schur index since the expression involves the $\vartheta _ { 4 }$ theta function. Hence to characterize the full modular orbit of Schur index, we should introduce certain new characters besides these boundary admissible weight.

# 1.4 Defect index and non-vacuum modules

# 1.5 AKM modules from spetral flow

# References

[1] Edward Frenkel, Victor Kac, and Minoru Wakimoto. Characters and fusion rules for W algebras via quantized Drinfeld-Sokolov reductions. Commun. Math. Phys., 147:295–328, 1992.

[2] Victor G. Kac, Shi-shyr Roan, and Minoru Wakimoto. Quantum Reduction for Affine Superalgebras. Commun. Math. Phys., 241(2-3):307–342, 2003.   
[3] Victor G Kac and Minoru Wakimoto. On rationality of w-algebras. Transformation Groups, 13(3):671–713, 2008.   
[4] Tomoyuki Arakawa and Jethro van Ekeren. Rationality and fusion rules of exceptional $\mathcal { W }$ -algebras. J. Eur. Math. Soc., 25(7):2763–2813, 2022.   
[5] Tomoyuki Arakawa, Igor Alarcon Blatt, Jethro van Ekeren, and Wenbin Yan. Characters and fusion rules of boundary W-algebras. 9 2025.   
[6] Victor G. Kac and Minoru Wakimoto. On Modular Invariance of Quantum Affine W-Algebras. Commun. Math. Phys., 406(2):44, 2025.