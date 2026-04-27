# Quantum Reduction for Affine Superalgebras

Victor Kac

Department of Mathematics, M. I. T.

Cambridge, MA 02139, USA

(e-mail: kac@math.mit.edu )

Shi-shyr Roan

Institute of Mathematics , Academia Sinica

Taipei , Taiwan

(e-mail: maroan@ccvax.sinica.edu.tw)

Minoru Wakimoto

Graduate School of Mathematics , Kyushu University

Fukuoka, 812-8581, Japan

(e-mail:wakimoto@math.kyushu-u.ac.jp)

# Abstract

We extend the homological method of quantization of generalized Drinfeld–Sokolov reductions to affine superalgebras. This leads, in particular, to a unified representation theory of superconformal algebras.

1991 MSC: 17B65, 17B67, 81R10

1990 PACS: 02.20, 03.65.F, 89

# 0 Introduction

A series of papers on $W$ -algebras written in the second half of the 1980’s and the early 1990’s (see [BS]) culminated in the work of Feigin and Frenkel [FF1, FF2] who showed that to a simple finitedimensional Lie algebra $\mathfrak { g }$ one canonically associates a $W$ -algebra $W _ { k } ( { \mathfrak { g } } )$ as a result of quantization of the classical Drinfeld–Sokolov reduction. Namely, $W _ { k } ( { \mathfrak { g } } )$ is realized as homology of a BRST complex involving the principal nilpotent element of $\mathfrak { g }$ (i.e., the nilpotent element the closure of whose orbit contains all other nilpotent elements), the universal enveloping algebra of the affine Kac–Moody algebra $\widehat { \mathfrak { g } }$ associated to $\mathfrak { g }$ , and the charged fermionic ghosts associated to the currents of a maximal nilpotent subalgebra $\mathfrak { n }$ of $\mathfrak { g }$ .

This approach allows one not only to define the $W$ -algebras, but also to construct a functor $H$ from the category of restricted $\widehat { \mathfrak { g } }$ -modules of level $k$ to the category of positive energy modules over $W _ { k } ( { \mathfrak { g } } )$ . Namely, the $W _ { k } ( { \mathfrak { g } } )$ -module corresponding to a $\widehat { \mathfrak { g } }$ -module is the homology $H ( M )$ of the BRST complex associated to $M$ . This functor was applied in [FKW] to the admissible $\widehat { \mathfrak { g } }$ - modules, classified in [KW1], [KW2], in order to compute the characters of $W _ { k } ( { \mathfrak { g } } )$ -modules. (In the

simplest case of ${ \mathfrak { g } } = s \ell _ { 2 }$ one recovers thereby the minimal series modules over the Virasoro algebra $= W _ { k } ( s \ell _ { 2 } )$ .)

It is straightforward to generalize this construction to the case when $f$ is an even nilpotent element, that is for the $s \ell _ { 2 }$ -triple $\langle e , x , f \rangle$ , such that $[ e , f ] = x$ , $[ x , e ] = e$ , $[ x , f ] = - f$ , all eigenvalues of ad $x$ are integers (for general $f$ they lie in $\scriptstyle { \frac { 1 } { 2 } } \mathbf { Z }$ ). One just takes instead of $\mathbf { n }$ the subalgebra $^ { 9 + }$ of $\mathfrak { g }$ spanned by eigenspaces with positive eigenvalues for ad $x$ . Unfortunately, most nilpotent elements are not even, but often one can replace $x$ by $x ^ { \prime }$ such that ad $x ^ { \prime }$ has integer eigenvalues, so that the construction gives the same homology (see e.g. [BT]). However, it remained unclear how to make it work for a general simple Lie algebra $\mathfrak { g }$ and a general nilpotent element $f$ . The situation gets worse if one tries to go to the Lie superalgebra case since already the simplest Lie superalgebra $s p o ( 2 | 1 )$ has no good $\mathbf { Z }$ - gradations.

In the present paper we show how to resolve this problem. It turns out that one needs only to add neutral fermionic ghosts associated to the currents of the eigenspace ${ \mathfrak { g } } _ { 1 / 2 }$ of ad $x$ .

This is done in Section 2, where to each quadruple $( { \mathfrak { g } } , x , f , k )$ , where $^ { 9 }$ is a simple finitedimensional Lie superalgebra with a fixed even invariant bilinear form (.|.), $x$ is an ad-diagonizable element of $\mathfrak { g }$ with eigenvalues in $\scriptstyle { \frac { 1 } { 2 } } \mathbf { Z }$ , $f$ is a nilpotent even element of $\mathfrak { g }$ such that $[ x , f ] = - f$ , and $k \in \mathbf { C }$ , we associate a BRST complex

$$
\left(\mathcal {C} (\mathfrak {g}, x, f, k) = V _ {k} (\mathfrak {g}) \otimes F ^ {\mathrm {c h}} \otimes F ^ {\mathrm {n e}}, d _ {0}\right).
$$

Here $V _ { k } ( { \mathfrak { g } } )$ is the universal affine vertex algebra of level $k$ associated to $\widehat { \mathfrak { g } }$ , $F ^ { \mathrm { c h } }$ is the vertex algebra of free charged fermions based on ${ \mathfrak { g } } _ { + } + { \mathfrak { g } } _ { + } ^ { * }$ with reversed parity, $F ^ { \mathrm { n e } }$ is the vertex algebra of free neutral fermions based on ${ \mathfrak { g } } _ { 1 / 2 }$ with the form $\langle a , b \rangle = ( f | [ a , b ] )$ , and $d _ { 0 }$ is an explicitly constructed odd derivation of the vertex algebra $\mathcal { C } ( \mathfrak { g } , x , f , k )$ whose square is 0 (see Section 2.1). The main object of our study is the $0 ^ { \mathrm { t h } }$ homology of this complex, which is a vertex algebra, denoted by $W _ { k } ( { \mathfrak { g } } , x , f )$ . In the case when the pair $( x , f )$ can be included in an $s \ell _ { 2 }$ -triple $( e , x , f )$ (then $x$ is determined by $f$ up to conjugation), we denote this vertex algebra by $W _ { k } ( { \mathfrak { g } } , f )$ . In this case the map $\operatorname { a d } f : { \mathfrak { g } } _ { 1 / 2 } \to { \mathfrak { g } } _ { - 1 / 2 }$ is an isomorphism, which suffices for the construction of the energy-momentum field $L ( z )$ of $W _ { k } ( { \mathfrak { g } } , x , f )$ (see Section 2.2); under the same assumption, we construct fields $J ^ { \{ v \} }$ in $W _ { k } ( { \mathfrak { g } } , x , f )$ of conformal weight 1, corresponding to each element $v \in { \mathfrak { g } } ^ { x , f }$ , the centralizer of $x$ and $f$ (see Section 2.4).

As in [FF2, FKW], given a restricted $\widehat { \mathfrak { g } }$ -module $M$ of level $k$ , hence a $V _ { k } ( { \mathfrak { g } } )$ -module, we extend it to a $\mathcal { C } ( \mathfrak { g } , x , f , k )$ -module $\mathcal { C } ( M ) = M \otimes F ^ { \mathrm { c h } } \otimes F ^ { \mathrm { n e } }$ , which gives rise to a complex $( { \mathcal { C } } ( M ) , d _ { 0 } )$ of $\mathcal { C } ( \mathfrak { g } , x , f , k )$ -modules. Its homology $H ( M )$ is a $W _ { k } ( { \mathfrak { g } } , x , f )$ -module. In Section 3.1 we compute the Euler–Poincar´e character of this module:

$$
\operatorname {c h} _ {H (M)} (h) = \sum_ {j \in {\bf Z}} (- 1) ^ {j} \mathrm {t r} _ {H _ {j} (M)} q ^ {L _ {0}} e ^ {2 \pi i J _ {0} ^ {\{h \}}},
$$

where $h$ is an element of a Cartan subalgebra of ${ \mathfrak { g } } ^ { x , f }$ and $J ^ { \{ h \} }$ is the corresponding field of $W _ { k } ( { \mathfrak { g } } , x , f )$ . Furthermore, in Section 3.2 we find necessary and sufficient conditions on the $\widehat { \mathfrak { g } }$ - module $M$ for the non-vanishing of $\operatorname { c h } _ { H ( M ) }$ . The $\widehat { \mathfrak { g } }$ -modules $M$ satisfying these conditions are called non-degenerate.

In Section 3.3 we recall the definition of admissible highest weight $\widehat { \mathfrak { g } }$ -modules $L ( \Lambda )$ in the Lie superalgebra case [KW4]. The characters of these modules in the Lie algebra case were computed in [KW1]. Unfortunately we do not know how to prove an analogous character formula even in its weaker form in the Lie superalgebra case. This character formula is our first fundamental conjecture (which is confirmed by many examples in [KW1], [KW2], [KW4]). The second fundamental conjecture states that the $W _ { k } ( { \mathfrak { g } } , x , f )$ -module $H ( M )$ is either zero or irreducible, provided that $( x , f )$ is a

“good” pair and $M$ is an admissible highest weight $\widehat { \mathfrak { g } }$ -module. Of course, these conjectures allow us to compute the characters of irreducible $W _ { k } ( { \mathfrak { g } } , x , f )$ -modules $H ( M )$ for non-degenerate admissible $\widehat { \mathfrak { g } }$ -modules, using the results of Section 3.1.

In Section 4 we study the vertex algebra $W _ { k } ( { \mathfrak { g } } , f )$ in the case a “minimal” nilpotent even element $f$ , namely when $f$ is a root vector corresponding to an even highest root of $\mathfrak { g }$ . These vertex algebras were considered from a quite different viewpoint in [FL], and they include all well known superconformal algebras, like the $N \leq 4$ superconformal algebras and the big $N = 4$ superconformal algebras.

In Section 5 we show (following [FKW]) that indeed all non-degenerate admissible $\hat { s } \hat { \ell } _ { 2 }$ -modules produce all minimal series Virasoro modules via the functor $M \to H ( M )$ . In Section 6 we show, in a similar fashion, that all non-degenerate admissible $s p o ( 2 | 1 )$ -modules (whose characters were computed in [KW1] as well) produce all characters of minimal series Neveu–Schwarz modules. Finally, in Section 7, using the conjectural character formulas for “boundary” admissible $s \ell ( 2 | 1 ) \cdot$ - modules, we recover the characters of all minimal series modules over the $N = 2$ superconformal algebra. Note that it was already established by Khovanova [Kh] that the classical reduction of $s \ell ( 2 | 1 )$ produces the $N = 2$ superconformal algebra.

Further examples and results are presented in [KW5], where, in particular, we give a proof of a stronger form of the fundamental Conjecture 2.1 of the present paper.

The results of this paper were reported at the ICM in Beijing [K5].

Throughout the paper all vector spaces, algebras and tensor products are considered over the field of complex numbers $\mathbf { C }$ , unless otherwise stated. We denote by $\mathbf { Z }$ , $\mathbf { Q }$ and $\mathbf { R }$ the rings of integers, rational and real numbers, respectively, and by $\mathbf { Z } _ { + }$ the set of non-negative integers..

# 1 An Overview of the Operator Product Expansion

In this section, we give a brief summary of some basic properties of the operator product expansion (OPE) which will be used in this paper (for the details, see [K4] or [W]).

Let $A$ be a Lie superalgebra with a central element $K$ and a $\mathbf { Z }$ -filtration by subspaces,

$$
\dots \supset A _ {(0)} \supset A _ {(1)} \supset A _ {(2)} \supset \dots
$$

where $\begin{array} { r } { \bigcup _ { j } A _ { ( j ) } = A , \bigcap _ { j } A _ { ( j ) } = 0 } \end{array}$ and $\left| A _ { ( i ) } , A _ { ( j ) } \right| \subset A _ { ( i + j ) }$ . Throughout this paper, we always write $[ , ]$ for the Lie superbracket. For a given complex number $k \in \mathbf { \mu } _ { \mathbf { C } }$ , we denote by $U _ { k } ( A )$ the quotient of the universal enveloping algebra of $A$ by the ideal generated by $K - k \cdot 1$ , and by $U _ { k } ( A ) ^ { \mathrm { c o m } }$ the completion of $U _ { k } ( A )$ , which consists of all series $\textstyle \sum _ { j } u _ { j }$ $( u _ { j } \in U _ { k } ( A ) )$ ), such that for each $N \in \mathbf { Z } _ { + }$ all but a finite number of the $u _ { j }$ ’s lie in $U _ { k } ( A ) A _ { ( N ) }$ . Then $U _ { k } ( A ) ^ { \mathrm { c o m } }$ is an associative algebra containing $U _ { k } ( A )$ . Any $A$ -module $M$ in which every element of $M$ is annihilated by some $A _ { ( N ) }$ , can be uniquely extended to a module over $U _ { k } ( A ) ^ { \mathrm { c o m } }$ . Such a module over $A$ is called a restricted $A$ -module.

A $U _ { k } ( A ) ^ { \mathrm { c o m } }$ -valued field is an expression of the form

$$
a (z) = \sum_ {n \in {\bf Z}} a _ {(n)} z ^ {- n - 1},
$$

where $a _ { ( n ) } \in U _ { k } ( A ) ^ { \mathrm { c o m } }$ satisfy the property that for each $N \in \mathbf { Z } _ { + }$ , $a _ { ( n ) } \in U _ { k } ( A ) ^ { \mathrm { c o m } } A _ { ( N ) }$ for $n \gg 0$ , and all $a _ { ( n ) }$ have the same parity, which will be denoted by $p ( a ) \in \mathbf { Z } / 2 \mathbf { Z }$ . Note that for a restricted $A$ -module $M$ , the image of a field in $\operatorname { E n d } ( M )$ gives rise to a usual $\operatorname { E n d } ( M )$ -valued field. It is easy to see that the derivative $\partial _ { z } a ( z )$ of a field $a ( z )$ is also a field. The normal ordered product of two

fields $a ( z )$ and $b ( z )$ is defined by

$$
: a (z) b (z) := a (z) _ {-} b (z) + (- 1) ^ {p (a) p (b)} b (z) a (z) _ {+},
$$

where $\begin{array} { r } { a ( z ) _ { + } \ = \ \sum _ { n < 0 } a _ { ( n ) } z ^ { - n - 1 } } \end{array}$ and $\begin{array} { r } { a ( z ) _ { - } ~ = ~ \sum _ { n \geq 0 } a _ { ( n ) } z ^ { - n - 1 } } \end{array}$ . For $\textbf { \textit { n } } \in \textbf { Z }$ , the $n$ -th product $a ( z ) _ { ( n ) } b ( z ) ( = ( a _ { ( n ) } b ) ( z ) )$ of $a ( z )$ and $b ( z )$ is defined as follows. For a non-negative integer $n$ ,

$$
a (z) _ {(n)} b (z) = \operatorname {R e s} _ {x} (x - z) ^ {n} [ a (x), b (z) ],
$$

and

$$
a (z) _ {(- n - 1)} b (z) = \frac {: \partial_ {z} ^ {n} a (z) b (z) :}{n !}.
$$

The $n ^ { \mathrm { t h } }$ products of fields $a ( z )$ and $b ( z )$ for $n \in \mathbf { Z } _ { + }$ are encoded in the λ-bracket defined by

$$
[ a _ {\lambda} b ] = \sum_ {n \in \mathbf {Z} _ {+}} \frac {\lambda^ {n}}{n !} a _ {(n)} b,
$$

which is in general a formal power series in $\lambda$ ( with coefficients in $U _ { k } ( A ) ^ { \mathrm { c o m } }$ ). Here and further on, we often drop the indeterminate $z$ , e.g., we shall write $\partial a$ in place of $\partial _ { z } a ( z )$ .

Proposition 1.1 [K4] The following properties hold for the λ-bracket:

(sesquilinearity)

$$
[ \partial a _ {\lambda} b ] = - \lambda [ a _ {\lambda} b ], \quad [ a _ {\lambda} \partial b ] = (\partial + \lambda) [ a _ {\lambda} b ];
$$

(Jacobi identity)

$$
\left[ a _ {\lambda} \left[ b _ {\mu} c \right] \right] = \left[ \left[ a _ {\lambda} b \right] _ {\lambda + \mu} c \right] + (- 1) ^ {p (a) p (b)} \left[ b _ {\mu} \left[ a _ {\lambda} c \right] \right];
$$

(noncommutative Wick formula)

$$
[ a _ {\lambda}: b c: ] =: [ a _ {\lambda} b ] c: + (- 1) ^ {p (a) p (b)}: b [ a _ {\lambda} c ]: + \int_ {0} ^ {\lambda} [ [ a _ {\lambda} b ] _ {\mu} c ] d \mu .
$$

Recall that a pair $( a ( z ) , b ( z ) )$ of fields is called local if

$$
\left. (z - w) ^ {N} [ a (z), b (w) ] = 0 \right., \text {f o r} N \gg 0 .
$$

Note that the $\lambda$ -bracket of two local fields is a polynomial in $\lambda$ .

Proposition 1.2 [K4] Let $( a ( z ) , b ( z ) )$ be a local pair of fields. Then

(a)

$$
\left[ a _ {(m)}, b _ {(n)} \right] = \sum_ {j \in \mathbf {Z} _ {+}} \binom {m} {j} \left(a _ {(j)} b\right) _ {(m + n - j)}. \tag {1.1}
$$

(b) The λ-bracket satisfies the properties:

(skewcommutativity)

$$
\left[ a _ {\lambda} b \right] = - (- 1) ^ {p (a) p (b)} \left[ b _ {- \lambda - \partial a} \right];
$$

(right noncommutative Wick formula)[BK]

$$
\begin{array}{l} [ \colon a b: _ {\lambda} c ] =: (e ^ {\partial \frac {d}{d \lambda}} a) [ b _ {\lambda} c ]: + (- 1) ^ {p (a) p (b)}: (e ^ {\partial \frac {d}{d \lambda}} b) [ a _ {\lambda} c ]: \\ + (- 1) ^ {p (a) p (b)} \int_ {0} ^ {\lambda} [ b _ {\mu} [ a _ {\lambda - \mu} c ] ] d \mu . \\ \end{array}
$$

(c) The normal order commutator of $a ( z )$ and $b ( z )$ is expressed via the λ-bracket:

$$
: a b: - (- 1) ^ {p (a) p (b)}: b a: = \int_ {- \partial} ^ {0} [ a _ {\lambda} b ] d \lambda . \tag {1.2}
$$

Note that formula (1.1) is nothing else but (the singular part of) the operator product expansion (OPE) for the local pair $( a ( z ) , b ( z ) )$ :

$$
[ a (z), b (w) ] = \sum_ {j = 0} ^ {N} \frac {\partial_ {w} ^ {j} \delta (z - w)}{j !} a (w) _ {(j)} b (w)
$$

where $\begin{array} { r } { \delta ( z - w ) = z ^ { - 1 } \sum _ { n \in \mathbf { Z } } ( \frac { w } { z } ) ^ { n } } \end{array}$ is the formal $\delta$ -function. Propositions 1.1 and 1.2 provide an efficient and convenient way of calculating the OPE of local pairs.

Of course, in the case of normal ordered products of any number of free fields, one can use the usual Wick formula (see e.g. [K4]). Note that (1.1) immediately implies the following corollary.

Corollary 1.1 If $( a ( z ) , b ( z ) )$ is a local pair with $a ( z ) _ { ( 0 ) } b ( z ) = 0$ , then $[ a _ { ( 0 ) } , b ( z ) ] = 0$ .

Given a collection $\nu$ of pairwise local ( $U _ { k } ( A ) ^ { \mathrm { c o m } }$ -valued) fields, we may consider its closure $V = V _ { k } ( A , \mathcal { V } )$ which is the minimal space of fields containing 1 and $\nu$ , closed under $\partial _ { z }$ and all $n$ -th products ( $\boldsymbol { n } \in \mathbf { Z }$ ). By Dong’s lemma [K4], $V$ consists of pairwise local fields, hence Propositions 1.1 and 1.2 also apply to fields in $V$ . Note that $V$ is a vertex algebra and any restricted $A$ -module $M$ extends uniquely to a $V$ -module.

Example 1.1 (energy-momentum field). Let V ir be the Virasoro algebra, i.e., the Lie algebra with the basis $L _ { j }$ ( $j \in \mathbf { Z }$ ) and a central element $C$ , with the commutation relations

$$
[ L _ {m}, L _ {n} ] = (m - n) L _ {m + n} + \delta_ {m, - n} \frac {(m ^ {3} - m) C}{1 2}.
$$

We take the filtration $\begin{array} { r } { V i r _ { ( j ) } = { \bf C } C + \sum _ { i \geq j } { \bf C } L _ { j } } \end{array}$ for $j \le 0$ , $\begin{array} { r } { V i r _ { ( j ) } = \sum _ { i \geq j } } \end{array}$ CLj for $j > 0$ . Let $L ( z ) =$ $\textstyle \sum _ { n \in \mathbf { Z } } L _ { n } z ^ { - n - 2 }$ (note that $L _ { n } = L _ { ( n + 1 ) }$ ). This field is local with itself, so that the commutation relations of $L _ { j }$ ’s are encoded by the $\lambda$ -bracket,

$$
\left[ L _ {\lambda} L \right] = (\partial + 2 \lambda) L + \frac {\lambda^ {3} c}{1 2}. \tag {1.3}
$$

Here $c \in \textbf { C }$ is the eigenvalue of $C$ .

A local field $L ( z )$ with the $\lambda$ -bracket (1.3) is called an energy-momentum field with central charge $c$ .

Fix an energy-momentum field $L = L ( z )$ . Let $a ( z )$ be a field such that $( L , a )$ is a local pair. One says that the field $a$ has conformal weight $\triangle \in \mathbf { C }$ (with respect to $L$ ) if the following relation holds:

$$
[ L _ {\lambda} a ] = (\partial + \triangle \lambda) a + o (\lambda).
$$

Note that in this case $\partial _ { z } a ( z ) ( = \partial a )$ has conformal weight $\triangle + 1$ . In the special case, when $[ L _ { \lambda } a ] = ( \partial + \triangle \lambda ) a$ , one calls $a$ a primary field. When $a ( z )$ is a field with conformal weight $\bigtriangleup$ , it is convenient to change the indexation of the modes of $a ( z )$ :

$$
a (z) = \sum_ {n \in \mathbf {Z}} a _ {(n)} z ^ {- n - 1} = \sum_ {n \in - \triangle + \mathbf {Z}} a _ {n} z ^ {- n - \triangle}, a _ {n} = a _ {(n + \triangle - 1)}.
$$

For example $L ( z )$ has the conformal weight 2, and we write $\begin{array} { r } { L ( z ) = \sum _ { n \in \mathbf { Z } } L _ { n } z ^ { - n - 2 } } \end{array}$ .

Proposition 1.3 Let $a ( z ) , b ( z )$ be fields of conformal weights $\triangle _ { a }$ and $\triangle _ { b }$ respectively. Then

(a) $\bigtriangleup _ { a ( n ) } b = \bigtriangleup _ { a } + \bigtriangleup _ { b } - n - 1$ ; in particular, $\triangle _ { : a b : } = \triangle _ { a } + \triangle _ { b }$ .   
(b) The commutator formula (1.1) takes the homogeneous form:

$$
[ a _ {m}, b _ {n} ] = \sum_ {j \in \mathbf {Z} _ {+}} \binom {\triangle_ {a} + m - 1} {j} (a _ {(j)} b) _ {m + n}.
$$

Recall that a vector superspace is a vector space $V$ decomposed into a direct sum of vector spaces $V _ { \bar { 0 } }$ and $V _ { \bar { 1 } }$ ( $0 , 1 \in \mathbf { Z } / 2 \mathbf { Z } )$ , called the even and odd part of $V$ , respectively. We write $p ( v ) = \alpha$ if $v \in V _ { \alpha }$ . Denoting by $\Gamma$ the endomorphism of $V$ that acts as $( - 1 ) ^ { \alpha }$ on $V _ { \alpha }$ , we may define the supertrace of $a \in \operatorname { E n d } V$ (provided that $\mathrm { d i m } V < \infty$ ) by [K1]

$$
\mathrm {s t r} _ {V} a = \mathrm {t r} _ {V} (\Gamma a).
$$

In particular, letting $\mathrm { s d i m } V = \mathrm { s t r } _ { V } I _ { V }$ , we have $\mathrm { s d i m } V = \mathrm { d i m } V _ { \bar { 0 } } - \mathrm { d i m } V _ { \bar { 1 } }$ .

Recall that a vertex algebra is called strongly generated by a collection of fields $\mathcal { F }$ if normally ordered products of fields from $\mathbf { C } [ \partial ] \mathcal { F }$ span the space of fields of this vertex algebra.

Example 1.2 (neutral free superfermions). Let $A = A _ { \bar { 0 } } \oplus A _ { \bar { 1 } }$ be a finite-dimensional superspace with a non-degenerate skew-supersymmetric bilinear form $\langle . , . \rangle$ , i.e., $\langle A _ { \bar { 0 } } , A _ { \bar { 1 } } \rangle = 0$ and $\langle . , . \rangle$ is skewsymmetric (resp. symmetric ) on $A _ { \bar { 0 } }$ ( resp. $A _ { \bar { 1 } }$ ). Let $\widehat { A }$ be the Clifford affinization of $A$ , which is the Lie superalgebra $\hat { { \cal A } } = { \cal A } \otimes { \bf C } [ t , t ^ { - 1 } ] + { \bf C } K$ with the commutation relations

$$
[ a t ^ {m}, b t ^ {n} ] = \langle a, b \rangle \delta_ {m, - n - 1} K, [ K, \widehat {A} ] = 0.
$$

We take the filtration $\begin{array} { r } { \widehat { A } _ { j } = \mathbf { C } K + \sum _ { i \geq j } A t ^ { i } } \end{array}$ for $j \le 0$ , $\begin{array} { r } { \widehat { A } _ { ( j ) } = \sum _ { i \geq j } A t ^ { i } } \end{array}$ for $j > 0$ , and let $k = 1$ . For $\Phi \in A$ , let $\begin{array} { r } { \Phi ( z ) = \sum _ { n \in \mathbf { Z } } ( \Phi t ^ { n } ) z ^ { - n - 1 } } \end{array}$ . Then $\{ \Phi ( z ) \} _ { \Phi \in A }$ , called a collection of neutral free superfermions, which consists of pairwise local fields with $\lambda$ -bracket

$$
\left[ \Phi_ {\lambda} \Psi \right] = \left\langle \Phi , \Psi \right\rangle 1, \quad \Phi , \Psi \in A.
$$

Let $\{ \Phi _ { i } \}$ and $\{ \Phi ^ { i } \}$ be a pair of dual bases of $A$ , i.e., $\langle \Phi _ { i } , \Phi ^ { j } \rangle = \delta _ { i , j }$ , and define

$$
L = \frac {1}{2} \sum_ {i}: (\partial \Phi^ {i}) \Phi_ {i}:. \tag {1.4}
$$

Then $L$ is an energy-momentum field with central charge $c = - { \frac { 1 } { 2 } } \mathrm { s d i m } A$ . Furthermore, the neutral free superfermions $\Phi ( z )$ are all primary (with respective to this $L$ ) of conformal weight $\begin{array} { l } { \displaystyle { \frac { 1 } { 2 } } } \end{array}$ . The vertex algebra $F ( A )$ strongly generated by these superfermions, with the above energy-momentum field $L$ , is called the vertex algebra of neutral free superfermions. Via the state-field correspondence, $F ( A )$ is identified with the space $U _ { 1 } ( \widehat { A } ) / U _ { 1 } ( \widehat { A } ) \widehat { A } _ { ( 0 ) }$ , and all fields of $F ( A )$ act on this space from the left.

Example 1.3 (charged free superfermions). Let $A _ { \mathrm { c h } }$ be a finite-dimensional superspace with a non-degenerate skew-supersymmetric bilinear form $\langle . , . \rangle$ , and suppose that $A _ { \mathrm { c h } } = A _ { + } \oplus A _ { - }$ , where both $A _ { \pm }$ are isotropic subspaces of $A _ { \mathrm { c h } }$ . We have the Clifford affinization $\widehat { A } _ { \mathrm { c h } }$ with the filtration $( \widehat { A } _ { \mathrm { c h } } ) _ { ( j ) }$ (j ∈ Z+), $k = 1$ , and the fields $\varphi ( z ) , \varphi ^ { * } ( z )$ for $\varphi \in A _ { + } , \varphi ^ { * } \in A _ { - }$ , as in Example 1.2. We define the charges of the fields by

$$
\operatorname {c h a r g e} \varphi (z) = - \operatorname {c h a r g e} \varphi^ {*} (z) = 1. \tag {1.5}
$$

Let $\left\{ \varphi _ { i } \right\}$ ( resp. $\{ \varphi _ { i } ^ { * } \}$ ) be a basis of $A _ { + }$ (resp. $A _ { - }$ ) such that $\langle \varphi _ { i } , \varphi _ { j } ^ { * } \rangle = \delta _ { i , j }$ . The set of pairwise local fields $\{ \varphi _ { i } ( z ) \} \cup \{ \varphi _ { i } ^ { * } ( z ) \}$ is called a collection of charged free superfermions. In this case, we can define a family of energy-momentum fields parametrized by ${ \vec { m } } = ( m _ { i } ) _ { i } , m _ { i } \in { \bf C }$ :

$$
L ^ {\vec {m}} = - \sum_ {i} m _ {i}: \varphi_ {i} ^ {*} \partial \varphi_ {i}: + \sum_ {i} (1 - m _ {i}): (\partial \varphi_ {i} ^ {*}) \varphi_ {i}: .
$$

The central charge of $L ^ { \vec { m } } ( z )$ is equal to

$$
- \sum_ {i} (- 1) ^ {p (\varphi^ {i})} (1 2 m _ {i} ^ {2} - 1 2 m _ {i} + 2).
$$

Furthermore, the fields $\varphi _ { i } ^ { * } ( z )$ and $\varphi _ { i } ( z )$ are primary ( with respect to $L ^ { \vec { m } }$ ) of conformal weights $m _ { i }$ and $1 - m _ { i }$ respectively. The vertex algebra $F ( A _ { \mathrm { c h } } )$ with one of the energy-momentum fields $L ^ { \vec { m } }$ is called the vertex algebra of charged free superfermions. The relations (1.5) give rise to the charge decomposition of $F ( A _ { \mathrm { c h } } )$ :

$$
F \left(A _ {\mathrm {c h}}\right) = \bigoplus_ {m \in \mathbf {Z}} F _ {m} \left(A _ {\mathrm {c h}}\right). \tag {1.6}
$$

Example 1.4 (currents and the Sugawara construction). Let $^ { 9 }$ be a simple finite-dimensional Lie superalgebra with an even non-degenerate supersymmetric invariant bilinear form (.|.). Let $\widehat { \mathfrak { g } }$ be the Kac-Moody affinization of $\mathfrak { g }$ , i.e., $\widehat { \mathfrak { g } } = \mathfrak { g } \otimes \mathbf { C } [ t , t ^ { - 1 } ] \oplus \mathbf { C } K \oplus \mathbf { C } D$ with the commutation relations:

$$
[ a t ^ {m}, b t ^ {n} ] = [ a, b ] t ^ {m + n} + m \delta_ {m, - n} (a | b) K, [ D, a t ^ {m} ] = m a t ^ {m}, [ K, \widehat {\mathfrak {g}} ] = 0.
$$

The filtration in this situation is defined as in Example 1.2, and we fix $k \in \mathbf { C }$ .

For an element $a \in { \mathfrak { g } }$ , one associates the current field $\begin{array} { r } { a ( z ) = \sum _ { n \in \mathbf { Z } } ( a t ^ { n } ) z ^ { - n - 1 } } \end{array}$ . The collection $\{ a ( z ) \} _ { a \in { \mathfrak { g } } }$ consists of pairwise local fields with the following $\lambda$ -bracket,

$$
[ a _ {\lambda} b ] = [ a, b ] + \lambda (a | b) k , a, b \in \mathfrak {g} .
$$

The vertex algebra $V _ { k } ( { \mathfrak { g } } )$ strongly generated by the current fields $a ( z )$ is called the universal affine vertex algebra. Via the state-field correspondence, $V _ { k } ( { \mathfrak { g } } )$ is identified with the space $U _ { k } ( \widehat { \mathfrak { g } } ) / U _ { k } ( \widehat { \mathfrak { g } } ) \widehat { \mathfrak { g } } _ { ( 0 ) }$ and all fields of $V _ { k } ( { \mathfrak { g } } )$ act on this space from the left.

Let $\left\{ { a } _ { i } \right\}$ and $\{ a ^ { i } \}$ be a pair of dual bases of $\mathfrak { g }$ : $( a _ { i } | a ^ { j } ) = \delta _ { i , j }$ . Then $\begin{array} { r } { \Omega = \sum _ { i } ( - 1 ) ^ { p ( a _ { i } ) } a _ { i } a ^ { i } } \end{array}$ is the Casimir operator of $\mathfrak { g }$ , and it lies in the center of $U ( { \mathfrak { g } } )$ . The one-half of the eigenvalue of $\Omega$ in the adjoint representation, denoted by $h ^ { \vee }$ , is called the dual Coxeter number of $\mathfrak { g }$ (it depends on the normalization of (.|.)).

Recall the following relation between the Killing form and the form (.|.) [KW3]:

$$
\operatorname {s t r} _ {\mathfrak {g}} (\operatorname {a d} a) (\operatorname {a d} b) = 2 h ^ {\vee} (a | b), \quad a, b \in \mathfrak {g}. \tag {1.7}
$$

(Since the LHS is the Killing form, it is equal to $\gamma ( a | b )$ for some $\gamma$ . Hence $\operatorname { s t r } _ { \mathfrak { g } } \Omega = \gamma \operatorname { s d i m } { \mathfrak { g } }$ . Since $\Omega = 2 h ^ { \vee } I _ { \mathfrak { g } }$ , we conclude that $\gamma = 2 h ^ { \vee }$ provided that $\sin { \mathfrak { g } } \neq 0$ . Hence (1.7) holds for all exceptional Lie superalgebras and also for all the series $s \ell ( m | n )$ , etc. apart for the values $( m , n )$ on a hyperplane. Hence (1.7) holds for all values $( m , n )$ .)

Assuming $k + h ^ { \vee } \neq 0$ , introduce the so called Sugawara construction:

$$
L (z) = \frac {1}{2 (k + h ^ {\vee})} \sum_ {i} (- 1) ^ {p (a _ {i})}: a _ {i} (z) a ^ {i} (z).
$$

This is an energy momentum field with the central charge

$$
c (k) = \frac {k \operatorname {s d i m} \mathfrak {g}}{k + h ^ {\vee}}. \tag {1.8}
$$

All currents are primary with respect to $L$ of conformal weight 1. We shall also use the following well known modification of the Sugawara construction. For a given $a \in { \mathfrak { g } } _ { \bar { 0 } }$ , let

$$
L ^ {(a)} = L + \partial a.
$$

This is again an energy momentum field, and its central charge becomes

$$
c (k, a) = c (k) - 1 2 k (a \mid a). \tag {1.9}
$$

With respect to $L ^ { ( a ) }$ , the currents are not primary anymore:

$$
\left[ L ^ {(a)} _ {\lambda} b \right] = \partial b + \lambda (b - [ a, b ]) - \lambda^ {2} k (a | b). \tag {1.10}
$$

However, one has

$$
[ L ^ {(a)} _ {\lambda} b ] = (\partial + (1 - m) \lambda) b, \quad \text {i f} [ a, b ] = m b, m \neq 0, \tag {1.11}
$$

since in this case $( a | b ) = 0$

# 2 The Quantum Reduction

# 2.1 The complex $\mathcal { C } ( { \mathfrak { g } } , x , f , k )$ and the associated vertex algebra $W _ { k } ( { \mathfrak { g } } , x , f )$

Here we describe a general construction of a vertex algebra via a differential complex, associated to a simple finite-dimensional Lie superalgebra and some additional data, by a quantum reduction procedure, generalizing that of [FF1], [FF2], [FKW], [BT].

Let $\mathfrak { g }$ be a simple finite-dimensional Lie superalgebra with a non-degenerate even supersymmetric invariant bilinear form (.|.). Fix a pair $x$ and $f$ of even elements of $\mathfrak { g }$ satisfying the following properties:

(A1) ad x is diagonizable with half-integer eigenvalues, i.e., we have the following eigenspace decomposition with respect to ad x:

$$
\mathfrak {g} = \oplus_ {j \in \frac {1}{2}} \mathbf {Z} \mathfrak {g} _ {j}. \tag {2.1}
$$

(A2) $f \in { \mathfrak { g } } _ { - 1 } , { \mathrm { ~ i . e . , ~ } } [ x , f ] = - f .$

It follows that $f$ is a nilpotent element of $\mathfrak { g }$ . We shall also assume

(A3) ad $\mid f : { \mathfrak { g } } _ { \frac { 1 } { 2 } } \to { \mathfrak { g } } _ { - { \frac { 1 } { 2 } } }$ is a vector space isomorphism.

The element $f$ defines a skew-supersymmetric even bilinear form on by the formula: $\mathfrak { g } _ { \frac { 1 } { 2 } }$

$$
\langle a, b \rangle = (f | [ a, b ]) . \tag {2.2}
$$

It follows from (A3) that this form is non-degenerate, since $\langle a , b \rangle = ( [ f , a ] | b )$ and (.|.) gives a non-degenerate pairing between and . Denote by $A _ { \mathrm { n e } }$ the vector superspace with the $\Theta _ { - \frac { 1 } { 2 } }$ $\mathfrak { g } _ { \frac { 1 } { 2 } }$ $\mathfrak { g } _ { \frac { 1 } { 2 } }$ non-degenerate supersymmetric bilinear form $\langle \cdot , \cdot \rangle$ .

Furthermore, let

$$
\mathfrak {g} _ {+} = \bigoplus_ {j > 0} \mathfrak {g} _ {j}, \mathfrak {g} _ {-} = \bigoplus_ {j <   0} \mathfrak {g} _ {j}, \tag {2.3}
$$

and let

$$
A = \sqcap \mathfrak {g} _ {+}, A ^ {*} = \sqcap \mathfrak {g} _ {+} ^ {*}, A _ {\mathrm {c h}} = A \oplus A ^ {*},
$$

where ⊓ stands the reversing the parity of a vector superspace. Let $\langle \cdot , \cdot \rangle$ be the skew-supersymmetric bilinear form on $A _ { \mathrm { c h } }$ defined by

$$
\langle A, A \rangle = \langle A ^ {*}, A ^ {*} \rangle = 0, \langle a, b ^ {*} \rangle = b ^ {*} (a) \quad \text {f o r} a \in A, b ^ {*} \in A ^ {*}.
$$

Define gradations of $A , A ^ { * }$ by (2.3):

$$
A = \bigoplus_ {j > 0} A _ {j} , \quad A ^ {*} = \bigoplus_ {j > 0} A _ {j} ^ {*}.
$$

Finally, fix a complex number $k$ such that $\boldsymbol { k } + \boldsymbol { h } ^ { \vee } \neq 0$ , where $h ^ { \vee }$ is the dual Coxeter number of $\mathfrak { g }$

We shall associate to the data $( { \mathfrak { g } } , x , f , k )$ a differential vertex algebra $( { \mathcal { C } } ( { \mathfrak { g } } , x , f , k ) , d _ { 0 } )$ (by this we mean that $\mathcal { C }$ is a vertex algebra and $d _ { 0 }$ is an odd derivation of all $n$ -th products of $\mathcal { C }$ , and $d _ { 0 } ^ { 2 } = 0$ ).

Let ${ \widehat { \mathfrak { g } } } , { \widehat { A } } _ { \mathrm { n e } }$ and $\widehat { A } _ { \mathrm { c h } }$ be the Kac–Moody and Clifford affinizations corresponding to ${ \mathfrak { g } } , A _ { \mathrm { n e } }$ and $A _ { \mathrm { c h } }$ respectively (see Examples 1.4, 1.2 and 1.3). Let $U _ { k } = U _ { k } ( \widehat { \mathfrak { g } } ) \otimes U _ { 1 } ( \widehat { A } _ { \mathrm { c h } } ) \otimes U _ { 1 } ( \widehat { A } _ { \mathrm { n e } } )$ , and let $U _ { k } ^ { \mathrm { c o m } }$ be the completion of $U _ { k }$ as defined in Section 1. Consider the corresponding vertex algebras $V _ { k } ( { \mathfrak { g } } )$ , $F ( A _ { \mathrm { c h } } )$ and $F ( A _ { \mathrm { n e } } )$ , generated by the currents (based on $\mathfrak { g }$ ), charged free fermions (based on $A _ { \mathrm { c h } }$ ), and neutral free fermions ( based on $A _ { \mathrm { n e } }$ ) respectively. Consider the vertex algebras

$$
F (\mathfrak {g}, x, f) = F (A _ {\mathrm {c h}}) \otimes F (A _ {\mathrm {n e}}), \quad \mathcal {C} (\mathfrak {g}, x, f, k) = V _ {k} (\mathfrak {g}) \otimes F (\mathfrak {g}, x, f).
$$

By letting $\mathrm { c h a r g e } ( V _ { k } ( { \mathfrak { g } } ) ) \ = \ \mathrm { c h a r g e } ( F ( A _ { \mathrm { n e } } ) ) \ = \ 0$ , and using (1.6), one has the induced charge decompositions of $F ( { \mathfrak { g } } , f )$ and $\mathcal { C } ( { \mathfrak { g } } , f , k )$ :

$$
F (\mathfrak {g}, x, f) = \bigoplus_ {m \in \mathbf {Z}} F _ {m} , \quad \mathcal {C} (\mathfrak {g}, x, f, k) = \bigoplus_ {m \in \mathbf {Z}} \mathcal {C} _ {m} .
$$

Next, we define a differential on $\mathcal { C } ( \mathfrak { g } , x , f , k )$ , which makes it a homology complex. For this purpose, choose a basis $\{ u _ { i } \} _ { i \in S ^ { \prime } }$ of ${ \mathfrak { g } } _ { \frac { 1 } { 2 } }$ , and extend it a basis $\{ u _ { i } \} _ { i \in S }$ of $^ { 9 + }$ compatible with the gradation (2.3). Furthermore, extend the latter basis to a basis $\{ u _ { i } \} _ { i \in \tilde { S } }$ of $^ { 9 }$ , compatible with this gradation, and define the structure constants $c _ { i j } ^ { \ell }$ by: $\begin{array} { r } { [ u _ { i } , u _ { j } ] = \sum _ { \ell } c _ { i j } ^ { \ell } u _ { \ell } } \end{array}$ . Denote by $\{ u ^ { i } \} _ { i \in S ^ { \prime } }$ the dual basis of ${ \mathfrak { g } } _ { \underline { { 1 } } }$ with respect to the form $\langle , \rangle$ , i.e., $\langle u _ { i } , u ^ { j } \rangle = \delta _ { i j }$ .

Denote by $\stackrel { 2 } { \{ \varphi _ { i } \} } _ { i \in { \cal S } } , \{ \varphi _ { i } ^ { * } \} _ { i \in { \cal S } }$ the corresponding bases of $A$ and $A ^ { * }$ , and by $\{ \Phi _ { i } \} _ { i \in S ^ { \prime } }$ the corresponding basis of $A _ { \mathrm { n e } }$ . The fields $\varphi _ { i } ( z ) , \varphi _ { i } ^ { * } ( z )$ $( i \in S )$ and $\Phi _ { i } ( z )$ $i \in S ^ { \prime }$ ) are called ghosts. Introduce the following field of the vertex algebra $\mathcal { C } ( \mathfrak { g } , x , f , k )$ :

$$
\begin{array}{l} d (z) = \sum_ {i \in S} (- 1) ^ {p (u _ {i})} u _ {i} (z) \otimes \varphi_ {i} ^ {*} (z) \otimes 1 - \frac {1}{2} \sum_ {i, j, \ell \in S} (- 1) ^ {p (u _ {i}) p (u _ {\ell})} c _ {i j} ^ {\ell} \otimes \varphi_ {\ell} (z) \varphi_ {i} ^ {*} (z) \varphi_ {j} ^ {*} (z) \otimes 1 \\ + \sum_ {i \in S} (f | u _ {i}) \otimes \varphi_ {i} ^ {*} (z) \otimes 1 + \sum_ {i \in S ^ {\prime}} 1 \otimes \varphi_ {i} ^ {*} (z) \otimes \Phi_ {i} (z). \\ \end{array}
$$

For simplicity of notation, we shall omit the tensor sign $\otimes$ in the expression of fields. Note that in the second term of the expression of $d ( z )$ , one has

$$
\varphi_ {\ell} (z) \varphi_ {i} ^ {*} (z) \varphi_ {j} ^ {*} (z) =: \varphi_ {\ell} (z) \varphi_ {i} ^ {*} (z) \varphi_ {j} ^ {*} (z): \quad \mathrm {i f} c _ {i j} ^ {\ell} \neq 0 ,
$$

hence $d ( z )$ is a vertex algebra field. Also, it is easy to see that $d ( z )$ is an odd field independent of the choice of the basis. By the right non-commutative Wick formula one has the following $\lambda$ -brackets of $d ( z )$ and the currents $u _ { j } ( z )$ $j \in \tilde { S }$ ), and the ghosts $\varphi _ { j } ( z ) , \varphi _ { j } ^ { * } ( z ) ( j \in S )$ and $\Phi _ { j } ( z )$ (j ∈ S ′ ):

$$
{[ d _ {\lambda} u _ {j} ] =} {\sum_ {i \in S \atop \ell \in \bar {S}} (- 1) ^ {p (u _ {j}) + p (u _ {\ell}) p (u _ {i})} c _ {i j} ^ {\ell} u _ {\ell} \varphi_ {i} ^ {*} + (\partial + \lambda) k \sum_ {i \in S} (u _ {j} | u _ {i}) \varphi_ {i} ^ {*};}
$$

$$
\left[ d _ {\lambda} \varphi_ {j} \right] = u _ {j} + (f | u _ {j}) + \sum_ {i, \ell \in S} (- 1) ^ {p (u _ {\ell})} c _ {j i} ^ {\ell} \varphi_ {\ell} \varphi_ {i} ^ {*} + \sum_ {i \in S ^ {\prime}} (- 1) ^ {p (u _ {i})} \delta_ {i, j} \Phi_ {i}; \tag {2.4}
$$

$$
{[ d _ {\lambda} \varphi_ {j} ^ {*} ] =} {- \frac {1}{2} \sum_ {i, s \in S} (- 1) ^ {p (u _ {i}) p (u _ {j})} c _ {i s} ^ {j} \varphi_ {i} ^ {*} \varphi_ {s} ^ {*};}
$$

$$
{[ d _ {\lambda} \Phi_ {j} ] =} {\sum_ {i \in S ^ {\prime}} (f | [ u _ {i}, u _ {j} ]) \varphi_ {i} ^ {*}, [ d _ {\lambda} \Phi^ {j} ] = \varphi_ {j} ^ {*}.}
$$

Theorem 2.1 One has: $[ d ( z ) _ { \lambda } d ( z ) ] = 0$ .

Proof. We express the field $d ( z )$ as

$$
d (z) = d (z) ^ {\mathrm {s t}} + d (z) ^ {(I I I)} + d (z) ^ {(I V)}, \quad d (z) ^ {\mathrm {s t}} := d (z) ^ {(I)} + d (z) ^ {(I I)},
$$

where

$$
d ^ {(I)} = \sum_ {i \in S} (- 1) ^ {p (u _ {i})} u _ {i} \varphi_ {i} ^ {*}, d ^ {(I I)} = \frac {- 1}{2} \sum_ {i, j, \ell \in S} (- 1) ^ {p (u _ {i}) p (u _ {l})} c _ {i j} ^ {l} \varphi_ {\ell} \varphi_ {i} ^ {*} \varphi_ {j} ^ {*},
$$

$$
d ^ {(I I I)} = \sum_ {i \in S} (f | u _ {i}) \varphi_ {i} ^ {*}, \quad d ^ {(I V)} = \sum_ {i \in S ^ {\prime}} \varphi_ {i} ^ {*} \Phi_ {i}.
$$

Then

$$
[ d _ {\lambda} d ] = [ d _ {\lambda} ^ {\mathrm {s t}} d ^ {\mathrm {s t}} ] + [ d _ {\lambda} ^ {(I I)} d ^ {(I I I)} ] + [ d _ {\lambda} ^ {(I I I)} d ^ {(I I)} ] + [ d ^ {(I I) _ {\lambda}} d ^ {(I V)} ] + [ d ^ {(I V) _ {\lambda}} d ^ {(I I)} ] + [ d ^ {(I V) _ {\lambda}} d ^ {(I V)} ].
$$

It is well known that $\left[ d _ { \lambda } ^ { \mathrm { s t } } d ^ { \mathrm { s t } } \right] = 0$ , which follows from $[ d ^ { ( I I ) } \lambda d ^ { ( I I ) } ] = 0$ by the Jacobi identity. In the expression of $d ^ { ( I I ) }$ , one has $\ell \notin S ^ { \prime }$ whenever $c _ { i j } ^ { \ell } \neq 0$ , hence $[ \varphi \ell \lambda \varphi _ { k } ^ { * } ] = 0$ for $k \in S ^ { \prime }$ . This implies $[ d ^ { ( I I ) } \lambda d ^ { ( I V ) } ] = [ d ^ { ( I V ) } \lambda d ^ { ( I I ) } ] = 0$ . Hence

$$
\left[ d _ {\lambda} d \right] = \left[ d ^ {(I I)} _ {\lambda} d ^ {(I I I)} \right] + \left[ d ^ {(I I I)} _ {\lambda} d ^ {(I I)} \right] + \left[ d ^ {(I V)} _ {\lambda} d ^ {(I V)} \right].
$$

Note that $p ( u _ { \ell } ) = p ( u _ { i } ) + p ( u _ { j } )$ whenever $c _ { i j } ^ { \ell } \neq 0$ . We have

$$
\begin{array}{l} {[ d ^ {(I V)} _ {\lambda} d ^ {(I V)} ]} {= \sum_ {i, j \in S ^ {\prime}} [ \varphi_ {i} ^ {*} \Phi_ {i \lambda} \varphi_ {j} ^ {*} \Phi_ {j} ] = \sum_ {i, j \in S ^ {\prime}} (- 1) ^ {p (u _ {i}) (p (u _ {j}) + 1)} \varphi_ {i} ^ {*} \varphi_ {j} ^ {*} (f | [ u _ {i}, u _ {j} ])} \\ = \sum_ {i, j \in S ^ {\prime}, \ell \in S} (- 1) ^ {p (u _ {i}) (p (u _ {j}) + 1)} c _ {i j} ^ {\ell} (f | u _ {\ell}) \varphi_ {i} ^ {*} \varphi_ {j} ^ {*}; \\ \end{array}
$$

$$
\begin{array}{l} {[ d ^ {(I I)} _ {\lambda} d ^ {(I I I)} ]} {= \frac {- 1}{2} \sum_ {i, j, \ell \in S} (- 1) ^ {p (u _ {i}) p (u _ {\ell})} c _ {i j} ^ {\ell} (f | u _ {\ell}) [ \varphi_ {\ell} \varphi_ {i} ^ {*} \varphi_ {j} ^ {*} _ {\lambda} \varphi_ {\ell} ^ {*} ]} \\ = \frac {- 1}{2} \sum_ {i, j, \ell \in S} (- 1) ^ {p (u _ {i}) (p (u _ {i}) + p (u _ {j}))} c _ {i j} ^ {\ell} (f | u _ {\ell}) \varphi_ {i} ^ {*} \varphi_ {j} ^ {*} \\ = \frac {- 1}{2} \sum_ {i, j \in S ^ {\prime}, \ell \in S} (- 1) ^ {p (u _ {i}) (p (u _ {j}) + 1)} c _ {i j} ^ {\ell} (f | u _ {\ell}) \varphi_ {i} ^ {*} \varphi_ {j} ^ {*} \text {(s i n c e} i, j \in S ^ {\prime} \text {i f} (f | [ u _ {i}, u _ {j} ]) \neq 0). \\ \end{array}
$$

Therefore $[ d ^ { ( I I ) } \lambda d ^ { ( I I I ) } ] = [ d ^ { ( I I I ) } \lambda d ^ { ( I I ) } ] = \textstyle { \frac { - 1 } { 2 } } [ d ^ { ( I V ) } \lambda d ^ { ( I V ) } ]$ , hence $[ d _ { \lambda } d ] = 0$ . ✷

Let $d _ { 0 } = \mathrm { R e s } _ { z } d ( z )$ . Note that $d _ { 0 }$ is an odd element of $U _ { k } ^ { \mathrm { c o m } }$ , and that $[ d _ { 0 } , { \mathcal { C } } _ { m } ] \subset { \mathcal { C } } _ { m - 1 }$ Theorem 2.1 implies that $[ d ( z ) , d ( w ) ] = 0$ , hence $[ d _ { 0 } , d _ { 0 } ] = 2 d _ { 0 } ^ { 2 } = 0$ . Thus $( { \mathcal { C } } ( { \mathfrak { g } } , x , f , k ) , d _ { 0 } )$ is a homology complex. We denote the 0-th homology of this complex by $W _ { k } ( { \mathfrak { g } } , x , f )$ . Since $\mathcal { C } _ { 0 }$ is a vertex subalgebra of $\mathcal { C } ( { \mathfrak { g } } , f , k )$ , and since $d _ { 0 }$ is a derivation of all of its $n ^ { \mathrm { t h } }$ products, we conclude that $W _ { k } ( { \mathfrak { g } } , x , f )$ is a vertex algebra. This vertex algebra is called the quantum reduction for the quadruple $( { \mathfrak { g } } , x , f , k )$ .

The most interesting pair $x , f$ satisfying properties (A1), (A2), (A3) comes from an $s \ell _ { 2 }$ -triple $\{ e , x , f \}$ , where $[ x , e ] = e$ , $[ x , f ] = - f$ , $[ e , f ] = x$ . The validity of these properties is immediate by the $s \ell _ { 2 }$ -representation theory. Since a nilpotent even element $f$ determines uniquely (up to conjugation) the element $x$ of an $s \ell _ { 2 }$ -triple (by a theorem of Dynkin), we shall use in this case the notation $W _ { k } ( { \mathfrak { g } } , f )$ for the quantum reduction.

The vertex algebra $W _ { k } ( { \mathfrak { g } } , f )$ is a generalization of the quantum Drinfeld–Sokolov reduction, studied in [FF1], [FF2], [FKW] and many other papers, when $\mathfrak { g }$ is a simple Lie algebra and $f$ is the principal nilpotent element. The case studied in [B] is when ${ \mathfrak { g } } = { \mathfrak { s } } l _ { 3 }$ and $f$ is a non-principal nilpotent element. Our construction is a development of the generalizations proposed in [FKW] and in [BT].

Remark 2.1. (a) The assumption (A3) is not used in the proof of Theorem 2.1. However, this condition is essential for the construction of the energy-momentum field $L ( z )$ in Section 2.2.

(b) One can take for $x$ a diagonalizable derivation of $\mathfrak { g }$   
(c) Let $\mathfrak { n }$ be an ad $x$ -invariant subalgebra of $^ { 9 + }$ . The above construction when applied to $\mathfrak { n }$ in place of $^ { 9 + }$ produces a complex $( { \mathcal { C } } ( { \mathfrak { g } } , { \mathfrak { n } } , x , f , k ) , d { \mathfrak { n } } )$ . The corresponding vertex algebra $W _ { k } ( { \mathfrak { g } } , { \mathfrak { n } } , x , f )$ is naturally a subalgebra of $W _ { k } ( { \mathfrak { g } } , x , f )$ .

# 2.2 The energy-momentum field of $W _ { k } ( { \mathfrak { g } } , x , f )$

Denote by $L ^ { 9 } ( z )$ the Sugawara energy momentum field of $\widehat { \mathfrak { g } }$ (see Example 1.4), by $L ^ { \mathrm { n e } }$ the energy momentum field for $F ( A _ { \mathrm { n e } } )$ (see Example 1.2), and by $L ^ { \mathrm { c h } }$ the energy momentum field $L ^ { \vec { m } }$ for $F ( A _ { \mathrm { c h } } )$ (see Example 1.3) with $m _ { i }$ ’s defined by

$$
[ x, u _ {i} ] = m _ {i} u _ {i}.
$$

Let

$$
L (z) = L ^ {\mathfrak {g}} (z) + \partial_ {z} x (z) + L ^ {\mathrm {c h}} (z) + L ^ {\mathrm {n e}} (z). \tag {2.5}
$$

The discussion in Section 1 immediately implies the following result.

Theorem 2.2 (a) The field $L ( z )$ is the energy-momentum field for the vertex algebra $\mathcal { C } ( \mathfrak { g } , x , f , k )$ , and its central charge equals to

$$
c (\mathfrak {g}, x, f, k) = \frac {k \operatorname {s d i m} \mathfrak {g}}{k + h ^ {\vee}} - 1 2 k (x | x) - \sum_ {i \in S} (- 1) ^ {p \left(u _ {i}\right)} \left(1 2 m _ {i} ^ {2} - 1 2 m _ {i} + 2\right) - \frac {1}{2} \operatorname {s d i m} \mathfrak {g} _ {\frac {1}{2}}. \tag {2.6}
$$

(b) With respect to $L ( z )$ , the fields $\varphi _ { i } ( z ) , \varphi _ { i } ^ { * } ( z )$ ( $\mathbf { \chi } _ { i } ~ \in ~ S$ ), are primary of conformal weights $1 - m _ { i } , m _ { i }$ respectively, and the fields $\Phi _ { i } ( z )$ (i ∈ S′), are primary of conformal weight $\frac { 1 } { 2 }$ . The fields $u ( z )$ for $u \in { \mathfrak { g } } _ { j }$ have conformal weight $1 - j$ , and are primary unless $j = 0$ and $( x | u ) \neq 0$ .

✷

Remark 2.2. In the same way as in [FKW], formula (2.6) can be rewritten as follows (see also [BT]):

$$
c (\mathfrak {g}, x, f, k) = \mathrm {s d i m} \mathfrak {g} _ {0} - \frac {1}{2} \mathrm {s d i m} \mathfrak {g} _ {1 / 2} - 1 2 | \frac {\rho}{(k + h ^ {\vee}) ^ {1 / 2}} - x (k + h ^ {\vee}) ^ {1 / 2} | ^ {2}.
$$

The next theorem says that the field $L ( z )$ defined by (2.5) is the energy-momentum field for the vertex algebra $W _ { k } ( { \mathfrak { g } } , x , f )$ .

Theorem 2.3 We have $[ d _ { 0 } , L ( z ) ] = 0$ .

Proof. We compute the $\lambda$ -bracket $[ L _ { \lambda } d ]$ . Using Theorem 2.2 (b) and the Wick formula (the ”noncommutative” terms vanish everywhere) , we have ( recall that $u _ { i } \in { \mathfrak { g } } _ { + }$ ):

$$
[ L _ {\lambda} (u _ {i} \varphi_ {i} ^ {*}) ] = \partial (u _ {i} \varphi_ {i} ^ {*}) + \lambda u _ {i} \varphi_ {i} ^ {*},
$$

$$
[ L _ {\lambda} (\varphi_ {i} ^ {*} \varphi_ {j} ^ {*}) ] = \partial (\varphi_ {i} ^ {*} \varphi_ {j} ^ {*}) + (m _ {i} + m _ {j}) \lambda \varphi_ {i} ^ {*} \varphi_ {j} ^ {*},
$$

$$
{[ L _ {\lambda} (\varphi_ {\ell} \varphi_ {i} ^ {*} \varphi_ {j} ^ {*}) ]} {= \partial (\varphi_ {\ell} \varphi_ {i} ^ {*} \varphi_ {j} ^ {*}) + (1 - m _ {\ell} + m _ {i} + m _ {j}) \lambda \varphi_ {\ell} \varphi_ {i} ^ {*} \varphi_ {j} ^ {*},}
$$

$$
{[ L _ {\lambda} (\varphi_ {i} ^ {*} \Phi_ {i}) ]} {= \partial (\varphi_ {i} ^ {*} \Phi_ {i}) + (\frac {1}{2} + m _ {i}) \lambda \varphi_ {i} ^ {*} \Phi_ {i},}
$$

therefore

$$
[ L _ {\lambda} d ] = (\partial + \lambda) d + \lambda (\sum_ {i \in S} (m _ {i} - 1) (f | u _ {i}) \varphi_ {i} ^ {*} + \sum_ {i \in S ^ {\prime}} (m _ {i} - \frac {1}{2}) \varphi_ {i} ^ {*} \Phi_ {i}).
$$

Since $( f | u _ { i } ) = 0$ unless $m _ { i } = 1$ , and $m _ { i } = \frac { 1 } { 2 }$ if $i \in S ^ { \prime }$ , we have $[ L _ { \lambda } d ] \ : = \ : ( \partial + \lambda ) d$ . Hence, by skew-commutativity, $[ d _ { \lambda } L ] = \lambda d$ , and therefore, by Corollary 1.1, $[ d _ { 0 } , L ( z ) ] = 0$ . ✷

# 2.3 The quasiclassical limit

Here we briefly discuss the standard construction of the quasiclassical limit for the complex $\mathcal { C } ( \mathfrak { g } , x , f , k )$ . Denote by $A _ { \hbar }$ the space of the Lie superalgebra $A$ with the new bracket.

$$
[ a, b ] _ {\hbar} = \hbar [ a, b ], \quad a, b \in A.
$$

Then $U _ { k } ( A _ { \hbar } )$ is the quotient of the tensor algebra over the vector space $A$ by the ideal generated by the elements $( K - k )$ and $a \otimes b - ( - 1 ) ^ { p ( a ) p ( b ) } b \otimes a - \hbar [ a , b ] ( a , b \in A )$ . Hence the limit of $U _ { k } ( A _ { \hbar } )$ as $\hbar  0$ is $S _ { k } ( A )$ , the symmetric superalgebra over $A$ quotiented by the ideal $( K - k )$ , with the Poisson bracket:

$$
\{u, v \} = \lim  _ {\hbar \rightarrow 0} \frac {1}{\hbar} [ u, v ] _ {\hbar}, \quad u, v \in S _ {k} (A) .
$$

In the same way as in Section 1, we construct the Poisson superalgebra $S _ { k } ( A ) ^ { \mathrm { c o m } } \supset S _ { k } ( A )$

As in Section 1, we consider $S _ { k } ( A ) ^ { \mathrm { c o m } }$ -valued fields, and define their $n$ -th product for $n \in \mathbf { Z } _ { + }$ , by $a ( z ) _ { ( n ) } b ( z ) = \mathrm { R e s } _ { x } ( x - z ) ^ { n } \{ a ( x ) , b ( z ) \}$ , and let $\begin{array} { r } { \left\{ a _ { \lambda } b \right\} = \sum _ { n \in \mathbf { Z } _ { + } } \frac { \lambda ^ { n } } { n ! } a _ { ( n ) } b } \end{array}$ . Since the product in $S _ { k } ( A ) ^ { \mathrm { c o m } }$ is (super)commutative, the normal ordered product becomes the usual product. Then Proposition 1.1 holds for $\{ a _ { \lambda } b \}$ , except that “non-commutative” Wick product formula turns into the Leibniz rule:

$$
\{a _ {\lambda} b c \} = \{a _ {\lambda} b \} c + (- 1) ^ {p (a) p (b)} b \{a _ {\lambda} c \}.
$$

Proposition 1.2 (a) and $( b )$ hold as well, while (c) turns into the supercommutativity of the product. The vertex algebra of free superfermions of Examples 1.1 and 1.2 in the quasiclassical limit turns into the Poisson vertex algebra generated by the fields $\{ a ( z ) \} _ { a \in A }$ with the $\lambda$ -bracket

$$
\{a _ {\lambda} b \} = \langle a, b \rangle 1.
$$

All formulas of Example 1.2—1.4 hold in the limit, except that the Virasoro central charge becomes 0 for Examples 1.2, 1.3 and in the formula (1.8) in Example 1.4, hence (1.9) becomes $- 1 2 k ( a | a )$ . Thus, the central charge of $L ( z )$ in the limit becomes $- 1 2 k ( x | x )$ (cf (2.6)).

The method of quantum reduction turns the Poisson structure of the quasiclassical limit into the generalized Drinfeld-Sokolov reduction. The complex $\mathcal { C } ( \mathfrak { g } , x , f , k )$ turns into the tensor product of the corresponding Poisson vertex algebras, and the differential $d$ is given by the same formula, (except that the commutator with $d$ is replaced by the Poisson bracket with $d$ ). Finally, the energy-momentum field is given by the same formula, but the central charge is $- 1 2 k ( x | x )$ .

# 2.4 The basic conjecture on the structure of $W _ { k } ( { \mathfrak { g } } , x , f )$

Let ${ \mathfrak { g } } ^ { f }$ be the centralizer of $f$ in $\mathfrak { g }$ . The gradation (2.1) induces a $\scriptstyle { \frac { 1 } { 2 } } \mathbf { Z }$ -gradation

$$
\mathfrak {g} ^ {f} = \bigoplus_ {j} \mathfrak {g} _ {j} ^ {f}. \tag {2.7}
$$

For a good description of the vertex algebra $W _ { k } ( { \mathfrak { g } } , x , f )$ the following additional condition is apparently necessary:

(A4) The operator ad $f$ maps ${ \mathfrak { g } } _ { \ v j }$ to ${ \mathfrak { g } } _ { j - 1 }$ injectively for $j \geq 1$ and surjectively for $j \le 0$ .

By the representation theory of $s \ell _ { 2 }$ , this condition holds if the pair $x , f$ can be embedded in an $s \ell _ { 2 }$ -triple (but there are many more examples).

We shall call a pair $( x , f )$ satisfying conditions (A1)–(A4) to be a good pair, and the one coming from an $s \ell _ { 2 }$ -triple a Dynkin pair. The corresponding $\scriptstyle { \frac { 1 } { 2 } } \mathbf { Z }$ -gradations are called good and Dynkin gradations, respectively. Note that these gradations uniquely determine $x$ (by definition) and also determine $f$ up to conjugation by $G _ { 0 } = \exp ( \mathfrak { g } _ { 0 , \bar { 0 } } )$ (preserving the gradation), since $[ { \mathfrak { g } } _ { 0 , \bar { 0 } } , f ] = { \mathfrak { g } } _ { - 1 , \bar { 0 } }$ and therefore $f$ lies in the open orbit of $G _ { 0 }$ .

Conjecture 2.1 Suppose that conditions (A1)–(A4) hold. Then for each $a \in { \mathfrak { g } } _ { - j } ^ { f }$ (j ≥ 0) there exists a field $F _ { a } ( z )$ of the vertex algebra $\mathcal { C } ( \mathfrak { g } , x , f , k )$ , such that the following properties hold:

(i) $[ d _ { 0 } , F _ { a } ( z ) ] = 0$   
(ii) $F _ { a } ( z )$ has conformal weight $1 + j$ with respect to $L ( z )$ ,   
(iii) $F _ { a } ( z ) - a ( z )$ is a linear combination of normally ordered products of the fields $b ( z )$ , where $b \in { \mathfrak { g } } _ { s }$ with $s > - j$ , the ghosts $\varphi _ { i } ( z ) , \varphi _ { i } ^ { * } ( z ) , \Phi _ { i } ( z )$ , and their derivatives.

Furthermore, the images of the fields $F _ { a _ { i } } ( z )$ in $W _ { k } ( { \mathfrak { g } } , x , f )$ , where $\left\{ { a } _ { i } \right\}$ is a basis of ${ \mathfrak { g } } ^ { f }$ compatible with the gradation (2.7), strongly generate the vertex algebra $W _ { k } ( { \mathfrak { g } } , x , f )$ .

Given $v \in { \mathfrak { g } }$ , introduce the fields (we assume here condition (A3)):

$$
v ^ {\mathrm {c h}} (z) = - \sum_ {i, j \in S} (- 1) ^ {p (\varphi_ {i})} c _ {i j} (v): \varphi_ {i} (z) \varphi_ {j} ^ {*} (z):,
$$

$$
{v ^ {\mathrm {n e}} (z)} {= - \frac {1}{2} \sum_ {i, j \in S ^ {\prime}} (- 1) ^ {p (\Phi_ {i})} c _ {i j} (v): \Phi_ {i} (z) \Phi^ {j} (z):,}
$$

where the $c _ { i j } ( \upsilon )$ are defined by $\begin{array} { r } { [ v , u _ { j } ] \ = \ \sum _ { i } c _ { i j } ( v ) u _ { i } } \end{array}$ and , as before, $\langle \Phi _ { i } , \Phi ^ { j } \rangle = \delta _ { i j }$ . Note that $v ^ { \mathrm { n e } } ( z ) = 0$ unless $v \in { \mathfrak { g } } _ { 0 }$ , and that all pairs of distinct fields from $\{ v , v ^ { \mathrm { c h } } , v ^ { \mathrm { n e } } \}$ have zero $\lambda$ -brackets. Let

$$
J ^ {\{v \}} (z) = v (z) + v ^ {\mathrm {c h}} (z) + v ^ {\mathrm {n e}} (z), J ^ {(v)} (z) = v (z) + v ^ {\mathrm {c h}} (z).
$$

The calculations with $v ^ { \mathrm { c h } }$ and $v ^ { \mathrm { n e } }$ will use the following lemma.

Lemma 2.1 (a) Let $v \in { \mathfrak { g } } _ { 0 }$ . Then

$$
[ v ^ {\mathrm {c h}} _ {\lambda} \varphi_ {k} ] = (- 1) ^ {p (v)} \sum_ {i \in S} c _ {i k} (v) \varphi_ {i},
$$

$$
[ v ^ {\mathrm {c h}} _ {\lambda} \varphi_ {k} ^ {*} ] = - (- 1) ^ {p (v) p (\varphi_ {k} ^ {*})} \sum_ {j \in S} c _ {k j} (v) \varphi_ {j} ^ {*}.
$$

(b) Let $v \in { \mathfrak { g } } _ { 0 } ^ { f }$ . Then

$$
[ v ^ {\mathrm {n e}} _ {\lambda} \Phi_ {k} ] = (- 1) ^ {p (v)} \sum_ {i \in S ^ {\prime}} c _ {i k} (v) \Phi_ {i}.
$$

Proof. The proof of (a) is straightforward, using the Wick formula and the observation that

$$
p \left(u _ {i}\right) + p \left(u _ {j}\right) = p (v) \text {i f} c _ {i j} (v) \neq 0. \tag {2.8}
$$

For the proof of (b) we choose a basis $\{ u ^ { i } \} _ { i \in S ^ { \prime } }$ of ${ \mathfrak { g } } _ { 1 / 2 }$ such that $\langle u _ { i } , u ^ { j } \rangle = \delta _ { i j }$ (recall that the skew-supersymmetric bilinear form $\langle . , . \rangle$ on ${ \mathfrak { g } } _ { 1 / 2 }$ defined by (2.2) is nondegenerate). Then we have:

$$
c _ {i j} (v) = \langle [ v, u _ {j} ], u ^ {i} \rangle . \tag {2.9}
$$

Furthermore, by the Jacobi identity, we have for $a , b \in { \mathfrak { g } } _ { 1 / 2 }$ and $v \in { \mathfrak { g } } _ { 0 } ^ { f }$

$$
\langle [ v, a ], b \rangle = (- 1) ^ {p (a) p (b)} \langle [ v, b ], a \rangle . \tag {2.10}
$$

The proof of (b) is straightforward, using the Wick formula and (2.8), (2.9), (2.10).

![](images/c4dcd17b70683b4578824a46d2cb6d2f005d22d39023442af9145b76433e22ea.jpg)

Let ${ \mathfrak { h } } ^ { f }$ be a maximal ad-diagonizable subalgebra of ${ \mathfrak { g } } _ { 0 } ^ { f }$ and let $\mathfrak { h }$ be a Cartan subalgebra of ${ \mathfrak { g } } _ { 0 }$ containing ${ \mathfrak { h } } ^ { f }$ (it contains $x$ ). We can choose a basis $\{ e _ { \alpha } \} _ { \alpha \in S ^ { \prime } }$ of ${ } ^ { 9 } { \frac { 1 } { 2 } }$ consisting of root vectors, and extend it to a basis $\{ e _ { \alpha } \} _ { \alpha \in S }$ of $^ { 9 + }$ consisting of root vectors. Thus we may think of $S ^ { \prime }$ and $S$ as subsets of the set of roots $\Delta \subset { \mathfrak { h } } ^ { * }$ of $\mathfrak { g }$ .

Lemma 2.1(a) implies that $\left[ h ^ { \mathrm { c h } } \lambda \varphi _ { \alpha } \right] = \alpha ( h ) \varphi _ { \alpha }$ , and $\left[ h ^ { \mathrm { c h } } \lambda \varphi _ { \alpha } ^ { * } \right] = - \alpha ( h ) \varphi _ { \alpha } ^ { * }$ for $h \in { \mathfrak { h } }$ and $\alpha \in S$ , hence

$$
\left[ J ^ {\{h \}} _ {\lambda} \varphi_ {\alpha} \right] = \alpha (h) \varphi_ {\alpha}, \left[ J ^ {\{h \}} _ {\lambda} \varphi_ {\alpha} ^ {*} \right] = - \alpha (h) \varphi_ {\alpha} ^ {*} \text {i f} h \in \mathfrak {h}, \alpha \in S. \tag {2.11}
$$

Likewise, Lemma 2.1(b) implies that $\left[ h ^ { \mathrm { n e } } \lambda \Phi _ { \alpha } \right] = \alpha ( h ) \Phi _ { \alpha }$ if $h \in \mathfrak { h } ^ { f }$ , hence

$$
\left[ J ^ {\{h \}} _ {\lambda} \Phi_ {\alpha} \right] = \alpha (h) \Phi_ {\alpha} \text {i f} h \in \mathfrak {h} ^ {f}, \alpha \in S ^ {\prime}. \tag {2.12}
$$

Part (a) of the following theorem confirms Conjecture 2.1 in the case $j = 0$ .

Theorem 2.4 (a) If $v \in { \mathfrak { g } } _ { 0 } ^ { f }$ , then $[ d _ { \lambda } J ^ { \{ v \} } ] = 0$ , hence the image of each $J ^ { \{ v \} }$ $( v \in { \mathfrak { g } } _ { 0 } ^ { f } )$ ) is a field of the vertex algebra $W _ { k } ( { \mathfrak { g } } , x , f )$ .

(b) $\begin{array} { r } { [ L _ { \lambda } J ^ { ( v ) } ] = ( \partial + ( 1 - j ) \lambda ) J ^ { ( v ) } + \delta _ { j 0 } \lambda ^ { 2 } ( \frac { 1 } { 2 } \operatorname { s t r } _ { \mathfrak { g } _ { + } } ( \operatorname { a d } v ) - ( k + h ^ { \vee } ) ( v | x ) ) \ i _ { 2 } } \end{array}$ f $v \in { \mathfrak { g } } _ { j }$

and the same formula holds for ${ \cal J } ^ { \{ v \} }$ if $v \in { \mathfrak { g } } _ { 0 }$

(c) $\begin{array} { r } { [ J ^ { \{ v \} } \lambda J ^ { \{ v ^ { \prime } \} } ] = J ^ { \{ [ v , v ^ { \prime } ] \} } + \lambda ( k ( v | v ^ { \prime } ) + \operatorname { s t r } _ { \mathfrak { g } _ { + } } ( \mathrm { a d } v ) ( \mathrm { a d } v ^ { \prime } ) - \frac { 1 } { 2 } \operatorname { s t r } _ { \mathfrak { g } _ { \frac { 1 } { 2 } } } ( \mathrm { a d } v ) ( \mathrm { a d } v ^ { \prime } ) ) i f v , v ^ { \prime } \in \mathfrak { g } _ { - } \mathfrak { g } _ { \frac { 1 } { 2 } } , } \end{array}$ gf0 ,

$$
[ J ^ {(v)} _ {\lambda} J ^ {(v ^ {\prime})} ] = J ^ {([ v, v ^ {\prime} ])} + \delta_ {i 0} \delta_ {j 0} \lambda (k (v | v ^ {\prime}) + \mathrm {s t r} _ {\mathfrak {g} _ {+}} (\mathrm {a d} v) (\mathrm {a d} v ^ {\prime})) i f v \in \mathfrak {g} _ {i}, v ^ {\prime} \in \mathfrak {g} _ {j} a n d i j \geq 0.
$$

Proof. Let $v \in { \mathfrak { g } } _ { 0 }$ . Due to (2.8) and $( { \mathfrak { g } } _ { 0 } | { \mathfrak { g } } _ { + } ) = 0$ , we obtain from (2.4):

$$
\left[ d _ {\lambda} v \right] = - \sum_ {i, j \in S} (- 1) ^ {p \left(u _ {i}\right)} c _ {i j} (v) u _ {i} \varphi_ {j} ^ {*}. \tag {2.13}
$$

Next, assuming that $c _ { i j } ( v ) \neq 0$ , we compute $[ d _ { \lambda } : \varphi _ { i } \varphi _ { j } ^ { * } : ]$ . Our assumption implies that the elements $u _ { i }$ and $u _ { j }$ have the same degree in the gradation (2.3), hence the degree of their commutator is larger. This implies that the integral term in the non-commutative Wick formula vanishes, i.e., $[ d _ { \lambda } : \varphi _ { i } \varphi _ { j } ^ { * } : ] = : [ d _ { \lambda } \varphi _ { i } ] \varphi _ { j } ^ { * } : + ( - 1 ) ^ { p ( \varphi _ { i } ) } : \varphi _ { i } [ d _ { \lambda } \varphi _ { j } ^ { * } ]$ :. Therefore, by (2.4) we obtain:

$$
[ d _ {\lambda} v ^ {\mathrm {c h}} ] = \sum_ {r = 1} ^ {5} [ d _ {\lambda} v ^ {\mathrm {c h}} ] _ {r},
$$

where

$$
[ d _ {\lambda} v ^ {\mathrm {c h}} ] _ {1} = - \sum_ {i, j \in S} (- 1) ^ {p (\varphi_ {i})} c _ {i j} (v): u _ {i} \varphi_ {j} ^ {*}:,
$$

$$
[ d _ {\lambda} v ^ {\mathrm {c h}} ] _ {2} = \sum_ {i, j \in S} c _ {i j} (v) (f | u _ {i}) \varphi_ {j} ^ {*} = \sum_ {j \in S} (f | [ v, u _ {j} ]) \varphi_ {j} ^ {*},
$$

$$
[ d _ {\lambda} v ^ {\mathrm {c h}} ] _ {3} = \sum_ {i, j, \ell , k \in S} (- 1) ^ {p (u _ {k})} c _ {i j} (v) c _ {i k} ^ {\ell}: \varphi_ {\ell} \varphi_ {k} ^ {*} \varphi_ {j} ^ {*}:,
$$

$$
[ d _ {\lambda} v ^ {\mathrm {c h}} ] _ {4} = \frac {1}{2} \sum_ {i, j, k, \ell \in S} (- 1) ^ {p (u _ {k}) p (u _ {j})} c _ {i j} (v) c _ {k \ell} ^ {j}: \varphi_ {i} \varphi_ {k} ^ {*} \varphi_ {\ell} ^ {*}:,
$$

$$
[ d _ {\lambda} v ^ {\mathrm {c h}} ] _ {5} = \sum_ {i, j \in S ^ {\prime}} c _ {i j} (v): \Phi_ {i} \varphi_ {j} ^ {*}:.
$$

It follows from (2.13) that

$$
\left[ d _ {\lambda} v \right] + \left[ d _ {\lambda} v ^ {\mathrm {c h}} \right] _ {1} = 0, \tag {2.14}
$$

and that

$$
\left[ d _ {\lambda} v ^ {\mathrm {c h}} \right] _ {2} = \sum_ {j \in S} \left(\left[ f, v \right] \mid u _ {j}\right) \varphi_ {j} ^ {*} (= 0 \text {i f} v \in \mathfrak {g} _ {0} ^ {f}). \tag {2.15}
$$

Furthermore, by relabeling the indices, one can write:

$$
[ d _ {\lambda} v ^ {\mathrm {c h}} ] _ {3} = - \sum_ {i, j, k, \ell \in S} (- 1) ^ {p (u _ {\ell}) + p (u _ {\ell}) p (u _ {k})} c _ {j \ell} (v) c _ {j k} ^ {i}: \varphi_ {i} \varphi_ {\ell} ^ {*} \varphi_ {k} ^ {*}:
$$

hence

$$
\left[ d _ {\lambda} v ^ {\mathrm {c h}} \right] _ {3} = \frac {1}{2} \sum_ {i, j, k, \ell \in S} (- 1) ^ {p (u _ {k})} \left(c _ {j \ell} (v) c _ {j k} ^ {i} - (- 1) ^ {p (u _ {k}) p (u _ {\ell})} c _ {j k} (v) c _ {j \ell} ^ {i}\right): \varphi_ {i} \varphi_ {k} ^ {*} \varphi_ {\ell} ^ {*}:.
$$

Therefore

$$
[ d _ {\lambda} v ^ {\mathrm {c h}} ] _ {3} + [ d _ {\lambda} v ^ {\mathrm {c h}} ] _ {4} = \frac {1}{2} \sum_ {i, \ell , k \in S} (- 1) ^ {p (u _ {k}) + p (u _ {k}) p (u _ {\ell})} A (i, \ell , k): \varphi_ {i} \varphi_ {k} ^ {*} \varphi_ {\ell} ^ {*}:,
$$

where $\begin{array} { r } { A ( i , \ell , k ) : = \sum _ { j \in S } ( ( - 1 ) ^ { p ( u _ { k } ) p ( u _ { \ell } ) } c _ { j \ell } ( v ) c _ { j k } ^ { i } - c _ { j k } ( v ) c _ { j \ell } ^ { i } + c _ { i j } ( v ) c _ { k \ell } ^ { j } ) } \end{array}$ . From the Jacobi identity: $\left[ v , [ u _ { k } , u _ { \ell } ] \right] = \left[ [ v , u _ { k } ] , u _ { \ell } \right] + ( - 1 ) ^ { p ( v ) p ( u _ { k } ) } [ u _ { k } , [ v , u _ { \ell } ] ] ,$ one has

$$
0 = \sum_ {i, j \in S} (c _ {i j} (v) c _ {k \ell} ^ {j} - c _ {j k} (v) c _ {j \ell} ^ {i} + (- 1) ^ {p (u _ {k}) p (u _ {\ell})} c _ {j \ell} (v) c _ {j k} ^ {i}) u _ {i} = \sum_ {i \in S} A (i, \ell , k) u _ {i},
$$

which implies $A ( i , \ell , k ) = 0$ for all $i , \ell , k$ . Thus we obtain

$$
\left[ d _ {\lambda} v ^ {\mathrm {c h}} \right] _ {3} + \left[ d _ {\lambda} v ^ {\mathrm {c h}} \right] _ {4} = 0. \tag {2.16}
$$

Next, we compute $[ d _ { \lambda } v ^ { \mathrm { n e } } ]$ for $v \in { \mathfrak { g } } _ { 0 } ^ { f }$ . For that recall the skew-supersymmetric bilinear form $\langle . , . \rangle$ on ${ \mathfrak { g } } _ { 1 / 2 }$ , given by (2.2), and formulas (2.9) and (2.10).

Using (2.4) and (2.9), we obtain:

$$
\left[ d _ {\lambda} v ^ {\mathrm {n e}} \right] = \left[ d _ {\lambda} v ^ {\mathrm {n e}} \right] _ {1} + \left[ d _ {\lambda} v ^ {\mathrm {n e}} \right] _ {2},
$$

where

$$
[ d _ {\lambda} v ^ {\mathrm {n e}} ] _ {1} = \frac {1}{2} \sum_ {i, j, k \in S ^ {\prime}} c _ {i j} (v) \langle u _ {i}, u _ {k} \rangle \varphi_ {k} ^ {*} \Phi_ {j},
$$

$$
\left[ d _ {\lambda} v ^ {\mathrm {n e}} \right] _ {2} = - \frac {1}{2} \sum_ {i, j \in S ^ {\prime}} c _ {i j} (v) \Phi_ {i} \varphi_ {j} ^ {*}. \tag {2.17}
$$

We have:

$$
[ d _ {\lambda} v ^ {\mathrm {n e}} ] _ {1} = \frac {1}{2} \sum_ {i, j, k} (- 1) ^ {p (u _ {j}) (p (u _ {k}) + 1)} c _ {i j} (v) \langle u _ {i}, u _ {k} \rangle \Phi^ {j} \varphi_ {k} ^ {*}.
$$

Using that $\begin{array} { r } { \Phi ^ { j } = \sum _ { r \in S ^ { \prime } } \langle u ^ { j } , u ^ { r } \rangle \Phi _ { r } } \end{array}$ and that $p ( u _ { i } ) = p ( u _ { k } )$ if $\langle u _ { i } , u _ { k } \rangle \neq 0$ , we obtain:

$$
[ d _ {\lambda} v ^ {\mathrm {n e}} ] _ {1} = - \frac {1}{2} \sum_ {i, j, k, r} (- 1) ^ {p (u _ {i}) p (u _ {j})} c _ {i j} (v) \langle u _ {i}, u _ {k} \rangle \langle u ^ {r}, u ^ {j} \rangle \Phi_ {r} \varphi_ {k} ^ {*}.
$$

Using (2.9) and (2.10), we obtain: $( - 1 ) ^ { p ( u _ { i } ) p ( u _ { j } ) } c _ { i j } ( v ) = \langle [ v , u ^ { i } ] , u _ { j } \rangle$ , hence:

$$
[ d _ {\lambda} v ^ {\mathrm {n e}} ] _ {1} = - \frac {1}{2} \sum_ {k, r} \langle [ v, \sum_ {i} \langle u _ {i}, u _ {k} \rangle u ^ {i} ], \sum_ {j} \langle u ^ {r}, u ^ {j} \rangle u _ {j} \rangle \Phi_ {r} \varphi_ {k} ^ {*} = - \frac {1}{2} \sum_ {k, r} c _ {r k} (v) \Phi_ {r} \varphi_ {k} ^ {*},
$$

and, by (2.17), we obtain:

$$
\left[ d _ {\lambda} v ^ {\mathrm {n e}} \right] = - \sum_ {i, j \in S ^ {\prime}} c _ {i j} (v) \Phi_ {i} \varphi_ {j} ^ {*}.
$$

Thus, we see that for $v \in { \mathfrak { g } } _ { 0 } ^ { f }$ one has

$$
\left[ d _ {\lambda} v ^ {\mathrm {c h}} \right] _ {5} + \left[ d _ {\lambda} v ^ {\mathrm {n e}} \right] = 0. \tag {2.18}
$$

Comparing (2.14), (2.15), (2.16) and (2.18) gives (a).

By (1.10), one has

$$
[ L _ {\lambda} v ] = (\partial + \lambda) v - \lambda^ {2} k (x | v),
$$

and by Theorem 2.2(b) and the noncommutative Wick formula the following relations hold:

$$
[ L _ {\lambda}: \varphi_ {i} \varphi_ {j} ^ {*}: ] = (\partial + \lambda): \varphi_ {i} \varphi_ {j} ^ {*}: + \left(\frac {1}{2} - m _ {i}\right) \delta_ {i, j} \lambda^ {2}, [ L _ {\lambda}: \Phi_ {i} \Phi^ {j}: ] = (\partial + \lambda): \Phi_ {i} \Phi^ {j}:.
$$

Hence:

$$
[ L _ {\lambda} v ^ {\mathrm {c h}} ] = (\partial + \lambda) v ^ {\mathrm {c h}} + \lambda^ {2} \left(\frac {1}{2} \operatorname {s t r} _ {\mathfrak {g} _ {+}} (\operatorname {a d} v) - \sum_ {i \in S} (- 1) ^ {p (u _ {i})} m _ {i} c _ {i i} (v)\right), \quad [ L _ {\lambda} v ^ {\mathrm {n e}} ] = (\partial + \lambda) v ^ {\mathrm {n e}}.
$$

As $\begin{array} { r } { \sum _ { i \in S } ( - 1 ) ^ { p ( u _ { i } ) } m _ { i } c _ { i i } ( v ) = h ^ { \vee } ( x | v ) } \end{array}$ by (1.7), (b) follows.

The proof of (c) is similar. It uses only the usual Wick formula. We omit the details.

✷

# 2.5 Construction of the $W _ { k } ( { \mathfrak { g } } , x , f )$ -modules

Let $M$ be a restricted $\widehat { \mathfrak { g } }$ -module of level $k$ (i.e. $K = k \ I _ { M }$ ). It extends to the $V _ { k } ( { \mathfrak { g } } )$ -module, and then to the $\mathcal { C } ( \mathfrak { g } , x , f , k )$ -module

$$
\mathcal {C} (M) = M \bigotimes F (\mathfrak {g}, x, f).
$$

One has the charge decomposition of $\mathcal { C } ( M )$ induced by that of $F ( { \mathfrak { g } } , x , f )$ by setting the charge of $M$ to be zero:

$$
\mathcal {C} (M) = \bigoplus_ {m \in \mathbf {Z}} \mathcal {C} (M) _ {m} .
$$

Furthermore, $( { \mathcal { C } } ( M ) , d _ { 0 } )$ form a $\mathcal { C } ( \mathfrak { g } , x , f , k )$ -module complex, hence its homology, $H ( M ) = \oplus _ { j \in \mathbf { Z } }$ $H _ { j } ( M )$ , is a direct sum of $W _ { k } ( { \mathfrak { g } } , x , f )$ -modules. We thus get a functor, which we denote by $H$ , from the category of restricted $\widehat { \mathfrak { g } }$ -modules to the category of $\mathbf { Z }$ -graded $W _ { k } ( { \mathfrak { g } } , x , f )$ -modules, that send $M$ to $H ( M )$ .

Remark 2.3. Let $| 0 \rangle$ be the vacuum vector of the vertex algebra $F ( { \mathfrak { g } } , x , f )$ and let $v \in M$ be such that $( { \mathfrak { g } } _ { + } t ^ { m } ) ( v ) = 0$ for all $m \geq 0$ . Then

$$
d _ {0} (v \otimes | 0 \rangle) = 0.
$$

In particular if $M$ is a highest weight $\widehat { \mathfrak { g } }$ -module with highest weight $\Lambda$ of level $k \neq - h ^ { \vee }$ , and $v _ { \Lambda }$ is the highest weight vector, then $d _ { 0 } ( v _ { \Lambda } \otimes | 0 \rangle ) = 0$ . So, if the vector $v _ { \Lambda } \otimes | 0 \rangle$ is not in the image of

$d _ { 0 }$ , its image in $H _ { 0 } ( M )$ , which we denote by ${ \tilde { v } } _ { \Lambda }$ , generates a non-zero $W _ { k } ( { \mathfrak { g } } , x , f )$ -submodule. Its central charge is given by formula (2.6). The eigenvalue of $L _ { 0 }$ on ${ \tilde { v } } _ { \Lambda }$ is equal to (cf. Section 3.1):

$$
\frac {(\Lambda | \Lambda + 2 \widehat {\rho})}{2 (k + h ^ {\vee})} - (x + D | \Lambda). \tag {2.19}
$$

The eigenvalue of ${ J _ { 0 } ^ { \{ h \} } }$ ( $h \in \mathfrak { h } ^ { f } )$ ) on ${ \tilde { v } } _ { \Lambda }$ is equal to $\Lambda ( h )$ .

# 3 Character Formulas

# 3.1 The Euler-Poincar´e character of $H ( M )$

Let $\mathfrak { g }$ be one of the basic simple finite-dimensional Lie superalgebras. Recall that, apart from the five exceptional Lie algebras, they are as follows: $s \ell ( m | n ) / \delta _ { m , n } \mathbf { C } I$ , $o s p ( m | n )$ , $D ( 2 , 1 ; a )$ , $F ( 4 )$ and $G ( 3 )$ [K1]. Recall that $\mathfrak { g }$ carries a unique (up to a constant factor) non-degenerate invariant bilinear form [K1], and it is automatically even supersymmetric. We choose one of them, and denote it by (.|.).

Given $h \in \mathfrak { h } ^ { f }$ , define, as before, the fields $h ^ { \mathrm { c h } } ( z )$ and $h ^ { \mathrm { n e } } ( z )$ . They are given by the following slightly simpler formulas:

$$
{h ^ {\mathrm {c h}} (z) = - \sum_ {\alpha \in S} \alpha (h): \varphi_ {\alpha} ^ {*} (z) \varphi_ {\alpha} (z):, h ^ {\mathrm {n e}} (z)} {= - \frac {1}{2} \sum_ {\alpha \in S ^ {\prime}} \alpha (h): \Phi^ {\alpha} (z) \Phi_ {\alpha} (z):.}
$$

Since these fields are of conformal weight 1, we write: $\begin{array} { r } { h ^ { \mathrm { c h } } ( z ) = \sum _ { n \in \mathbf { Z } } h _ { n } ^ { \mathrm { c h } } z ^ { - n - 1 } } \end{array}$ , and $h ^ { \mathrm { n e } } ( z ) =$ $\begin{array} { r } { \sum _ { n \in \mathbf { Z } } h _ { n } ^ { \mathrm { n e } } z ^ { - n - 1 } } \end{array}$ . Likewise, we write $\begin{array} { r } { J ^ { \{ h \} } ( z ) = \sum _ { n \in \mathbf { Z } } J _ { n } ^ { \{ h \} } z ^ { - n - 1 } } \end{array}$ .

Let $\widehat { \mathfrak { h } } = \mathfrak { h } + \mathbf { C } K + \mathbf { C } D$ be the Cartan subalgebra of the affine Lie superalgebra $\widehat { \mathfrak { g } }$ . As usual, we extend a root $\alpha \in \triangle$ to $\mathfrak { h }$ by letting $\alpha ( K ) = \alpha ( D ) = 0$ . We extend the bilinear form (.|.) from $\mathfrak { h }$ (on which it is non-degenerate) to $\widehat { \mathfrak { h } }$ by letting:

$$
\left(\mathfrak {h} \right| \mathbf {C} K + \mathbf {C} D) = 0, \quad \left(K \right| K) = \left(D \right| D) = 0, \quad \left(K \right| D) = 1.
$$

We shall identify $\widehat { \mathfrak { h } }$ with $\widehat { \mathfrak { h } } ^ { * }$ via this form. The bilinear form (.|.) extends further to the whole $\widehat { \mathfrak { g } }$ by letting $( t ^ { m } a | t ^ { n } b ) = \delta _ { m , - n } ( a | b )$ . Let $\widehat \Omega$ be the Casimir operator for $\widehat { \mathfrak { g } }$ and this bilinear form. Recall that its eigenvalue for a $\widehat { \mathfrak { g } }$ -module with the highest weight $\Lambda$ is equal to $( \Lambda | \Lambda + 2 \widehat { \rho } )$ [K3]. Denote by $\widehat { \triangle } \subset \widehat { \mathfrak { h } } ^ { * } = \widehat { \mathfrak { h } }$ the set of roots of $\widehat { \mathfrak { g } }$ with respect to $\widehat { \mathfrak { h } }$ .

Recall that $\widehat { \triangle } = \widehat { \triangle } ^ { \mathrm { { r e } } } \cup \widehat { \triangle } ^ { \mathrm { { i m } } }$ b, where

$$
\widehat {\triangle} ^ {\mathrm {r e}} = \left\{\alpha + n K \mid \alpha \in \bigtriangleup , n \in \mathbf {Z} \right\}, \quad \widehat {\triangle} ^ {\mathrm {i m}} = \left\{n K \mid n \in \mathbf {Z} \setminus \left\{0 \right\} \right\}
$$

are the sets of real and imaginary roots respectively. Choosing a set of positive roots $\triangle _ { 0 + }$ of the set of roots $\bigtriangleup _ { 0 } = \{ \alpha \in \triangle \ | \ ( \alpha | x ) = 0 \}$ , we get a set of positive roots $\triangle _ { + } = \{ \alpha \in \triangle \mid ( \alpha | x ) > 0 \} \cup \triangle _ { 0 + }$ of $\mathfrak { g }$ and the set of positive roots

$$
\widehat {\triangle} _ {+} = \triangle_ {+} \bigcup \{\alpha + n K \mid \alpha \in \triangle \cup \{0 \}, n > 0 \}
$$

of $\widehat { \mathfrak { g } }$ . We shall denote by $\widehat { \triangle } _ { \mathrm { e v e n } }$ and $\widehat { \triangle } _ { \mathrm { o d d } }$ , $\widehat { \Delta } _ { + \mathrm { e v e n } }$ and $\widehat { \Delta } _ { \mathrm { + o d d } }$ , etc. the sets of even and odd roots respectively.

Introduce the following subsets of $\widehat { \triangle }$ (where $n \in \mathbf { Z }$ ):

$$
\widehat {S} = \left\{\alpha + n K \mid \alpha \in S, n \geq 0 \right\} \bigcup \left\{- \alpha + n K \mid \alpha \in S, n > 0 \right\},
$$

$$
\widehat {S} ^ {\prime} = \left\{- \alpha + n K \mid \alpha \in S ^ {\prime}, n > 0 \right\}.
$$

As usual, we write $\begin{array} { r } { L ( z ) \ = \ \sum _ { n \in \mathbf { Z } } L _ { n } z ^ { - n - 2 } } \end{array}$ , $\begin{array} { r } { L ^ { \mathfrak { g } } ( z ) = \sum _ { n \in \mathbf { Z } } L _ { n } ^ { \mathfrak { g } } z ^ { - n - 2 } } \end{array}$ , etc (see Section 2.2). Recall that we have [K3]:

$$
L _ {0} ^ {\mathfrak {g}} = \frac {\widehat {\Omega}}{2 (k + h ^ {\vee})} - D
$$

for any highest weight $\widehat { \mathfrak { g } }$ -module $M$ of level $k$ , $k \neq - h ^ { \vee }$

We shall coordinatize $\widehat { \mathfrak { h } }$ by letting

$$
(\tau , z, u) = 2 \pi i (z - \tau D + u K),
$$

where $z \in { \mathfrak { h } }$ , $\tau , u \in \mathbf { C }$ . We shall assume that Im $\tau > 0$ in order to guarantee the convergence of characters, and set $q = e ^ { 2 \pi \imath \tau }$ . Define the character of a $\widehat { \mathfrak { g } }$ -module $M$ by $\mathrm { c h } _ { M } : = \mathrm { t r } _ { M } e ^ { 2 \pi i ( z - \tau D + u K ) }$ . For any highest weight $\widehat { \mathfrak { g } }$ -module $M$ of level $k \ne - h ^ { \vee }$ the series $\mathrm { c h } _ { M }$ converges to an analytic function in the interior of the domain $Y _ { > } : = \{ h \in \widehat { \mathfrak { h } } \ | \ ( \alpha | h ) > 0$ for all $\alpha \in \widehat { \triangle } _ { + } \}$ ; moreover, the domain of convergence $Y ( M )$ is a convex domain contained in the upper half space $Y = \{ h \in$ ${ \mathfrak { h } } \mid \operatorname { R e } ( h | K ) > 0 \} = \{ ( \tau , z , u ) \mid \operatorname { I m } \tau > 0 \}$ ([K3], Lemma 10.6).

Lemma 3.1 For a regular element $b \in { \mathfrak { h } }$ (i.e., $\alpha ( b ) \neq 0$ for all $\alpha \in \triangle$ ), $h \in \mathfrak { h } ^ { f }$ and any sufficiently small $\epsilon \in { \bf C } \backslash \{ 0 \}$ , one has the following formula for the Euler-Poincar´e character of $F ( { \mathfrak { g } } , x , f )$ :

$$
\sum_ {j \in \mathbf {Z}} (- 1) ^ {j} \mathrm {t r} _ {F _ {j}} \left(q ^ {L _ {0} ^ {\mathrm {c h}} + \epsilon b _ {0} ^ {\mathrm {c h}} + L _ {0} ^ {\mathrm {n e}}} e ^ {2 \pi i (h _ {0} ^ {\mathrm {c h}} + h _ {0} ^ {\mathrm {n e}})}\right) = \prod_ {\alpha \in \widehat {S} \backslash \widehat {S} ^ {\prime}} (1 - s (\alpha) e ^ {- \alpha}) ^ {s (\alpha)} (\tau , \tau (\epsilon b - x) + h, 0), \quad (3. 1)
$$

where $s ( \alpha ) : = ( - 1 ) ^ { p ( \alpha ) } , \alpha \in \widehat { \triangle }$

Proof. Since the fields $\varphi _ { \alpha } ^ { * }$ and $\varphi _ { \alpha }$ (resp. $\Phi _ { \alpha }$ ) are primary with respect to $L ^ { \mathrm { c h } }$ (resp. $L ^ { \mathrm { n e } }$ ) of conformal weights $( \alpha | x )$ and $1 - \left( \alpha | x \right)$ (resp. $\begin{array} { l } { { \frac { 1 } { 2 } } } \end{array}$ ), we have:

$$
\left[ L _ {0} ^ {\mathrm {c h}}, \varphi_ {\alpha (- n)} \right] = (n - (\alpha | x)) \varphi_ {\alpha (- n)},
$$

$$
{[ L _ {0} ^ {\mathrm {c h}}, \varphi_ {\alpha (- n)} ^ {*} ]} = {(n - 1 + (\alpha | x)) \varphi_ {\alpha (- n)} ^ {*},}
$$

$$
\left[ L _ {0} ^ {\mathrm {n e}}, \Phi_ {\alpha (- n)} \right] = (n - \frac {1}{2}) \Phi_ {\alpha (- n)}.
$$

Using this, we get:

$$
\sum_ {j \in \mathbf {Z}} (- 1) ^ {j} \mathrm {t r} _ {F (A) _ {j}} q ^ {L _ {0} ^ {\mathrm {c h}}} = \prod_ {\alpha \in S} \prod_ {n = 1} ^ {\infty} \left(1 - s (\alpha) q ^ {(n K - \alpha | D + x)}\right) ^ {s (\alpha)},
$$

$$
\sum_ {j \in {\bf Z}} (- 1) ^ {j} \mathrm {t r} _ {F (A ^ {*}) _ {j}} q ^ {L _ {0} ^ {\mathrm {c h}}} = \prod_ {\alpha \in S} \prod_ {n = 1} ^ {\infty} \left(1 - s (\alpha) q ^ {((n - 1) K + \alpha | D + x)}\right) ^ {s (\alpha)},
$$

$$
\mathrm {t r} _ {F (A _ {\mathrm {n e}})} q ^ {L _ {0} ^ {\mathrm {n e}}} = \prod_ {\alpha \in S ^ {\prime}} \prod_ {n = 1} ^ {\infty} \Big (1 - s (\alpha) q ^ {(n K - \alpha | D + x)} \Big) ^ {- s (\alpha)}.
$$

Using these formulas along with (2.11) and (2.12), we get for $h \in { \mathfrak { h } }$

$$
\sum_ {j \in \mathbf {Z}} (- 1) ^ {j} \mathrm {t r} _ {F (A) _ {j}} q ^ {L _ {0} ^ {\mathrm {c h}}} e ^ {2 \pi i J _ {0} ^ {\{h \}}} = \prod_ {\alpha \in S} \prod_ {n = 1} ^ {\infty} \left(1 - s (\alpha) e ^ {2 \pi i (- n K + \alpha | - \tau (D + x) + h)}\right) ^ {s (\alpha)}
$$

$$
\sum_ {j \in {\bf Z}} (- 1) ^ {j} \mathrm {t r} _ {F (A ^ {*}) _ {j}} q ^ {L _ {0} ^ {\mathrm {c h}}} e ^ {2 \pi i J _ {0} ^ {\{h \}}} = \prod_ {\alpha \in S} \prod_ {n = 1} ^ {\infty} \Big (1 - s (\alpha) e ^ {2 \pi i (- (n - 1) K - \alpha | - \tau (D + x) + h)} \Big) ^ {s (\alpha)},
$$

and for $h \in \mathfrak { h } ^ { f }$ :

$$
\mathrm {t r} _ {F (A _ {\mathrm {n e}})} q ^ {L _ {0} ^ {\mathrm {n e}}} e ^ {2 \pi i J _ {0} ^ {\{h \}}} = \prod_ {\alpha \in S ^ {\prime}} \prod_ {n = 1} ^ {\infty} \left(1 - s (\alpha) e ^ {2 \pi i (- n K + \alpha | - \tau (D + x) + h)}\right) ^ {- s (\alpha)}.
$$

The lemma follows immediately from the last three identities.

✷

Note that the right hand side of (3.1) defines a meromorphic function on $Y$ with simple poles on the hyperplanes $T _ { \alpha } : = \{ h \in \widehat { { \mathfrak { h } } } | \alpha ( h ) = 0 \}$ , $\alpha \in \widehat { \triangle } _ { \mathrm { e v e n } } ^ { \mathrm { r e } }$ .

Let $M$ be a highest weight ${ \widehat { \mathfrak { g } } } .$ -module of level $k \neq - h ^ { \vee }$ . We shall assume that its character $\mathrm { c h } _ { M }$ extends to a meromorphic function in the whole upper half space $Y$ with at most simple poles at the hyperplanes $T _ { \alpha }$ , where $\alpha \in \widehat { \triangle } _ { \mathrm { e v e n } } ^ { \mathrm { r e } }$ . (We conjecture that this is always the case.)

Let $H ( M )$ be the $W _ { k } ( { \mathfrak { g } } , x , f )$ -module defined in Section 2.5. Define the Euler-Poincar´e character of $H ( M )$ :

$$
\mathrm {c h} _ {H (M)} (h) = \sum_ {j \in {\bf Z}} (- 1) ^ {j} \mathrm {t r} _ {H _ {j} (M)} q ^ {L _ {0}} e ^ {2 \pi i J _ {0} ^ {\{h \}}},
$$

where $h \in \mathfrak { h } ^ { f }$ (see Theorem 2.4). We have the following formula for this character:

$$
\operatorname {c h} _ {H (M)} (h) = q ^ {\frac {\widehat {\Omega} | M}{2 (k + h ^ {\vee})}} \lim  _ {\epsilon \rightarrow 0} \left(\operatorname {c h} _ {M} \prod_ {\alpha \in \widehat {S} \backslash \widehat {S} ^ {\prime}} (1 - s (\alpha) e ^ {- \alpha}) ^ {s (\alpha)}\right) (\tau , \tau (\epsilon b - x) + h, 0). \tag {3.2}
$$

Indeed, by the Euler-Poincar´e principle we have

$$
\begin{array}{l} \operatorname {c h} _ {H (M)} (h) = \lim _ {\epsilon \to 0} \sum_ {j \in \mathbf {Z}} (- 1) ^ {j} \operatorname {t r} _ {\mathcal {C} _ {j} (M)} q ^ {L _ {0} + \epsilon (b + b _ {0} ^ {\mathrm {c h}})} e ^ {2 \pi i J _ {0} ^ {\{h \}}} \\ = q ^ {\frac {\widehat {\Omega} | M}{2 (k + h \vee)}} \lim _ {\epsilon \to 0} \{\mathrm {t r} _ {M} q ^ {- D + \epsilon b - x} e ^ {2 \pi i h} \sum_ {j \in \mathbf {Z}} (- 1) ^ {j} \mathrm {t r} _ {F _ {j}} q ^ {L _ {0} ^ {\mathrm {c h}} + \epsilon b _ {0} ^ {\mathrm {c h}} + L _ {0} ^ {\mathrm {n e}}} e ^ {2 \pi i (h _ {0} ^ {\mathrm {n e}} + h _ {0} ^ {\mathrm {c h}})} \}. \\ \end{array}
$$

Now (3.2) follows from Lemma 3.1.

Introduce the Weyl denominator

$$
\widehat {R} = \prod_ {\alpha \in \widehat {\Delta} _ {+}} (1 - s (\alpha) e ^ {- \alpha}) ^ {s (\alpha) \mathrm {m u l t} \alpha}.
$$

Rewriting the RHS of (3.2) using $\widehat { R }$ , we arrive at the following result.

Theorem 3.1 Let $M$ be the highest weight $\widehat { \mathfrak { g } }$ -module with the highest weight Λ of level $k \neq - h ^ { \vee }$ , and suppose that $\mathrm { c h } _ { M }$ extends to a meromorphic function on $Y$ with at most simple poles at the hyperplanes $T _ { \alpha }$ , where $\alpha \in \widehat { \triangle } _ { \mathrm { e v e n } } ^ { \mathrm { r e } }$ . Then

$$
\begin{array}{l} \operatorname {c h} _ {H (M)} (h) = \frac {q ^ {\frac {(\Lambda | \Lambda + 2 \widehat {\rho})}{2 (k + h ^ {\vee})}}}{\prod_ {j = 1} ^ {\infty} (1 - q ^ {j}) ^ {\dim \mathfrak {h}}} (\widehat {R} \operatorname {c h} _ {M}) (H) \\ \times \prod_ {n = 1} ^ {\infty} \prod_ {\alpha \in \triangle_ {+}, (\alpha | x) = 0} \left((1 - s (\alpha) e ^ {- (n - 1) K - \alpha}) ^ {- s (\alpha)} (1 - s (\alpha) e ^ {- n K + \alpha}) ^ {- s (\alpha)}\right) (H), \\ \times \prod_ {n = 1} ^ {\infty} \prod_ {\alpha \in \Delta_ {+}, (\alpha | x) = \frac {1}{2}} (1 - s (\alpha) e ^ {- n K + \alpha}) ^ {- s (\alpha)} (H), \tag {3.3} \\ \end{array}
$$

where, as before, $s ( \alpha ) = ( - 1 ) ^ { p ( \alpha ) }$ and $H : = ( \tau , - \tau x + h , 0 ) = 2 \pi i ( - \tau D - \tau x + h ) , h \in \mathfrak { h } ^ { f } .$

Remark 3.1. Here is a slightly more explicit expression for $\operatorname { c h } _ { H ( M ) }$

$$
\begin{array}{l} \operatorname {c h} _ {H (M)} (h) = \frac {q ^ {\frac {(\Lambda | \Lambda + 2 \hat {\rho})}{2 (k + h ^ {\vee})}}}{\prod_ {j = 1} ^ {\infty} (1 - q ^ {j}) ^ {\dim \mathfrak {h}}} (\widehat {R} \operatorname {c h} _ {M}) (\tau , - \tau x + h, 0) \\ \times \prod_ {n = 1} ^ {\infty} \prod_ {\alpha \in \triangle_ {+}, (\alpha | x) = \frac {1}{2}} (1 - s (\alpha) q ^ {n - \frac {1}{2}} e ^ {2 \pi i (\alpha | h)}) ^ {- s (\alpha)} \\ \times \prod_ {n = 1} ^ {\infty} \prod_ {\alpha \in \triangle_ {+}, (\alpha | x) = 0} (1 - s (\alpha) q ^ {n - 1} e ^ {- 2 \pi i (\alpha | h)}) ^ {- s (\alpha)} (1 - s (\alpha) q ^ {n} e ^ {2 \pi i (\alpha | h)}) ^ {- s (\alpha)}. \\ \end{array}
$$

Since we may assume that $( \gamma _ { i } | x ) \geq 0$ , for a set of simple roots $\{ \gamma _ { i } \}$ of $\Delta _ { + \mathrm { e v e n } }$ , is easy to show that if the set $\{ \alpha \in \Delta _ { + \mathrm { e v e n } } | ( \alpha | x ) = 0 \}$ is non-empty, then the restriction of each $\alpha$ from this set to ${ \mathfrak { h } } ^ { f }$ is a non-zero linear function.

# 3.2 Conditions of non-vanishing of $H ( M )$

Using Theorem 3.1, we can establish a necessary and sufficient condition for $\operatorname { c h } _ { H ( M ) }$ to be not identically zero, hence a sufficient condition for the non-vanishing of $H ( M )$ .

Theorem 3.2 Let M be as in Theorem 3.1. Then $\operatorname { c h } _ { H ( M ) }$ is not identically zero if and only if the $\widehat { \mathfrak { g } }$ -module M is not locally nilpotent with respect to all root spaces ${ \mathfrak { g } } _ { - \alpha }$ , where $\alpha$ are positive even real roots satisfying the following three properties:

$$
(i) (\alpha | D + x) = 0, \quad (i i) (\alpha | \mathfrak {h} ^ {f}) = 0, \quad (i i i) | (\alpha | x) | \geq 1.
$$

In particular, these conditions guarantee that $H ( M ) \neq 0$ .

Lemma 3.2 Let $\alpha \in \widehat { \triangle } _ { + \mathrm { e v e n } } ^ { \mathrm { r e } }$ . Then the function ch $M$ is analytic on a non-empty open subset of the hyperplane $T _ { \alpha }$ if and only if ${ \widehat { \mathfrak { g } } } _ { - \alpha }$ is locally nilpotent on $M$ .

Proof. If ${ \widehat { \mathfrak { g } } } _ { - \alpha }$ is locally nilpotent on $M$ , then $r _ { \alpha } \mathrm { c h } _ { M } = \mathrm { c h } _ { M }$ , where $r _ { \alpha }$ is a reflection with respect to the hyperplane $T _ { \alpha }$ [K3]. Hence $Y ( M )$ is an $r _ { \alpha }$ -invariant convex domain and therefore, $Y ( M ) \cap T _ { \alpha }$ contains a non-empty open set (it is because any segment connecting $a$ and $r _ { \alpha } a$ , where $a \in Y _ { > }$ , has a non-empty intersection with the hyperplane $T _ { \alpha }$ ).

Conversely, suppose that ${ \widehat { \mathfrak { g } } } _ { - \alpha }$ is not locally nilpotent on $M$ . Consider $s \ell _ { 2 } \subset { \mathfrak { g } }$ generated by ${ \widehat { \mathfrak { g } } } _ { - \alpha }$ and ${ \widehat { \mathfrak { g } } } _ { \alpha }$ , and let $M _ { \mathrm { { i n t } } }$ denote the subspace of $M$ consisting of locally finite vectors with respect to this sl2. Then $\mathrm { c h } _ { M _ { \mathrm { i n t } } }$ is $r _ { \alpha }$ -invariant, hence (as above) it is analytic on an open subset of $T _ { \alpha }$ . On the other hand, $\mathrm { c h } _ { M / M _ { \mathrm { i n t } } }$ is a sum of functions of the form $\frac { e ^ { \lambda } } { 1 - e ^ { - \alpha } }$ eλ , where $\lambda$ is a weight of $M$ . Hence $\begin{array} { r } { \mathrm { c h } _ { M } = \mathrm { c h } _ { M _ { \mathrm { i n t } } } + \frac { f } { 1 - e ^ { - \alpha } } } \end{array}$ , where $f$ is a meromorphic function on $Y$ , which is analytic and non-zero on a non-empty open subset of $T _ { \alpha }$ . ✷

Proof of Theorem 3.2. It follows from Theorem 3.1 and Lemma 3.2 that $\operatorname { c h } _ { H ( M ) }$ is not identically zero if and only if $\widehat { R } \mathrm { c h } _ { M }$ cannot be decomposed as the product of $1 - e ^ { - \alpha }$ and a meromorphic function which is analytic in a non-zero open subset of $T _ { \alpha }$ for each positive even real root $\alpha$ such that $( \alpha | H ) =$ 0 and $\alpha \not \in \{ n K - \gamma | ( \gamma | x ) = 0$ or ${ \scriptstyle { \frac { 1 } { 2 } } } \cup \{ n K + \gamma | ( \gamma | x ) = 0 \}$ . But $( \alpha | H ) = 2 \pi i ( - \tau ( \alpha | D + x ) + ( \alpha | h ) )$ , hence $( \alpha | H ) = 0$ is equivalent to $( i )$ and $( i i )$ . The second condition on $\alpha$ is equivalent to $( i i i )$ . Hence $\operatorname { c h } _ { H ( M ) }$ is not identically zero if and only if conditions $( i ) - ( i i i )$ hold. ✷

A $\widehat { \mathfrak { g } }$ -module $M$ is called non-degenerate if each ${ \widehat { \mathfrak { g } } } _ { - \alpha }$ , where $\alpha$ is a positive real even root satisfying properties (i)—(iii) (in Theorem 3.2), is not locally nilpotent on $M$ . Otherwise $M$ is called degenerate.

# 3.3 Admissible highest weight $\widehat { \mathfrak { g } }$ -modules

Fix a non-degenerate invariant bilinear form (.|.) on $\mathfrak { g }$ such that all $( { \boldsymbol { \alpha } } | { \boldsymbol { \alpha } } ) \in \mathbf { R }$ for $\alpha \in \Delta$ . Then we have a decomposition of the set of even roots $\Delta _ { \overline { { 0 } } }$ into a disjoint union of $\Delta _ { \overline { { 0 } } } ^ { > }$ and $\Delta _ { \overline { { 0 } } } ^ { < }$ , where $\Delta _ { \overline { { 0 } } } ^ { > }$ (resp. $\Delta _ { \overline { { 0 } } } ^ { < }$ <) is the set of $\alpha \in \Delta _ { \overline { { 0 } } }$ such that $( \alpha | \alpha ) > 0$ (resp. $< 0$ ). Let ${ \mathfrak { g } } _ { \overline { { 0 } } } ^ { > }$ be the semisimple subalgebra of the reductive Lie algebra ${ \mathfrak { g } } _ { \overline { { 0 } } }$ with root system $\Delta _ { \overline { { 0 } } } ^ { > }$ , and let $\widehat { \mathfrak { g } } _ { 0 } ^ { > }$ be the affine subalgebra of g associated to g>. $\widehat { \mathfrak { g } }$ ${ \mathfrak { g } } _ { \overline { { 0 } } } ^ { > }$

Recall that a $\widehat { \mathfrak { g } }$ -module $L ( \Lambda )$ is called integrable if it is integrable with respect to $\widehat { \mathfrak { g } } _ { 0 } ^ { > }$ and is locally finite with respect to $\mathfrak { g }$ . In [KW4] a complete classification of integrable ${ \widehat { \mathfrak { g } } } .$ -modules was obtained.

Definition. ([KW1], [KW2], [KW4]). Let $\widehat { \Delta } ^ { \prime } \subset \widehat { \Delta }$ be a subset such that $\mathbf { Q } \widehat { \Delta } ^ { \prime } = \mathbf { Q } \widehat { \Delta }$ and $\widehat { \Delta } ^ { \prime }$ is isomorphic to a set of roots of an affine superalgebra ${ \widehat { \mathfrak { g } } } ^ { \prime }$ (which is not necessarily a subalgebra of $\widehat { \mathfrak { g } }$ ). Let $\widehat { \Pi } ^ { \prime } \subset \widehat { \Delta } _ { + }$ be the set of simple roots of $\widehat { \Delta } ^ { \prime }$ (for the subset of positive roots $\widehat { \Delta } \cap \widehat { \Delta } _ { + }$ ). Let $\widehat { \rho } ^ { \prime } \in { \mathfrak { h } } ^ { \prime }$ be the Weyl vector, i.e., $2 ( \widetilde { \rho } ^ { \prime } | \alpha ^ { \prime } ) = ( \alpha ^ { \prime } | \alpha ^ { \prime } )$ for all $\alpha ^ { \prime } \in \widehat { \Pi } ^ { \prime }$ . A $\widehat { \mathfrak { g } }$ -module $L ( \Lambda )$ (and the weight $\Lambda$ ) is called admissible for $\Delta ^ { \prime }$ if the ${ \widehat { \mathfrak { g } } } ^ { \prime }$ -module $L ^ { \prime } ( \Lambda + \widehat { \rho } - \widehat { \rho } ^ { \prime } )$ is integrable and this condition does not hold for any $\widehat { \Delta ^ { \prime \prime } } \frac { \supset } { \neq } \widehat { \Delta ^ { \prime } }$ . It is called principal admissible if $\widehat { \Delta } ^ { \prime }$ is isomorphic to $\widehat { \Delta }$ .

Conjecture 3.3A.[KW2], [KW4] The character of an admissible $\widehat { \mathfrak { g } }$ -module $L ( \Lambda )$ is related to the character of an integrable ${ \widehat { \mathfrak { g } } } ^ { \prime }$ -module by the formula:

$$
e ^ {\widehat {\rho}} \widehat {R} \mathrm {c h} _ {L (\Lambda)} = e ^ {\widehat {\rho} ^ {\prime}} \widehat {R} ^ {\prime} \mathrm {c h} _ {L ^ {\prime} (\Lambda + \widehat {\rho} - \widehat {\rho} ^ {\prime})}. (3. 4)
$$

Remark 3.3. Formula (3.4) holds for general symmetrizable Kac–Moody Lie algebra. It is immediate from the character formula for admissible modules [KW1], [KW2]. In fact (3.4) holds for these Lie algebras in the much more difficult case when “integrable” is replaced by “integral” [F].

Definition. [KW4]. A ${ \widehat { \mathfrak { g } } } .$ -module $L ( \Lambda )$ is called boundary admissible for $\widehat { \Delta } ^ { \prime }$ if $\Lambda + \widehat { \rho } - \widetilde { \rho } = 0$ (i.e., $\mathrm { d i m } L ^ { \prime } ( \Lambda + \widehat { \rho } - \widehat { \rho } ^ { \prime } ) = 1$ ).

Of course, (3.4) provides an explicit product formula for the boundary admissible $\widehat { \mathfrak { g } }$ -modules:

$$
\mathrm {c h} _ {\Lambda} = e ^ {\Lambda} \widehat {R} ^ {\prime} / \widehat {R}. (3. 5)
$$

Conjecture 3.3B. If $L ( \Lambda )$ is an admissible $\widehat { \mathfrak { g } }$ -module, then the $W _ { k } ( { \mathfrak { g } } , x , f )$ -module $H ( L ( \Lambda ) )$ is either zero or irreducible.

If Conjecture 3.3B holds, then Theorem 3.2 gives necessary and sufficient conditions for the vanishing of $H ( L ( \Lambda ) )$ .

# 4 Vertex Algebras $W _ { k } ( { \mathfrak { g } } , e _ { - \theta } )$ , where $\theta$ is a Highest Root

We now choose a subset of positive roots in the set of roots $\bigtriangleup$ such that the highest root $\theta$ (i.e., $\theta { + } \alpha$ is not a root for any positive root $\alpha$ ) is even. In this section, we shall classify all the examples of vertex algebras $W _ { k } ( { \mathfrak { g } } , f )$ where $f = e _ { - \theta }$ . Denote by $e = e _ { \theta }$ the root vector such that $( e | f ) = ( \theta | \theta ) ^ { - 1 }$ . Let $\begin{array} { r } { x = \frac { \theta } { ( \theta | \theta ) } } \end{array}$ , so that $\theta ( x ) = 1$ . Then $\langle e , x , f \rangle$ is an $s \ell _ { 2 }$ -triple. Furthermore, we have :

$$
S = S ^ {\prime} \cup \{\theta \}. \tag {4.1}
$$

Indeed, otherwise there exists an element $\alpha \in \triangle \backslash \{ \theta \}$ such that $\begin{array} { r } { \frac { 2 ( \alpha | \theta ) } { ( \theta | \theta ) } \ge 2 } \end{array}$ , hence $\alpha - 2 \theta \in \triangle$ . This is impossible since $\alpha - 2 \theta < - \theta$ .

Thus, the $\scriptstyle { \frac { 1 } { 2 } } \mathbf { Z }$ -gradation (2.1) of $\mathfrak { g }$ has the form:

$$
\mathfrak {g} = \mathfrak {g} _ {- 1} + \mathfrak {g} _ {- \frac {1}{2}} + \mathfrak {g} _ {0} + \mathfrak {g} _ {\frac {1}{2}} + \mathfrak {g} _ {1}, \text {w h e r e} \mathfrak {g} _ {- 1} = \mathbf {C} f, \mathfrak {g} _ {1} = \mathbf {C} e. \tag {4.2}
$$

One also has:

$$
\mathfrak {g} ^ {f} = \mathfrak {g} _ {- 1} + \mathfrak {g} _ {- \frac {1}{2}} + \mathfrak {g} _ {0} ^ {f}, \mathfrak {g} _ {0} = \mathfrak {g} _ {0} ^ {f} \oplus \mathbf {C} x, \tag {4.3}
$$

where

$$
\mathfrak {g} _ {0} ^ {f} = \{a \in \mathfrak {g} _ {0} | (a | x) = 0 \} = \mathfrak {h} ^ {f} \oplus (\oplus_ {\alpha \in \Delta_ {0}} \mathbf {C} e _ {\alpha}), \mathfrak {h} ^ {f} = \{h \in \mathfrak {h} | (h | x) = 0 \}, \Delta_ {0} = \{\alpha \in \Delta | (\alpha | x) = 0 \}.
$$

It is easy to see now that formula (2.5) for the central charge of the Virasoro algebra of $W _ { k } ( { \mathfrak { g } } , e _ { - \theta } )$ becomes:

$$
c = \frac {k}{k + h ^ {\vee}} \operatorname {s d i m} \mathfrak {g} - \frac {1 2 k}{(\theta | \theta)} + \frac {1}{4} (\operatorname {s d i m} \mathfrak {g} - \operatorname {s d i m} \mathfrak {g} _ {0} ^ {f}) - \frac {1 1}{4}. \tag {4.4}
$$

Furthermore, it is easy to see from Theorem 2.4(b) that all fields $J ^ { \{ v \} }$ , $v \in { \mathfrak { g } } _ { 0 } ^ { f }$ , of the vertex algebra $W _ { k } ( { \mathfrak { g } } , e _ { - \theta } )$ are primary (of conformal weight 1). Indeed, we have for $v \in { \mathfrak { g } } _ { 0 } ^ { f }$ :

$$
\operatorname {s t r} _ {\mathfrak {g} _ {+}} \operatorname {a d} v = \frac {1}{2} \operatorname {s t r} _ {\mathfrak {g}} (\operatorname {a d} v) (\operatorname {a d} x) = h ^ {\vee} (v | x)
$$

by (1.7). But $( v | x ) = - ( v | [ f , e ] ) = - ( [ v , f ] | e ) = 0$ . (Hence $( { \mathfrak { g } } ^ { f } | x ) = 0$ in any Dynkin gradation.)

Likewise, by Theorem 2.4(c), the 2-cocycle of the affine subalgebra $\big ( \mathfrak { g } _ { 0 } ^ { f } \big ) ^ { }$ of $W _ { k } ( { \mathfrak { g } } , e _ { - \theta } )$ equals:

$$
\alpha (v, v ^ {\prime}) = k (v | v ^ {\prime}) + \frac {1}{2} h ^ {\vee} (v | v ^ {\prime}) - \frac {1}{4} \operatorname {s t r} _ {\mathfrak {g} _ {0}} (\operatorname {a d} v) (\operatorname {a d} v ^ {\prime}). \tag {4.5}
$$

In the case when ${ \mathfrak { g } } _ { 0 } ^ { f }$ simple, denoting by $h _ { 0 } ^ { \vee }$ its dual Coxeter number for (.|.) restricted to ${ \mathfrak { g } } _ { 0 } ^ { f }$ , we can rewrite (4.5):

$$
\alpha (v, v ^ {\prime}) = \left(v \mid v ^ {\prime}\right) \left(k + \frac {1}{2} \left(h ^ {\vee} - h _ {0} ^ {\vee}\right)\right). \tag {4.6}
$$

The following proposition lists all vertex algebras $W _ { k } ( { \mathfrak { g } } , e _ { - \theta } )$ .

Proposition 4.1 All cases of $( { \mathfrak { g } } , \theta )$ along with the description of the ${ \mathfrak { g } } _ { 0 } ^ { f }$ -module ${ \mathfrak { g } } _ { \frac { 1 } { 2 } }$ are as follows: I . is a simple Lie algebra, and $\theta$ is the highest root. $\mathfrak { g }$

<table><tr><td>g</td><td>g0f</td><td>g1/2</td><td>g</td><td>g0f</td><td>g1/2</td></tr><tr><td>sln (n ≥ 3)</td><td>gln-2</td><td>Cn-2 ⊕ Cn-2*</td><td>F4</td><td>sp6</td><td>Λ03C6</td></tr><tr><td>so(n ≥ 5)</td><td>sl2 ⊕ so(n-4)</td><td>C2 ⊗ Cn-4</td><td>E6</td><td>sl6</td><td>Λ3C6</td></tr><tr><td>sp(n ≥ 2)</td><td>sp(n-2)</td><td>Cn-2</td><td>E7</td><td>so12</td><td>spin12</td></tr><tr><td>G2</td><td>sl2</td><td>S4C2</td><td>E8</td><td>E7</td><td>56-dim</td></tr></table>

II. g is a simple Lie superalgebra but not a Lie algebra, $s \ell _ { 2 }$ is a simple component of ${ \mathfrak { g } } _ { 0 }$ and $\theta$ is the highest root of this component. Below are all cases when ${ \mathfrak { g } } _ { 0 } ^ { f }$ is a Lie algebra ( $\stackrel { \prime } { m } \geq 1$ and g 1 ${ ^ { 9 } \mathrm { \frac { 1 } { 2 } } }$ is odd):

<table><tr><td>g</td><td>g0f</td><td>g1/2</td><td>g</td><td>g0f</td><td>g1/2</td></tr><tr><td>sl(2|m) (m≠2)</td><td>glm</td><td>Cm ⊕ Cm*</td><td>D(2,1;a)</td><td>sl2 ⊕ sl2</td><td>C2 ⊗ C2</td></tr><tr><td>sl(2|2)/CI</td><td>sl2</td><td>C2 ⊕ C2</td><td>F(4)</td><td>so7</td><td>spin7</td></tr><tr><td>spo(2|m)</td><td>so m</td><td>Cm</td><td>G(3)</td><td>G2</td><td>7-dim</td></tr><tr><td>osp(4|m)</td><td>sl2 ⊕ spm</td><td>C2 ⊗ Cm</td><td></td><td></td><td></td></tr></table>

III. g is a simple Lie superalgebra but not a a Lie algebra. The remaining possibilities are:

<table><tr><td>g</td><td>g0f</td><td>g1/2</td></tr><tr><td>sl(m|n) (m ≠ n, m &gt; 2)</td><td>gl(m - 2|n)</td><td>Cm-2|n ⊕ Cm-2|n*</td></tr><tr><td>sl(m|m)/CI (m &gt; 2)</td><td>sl(m - 2|m)</td><td>Cm-2|m ⊕ Cm-2|m*</td></tr><tr><td>spo(n|m) (n ≥ 4)</td><td>spo(n - 2|m)</td><td>Cn-2|m</td></tr><tr><td>osp(m|n) (m ≥ 5)</td><td>osp(m - 4|n) ⊕ sl2</td><td>Cm-4|n ⊗ C2</td></tr><tr><td>F(4)</td><td>D(2,1;2)</td><td>1
○ - ∩ - ○ ((6|4)-dim)</td></tr><tr><td>G(3)</td><td>osp(3|2)</td><td>-3
⊗→0
((4|4)-dim)</td></tr></table>

Proof. The proof of this proposition is straightforward by looking at all highest roots $\theta$ of simple components of ${ \mathfrak { g } } _ { 0 }$ and choosing an ordering for which this $\theta$ is the highest root of $\mathfrak { g }$ . ✷

All examples from Table I (resp. Table II) of Proposition 4.1 occur in the Fradkin–Linetsky list of quasisuperconformal (resp. superconformal) algebras [FL], but the last two examples from Table III are missing there.

One can check that in all cases of Table II when ${ \mathfrak { g } } _ { 0 } ^ { f }$ is simple one has:

$$
\frac {1}{2} (h ^ {\vee} - h _ {0} ^ {\vee}) = - 1,
$$

if we consider the normalization of the form (.|.) which restricts to the standard one on ${ \mathfrak { g } } _ { 0 } ^ { f }$ (i.e., $( \alpha | \alpha )$ =2 for a long root of ${ \mathfrak { g } } _ { 0 } ^ { f }$ ). Then $( \theta | \theta ) = 4$ for $o s p ( m | n ) , \ : = \ : - 3$ for $F ( 4 )$ , $\mathrm { a n d } = - \frac { 8 } { 3 }$ for $G ( 3 )$ .

It follows from (4.6) that the affine central charge in [FL] equals $k - 1$ in all these cases, and this leads to a perfect agreement of (4.4) with the Virasoro central charges in [FL].

In the case of Table I we take the usual normalization $( \theta | \theta ) = 2$ . Then (4.4) becomes

$$
c = \frac {k}{k + h ^ {\vee}} \dim \mathfrak {g} - 6 k + \frac {1}{4} (\dim \mathfrak {g} - \dim \mathfrak {g} ^ {f}) - \frac {1 1}{4}.
$$

If, in addition, ${ \mathfrak { g } } _ { 0 } ^ { f }$ is simple, then $h ^ { \vee } - h _ { 0 } ^ { \vee } = 1 , 6 , 8 , 1 2 , 5$ and $\textstyle { \frac { 1 0 } { 3 } }$ for $\mathfrak { g }$ of type $C _ { n } , E _ { 6 } , E _ { 7 } , E _ { 8 } , F _ { 4 }$ and $G _ { 2 }$ , respectively, and again we are in agreement with the Virasoro central charge of [FL].

Remark 4.1. Many examples of vertex algebras from Proposition 4.1 are well known:

$W _ { k } ( s \ell _ { 2 } , e _ { - \theta } )$ is the Virasoro vertex algebra,

$W _ { k } ( s \ell _ { 3 } , e _ { - \theta } )$ is the Bershadsky–Polyakov algebra [B],

$W _ { k } ( s p o ( 2 | 1 ) , e _ { - \theta } )$ is the Neveu–Schwarz algebra,

$W _ { k } ( s p o ( 2 | m ) , e _ { - \theta } )$ for $m \geq 3$ are the Bershadsky–Knizhnik algebras [BeK],

$W _ { k } ( s \ell ( 2 | 1 ) = s p o ( 2 | 2 ) , e _ { - \theta } )$ is the $N = 2$ superconformal algebra,

$W _ { k } ( s \ell ( 2 | 2 ) / { \bf C } I , e _ { - \theta } )$ is the $N = 4$ superconformal algebra,

$W _ { k } ( s p o ( 2 | 3 ) , e _ { - \theta } )$ tensored with one fermion is the $N = 3$ superconformal algebra (cf. [GS]),

$W _ { k } ( D ( 2 , 1 ; a ) , e _ { - \theta } )$ tensored with four fermions and one boson is the big $N = 4$ superconformal algebra (cf. [GS]).

# 5 The example of $\hat { s } \hat { \ell _ { 2 } }$ and Virasoro algebra

(See [KW1, FKW] for details).

Let ${ \mathfrak { g } } = s \ell _ { 2 }$ with the invariant bilinear form $( a | b ) = \operatorname { t r } a b$ ; then $\Delta _ { + } = \{ \alpha \}$ . All possibilities for Πb′ are as follows: $\widehat { \Pi } ^ { \prime }$

$$
\widehat {\Pi} _ {u, j} = \left\{(u - j) K - \alpha , j K + \alpha \right\} \text {w h e r e} 0 \leq j \leq u - 1, u \geq 1.
$$

All possible levels of the admissible weights for $\widehat { \Pi } _ { u , j }$ are rational numbers $k$ with a positive denominator $u$ (relatively prime to the numerator) such that $u ( k + 2 ) \geq 2$ . The set of all admissible weights of such a level $k$ is

$$
\left\{\Lambda_ {k, j, n} = k D + \frac {1}{2} (n - j (k + 2)) \alpha | 0 \leq j \leq u - 1, 0 \leq n \leq u (k + 2) - 2 \right\}.
$$

Such a weight is degenerate iff it is integrable with respect to the root $K - \alpha$ , which happens iff $j = u - 1$ . In particular, all such weights corresponding to $u = 1$ are degenerate.

We have: $W _ { k } ( s \ell _ { 2 } , e _ { - \alpha } )$ is generated by the Virasoro field $L ( z )$ . Furthermore, by Theorem 3.2, $H ( L ( \Lambda _ { k , j , n } ) )$ is zero iff $j = u - 1$ . Otherwise, $H ( L ( \Lambda _ { k , j , n } ) ) = H _ { 0 } ( L ( \Lambda _ { k , j , n } ) )$ is an irreducible highest weight module over the Virasoro algebra defined by $L ( z )$ (given by (2.5)), corresponding to the parameters $p = u ( k + 2 )$ , $p ^ { \prime } = u$ of the so-called minimal series:

$$
c ^ {(p, p ^ {\prime})} = 1 - 6 \frac {(p - p ^ {\prime}) ^ {2}}{p p ^ {\prime}}, h _ {j + 1, n + 1} ^ {(p, p ^ {\prime})} = \frac {(p (j + 1) - p ^ {\prime} (n + 1)) ^ {2} - (p - p ^ {\prime}) ^ {2}}{4 p p ^ {\prime}}.
$$

Here $p , p ^ { \prime } \in \mathbf { Z }$ , $2 \leq p ^ { \prime } < p$ , $g c d ( p , p ^ { \prime } ) = 1$ , $1 \leq j + 1 \leq p ^ { \prime } - 1$ , $1 \leq n + 1 \leq p - 1$ , which are precisely all minimal series Virasoro modules. The character formula for ${ \cal M } = { \cal L } ( \Lambda _ { k , j , n } )$ plugged in (3.3) gvector $\tilde { v } _ { \Lambda _ { k , j , n } }$ e well-known characters of the minimal series modu(see Remark 2.3) is the eigenvector with the lowest $L _ { 0 }$ over the Virasoro alg-eigenvalue (equal to $h _ { j + 1 , n + 1 } ^ { ( p , p ^ { \prime } ) }$

# 6 The Example of spo(2|1)b and Neveu–Schwarz algebra

In this section, ${ \mathfrak { g } } = { \mathrm { s p o } } ( 2 | 1 )$ with the invariant bilinear for $( a | b ) = { \textstyle { \frac { 1 } { 2 } } } \mathrm { s t r } a b$ . This is a 5-dimensional Lie superalgebra with the basis consisting of odd elements $e _ { \alpha } , e _ { - \alpha }$ and even elements $\begin{array} { r l } { e _ { 2 \alpha } } & { { } = } \end{array}$ $[ e _ { \alpha } , e _ { \alpha } ] , e _ { - 2 \alpha } = - [ e _ { - \alpha } , e _ { - \alpha } ]$ and $h = 2 [ e _ { \alpha } , e _ { - \alpha } ]$ such that $[ h , e _ { \alpha } ] = e _ { \alpha }$ , $[ h , e _ { - \alpha } ] = - e _ { - \alpha }$ . Then $[ h , e _ { 2 \alpha } ] = 2 e _ { 2 \alpha }$ , $| h , e _ { - 2 \alpha } | = - 2 e _ { - 2 \alpha }$ , $[ e _ { 2 \alpha } , e _ { - 2 \alpha } ] = h$ , $[ e _ { \alpha } , e _ { - 2 \alpha } ] = e _ { - \alpha }$ , $[ e _ { - \alpha } , e _ { 2 \alpha } ] = e _ { \alpha }$ ; $( e _ { \alpha } | e _ { - \alpha } ) =$ $\begin{array} { r } { ( e _ { 2 \alpha } | e _ { - 2 \alpha } ) = \frac { 1 } { 2 } } \end{array}$ , $( h | h ) = 1$ . We have: $h ^ { \vee } = 3$ and $\Delta _ { + } = \{ \alpha , 2 \alpha \}$ . The element $f = 2 e _ { - 2 \alpha }$ is the only, up to conjugacy, nilpotent even element, and then $x = { \textstyle { \frac { 1 } { 2 } } } h$ .

We have the charged free superfermions $\varphi _ { j \alpha } = \varphi _ { j \alpha } ( z )$ and $\varphi _ { j \alpha } ^ { * } = \varphi _ { j \alpha } ^ { * } ( z )$ , $j = 1 , 2$ , and the neutral free fermion $\Phi = \Phi ( z )$ such that $[ \Phi _ { \lambda } \Phi ] = 1$ (since $( f | [ e _ { \alpha } , e _ { \alpha } ] ) = 1 _ { , }$ ). Hence we have:

$$
d = d (z) = - e _ {\alpha} \varphi_ {\alpha} ^ {*} + e _ {2 \alpha} \varphi_ {2 \alpha} ^ {*} - \frac {1}{2}: \varphi_ {2 \alpha} (\varphi_ {\alpha} ^ {*}) ^ {2}: + \varphi_ {2 \alpha} ^ {*} + \varphi_ {\alpha} ^ {*} \Phi ,
$$

and the $\lambda$ -brackets of $d$ with all generators of the complex $\mathcal { C } ( { \mathfrak { g } } , f , k )$ are:

$$
[ d _ {\lambda} e _ {2 \alpha} ] = 0, [ d _ {\lambda} e _ {\alpha} ] = - e _ {2 \alpha} \varphi_ {\alpha} ^ {*}, [ d _ {\lambda} h ] = e _ {\alpha} \varphi_ {\alpha} ^ {*} - 2 e _ {2 \alpha} \varphi_ {2 \alpha} ^ {*},
$$

$$
[ d _ {\lambda} e _ {- \alpha} ] = - \frac {1}{2} h \varphi_ {\alpha} ^ {*} + e _ {\alpha} \varphi_ {2 \alpha} ^ {*} - \frac {k}{2} (\partial + \lambda) \varphi_ {\alpha} ^ {*}, [ d _ {\lambda} e _ {- 2 \alpha} ] = - e _ {- \alpha} \varphi_ {\alpha} ^ {*} + h \varphi_ {2 \alpha} ^ {*} + \frac {k}{2} (\partial + \lambda) \varphi_ {2 \alpha} ^ {*},
$$

$$
[ d _ {\lambda} \varphi_ {2 \alpha} ] = e _ {2 \alpha} + 1, [ d _ {\lambda} \varphi_ {\alpha} ] = e _ {\alpha} + e _ {2 \alpha} \varphi_ {\alpha} ^ {*} - \Phi , [ d _ {\lambda} \varphi_ {2 \alpha} ^ {*} ] = - \frac {1}{2} (\varphi_ {\alpha} ^ {*}) ^ {2}, [ d _ {\lambda} \varphi_ {\alpha} ^ {*} ] = 0, [ d _ {\lambda} \Phi ] = \varphi_ {\alpha} ^ {*}.
$$

Since : $\Phi \Phi : = 0$ , we have:

$$
J ^ {(h)} (z) = h (z) -: \varphi_ {\alpha} \varphi_ {\alpha} ^ {*}: + 2: \varphi_ {2 \alpha} \varphi_ {2 \alpha} ^ {*}:, J ^ {(e _ {- \alpha})} (z) = e _ {- \alpha} (z) - \varphi_ {\alpha} \varphi_ {2 \alpha} ^ {*}, J ^ {(e _ {- 2 \alpha})} (z) = e _ {- 2 \alpha} (z).
$$

It is not difficult to check that the following fields are closed under $d _ { 0 }$ :

$$
{ G } { = } { \frac { 2 } { ( k + 3 ) ^ { 1 / 2 } } ( J ^ { ( e _ { - \alpha } ) } + \frac { 1 } { 2 } \Phi J ^ { ( h ) } + \frac { k + 2 } { 2 } \partial \Phi ) , }
$$

$$
L = \frac {2}{k + 3} (- J ^ {(e _ {- 2 \alpha})} - \Phi J ^ {(e _ {- \alpha})} + \frac {1}{4}: J ^ {(h)} J ^ {(h)}: + \frac {k + 2}{4} \partial J ^ {(h)}) - \frac {1}{2}: \Phi \partial \Phi:,
$$

and that the field $L$ is equal to the Virasoro field, defined by (2.5), modulo the image of $d _ { 0 }$ so that they define the same field of $W _ { k } ( { \mathfrak { g } } , f )$ .

Furthermore, a direct calculation with $\lambda$ -brackets in $W _ { k } ( { \mathfrak { g } } , f )$ shows that $L$ and $G$ form the Neveu–Schwarz algebra with central charge $c$ :

$$
[ L _ {\lambda} L ] = (\partial + 2 \lambda) L + \frac {\lambda^ {3}}{1 2} c, [ L _ {\lambda} G ] = (\partial + \frac {3}{2} \lambda) G, [ G _ {\lambda} G ] = 2 L + \frac {\lambda^ {2}}{3} c, \tag {6.1}
$$

$$
c = \frac {3}{2} \left(1 - \frac {2 (k + 2) ^ {2}}{k + 3}\right). \tag {6.2}
$$

The set of positive roots of $\widehat { \mathfrak { g } }$ is $( n \in \mathbf { Z } )$ ):

$$
\widehat {\Delta} _ {+} = \{n K | n > 0 \} \cup \{j \alpha + n K | n \geq 0, j = 1, 2 \}, \cup \{- j \alpha + n K | n > 0, j = 1, 2 \},
$$

and the set of simple roots is

$$
\widehat {\Pi} = \left\{\alpha_ {0} = K - \alpha , \alpha_ {1} = 2 \alpha \right\}.
$$

All possibilities for the sets $\widehat { \Pi ^ { \prime } }$ of simple roots of subsets $\widehat { \Delta } _ { + } ^ { \prime }$ of $\widehat { \Delta } _ { + }$ that are isomorphic to a set of positive roots of an affine superalgebra, are of three types: the principal ones (isomorphic to $\widehat { \Pi }$ ), the even type ones, isomorphic to the set of simple roots of type ${ A } _ { 1 } ^ { ( 1 ) }$ , and the subprincipal ones, isomorphic to the set of simple roots of the twisted affine superalgebra $C ^ { ( 2 ) } ( 2 )$ [K2].

All admissible weights for $\widehat { \mathfrak { g } }$ are of the form:

$$
\Lambda_ {k, j, n} = 2 (n - j (k + 3)) \Lambda_ {0} + (\frac {1}{2} k - n + j (k + 3)) \Lambda_ {1},
$$

where $\Lambda _ { 0 } , \Lambda _ { 1 }$ are the fundamental weights, $\textstyle k = { \frac { v } { u } } \in \mathbf { Q }$ is the level (u, v ∈ Z , $u \geq 1$ , $\operatorname* { g c d } ( u , v ) = 1$ ), and $j , n \in \frac 1 2 \mathbf { Z } _ { + }$ . The ranges of $k$ and $j , n$ are described below.

All principal admissible weights have level $k$ such that its denominator $u$ is a (positive) odd integer, $v$ is an even integer, and $u ( k + 3 ) \geq 3$ . Both $j , n$ are integers satisfying the following conditions:

(i)   
(ii) u+1 ≤ j ≤ u − 1 and u(k+3)+1

In case (i), $\widehat { \Pi } ^ { \prime } = \{ j K + \alpha _ { 0 }$ , $( u - 1 - 2 j ) K + \alpha _ { 1 } \}$ . Hence, by Theorem 3.2, the principal admissible weight $\Lambda _ { k , j , n }$ is degenerate iff $\begin{array} { r } { j = \frac { u - 1 } { 2 } } \end{array}$ .

In case (ii), $\widehat { \Pi } ^ { \prime } = \{ ( u - j ) K - \alpha _ { 0 }$ , $( 2 j + 1 - u ) K - \alpha _ { 1 } \big \}$ and all the admissible weights are non-degenerate.

For the even type admissible weights, $u$ is even and $v$ is odd, and $u ( k { + } 3 ) \geq 2$ . Both $j$ , $n \in { \frac { 1 } { 2 } } + \mathbf { Z }$ and satisfy the inequalities: $\begin{array} { r } { 0 < j \le u - \frac { 1 } { 2 } } \end{array}$ , $0 < n < u ( k + 3 ) - 1$ . In this case $\widehat { \Pi } = \{ ( 2 j + 1 ) K -$ $\alpha _ { 1 }$ , $( 2 ( u - j ) - 1 ) K + \alpha _ { 1 } \}$ , and $\Lambda _ { k , j , n }$ is degenerate iff $\begin{array} { r } { j = u - \frac { 1 } { 2 } } \end{array}$ .

For the subprincipal admissible weights, both $u$ and $v$ are odd integers, and $u ( k + 3 ) \geq 1$ . Both $j , n$ are integers, satisfying the inequalities: $0 \leq j \leq u - 1$ , $0 \leq n \leq u ( k + 3 ) - 1$ . In this case $\widehat { \Pi } ^ { \prime } = \{ j K + \alpha _ { 0 }$ , $( u - j ) K - \alpha _ { 0 } \}$ and all the admissible weights are non-degenerate.

The characters of all admissible $s p o ( 2 | 1 )$ -modules are known [KW1]. Applying to them Theorem 3.1 we obtain the well known characters of all minimal series modules of the Neveu-Schwarz algebra (see e.g. [KW1]).

Recall that these minimal series correspond to central charges equal

$$
c ^ {(p, p ^ {\prime})} = \frac {3}{2} \left(1 - \frac {2 (p - p ^ {\prime}) ^ {2}}{p p ^ {\prime}}\right), \tag {6.3}
$$

where $p , p ^ { \prime } \in \mathbf { Z }$ , $2 \leq p ^ { \prime } < p$ , $p - p ^ { \prime } \in 2 \mathbf { Z }$ , $\begin{array} { r } { g c d \left( \frac { p - p ^ { \prime } } { 2 } , p ^ { \prime } \right) = 1 } \end{array}$ , and the minimal eigenvalue of $L _ { 0 }$ equals

$$
h _ {r, s} ^ {(p, p ^ {\prime})} = \frac {(p r - p ^ {\prime} s) ^ {2} - (p - p ^ {\prime}) ^ {2}}{8 p p ^ {\prime}}, \tag {6.4}
$$

where $r , s \in \mathbf { Z }$ , $1 \leq r \leq p ^ { \prime } - 1$ , $1 \leq s \leq p - 1$ , $r - s \in 2 \mathbf { Z }$ . The corresponding normalized character is as follows:

$$
\chi_ {r, s} ^ {(p, p ^ {\prime})} (\tau) = \frac {1}{\eta_ {1 / 2} (\tau)} \left(\theta_ {\frac {p r - p ^ {\prime} s}{2}, \frac {p p ^ {\prime}}{2}} (\tau) - \theta_ {\frac {p r + p ^ {\prime} s}{2}, \frac {p p ^ {\prime}}{2}} (\tau)\right), \tag {6.5}
$$

where η1/2(τ ) = η(τ /2)η(2τ)η(τ) and θn,m(τ ) = Pk∈Z+ n2m e $\begin{array} { r } { \eta _ { 1 / 2 } ( \tau ) = \frac { \eta ( \tau / 2 ) \eta ( 2 \tau ) } { \eta ( \tau ) } } \end{array}$ $\begin{array} { r } { \theta _ { n , m } ( \tau ) = \sum _ { k \in \mathbf { Z } + \frac { n } { 2 m } } e ^ { 2 \pi i m k ^ { 2 } \tau } } \end{array}$ . Another way of writing these characters, via the Weyl group $\widehat { W }$ of $\widehat { \mathfrak { g } }$ , is as follows:

$$
\chi_ {r, s} ^ {(p, p ^ {\prime})} (\tau) = \frac {1}{\eta_ {1 / 2} (\tau)} \sum_ {w \in \widehat {W}} \epsilon (w) q ^ {\frac {p p ^ {\prime}}{4} | \frac {w (\Lambda + \widehat {\rho})}{p} - \frac {\Lambda^ {\prime} + \widehat {\rho}}{p ^ {\prime}} | ^ {2}}, \tag {6.6}
$$

where $\Lambda + \widehat { \rho } = p \Lambda _ { 0 } + s \frac { \alpha _ { 1 } } { 2 }$ , $\Lambda ^ { \prime } + \widehat { \rho } = p ^ { \prime } \Lambda _ { 0 } + r \frac { \alpha _ { 1 } } { 2 }$ , $1 \leq s \leq p - 1$ , $1 \leq r \leq p ^ { \prime } - 1$

In the principal case we let $p = u ( k + 3 )$ , $p ^ { \prime } = u$ . Then (6.2) becomes $c = c ^ { ( p , p ^ { \prime } ) }$ , given by (6.3). Using Theorem 3.1 and (6.6) we obtain:

$$
q ^ {- c ^ {(p, p ^ {\prime})} / 2 4} \mathrm {c h} _ {H (L (\Lambda_ {k, j, n}))} = \left\{ \begin{array}{l l} \chi_ {p ^ {\prime} - 2 j - 2, p - 2 n - 2} ^ {(p, p ^ {\prime})} (\tau) \text {i n c a s e (i)} \\ \chi_ {2 j - p ^ {\prime}, 2 n - p} ^ {(p, p ^ {\prime})} (\tau) \text {i n c a s e (i i)} \end{array} \right.,
$$

so we get all characters of minimal series for which both $p$ and $p ^ { \prime }$ are odd.

In the even type cases and subprincipal cases we let $p = 2 u ( k + 3 )$ , $p ^ { \prime } = 2 u$ . Then again (6.2) becomes $c = c ^ { ( p , p ^ { \prime } ) }$ , and we obtain

$$
q ^ {- c ^ {(p, p ^ {\prime})} / 2 4} \mathrm {c h} _ {H (L (\Lambda_ {k, j, n}))} = \chi_ {2 p ^ {\prime} - 2 j - 1, 2 p - 2 n - 1} ^ {(p, p ^ {\prime})},
$$

so we get all characters of minimal series for which both $p$ and $p ^ { \prime }$ are even. Both $r$ and $s$ are either even (in the even type case) or odd (in the subprincipal case).

# 7 The example of $s \ell ( 2 | 1 )$ and $N = 2$ superconformal algebra

In this section, ${ \mathfrak { g } } = s \ell ( 2 | 1 )$ with the invariant bilinear form $( a | b ) = \operatorname { s t r } a b$ . This is the Lie superalgebra of traceless matrices in the superspace $\mathbf { C } ^ { 2 | 1 }$ whose even part is $\mathbf { C } \epsilon _ { 1 } + \mathbf { C } \epsilon _ { 3 }$ and odd part is $\mathbf { C } \epsilon _ { 2 }$ , where $\epsilon _ { 1 } , \epsilon _ { 2 } , \epsilon _ { 3 }$ is the standard basis. We shall denote by $E _ { i j }$ the standard basis of the space of matrices. We shall work in the following basis of $^ { 9 }$ :

$$
e _ {1} = E _ {1 2}, e _ {2} = E _ {2 3}, e _ {1 2} = - E _ {1 3}, f _ {1} = E _ {2 1}, f _ {2} = - E _ {3 2},
$$

$$
f _ {1 2} = - E _ {3 1}, h _ {1} = E _ {1 1} + E _ {2 2}, h _ {2} = - E _ {2 2} - E _ {3 3}.
$$

The elements $e _ { i } , f _ { i } , h _ { i }$ ( $i = 1 , 2$ ) are the Chevalley generators of $\mathfrak { g }$ [K1]. Elements $e _ { i } , f _ { i }$ ( $i = 1 , 2$ ) are all odd elements of $\mathfrak { g }$ . We pick the Cartan subalgebra $\mathfrak { h } = \mathbf { C } h _ { 1 } + \mathbf { C } h _ { 2 }$ . The simple roots $\alpha _ { 1 }$ and $\alpha _ { 2 }$ are the roots attached to $e _ { 1 }$ and $e _ { 2 }$ , and $\Delta _ { + } = \{ \alpha _ { 1 } , \alpha _ { 2 } , \alpha _ { 1 } + \alpha _ { 2 } \}$ . We have: $\alpha _ { i } = h _ { i }$ ( $i = 1 , 2$ ) (under the identification of $\mathfrak { h }$ with ${ \mathfrak { h } } ^ { * }$ ).

Since ${ \mathfrak { g } } _ { \bar { 0 } } = \mathbf { C } e _ { 1 2 } + \mathbf { C } f _ { 1 2 } + { \mathfrak { h } } \left( \simeq { g } \ell _ { 2 } \right)$ , there is only one, up to conjugacy, nilpotent element $f = f _ { 1 2 }$ , which embeds in the following $s \ell _ { 2 }$ -triple $\begin{array} { r } { \langle e = e _ { 1 2 } , x = \frac { 1 } { 2 } ( h _ { 1 } + h _ { 2 } ) , f \rangle } \end{array}$ . The corresponding $\scriptstyle { \frac { 1 } { 2 } } \mathbf { Z }$ -gradation looks as follows:

$$
\mathfrak {g} = \mathbf {C} f \oplus (\mathbf {C} f _ {1} + \mathbf {C} f _ {2}) \oplus \mathfrak {h} \oplus (\mathbf {C} e _ {1} + \mathbf {C} e _ {2}) \oplus \mathbf {C} e.
$$

We have: ${ \mathfrak { g } } ^ { f } = \mathbf { C } f + \mathbf { C } f _ { 1 } + \mathbf { C } f _ { 2 } + \mathbf { C } ( h _ { 1 } - h _ { 2 } )$ . There is only one other good $\scriptstyle { \frac { 1 } { 2 } } \mathbf { Z }$ -gradation (which is non-Dynkin). It will be considered after the discussion related to the Dynkin gradation is completed.

We have three pairs of charged free fermions: $\varphi _ { 1 } = \varphi _ { 1 } ( z )$ , $\varphi _ { 1 } ^ { * } = \varphi _ { 1 } ^ { * } ( z ) ,$ , $\varphi _ { 2 } = \varphi _ { 2 } ( z )$ , $\varphi _ { 2 } ^ { * } = \varphi _ { 2 } ^ { * } ( z )$ (which are even fields), and $\varphi _ { 1 2 } = \varphi _ { 1 2 } ( z )$ , $\varphi _ { 1 2 } ^ { * } = \varphi _ { 1 2 } ^ { * } ( z )$ (which are odd fields). There are two neutral free fermions: $\Phi _ { i } = \Phi _ { i } ( z )$ ( $i = 1 , 2$ ), they are odd, and their $\lambda$ -bracket is easily seen to be:

$$
\left[ \Phi_ {i \lambda} \Phi_ {j} \right] = - 1 \text {i f} i \neq j, = 0 \text {o t h e r w i s e}.
$$

Hence the field $d = d ( z )$ is as follows:

$$
d = - e _ {1} \varphi_ {1} ^ {*} - e _ {2} \varphi_ {2} ^ {*} + e _ {1 2} \varphi_ {1 2} ^ {*} + \varphi_ {1 2} \varphi_ {1} ^ {*} \varphi_ {2} ^ {*} + \varphi_ {1 2} ^ {*} + \varphi_ {1} ^ {*} \Phi_ {1} + \varphi_ {2} ^ {*} \Phi_ {2}.
$$

Its $\lambda$ -brackets with the generators of the complex $C ( { \mathfrak { g } } , f , k )$ are as follows:

$$
\left[ d _ {\lambda} e _ {1} \right] = e _ {1 2} \varphi_ {2} ^ {*}, \left[ d _ {\lambda} e _ {2} \right] = e _ {1 2} \varphi_ {1} ^ {*}, \left[ d _ {\lambda} e _ {1 2} \right] = 0,
$$

$$
\left[ d _ {\lambda} f _ {1} \right] = - h _ {1} \varphi_ {1} ^ {*} - e _ {1} \varphi_ {1 2} ^ {*} - (\partial + \lambda) k \varphi_ {1} ^ {*},
$$

$$
\left[ d _ {\lambda} f _ {2} \right] = h _ {2} \varphi_ {2} ^ {*} - e _ {1} \varphi_ {1 2} ^ {*} - (\partial + \lambda) k \varphi_ {2} ^ {*},
$$

$$
\left[ d _ {\lambda} f _ {1 2} \right] = f _ {2} \varphi_ {1} ^ {*} + f _ {1} \varphi_ {2} ^ {*} + \left(h _ {1} + h _ {2}\right) \varphi_ {1 2} ^ {*} + (\partial + \lambda) k \varphi_ {1 2} ^ {*},
$$

$$
\left[ d _ {\lambda} h _ {1} \right] = e _ {2} \varphi_ {2} ^ {*} - e _ {1 2} \varphi_ {1 2} ^ {*}, \left[ d _ {\lambda} h _ {2} \right] = e _ {1} \varphi_ {1} ^ {*} - e _ {1 2} \varphi_ {1 2} ^ {*},
$$

$$
\left[ d _ {\lambda} \varphi_ {1} \right] = e _ {1} - \varphi_ {1 2} \varphi_ {2} ^ {*} - \Phi_ {1}, \left[ d _ {\lambda} \varphi_ {2} \right] = e _ {2} - \varphi_ {1 2} \varphi_ {1} ^ {*} - \Phi_ {2},
$$

$$
\left[ d _ {\lambda} \varphi_ {1 2} \right] = e _ {1 2} + 1, \left[ d _ {\lambda} \varphi_ {j} ^ {*} \right] = 0, \left[ d _ {\lambda} \varphi_ {1 2} ^ {*} \right] = \varphi_ {1} ^ {*} \varphi_ {2} ^ {*},
$$

$$
\left[ d _ {\lambda} \Phi_ {1} \right] = - \varphi_ {2} ^ {*}, \left[ d _ {\lambda} \Phi_ {2} \right] = - \varphi_ {1} ^ {*}.
$$

Furthermore, we have the fields:

$$
J ^ {(h _ {1})} (z) = h _ {1} (z) -: \varphi_ {2} \varphi_ {2} ^ {*}: +: \varphi_ {1 2} \varphi_ {1 2} ^ {*}:,
$$

$$
J ^ {(h _ {2})} (z) = h _ {2} (z) -: \varphi_ {1} \varphi_ {1} ^ {*}: +: \varphi_ {1 2} \varphi_ {1 2} ^ {*}:,
$$

$$
J ^ {(f _ {1})} (z) = f _ {1} (z) +: \varphi_ {2} \varphi_ {1 2} ^ {*}:, J ^ {(f _ {2})} (z) = f _ {2} (z) +: \varphi_ {1} \varphi_ {1 2} ^ {*}:, J ^ {(f _ {1 2})} (z) = f _ {1 2} (z).
$$

One easily calculates the $\lambda$ -brackets of $d$ with these fields, using (2.4):

$$
\left[ d _ {\lambda} J ^ {(h _ {1})} \right] = \varphi_ {1 2} ^ {*} + \varphi_ {2} ^ {*} \Phi_ {2}, [ d _ {\lambda} J ^ {(h _ {2})} ] = \varphi_ {1 2} ^ {*} + \varphi_ {1} ^ {*} \Phi_ {1},
$$

$$
\left[ d _ {\lambda} J ^ {(f _ {1})} \right] = -: \varphi_ {1} ^ {*} J ^ {(h _ {1})}: + \varphi_ {1 2} ^ {*} \Phi_ {2} - (k + 1) (\partial + \lambda) \varphi_ {1} ^ {*},
$$

$$
\left[ d _ {\lambda} J ^ {(f _ {2})} \right] = -: \varphi_ {2} ^ {*} J ^ {(h _ {2})}: +: \varphi_ {1 2} ^ {*} \Phi_ {1}: - (k + 1) (\partial + \lambda) \varphi_ {2} ^ {*},
$$

$$
\left[ d _ {\lambda} J ^ {(f _ {1 2})} \right] =: \varphi_ {1} ^ {*} J ^ {(f _ {2})}: +: \varphi_ {2} ^ {*} J ^ {(f _ {1})}: +: \varphi_ {1 2} ^ {*} J ^ {(h _ {1} + h _ {2})}: + k (\partial + \lambda) \varphi_ {1 2} ^ {*}.
$$

Using this, one checks directly that the following fields are closed under $d _ { 0 }$ :

$$
\begin{array}{l} J = J ^ {\left(h _ {1} - h _ {2}\right)} +: \Phi_ {1} \Phi_ {2}:, \\ L = - \frac {1}{k + 1} \left(J ^ {\left(f _ {1 2}\right)} +: \Phi_ {1} J ^ {\left(f _ {1}\right)}: +: \Phi_ {2} J ^ {\left(f _ {2}\right)}: -: J ^ {\left(h _ {1}\right)} J ^ {\left(h _ {2}\right)}:\right) \\ + \frac {1}{2} \left(\partial J ^ {\left(h _ {1} + h _ {2}\right)} +: \Phi_ {1} \partial \Phi_ {2}: +: \Phi_ {2} \partial \Phi_ {1}:\right), \\ \end{array}
$$

$$
G ^ {+} = - \frac {1}{k + 1} \left(J ^ {\left(f _ {1}\right)} -: \Phi_ {2} J ^ {\left(h _ {1}\right)}\right) + \partial \Phi_ {2},
$$

$$
G ^ {-} = J ^ {(f _ {2})} -: \Phi_ {1} J ^ {(h _ {2})}: - (k + 1) \partial \Phi_ {1}.
$$

Moreover, one can show that the field $L$ coincides with the Virasoro field defined by (2.5), modulo to the image of $d _ { 0 }$ , and therefore they give the same field of $W _ { k } ( { \mathfrak { g } } , f )$ .

A direct calculation with $\lambda$ -brackets shows that $J$ , $L$ , $G ^ { + }$ and $G ^ { - }$ form the $N = 2$ superconformal algebra with central charge $c = - 3 ( 2 k + 1 )$ :

$$
[ L _ {\lambda} L ] = (\partial + 2 \lambda) L + \lambda^ {2} \frac {c}{1 2}, [ J _ {\lambda} J ] = \lambda \frac {c}{3},
$$

$$
\left[ G ^ {\pm} _ {\lambda} G ^ {\pm} \right] = 0, \left[ J _ {\lambda} G ^ {\pm} \right] = \pm G ^ {\pm}, \left[ G ^ {+} _ {\lambda} G ^ {-} \right] = L + \frac {1}{2} (\partial + 2 \lambda) J + \lambda^ {2} \frac {c}{6}, \tag {7.1}
$$

$$
[ L _ {\lambda} J ] = (\partial + \lambda) J, [ L _ {\lambda} G ^ {\pm} ] = (\partial + \frac {3}{2} \lambda) G ^ {\pm}.
$$

The good non-Dynkin $\scriptstyle { \frac { 1 } { 2 } } \mathbf { Z }$ -gradation looks as follows:

$$
\mathfrak {g} = \left(\mathbf {C} f _ {1 2} + \mathbf {C} f _ {2}\right) \oplus 0 \oplus \left(\mathbf {C} e _ {1} + \mathbf {C} f _ {1} + \mathfrak {h}\right) \oplus 0 \oplus \left(\mathbf {C} e _ {1 2} + \mathbf {C} e _ {2}\right).
$$

It corresponds to $x = h _ { 1 }$ . As before, we take $f = f _ { 1 2 }$

In this case we have two pairs of charged free fermions. Hence the field $d = d ( z )$ is as follows:

$$
d = - e _ {2} \varphi_ {2} ^ {*} + e _ {1 2} \varphi_ {1 2} ^ {*} + \varphi_ {1 2} ^ {*},
$$

and its $\lambda$ -brackets with the generators of the complex are as follows:

$$
\left[ d _ {\lambda} e _ {1} \right] = e _ {1 2} \varphi_ {2} ^ {*}, \left[ d _ {\lambda} e _ {2} \right] = 0, \left[ d _ {\lambda} e _ {1 2} \right] = 0,
$$

$$
\left[ d _ {\lambda} h _ {1} \right] = e _ {2} \varphi_ {2} ^ {*} - e _ {1 2} \varphi_ {1 2} ^ {*}, \left[ d _ {\lambda} h _ {2} \right] = - e _ {1 2} \varphi_ {1 2} ^ {*},
$$

$$
\left[ d _ {\lambda} f _ {1} \right] = - e _ {2} \varphi_ {1 2} ^ {*}, \left[ d _ {\lambda} f _ {2} \right] = - h _ {2} \varphi_ {2} ^ {*} + \left(h _ {1} + h _ {2}\right) \varphi_ {1 2} ^ {*} - k (\partial + \lambda) \varphi_ {2} ^ {*},
$$

$$
\left[ d _ {\lambda} f _ {1 2} \right] = f _ {1} \varphi_ {2} ^ {*} + \left(h _ {1} + h _ {2}\right) \varphi_ {1 2} ^ {*} + k (\partial + \lambda) \varphi_ {1 2} ^ {*},
$$

$$
[ d _ {\lambda} \varphi_ {2} ] = e _ {2}, [ d _ {\lambda} \varphi_ {2} ^ {*} ] = 0, [ d _ {\lambda} \varphi_ {1 2} ] = e _ {1 2} + 1, [ d _ {\lambda} \varphi_ {1 2} ^ {*} ] = 0.
$$

Furthermore, we have the fields:

$$
J ^ {(e _ {1})} (z) = e _ {1} (z) - \varphi_ {1 2} \varphi_ {2} ^ {*}, J ^ {(h _ {1})} (z) = h _ {1} (z) -: \varphi_ {2} \varphi_ {2} ^ {*}: +: \varphi_ {1 2} \varphi_ {1 2} ^ {*}:,
$$

$$
J ^ {(h _ {2})} (z) = h _ {2} (z) +: \varphi_ {1 2} \varphi_ {1 2} ^ {*}:,
$$

$$
J ^ {(f _ {1})} (z) = f _ {1} (z) +: \varphi_ {2} \varphi_ {1 2} ^ {*}:, J ^ {(f _ {2})} (z) = f _ {2} (z), J ^ {(f _ {1 2})} (z) = f _ {1 2} (z).
$$

One easily calculates the $\lambda$ -brackets of $d$ with these fields, using (2.4):

$$
\left[ d _ {\lambda} J ^ {\left(h _ {1}\right)} \right] = \left[ d _ {\lambda}, J ^ {\left(h _ {2}\right)} \right] = \varphi_ {1 2} ^ {*},
$$

$$
\left[ d _ {\lambda} J ^ {(e _ {1})} \right] = - \varphi_ {2} ^ {*}, \left[ d _ {\lambda} J ^ {(f _ {1})} \right] = 0,
$$

$$
\left[ d _ {\lambda} J ^ {(f _ {2})} \right] = -: \varphi_ {2} ^ {*} J ^ {(h _ {2})}: +: \varphi_ {1 2} ^ {*} J ^ {(e _ {1})}: - k (\partial + \lambda) \varphi_ {2} ^ {*},
$$

$$
\left[ d _ {\lambda} J ^ {\left(f _ {1 2}\right)} \right] =: \varphi_ {1 2} ^ {*} J ^ {\left(h _ {1} + h _ {2}\right)}: +: \varphi_ {2} ^ {*} J ^ {\left(f _ {1}\right)}: + k (\partial + \lambda) \varphi_ {1 2} ^ {*}.
$$

Using this one checks directly that the following fields are closed under $d _ { 0 }$ :

$$
J = J ^ {\left(h _ {1} - h _ {2}\right)},
$$

$$
L ^ {\prime} = - \frac {1}{k + 1} \left(J ^ {(f _ {1 2})} +: J ^ {(e _ {1})} J ^ {(f _ {1})}: -: J ^ {(h _ {1})} J ^ {(h _ {2})}:\right) + \frac {1}{2} \partial J ^ {(h _ {1} + h _ {2})},
$$

$$
G ^ {+} = - \frac {1}{k + 1} J ^ {(f _ {1})}, G ^ {-} = J ^ {(f _ {2})} -: J ^ {(h _ {2})} J ^ {(e _ {1})}: - k \partial J ^ {(e _ {1})}.
$$

A direct calculation with $\lambda$ -brackets shows that $J$ , $L ^ { \prime }$ , $G ^ { + }$ and $G ^ { - }$ form the $N = 2$ superconformal algebra with central charge $c = - 3 ( 2 k + 1 )$ . However, in this case the relation between $L ^ { \prime }$ and the field $L$ , defined by (2.5), is more complicated. One can show that in $W _ { k } ( { \mathfrak { g } } , x , f )$ one has:

$$
L = L ^ {\prime} + \frac {1}{2} \partial J. \tag {7.2}
$$

The four fields $J$ , $L$ , $G ^ { + }$ and $G ^ { - }$ form the Ramond type basis of $N = 2$ superconformal algebra ([RY, R]):

$$
[ L _ {\lambda} L ] = (\partial + 2 \lambda) L, [ J _ {\lambda} J ] = \lambda \frac {c}{3},
$$

$$
[ G ^ {\pm} _ {\lambda} G ^ {\pm} ] = 0, [ J _ {\lambda} G ^ {\pm} ] = \pm G ^ {\pm}, [ G ^ {+} _ {\lambda} G ^ {-} ] = L + \lambda J - \lambda^ {2} \frac {c}{6}, (7. 3)
$$

$$
[ L _ {\lambda} J ] = (\partial + \lambda) J - \lambda^ {2} \frac {c}{6}, [ L _ {\lambda} G ^ {+} ] = (\partial + \lambda) G ^ {+}, [ L _ {\lambda} G ^ {-} ] = (\partial + 2 \lambda) G ^ {-}.
$$

The set of positive roots of the affine superalgebra $\widehat { \mathfrak { g } } = s \ell ( 2 | 1 )$ is $( n \in \mathbf { Z } )$ :

$$
\widehat {\Delta} _ {+} = \{n K \mathrm {o f m u l t i p l i c i t y} 2 | n > 0 \}
$$

$$
\cup \left\{\alpha + n K \mid \alpha \in \Delta_ {+}, n \geq 0 \right\} \cup \left\{- \alpha + n K \mid \alpha \in \Delta_ {+}, n > 0 \right\},
$$

and the set of simple roots is

$$
\widehat {\Pi} = \left\{\alpha_ {0} = K - \alpha_ {1} - \alpha_ {2}, \alpha_ {1}, \alpha_ {2} \right\}.
$$

All admissible subsets $\widehat { \Delta } _ { + } ^ { \prime }$ of $\widehat { \Delta } _ { + }$ are principal, and the corresponding sets of simple roots are as follows:

$$
\widehat {\Pi} _ {b} = \left\{b _ {0} K + \alpha_ {0}, b _ {1} K + \alpha_ {1}, b _ {2} K + \alpha_ {2} \right\}, \text {w h e r e} b = (b _ {0}, b _ {1}, b _ {2}) \in \mathbf {Z} _ {+} ^ {3},
$$

$$
\widehat {\Pi} _ {b} ^ {-} = \left\{b _ {0} K - \alpha_ {0}, b _ {1} K - \alpha_ {1}, b _ {2} K - \alpha_ {2} \right\}, \mathrm {w h e r e} b = (b _ {0}, b _ {1}, b _ {2}) \in (1 + \mathbf {Z} _ {+}) ^ {3}.
$$

For the set $\widehat { \Pi _ { b } }$ , the boundary admissible weights $\Lambda$ are determined from the equation

$$
\left(\Lambda + \widehat {\rho} | b _ {0} K + \alpha_ {0}\right) = 1, \left(\Lambda + \widehat {\rho} | b _ {1} K + \alpha_ {1}\right) = \left(\Lambda + \widehat {\rho} | b _ {2} K + \alpha_ {2}\right) = 0. \tag {7.4}
$$

Adding these equations, we get $( \Lambda + \widehat { \rho } | u K ) = 1$ , where $u = b _ { 0 } + b _ { 1 } + b _ { 2 } + 1$ . Since $( \widehat { \rho } | K ) = 1$ , we obtain that the level of $\Lambda$ is given by

$$
k = \frac {1}{u} - 1, \text {w h e r e} u = b _ {0} + b _ {1} + b _ {2} + 1, \tag {7.5}
$$

and from (7.4) we obtain: $\left( \Lambda | \alpha _ { i } \right) \ = \ - \frac { b _ { i } } u$ u , $i = 0 , 1 , 2$ . Hence, denoting by $\Lambda _ { i }$ $\mathit { i } \ = \ 0 , 1 , 2$ ) the fundamental weights, i.e., $( \Lambda _ { i } | \alpha _ { j } ) = \delta _ { i j }$ , $( \Lambda _ { i } | D ) = 0$ , we obtain the unique boundary admissible weight corresponding to $\widehat { \Pi } _ { b }$ :

$$
\Lambda_ {b} = - \frac {1}{u} \left(b _ {0} \Lambda_ {0} + b _ {1} \Lambda_ {1} + b _ {2} \Lambda_ {2}\right), u = b _ {0} + b _ {1} + b _ {2} + 1.
$$

It is easy to see that this weight is nondegenerate iff $b _ { 0 } \geq 1$ , which we will assume.

Recall that ${ \mathfrak { h } } ^ { f } = \mathbf { C } ( h _ { 1 } - h _ { 2 } )$ . We let in (3.3) $h = z ( h _ { 1 } - h _ { 2 } )$ , $z \in \mathbf { C }$ , and let $y = e ^ { 2 \pi \imath z }$ . We shall calculate the normalized Euler–Poincar´e character

$$
\chi_ {H (M)} (\tau , z) := q ^ {- c / 2 4} \mathrm {c h} _ {H (M)} (z (h _ {1} - h _ {2})), \tag {7.6}
$$

where $c$ is the central charge (given by formula (7.9) below).

The conjectural character formula (3.5) gives in this case:

$$
\widehat {R} \operatorname {c h} _ {L \left(\Lambda_ {b}\right)} = e ^ {\Lambda_ {b}} \Pi_ {j = 1} ^ {\infty} \frac {\left(1 - q ^ {u (j - 1) + b _ {0}} e ^ {- \alpha_ {0}}\right) \left(1 - q ^ {u j - b _ {0}} e ^ {\alpha_ {0}}\right) \left(1 - q ^ {j}\right) ^ {2}}{\left(1 + q ^ {u (j - 1) + b _ {1}} e ^ {- \alpha_ {1}}\right) \left(1 + q ^ {u _ {j} - b _ {1}} e ^ {\alpha_ {1}}\right) \left(1 + q ^ {u (j - 1) + b _ {2}} e ^ {- \alpha_ {2}}\right) \left(1 + q ^ {u j - b _ {2}} e ^ {- \alpha_ {2}}\right)}. \tag {7.7}
$$

Due to (3.3), $\chi _ { H ( L ( \Lambda _ { b } ) ) }$ is obtained from this formula in the case of the Dynkin gradation by the specialization

$$
e ^ {- \alpha_ {0}} = 1, e ^ {- \alpha_ {1}} = y q ^ {\frac {1}{2}}, e ^ {- \alpha_ {2}} = y ^ {- 1} q ^ {\frac {1}{2}} \tag {7.8}
$$

(and multiplication by the specialized product). In order to write down the explicit formula, it is convenient to introduce the following important function:

$$
F (\tau , z _ {1}, z _ {2}) = \Pi_ {n = 1} ^ {\infty} \frac {(1 - q ^ {n}) ^ {2} (1 - e ^ {- 2 \pi i (z _ {1} + z _ {2})} q ^ {n}) (1 - e ^ {2 \pi i (z _ {1} + z _ {2})} q ^ {n - 1})}{(1 - e ^ {- 2 \pi i z _ {1}} q ^ {n}) (1 - e ^ {- 2 \pi i z _ {1}} q ^ {n - 1}) (1 - e ^ {- 2 \pi i z _ {2}} q ^ {n}) (1 - e ^ {2 \pi i z _ {2}} q ^ {n - 1})}
$$

and the following its specializations:

$$
F _ {j, \ell} ^ {(u)} (\tau , z) = q ^ {\frac {j \ell}{u}} e ^ {\frac {2 \pi i (j - \ell) z}{u}} F \left(u \tau , j \tau - z - \frac {1}{2}, \ell \tau + z + \frac {1}{2}\right).
$$

Note that plugging (7.5) in the formula for the central charge $c = - 3 ( 2 k + 1 )$ , we obtain:

$$
c = 3 - \frac {6}{u}, u = 2, 3, \dots . \tag {7.9}
$$

This is precisely the central charge of the minimal series representations of the $N = 2$ superconformal algebra. Recall that all these representations with given central charge (7.9) are parameterized by a pair of numbers $j , \ell \in \frac { 1 } { 2 } \mathbf { Z }$ satisfying inequalities $0 < j , \ell , j + \ell < u$ , the minimal eigenvalue of $L _ { 0 }$ being $\textstyle { \frac { j \ell - 1 / 4 } { u } }$ and the corresponding eigenvalue of $J _ { 0 }$ being $\textstyle { \frac { j - \ell } { u } }$ .

The specialization (7.8) of the right-hand side of (7.7) gives $F _ { b _ { 1 } + \frac { 1 } { 2 } , b _ { 2 } + \frac { 1 } { 2 } } \ ( \tau , z )$ , and the specialization in (7.8) of the product in (3.3) gives F (2)1 , 1 ( $F _ { \frac { 1 } { 2 } , \frac { 1 } { 2 } } ^ { ( 2 ) } ( \tau , z ) ^ { - 1 }$ . Hence, letting $j = b _ { 2 } + \frac { 1 } { 2 }$ and $\begin{array} { r } { \ell = b _ { 1 } + \frac { 1 } { 2 } } \end{array}$ formula (3.3) gives the well known (normalized) characters of the minimal series of $N = 2$ superconformal algebra (cf. [D, Ki, M]):

$$
\chi_ {H \left(L \left(\Lambda_ {b}\right)\right)} (\tau , z) = \chi_ {j, \ell} ^ {(u)} (\tau , z) := F _ {j, \ell} ^ {(u)} (\tau , z) / F _ {\frac {1}{2}, \frac {1}{2}} ^ {(2)} (\tau , z). \tag {7.10}
$$

Note that, given $u \geq 2$ , the range of $j$ and $\ell$ exactly corresponds to the range of $b _ { 1 }$ and $b _ { 2 }$ (defined by (7.5)), since $b _ { 0 } \geq 1$ . It is also easy to see that (2.19) for $\Lambda = \Lambda _ { b }$ gives the minimal eigenvalue of $L _ { 0 }$ , and the corresponding eigenvalue of $J _ { 0 }$ is indeed $\Lambda _ { b } ( h _ { 1 } - h _ { 2 } )$ . Using Remark 2.3, one can conclude that $H _ { 0 } ( L ( \Lambda _ { b } ) ) \ne 0$ (if $\Lambda _ { b }$ is non-degenerate). Hence, by Conjecture 3.3B, $H _ { j } ( L ( \Lambda _ { b } ) ) = 0$ for $j \neq 0$ , and therefore $H _ { 0 } ( L ( \Lambda _ { b } ) )$ is the irreducible module of minimal series corresponding to the parameters $u , j , \ell$ .

In a similar fashion, for $\Pi _ { b } ^ { - }$ the only boundary admissible weight is

$$
\Lambda_ {b} ^ {-} = \left(\frac {b _ {0}}{u} - 2\right) \Lambda_ {0} + \frac {b _ {1}}{u} \Lambda_ {1} + \frac {b _ {2}}{u} \Lambda_ {2}, u = b _ {0} + b _ {1} + b _ {2} - 1.
$$

All these weights are non-degenerate.

In a similar fashion, $\chi _ { H ( L ( \Lambda _ { b } ^ { - } ) ) }$ is obtained from (3.3) by using (7.7) and the specialization (7.8). It turns out that we again recover all characters of the $N = 2$ minimal series (7.10), where we set $\begin{array} { r } { j = b _ { 1 } - \frac { 1 } { 2 } } \end{array}$ , $\begin{array} { r } { \ell = b _ { 2 } - \frac { 1 } { 2 } } \end{array}$ . All other statements made about $\Lambda _ { b }$ hold for $\Lambda _ { b } ^ { - }$ as well.

We proceed in exactly the same way in the case of a non-Dynkin gradation. In this case the specialization (7.8) is replaced by

$$
e ^ {- \alpha_ {0}} = 1, e ^ {- \alpha_ {1}} = y, e ^ {- \alpha_ {2}} = q y ^ {- 1}.
$$

In a similar fashion we recover all Ramond type characters of the $N = 2$ superconformal algebra (meaning that we use the Virasoro field from the Ramond type basis (7.2), cf. [RY], [R]):

$$
e ^ {- \pi i c z} \mathrm {c h} _ {H (L (\Lambda_ {b}))} = \chi_ {j, \ell} ^ {(u) R} (\tau , z) := F _ {j, \ell} ^ {(u)} (\tau , z) / F _ {1, 0} ^ {(2)} (\tau , z), \tag {7.11}
$$

where $j = b _ { 2 } + 1$ and $\ell = b _ { 1 }$ so that the range of $j , \ell$ is exactly right:

$$
0 <   j, j + \ell <   u, 0 \leq \ell <   u.
$$

Likewise, the same result holds for $\Lambda _ { b } ^ { - }$ if we let $j = b _ { 1 }$ , $\ell = b _ { 2 } - 1$ . (Incidentally, using $L ^ { \prime }$ instead of $L$ , see (7.2), we get again $\chi _ { j , \ell } ^ { ( u ) }$ .)

Note that for the Ramond type basis (7.3) the fields $G ^ { + }$ and $G ^ { - }$ have conformal weights 1 and 2, respectively. Letting $\begin{array} { r } { G ^ { + } ( z ) = \sum _ { n \in \mathbf { Z } } G _ { n } ^ { + } z ^ { - n - 1 } } \end{array}$ , $\begin{array} { r } { G ^ { - } ( z ) = \sum _ { n \in \mathbf { Z } } G _ { n } ^ { - } z ^ { - n - 2 } } \end{array}$ , and introducing the constant term corrections: $\begin{array} { r } { \ddot { L } ( z ) = L ( z ) + \frac { c } { 2 4 z ^ { 2 } } } \end{array}$ , $\begin{array} { r } { \bar { J } ( z ) = J ( z ) - \frac { c } { 6 z } } \end{array}$ , formula (7.3) gives us exactly the commutation relation of the Ramond type $N = 2$ superalgebra. Using ${ \tilde { L } } _ { 0 }$ and ${ \tilde { J } } _ { 0 }$ in place of $L _ { 0 }$ and $J _ { 0 }$ in the definition of the normalized Euler–Poincar´e character, the definition (7.11) turns into the standard definition (7.6).

Recall [RY], [KW3], that, given $u$ , the span of all $N = 2$ characters, Ramond type characters and the corresponding supercharacters (obtained, up to a constant factor, by replacing $\tau$ by $\tau + 1$ in the character) form the minimal $S L _ { 2 } ( \mathbf { Z } )$ -invariant subspace containing the “vacuum” character $\chi _ { \frac { 1 } { 2 } , \frac { 1 } { 2 } } ^ { ( u ) }$ . Thus, taking quantum reduction for all good gradations of $s \ell ( 2 | 1 )$ of all boundary admissible highest weight $s \ell ( 2 | 1 )$ -modules, we get an $S L _ { 2 } ( \mathbf { Z } )$ -invariant space spanned by all characters and supercharacters.

# Acknowledgments.

We would like to thank ESI, Vienna, where we began this work in the summer of 2000, MSRI, Berkeley, where this work was continued in the spring of 2002, and M.I.T., where this paper was completed in the fall of 2002, for their hospitality. This paper was partially supported by NSF grants DMS9970007 and DMS0201017, NSC grant 902115M001020 of Taiwan, and grant in aid 13440012 for scientific research Japan.

# References

[BK] B. Bakalov and V. G. Kac, Field algebras, IMRN 2003 (3), 123-159. QA/ 0204282 .   
[B] M. Bershadsky, Conformal field theory via Hamiltonian reduction, Comm. Math. Phys. 139 (1991) 71-82 .

[BeK] M. Bershadsky, Phys. Lett. 174B (1986) 285; V.G. Knizhnik, Theor. Math. Phys. 66 (1986) 68.   
[BT] J.de Boer and T. Tjin, The relation between quantum W-algebras and Lie algebras, Comm. Math. Phys. 160 (1994) 317-332 .   
[BS] P. Bouwknegt and K. Schoutens,W-symmetry, Advanced ser Math. Phys, vol 22, World Sci., 1995.   
[D] V. K. Dobrev, Characters of the unitarizable highest weight modules over $N = 2$ superconformal algebras, Phys. Lett. B 186 (1987) 43-51.   
[FF1] B.L. Feigin and E. Frenkel, Quantization of Drinfeld-Sokolov reduction, Phys. Lett. B 246 (1990) 75-81.   
[FF2] B.L. Feigin and E. Frenkel, Representations of affine Kac-Moody algebras, bozonization and resolutions, Lett. Math. Phys. 19 (1990) 307-317.   
[F] P. Fiebig, The combinatorics of category $\mathcal { O }$ for symmetrizable Kac–Moody algebras, 2002 preprint.   
[FL] E.S. Fradkin and V. Ya. Linetsky, Classification of superconformal and quasisuperconformal algebras in two dimensions, Phys. Lett. B 291 (1992), 71-76.   
[FKW] E. Frenkel, V. Kac and M. Wakimoto, Characters and fusion rules for W-algebras via quantized Drinfeld-Sokolov reduction, Comm. Math. Phys. 147 (1992) 295-328.   
[GS] P. Goddard and A. Schwimmer, Factoring out free fermions and superconformal algebras, Phys. Lett. 214B (1988) 209-214.   
[K1] V. G. Kac, Lie superalgebras, Adv. Math. 26 (1977) 8-96.   
[K2] V. G. Kac, Infinite-dimensional algebras, Dedekind’s $\eta$ -function, classical M¨obius function and the very strange formula, Adv. in Math., 30 (1978) 85-136.   
[K3] V. G. Kac, Infinite-dimensional Lie algebras, 3rd edition, Cambridge University Press, 1990.   
[K4] V. G. Kac, Vertex algebras for beginners, Providence: AMS, University Lecture Series, Vol. 10, 1996. Second edition, 1998.   
[K5] V. G. Kac, Classification of supersymmetries, ICM talk, August 2002.   
[KW1] V. G. Kac and M. Wakimoto, Modular invariant representations of infinite-dimensional Lie algebras and superalgebras, Proc. Natl. Acad. Sci. USA 85 (1988) 4956-4960.   
[KW2] V. G. Kac and M. Wakimoto, Classification of modular invariant representations of affine algebras, in Infinite-dimensional Lie algebras and groups, Advanced Ser. Math. Phys. vol. 7, World Scientific, 1989, 138-177.   
[KW3] V. G. Kac and M. Wakimoto, Integrable highest weight modules over affine superalgebras and number theory, Progress in Math., 123, 1994, pp. 415-456. Birkh¨auser, Boston.   
[KW4] V. G. Kac and M. Wakimoto, Integrable highest weight modules over affine superalgebras and Appell’s function, Commun. Math, Phys. 215 (2001) 631-682.

[KW5] V.G. Kac and M. Wakimoto, Quantum reduction and representation theory of superconformal algebras.   
[Kh] T. Khovanova, Super KdV equation related to the Neveu–Schwarz-2 Lie superalgebra of string theory, Teor. Mat. Phys. 72 (1987) 306-312.   
[Ki] E. B. Kiritsis, Character formulae and the structure of the presentations of the $N = 1$ , $N = 2$ superconformal algebras, Int. J. Mod. Phys. A , 3 (1988) 1871-1906.   
[M] Y. Matsuo, Character formula of $C < 1$ unitary representation of $N = 2$ superconformal algebra, Prog. Theor. Phys. 77 (1987) 793-797.   
[RY] F. Ravanini and S-K. Yang, Modular invariance in $N = 2$ superconformal field theories, Phys. Lett. B195 (1987) 202-208.   
[R] S. S. Roan, Heisenberg and modular invariance of N=2 conformal field theory, Intern. J. Mod. Phys. A, 15 (2000) 3065-3094, hep-th/9902198.   
[W] M. Wakimoto, Lectures on infinite-dimensional Lie algebra, World Scientific, 2001.