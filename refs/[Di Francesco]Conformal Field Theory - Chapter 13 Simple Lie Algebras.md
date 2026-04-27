# Simple Lie Algebras

This chapter presents a survey of the theory of Lie algebras. This might appear somewhat remote from our main subject of interest: affine Lie algebras and their applications to conformal field theory. However, it turns out that in many respects the theory of affine Lie algebras is a natural extension of the theory of simple Lie algebras, and as such cannot be studied efficiently in isolation. This is an immediate motivation for devoted a complete chapter to Lie algebras. But as subsequent developments will show, conformal field theories with nonaffine additional symmetries, such as $W$ algebras, parafermions, and son on, as well as related exactly solvable statistical models, also have a deep Lie-algebraic underlying structure, which can only be appreciated with a minimal background on simple Lie algebras.

No previous knowledge on Lie algebras is assumed except for a first encounter with $su(2)$ and the theory of angular momentum in quantum mechanics. Admittedly, for those readers unfamiliar with the subject, this chapter will appear to be somewhat dense. Nevertheless the presentation is conceptually self-contained. This is not so at the technical level, since some statements and constructions are given without proofs. Furthermore, the choice of material is not completely standard, being dictated by our subsequent applications.

Section 13.1 covers the basic elements of the theory of simple Lie algebras: roots, weights, Cartan matrices, Dynkin diagrams, and the Weyl group. The subsequent section is devoted to the study of highest-weight representations. This is followed by an explicit description of states in $su(N)$ highest-weight representations, in terms of tableaux and patterns. Characters of irreducible representations are introduced in Sect. 13.4. From the results of Part B, it should already be clear that characters play a central role in some aspects of conformal field theory, such as modular invariance.

One of the central problems in conformal field theory is the calculation of fusion rules. Given that most conformal field theories have a Lie-algebraic core, the fusion rules are, to a large extent, determined by the tensor-product coefficients of this Lie algebra. It is thus mandatory to review in detail some methods for calculating tensor products. This is the subject of two sections, Sects 13.5 and 13.6. In the first

we present efficient techniques for tensor-product calculations, and in the second one we reconsider the problem from a conformal-field-theoretical angle.

Quotienting two affine Lie algebras will prove to be one the key tools in constructing conformal field theories. At the heart of this construction, there is a finite Lie algebra embedding, to which Sect. 13.7 is dedicated.

The basic properties of simple Lie algebras are displayed in App. 13.A, in a form that should facilitate later consultation. Finally, all symbols used in this chapter are collected in App. 13.B.

Readers familiar with Lie algebras may skip this chapter and use it only as a reference—except for a glance at App. 13.B in order to fix the notations. To those readers, we indicate that only Sects. 13.5.4 and 13.6 do not present standard material. On the other hand, those who wish to learn the basics of Lie algebra in this chapter should not necessarily read it linearly. Sections 13.1, 13.2, and 13.4.1 are essential and must be read sequentially. But the rest can be consulted when needed. Furthermore, it is not essential to master all the techniques for calculating tensor products in order to proceed. The description in the main text of tools particular to $su(N)$ (tableaux, the Littlewood-Richardson rule, Berenstein-Zelevinsky triangles) has the main purpose of lightening the presentation.

# §13.1. The Structure of Simple Lie Algebras

# 13.1.1. The Cartan-Weyl Basis

A Lie algebra $\mathbf{g}$ is a vector space equipped with an antisymmetric binary operation $[,]$ , called a commutator, mapping $\mathbf{g} \times \mathbf{g}$ into $\mathbf{g}$ , and further constrained to satisfy the Jacobi identity

$$
[ X, [ Y, Z ] ] + [ Z, [ X, Y ] ] + [ Y, [ Z, X ] ] = 0 \quad \text {f o r} \quad X, Y, Z \in \mathrm {g} \tag {13.1}
$$

Roughly speaking, the exponential of $\mathbf{g}$ is the Lie group $G$ (more precisely, its connected component containing the unit element): to $X \in \mathbf{g}$ , there corresponds the group elements $e^{iaX}$ where $a$ is some parameter and the exponential is defined from its power expansion. Hence, the algebra describes the group in the vicinity of the identity.

A representation refers to the association of every element of $\mathbf{g}$ to a linear operator acting on some vector space $V$ , which respects the commutation relations of the algebra. The maximal number of linearly independent states that generate $V$ is the dimension of the representation. Relative to a given basis, each element of $\mathbf{g}$ can thus be represented in terms of a square matrix and the basis vectors are represented by column matrices. (In the representation, the commutator corresponds to the usual matrix commutation.) A representation is said to be irreducible if the matrices representing the elements of $\mathbf{g}$ cannot all be brought in a block-diagonal form by a change of basis.

These elementary notions are sufficient to start analyzing the structure of Lie algebras. A Lie algebra can be specified by a set of generators $\{J^a\}$ and their

commutation relations

$$
\left[ J ^ {a}, J ^ {b} \right] = \sum_ {c} i f ^ {a b} c J ^ {c} \tag {13.2}
$$

The number of generators is the dimension of the algebra. The constants $f^{ab}_c$ are the structure constants, real parameters when $(J^a)^\dagger = J^a$ .<sup>1</sup> We are concerned with simple Lie algebras, that is, Lie algebras that contain no proper ideal (meaning no proper subset of generators $\{L^a\}$ such that $[L^a, J^b] \in \{L^a\}$ for any $J^b$ ). A direct sum of simple algebras is said to be semisimple.

In the standard Cartan-Weyl basis, the generators are constructed as follows. We first find the maximal set of commuting Hermitian generators $H^i$ , $i = 1, \dots, r$ ( $r$ is the rank of the algebra):

$$
\left[ H ^ {i}, H ^ {j} \right] = 0 \tag {13.3}
$$

This set of generators form the Cartan subalgebra $\mathbf{h}$ . The generators of the Cartan subalgebra can all be diagonalized simultaneously. The remaining generators are chosen to be those particular combinations of the $J^a$ 's that satisfy the following eigenvalue equation:

$$
\left[ H ^ {i}, E ^ {\alpha} \right] = \alpha^ {i} E ^ {\alpha} \tag {13.4}
$$

The vector $\alpha = (\alpha^1, \dots, \alpha^r)$ is called a root and $E^\alpha$ is the corresponding ladder operator. Because $\mathfrak{h}$ is the maximal Abelian subalgebra of $\mathfrak{g}$ , the roots are nondegenerate. The root $\alpha$ naturally maps an element $H^i \in \mathfrak{h}$ to the number $\alpha^i$ by $\alpha(H^i) = \alpha^i$ . Hence, the roots are elements of the dual of the Cartan subalgebra: $\alpha \in \mathfrak{h}^*$ .

Equation (13.4), through its Hermitian conjugate, shows that $-\alpha$ is necessarily a root whenever $\alpha$ is, with

$$
E ^ {- \alpha} = \left(E ^ {\alpha}\right) ^ {\dagger} \tag {13.5}
$$

In the following, $\Delta$ will denote the set of all roots.

Root components can be regarded as the nonzero eigenvalues of the $H^{i}$ in the particular representation, called the adjoint, for which the Lie algebra itself serves as the vector space on which the generators act. In this representation, we have an identification

$$
\begin{array}{c c c}E ^ {\alpha}&\longmapsto&| E ^ {\alpha} \rangle \equiv | \alpha \rangle\\\rightarrow&&\end{array}\tag {13.6}
$$

$$
\begin{array}{c c c} H ^ {i} & \longmapsto & | H ^ {i} \rangle \end{array}
$$

between the generators and the states of the representation. It follows from Eq. (13.4) that in the adjoint representation the action of a generator $X$ is represented by $\operatorname{ad}(X)$ , defined as

$$
\operatorname {a d} (X) Y = [ X, Y ] \tag {13.7}
$$

so that

$$
\operatorname {a d} \left(H ^ {i}\right) E ^ {\alpha} = \alpha^ {i} E ^ {\alpha} \quad \mapsto \quad H ^ {i} | \alpha \rangle = \alpha^ {i} | \alpha \rangle \tag {13.8}
$$

The one-to-one correspondence between the states $|\alpha \rangle$ and the ladder operators $E^{\alpha}$ reflects the nondegenerate character of roots. In this representation, the zero eigenvalue has degeneracy $r$ (associated with the different states $|H^{i}\rangle$ ). By construction, the dimension of the adjoint is equal to the dimension of the algebra, itself equal to the total number of roots plus $r$ .

In view of specifying the remaining commutators, we first observe that the Jacobi identity implies

$$
\left[ H ^ {i}, \left[ E ^ {\alpha}, E ^ {\beta} \right] \right] = \left(\alpha^ {i} + \beta^ {i}\right) \left[ E ^ {\alpha}, E ^ {\beta} \right] \tag {13.9}
$$

If $\alpha + \beta \in \Delta$ , the commutator $[E^{\alpha}, E^{\beta}]$ must be proportional to $E^{\alpha + \beta}$ , and it must vanish if $\alpha + \beta \notin \Delta$ . When $\alpha = -\beta$ , $[E^{\alpha}, E^{-\alpha}]$ commutes with all $H^{i}$ , which is possible only if it is a linear combination of the generators of the Cartan subalgebra. The normalization of the ladder operators is fixed by setting this commutator equal to $2\alpha \cdot H / |\alpha|^{2}$ , where

$$
\alpha \cdot H = \sum_ {i = 1} ^ {r} \alpha^ {i} H ^ {i} \quad | \alpha | ^ {2} = \sum_ {i = 1} ^ {r} \alpha^ {i} \alpha^ {i} \tag {13.10}
$$

Summarizing, the full set of commutation relations in the Cartan-Weyl basis is

$$
\begin{array}{l l} \left[ H ^ {i}, H ^ {j} \right] = 0 \\ \left[ H ^ {i}, E ^ {\alpha} \right] = \alpha^ {i} E ^ {\alpha} \\ \left[ E ^ {\alpha}, E ^ {\beta} \right] = N _ {\alpha , \beta} E ^ {\alpha + \beta} & \text {i f} \quad \alpha + \beta \in \Delta \\ = \frac {2}{| \alpha | ^ {2}} \alpha \cdot H & \text {i f}, \alpha = - \beta \\ = 0 & \text {o t h e r w i s e} \end{array} \tag {13.11}
$$

where $N_{\alpha,\beta}$ is a constant.

# 13.1.2. The Killing Form

The normalization used to fix the commutators is usually introduced by means of the Killing form

$$
\tilde {K} (X, Y) = \operatorname {T r} (\operatorname {a d} X \operatorname {a d} Y) \tag {13.12}
$$

which gives a sort of scalar product for the Lie algebra. To calculate this trace in some basis of generators $\{T^a\}$ , we first evaluate $[X, [Y, T^b]]$ in terms of the elements of this basis; the coefficient of $T^b$ in the result gives the contribution of this term to the trace. For semisimple Lie algebras, the Killing form is nondegenerate: $\tilde{K}(X, Y) = 0$ for all $Y$ implies that $X = 0$ . This is in fact an alternate way of defining semisimplicity.

In the following, we will mainly use a renormalized version of the Killing form defined as

$$
K (X, Y) \equiv \frac {1}{2 g} \operatorname {T r} (\operatorname {a d} X \operatorname {a d} Y) \tag {13.13}
$$

where $g$ is a constant that will be defined later (it is the dual Coxeter number of the algebra $g$ ). The standard basis $\{J^a\}$ is understood to be orthonormal with respect to $K$ :

$$
K \left(J ^ {a}, J ^ {b}\right) = \delta^ {a, b} \tag {13:14}
$$

The same normalization holds for the generators of the Cartan subalgebra

$$
K \left(H ^ {i}, H ^ {j}\right) = \delta^ {i, j} \tag {13.15}
$$

Since the Killing form defines a scalar product, it can be used to lower or raise the indices, e.g.,

$$
\hat {f} _ {a b c} = \sum_ {d} f ^ {a d} _ {c} [ K (J ^ {d}, J ^ {b}) ] ^ {- 1} \tag {13.16}
$$

We note that $f_{abc}$ is antisymmetric in all three indices. In the $\{J^a\}$ orthonormal basis, the position of the indices (up or down) is thus irrelevant.

The cyclic property of the trace yields the identity

$$
K ([ Z, X ], Y) + K (X, [ Z, Y ]) = 0 \tag {13.17}
$$

The Killing form is actually uniquely characterized by this property. With appropriate choices for $X, Y, Z \in \mathbf{g}$ , it follows that

$$
\left[ E ^ {\alpha}, E ^ {- \alpha} \right] = K \left(E ^ {\alpha}, E ^ {- \alpha}\right) \alpha \cdot H \tag {13.18}
$$

(all other pairs involving a ladder operator have zero Killing form). Hence, the previously introduced normalization corresponds to

$$
K \left(E ^ {\alpha}, E ^ {- \alpha}\right) = \frac {2}{| \alpha | ^ {2}} \tag {13.19}
$$

However, the fundamental role of the Killing form is to establish an isomorphism between the Cartan subalgebra $\mathfrak{h}$ and its dual $\mathfrak{h}^*$ : the form $K(H^i,\cdot)$ ( $i$ fixed) maps every element of the Cartan subalgebra onto a number. Hence, to every element $\gamma \in \mathfrak{h}^*$ , there corresponds a $H^\gamma \in \mathfrak{h}$ through

$$
\gamma \left(H ^ {i}\right) = K \left(H ^ {i}, H ^ {\gamma}\right) \tag {13.20}
$$

(in particular for a root $\alpha$ , $H^{\alpha} = \alpha \cdot H = \sum_{i} \alpha^{i} H^{i}$ ). With this isomorphism, the Killing form can be transferred into a positive definite scalar product in the dual space

$$
(\gamma , \beta) = K \left(H ^ {\beta}, H ^ {\gamma}\right) \tag {13.21}
$$

Since roots are elements of $\mathbf{h}^*$ , this defines a scalar product in the root space. From now on, the scalar product between roots will be denoted as above, with the understanding that $|\alpha|^2 = (\alpha, \alpha)$ .

# 13.1.3. Weights

Up to this point, we have analyzed the structure of the algebra from the point of view of a particular representation (the adjoint), that for which the algebra itself plays the role of the vector space. In this representation, the eigenvalues of the Cartan generators are called the roots and the scalar product between roots is induced by the Killing form. Since the essential structure of the algebra is coded in this representation, it needs to be studied in more detail. For this, it is useful to first recast the problem in the general context of a finite-dimensional representation.

For an arbitrary representation, a basis $\{|\lambda\rangle\}$ can always be found such that

$$
H ^ {i} | \lambda \rangle = \lambda^ {i} | \lambda \rangle \tag {13.22}
$$

The eigenvalues $\lambda^i$ build the vector $\lambda = (\lambda^1, \dots, \lambda^r)$ , called a weight. Weights live in the space $\mathbf{h}^* \colon \lambda(H^i) = \lambda^i$ . Hence, the scalar product between weights is also fixed by the Killing form. In the adjoint representation, the weights deserve the special name of roots. The commutator (13.4) shows that $E^\alpha$ changes the eigenvalue of a state by $\alpha$ :

$$
H ^ {i} E ^ {\alpha} | \lambda \rangle = \left[ H ^ {i}, E ^ {\alpha} \right] | \lambda \rangle + E ^ {\alpha} H ^ {i} | \lambda \rangle = \left(\lambda^ {i} + \alpha^ {i}\right) E ^ {\alpha} | \lambda \rangle \tag {13.23}
$$

so that $E^{\alpha}|\lambda \rangle$ , if nonzero, must be proportional to a state $|\lambda + \alpha \rangle$ . This justifies the name ladder (or step) operator for $E^{\alpha}$ .

Representations of interest are the finite-dimensional ones. For these, we will derive an important relation, to be used shortly for the adjoint representation. For any state $|\lambda \rangle$ in a finite-dimensional representation, there are necessarily two positive integers $p$ and $q$ , such that

$$
\begin{array}{l} (E ^ {\alpha}) ^ {p + 1} | \lambda \rangle \sim E ^ {\alpha} | \lambda + p \alpha \rangle = 0 \tag {13.24} \\ (E ^ {- \alpha}) ^ {q + 1} | \lambda) \sim E ^ {\alpha} | \lambda - q \alpha) = 0 \\ \end{array}
$$

for any root $\alpha$ . Indeed, notice that the triplet of generators $E^{\alpha}, E^{-\alpha}$ , and $\alpha \cdot H / |\alpha|^2$ forms an $su(2)$ subalgebra analogue to the set $\{J^{+}, J^{-}, J^{3}\}$ , with commutation relations

$$
\left[ J ^ {+}, J ^ {-} \right] = 2 J ^ {3}, \quad \left[ J ^ {3}, J ^ {\pm} \right] = \pm J ^ {\pm} \tag {13.25}
$$

Therefore, if $|\lambda \rangle$ belongs to a finite-dimensional representation, its projection onto the $su(2)$ subalgebra associated with the root $\alpha$ must also be finite dimensional. Let the dimension of the latter be $2j + 1$ ; then from the state $|\lambda \rangle$ , the state with highest $J^3 = \alpha \cdot H / |\alpha|^2$ projection $(m = j)$ can be reached by a finite number, say $p$ , applications of $J^{+} = E^{\alpha}$ , whereas, say, $q$ applications of $J^{-} = E^{-\alpha}$ lead to the

state with $m = -j$ :

$$
j = \frac {(\alpha , \lambda)}{| \alpha | ^ {2}} + p, \quad - j = \frac {(\alpha , \lambda)}{| \alpha | ^ {2}} - q \tag {13.26}
$$

Eliminating $j$ from the above two equations yields

$$
2 \frac {(\alpha , \lambda)}{| \alpha | ^ {2}} = - (p - q) \tag {13.27}
$$

This is the relation we were looking for: any weight $\lambda$ in a finite-dimensional representation is such that $(\alpha, \lambda) / |\alpha|^2$ is an integer. This is true in particular for $\lambda = \beta$ , where $\beta$ is any root of the algebra. We now return to the analysis of the root properties.

# 13.1.4. Simple Roots and the Cartan Matrix

As already mentioned, the number of roots is equal to the dimension of the algebra minus its rank, and this number is in general much larger than the rank itself. This means that the roots are linearly dependent. We then fix a basis $\{\beta_1,\beta_2,\dots ,\beta_r\}$ in the space $\mathbf{h}^*$ , so that any root can be expanded as

$$
\alpha = \sum_ {i = 1} ^ {r} n _ {i} \beta_ {i} \tag {13.28}
$$

In this basis, an ordering can be defined as follows: $\alpha$ is said to be positive if the first nonzero number in the sequence $(n_{1}, n_{2}, \dots, n_{r})$ is positive. Denote by $\Delta_{+}$ the set of positive roots. The set of negative roots $\Delta_{-}$ is defined in the obvious way. We have already observed that whenever $\alpha$ is a root, $-\alpha$ is also a root; hence $\Delta_{-} = -\Delta_{+}$ .

A simple root $\alpha_{i}$ is defined to be a root that cannot be written as the sum of two positive roots. There are necessarily $r$ simple roots, and their set $\{\alpha_{1},\dots ,\alpha_{r}\}$ provides the most convenient basis for the $r$ -dimensional space of roots. Notice that the subindex is a labeling index: it does not refer to a root component. Two immediate consequences of the definition of simple roots are: (i) $\alpha_{i} - \alpha_{j}\notin \Delta$ (otherwise, if $\alpha_{i} - \alpha_{j} > 0$ , say, we would conclude that $\alpha_{i} = \alpha_{j} + (\alpha_{i} - \alpha_{j})$ , a contradiction); (ii) any positive root is a sum of positive roots (indeed, if a positive root is not simple, it can be written as a sum of two positive roots, which, if not simple, can also be written as the sum of two positive roots, and so on).

The scalar products of simple roots define the Cartan matrix

$$
A _ {i j} = \frac {2 \left(\alpha_ {i} , \alpha_ {j}\right)}{\alpha_ {j} ^ {2}} \tag {13.29}
$$

In view of Eq. (13.27), the entries of this matrix are necessarily integers. Its diagonal elements are all equal to 2 and it is not symmetric in general. The Schwarz inequality implies that $A_{ij}A_{ji} < 4$ for $i \neq j$ . Since $\alpha_i - \alpha_j$ is not a root, $E^{-\alpha_j}|\alpha_i\rangle = 0$ and $q = 0$ in Eq. (13.24) for $\lambda = \alpha_i$ and $\alpha = \alpha_j$ . Hence, from Eq. (13.27) it follows that

$$
\left(\alpha_ {i}, \alpha_ {j}\right) \leq 0, \quad i \neq j \tag {13.30}
$$

Thus for $i \neq j$ , $A_{ij}$ is a nonpositive integer, and in view of the above inequality, it can only be $0, -1, -2$ , or $-3$ . If $A_{ij} \neq 0$ , the inequality forces at least one of $A_{ij}$ or $A_{ji}$ to be $-1$ .

It can be shown that in the set of roots of a simple Lie algebra, at most two different lengths (long and short) are possible. The ratio of the length of the long roots over the short roots is bound to be 2 or 3, if different from 1. When all the roots have the same length, the algebra is said to be simply laced.

It is convenient for us to introduce a special notation for the quantity $2\alpha_{i} / |\alpha_{i}|^{2}$ :

$$
\alpha_ {i} ^ {\vee} = \frac {2 \alpha_ {i}}{\left| \alpha_ {i} \right| ^ {2}} \tag {13.31}
$$

$\alpha_{i}^{\vee}$ is called the coroot associated with the root $\alpha_{i}$ . The scalar product between roots and coroots is thus always an integer. The Cartan matrix now takes the compact form

$$
A _ {i j} = \left(\alpha_ {i}, \alpha_ {j} ^ {\vee}\right) \tag {13.32}
$$

A distinguished element of $\Delta$ is the highest root $\theta$ . It is the unique root for which, in the expansion $\sum m_i \alpha_i$ , the sum $\sum m_i$ is maximized. All elements of $\Delta$ can be obtained by repeated subtraction of simple roots from $\theta$ . The coefficients of the decomposition of $\theta$ in the bases $\{\alpha_i\}$ and $\{\alpha_i^\vee\}$ bear special names, being called, respectively, the marks $(a_i)$ and the comarks $(a_i^\vee)$ :

$$
\theta = \sum_ {i = 1} ^ {r} a _ {i} \alpha_ {i} = \sum_ {i = 1} ^ {r} a _ {i} ^ {\vee} \alpha_ {i} ^ {\vee}, \quad a _ {i}, a _ {i} ^ {\vee} \in \mathbb {N} \tag {13.33}
$$

Marks and comarks are related by

$$
a _ {i} = a _ {i} ^ {\vee} \frac {2}{| \alpha_ {i} | ^ {2}} \tag {13.34}
$$

The dual Coxeter number is defined as

$$
\boxed {g = \sum_ {i = 1} ^ {r} a _ {i} ^ {\vee} + 1} \tag {13.35}
$$

(The Coxeter number can be defined similarly, but it will not be used here. The superscript $\vee$ , which would naturally appear in the notation for the dual Coxeter number, is thus omitted.)

# 13.1.5. The Chevalley Basis

As will be shown below, the full set of roots can be reconstructed from the set of simple roots, and the latter can be extracted from the Cartan matrix in a very simple way. Moreover, the Cartan matrix fixes completely the commutation relations of the algebra. This point is made fully manifest in the Chevalley basis where to each simple root $\alpha_{i}$ there corresponds the three generators

$$
e ^ {i} = E ^ {\alpha_ {i}} \quad f ^ {i} = E ^ {- \alpha_ {i}} \quad h ^ {i} = \frac {2 \alpha_ {i} \cdot H}{\left| \alpha_ {i} \right| ^ {2}} \tag {13.36}
$$

whose commutation relations are

$$
\begin{array}{l} [ h ^ {i}, h ^ {j} ] = 0 \\ [ h ^ {i}, e ^ {j} ] = A _ {j i} e ^ {j} \\ \left[ h ^ {i}, f ^ {j} \right] = - A _ {j i} f ^ {j} \tag {13.37} \\ [ e ^ {i}, f ^ {j} ] = \delta_ {i j} h ^ {j} \\ \end{array}
$$

The remaining step operators are obtained by repeated commutations of these basic generators, subject to the Serre relations

$$
\begin{array}{l} [ \operatorname {a d} (e ^ {i}) ] ^ {1 - A _ {j i}} e ^ {j} = 0 \\ [ \operatorname {a d} (f ^ {i}) ] ^ {1 - A _ {j i}} f ^ {j} = 0 \tag {13.38} \\ \end{array}
$$

For instance, $[\mathrm{ad}(e^i)]^2 e^j = [e^i, [e^i, e^j]]$ . These constraints—the analogues of relations (13.24) for the adjoint representation—encode the rules for reconstructing the full root system from the simple roots. (For this specific problem, still another approach will be presented later.) The Serre relations do not mix the $e^i$ 's and the $f^i$ 's and this reflects the separation of the roots into two disjoint sets $\Delta_{\pm}$ . That the Serre relations and the basic commutation relations can be expressed in terms of the Cartan matrix shows that $A$ contains all the information on the structure of $g$ . Actually, the abstract formulation of Lie algebras in terms of Cartan matrices is the most efficient starting point for generalizations.

The Killing form of the generators of the Cartan subalgebra is easily transcribed from the Cartan-Weyl to the Chevalley basis:

$$
K \left(h ^ {i}, h ^ {j}\right) = \left(\alpha_ {i} ^ {\vee}, \alpha_ {j} ^ {\vee}\right) \tag {13.39}
$$

# 13.1.6. Dynkin Diagrams

All the information contained in the Cartan matrix can be encapsulated in a simple planar diagram: the Dynkin diagram. To every simple root $\alpha_{i}$ , we associate a

node (white for a long root and black for a short one) and join the nodes $i$ and $j$ with $A_{ij}A_{ji}$ lines. Hence orthogonal simple roots are disconnected, and those sustaining an angle of 120, 135, or 150 degrees are linked by one, two, or three lines, respectively.

The classification of simple Lie algebras boils down to a classification of Dynkin diagrams. The complete list contains four infinite families, the algebras $A_r, B_r, C_r$ and $D_r$ (the classical algebras, whose compact real forms are respectively $su(r + 1), so(2r + 1), sp(2r)$ , and $so(2r)$ ), and five exceptional cases: $E_6, E_7, E_8, F_4$ , and $G_2$ . The subscript gives the rank of the algebra. The Dynkin diagrams as well as basic properties of these Lie algebras are displayed in App. 13.A. Note that the $A, D, E$ algebras are simply laced. (The classification of simply-laced algebras has already been considered in Ex. 10.10.)

# 13.1.7. Fundamental Weights

As already pointed out, weights and roots live in the same $r$ -dimensional vector space. The weights can thus be expanded in the basis of simple roots. However, this expansion is not very useful since for irreducible finite-dimensional representations—the representations of interest—its coefficients are not integers. The convenient basis for weights is in fact the one dual to the simple coroot basis. It is denoted by $\{\omega_i\}$ and defined by

$$
\left(\omega_ {i}, \alpha_ {j} ^ {\vee}\right) = \delta_ {i j} \tag {13.40}
$$

The $\omega_{i}$ are called the fundamental weights.

The expansion coefficients $\lambda_{i}$ of a weight $\lambda$ in the fundamental weight basis are called Dynkin labels. Hence,

$$
\lambda = \sum_ {i = 1} ^ {r} \lambda_ {i} \omega_ {i} \quad \Longleftrightarrow \quad \lambda_ {i} = \left(\lambda , \alpha_ {i} ^ {\vee}\right) \tag {13.41}
$$

The Dynkin labels of weights in finite-dimensional irreducible representations are always integers (this follows from Eq. (13.27) and it will be made explicit in the next section); such weights are said to be integral. From now on, whenever a weight is written in component form

$$
\lambda = \left(\lambda_ {1}, \dots , \lambda_ {r}\right) \tag {13.42}
$$

(with entries separated by commas) it is understood that these components are the Dynkin labels. Note that the elements of the Cartan matrix are the Dynkin labels

of the simple roots

$$
\alpha_ {i} = \sum_ {j} A _ {i j} \omega_ {j} \tag {13.43}
$$

that is, the $i$ -th row of $\pmb{A}$ is the set of Dynkin labels for the simple root $\alpha_{i}$ .

The Dynkin labels are the eigenvalues of the Chevalley generators of the Cartan subalgebra:

$$
h ^ {i} | \lambda \rangle = \lambda \left(h ^ {i}\right) | \lambda \rangle = \left(\lambda , \alpha_ {i} ^ {\vee}\right) | \lambda \rangle \tag {13.44}
$$

that is

$$
\boxed {h ^ {i} | \lambda \rangle = \lambda_ {i} | \lambda \rangle} \tag {13.45}
$$

The position of the index has the following meaning: $\lambda_{i}$ refers to an eigenvalue of $h^{i}$ (a Dynkin label), whereas $\lambda^i$ is an eigenvalue of $H^{i}$ .

A weight of special importance, thus deserving a special notation, is the one for which all Dynkin labels are unity:

$$
\rho = \sum_ {i} \omega_ {i} = (1, 1, \dots , 1) \tag {13.46}
$$

This is called the Weyl vector (or principal vector) and has the following alternate definition (to be proved later):

$$
\rho = \frac {1}{2} \sum_ {\alpha \in \Delta_ {+}} \alpha . \tag {13.47}
$$

The scalar product of weights can be expressed in terms of a symmetric quadratic form matrix $F_{ij}$

$$
\left(\omega_ {i}, \omega_ {j}\right) = F _ {i j} \tag {13.48}
$$

The definition implies that $F_{ij}$ is the transformation matrix relating the two bases $\{\omega_i\}$ and $\{\alpha_i^\vee\}$

$$
\omega_ {i} = \sum_ {j} F _ {i j} \alpha_ {j} ^ {\vee} \tag {13.49}
$$

Indeed, the product of this equation with $\alpha_{j}^{\vee}$ reproduces (13.48). Hence $F_{ij}$ is the inverse of the matrix whose rows are the Dynkin labels of the simple coroots, and these can be read off the following rescaled version of (13.43).

$$
\alpha_ {i} ^ {\vee} = \sum_ {j} \frac {2}{| \alpha_ {i} | ^ {2}} A _ {i j} \omega_ {j} \tag {13.50}
$$

This leads to an explicit relation between the quadratic form and the Cartan matrix:

$$
F _ {i j} = \left(A ^ {- 1}\right) _ {i j} \frac {\alpha_ {j} ^ {2}}{2} \tag {13.51}
$$

The scalar product of the two weights $\lambda = \sum \lambda_{i}\omega_{i}$ and $\mu = \sum \mu_{i}\omega_{i}$ reads

$$
(\lambda , \mu) = \sum_ {i, j} \lambda_ {i} \mu_ {j} \left(\omega_ {i}, \omega_ {j}\right) = \sum_ {i, j} \lambda_ {i} \mu_ {j} F _ {i j} \tag {13.52}
$$

The quadratic form matrices of all the simple Lie algebras are tabulated in App. 13.A, with the normalization convention defined in Sect. 13.1.10.

# 13.1.8. The Weyl Group

We return for a moment to the projection of the adjoint representation onto the $su(2)$ subalgebra associated with the root $\alpha$ . Let $m$ be the eigenvalue of the $J^3$ operator $\alpha \cdot H / |\alpha|^2$ on the state $|\beta\rangle$ ; that is,

$$
2 m = \left(\alpha^ {\vee}, \beta\right) \tag {13.53}
$$

If $m \neq 0$ , this state must be paired with another one with $J^3$ eigenvalue $-m$ . Therefore, there must exist another state in the multiplet, say $|\beta + \ell \alpha\rangle$ ; whose projection on the $J^3$ axis is equal to

$$
\left(\alpha^ {\vee}, \beta + \ell \alpha\right) = \left(\alpha^ {\vee}, \beta\right) + 2 \ell = - \left(\alpha^ {\vee}, \beta\right) \tag {13.54}
$$

This shows that if $\beta$ is a root, $\beta - (\alpha^{\vee}, \beta)\alpha$ is also a root.

The operation $s_{\alpha}$ defined by

$$
s _ {\alpha} \beta = \beta - \left(\alpha^ {\vee}, \beta\right) \alpha \tag {13.55}
$$

is a reflection with respect to the hyperplane perpendicular to $\alpha$ . The set of all such reflections with respect to roots forms a group, called the Weyl group of the algebra, denoted $W$ . It is generated by the $r$ elements $s_i$ , the simple Weyl reflections,

$$
s _ {i} \equiv s _ {\alpha_ {i}} \tag {13.56}
$$

in the sense that every element $w \in W$ can be decomposed as

$$
w = s _ {i} s _ {j} \dots s _ {k} \tag {13.57}
$$

For the simple Weyl reflections, the following relations are easily checked

$$
s _ {i} ^ {2} = 1, \quad s _ {i} s _ {j} = s _ {j} s _ {i} \quad \text {i f} \quad A _ {i j} = 0 \tag {13.58}
$$

These generalize to

$$
(s _ {i} s _ {j}) ^ {m _ {i j}} = 1 \quad \text {w h e r e} \quad m _ {i j} = \left\{ \begin{array}{c c} 2 & \text {i f} i = j \\ \frac {\pi}{\pi - \theta_ {i j}} & \text {i f} i \neq j \end{array} \right. \tag {13.59}
$$

with $\theta_{ij}$ the angle between the simple root $\alpha_{i}$ and $\alpha_{j}$ .<sup>7</sup> Eq. (13.59) can be regarded as the defining relation of the Weyl group. We note again that it is expressed in

terms of data directly related to the Cartan matrix. On the simple roots, the action of $s_i$ takes the simple form

$$
s _ {i} \alpha_ {j} = \alpha_ {j} - A _ {j i} \alpha_ {i} \tag {13.60}
$$

It has just been shown that $W$ maps $\Delta$ into itself. In fact, it provides a simple way to generate the complete set $\Delta$ from the simple roots by acting with all the elements of $W$ on the set $\{\alpha_i\}$ :

$$
\Delta = \{w \alpha_ {1}, \dots , w \alpha_ {r} | w \in W \} \tag {13.61}
$$

From this construction, it is clear that any set $\{w' \alpha_i\}$ with $w'$ fixed, could serve as a basis of simple roots. (This gives the announced relation between the different bases of simple roots.)

As a short digression; we now prove, using the Weyl group, the equivalence between (13.46) and (13.47). From (13.46) it follows that $(\rho, \alpha_i^\vee) = 1$ for all $i$ . We want to show that the same result follows from the second definition. We set $\sigma = \sum_{\alpha > 0} \alpha / 2$ and consider $s_i \sigma$ . Since $s_i$ permutes all the positive roots—that is, $A_{ij} \leq 0$ if $i \neq j$ (except $\alpha_i$ which is mapped to $-\alpha_i$ ), we can write

$$
s _ {i} \sigma = \frac {1}{2} \sum_ {\substack {\alpha > 0 \\ \alpha \neq \alpha_ {i}}} \alpha - \frac {1}{2} \alpha_ {i} = \frac {1}{2} \sum_ {\alpha > 0} \alpha - \alpha_ {i} \tag{13.62}
$$

implying that

$$
\left(s _ {i} \sigma , \alpha_ {i} ^ {\vee}\right) = \left(\sigma - \alpha_ {i}, \alpha_ {i} ^ {\vee}\right) = \left(\sigma , \alpha_ {i} ^ {\vee}\right) - 2 \tag {13.63}
$$

On the other hand, from the invariance of the scalar product with respect to Weyl transformations, the same product can be written as

$$
\left(s _ {i} \sigma , \alpha_ {i} ^ {\vee}\right) = \left(\sigma , s _ {i} \alpha_ {i} ^ {\vee}\right) = - \left(\sigma , \alpha_ {i} ^ {\vee}\right) \tag {13.64}
$$

The compatibility of these two equations gives the desired result, namely $(\sigma, \alpha_{i}^{\vee}) = 1$ and thus $\sigma = \rho$ .

The action of the Weyl group, defined so far only for roots, extends naturally to weights:

$$
s _ {\alpha} \lambda = \lambda - (\alpha^ {\vee}, \lambda) \alpha \mid \tag {13.65}
$$

It is straightforward to verify from the above relation that the Weyl group leaves the scalar product invariant

$$
\left(s _ {\alpha} \lambda , s _ {\alpha} \mu\right) = (\lambda , \mu) \tag {13.66}
$$

or more generally

$$
(w \lambda , \mu) = (\lambda , w ^ {- 1} \mu) \tag {13.67}
$$

The Weyl group induces a natural splitting of the $r$ -dimensional weight vector space into chambers, whose number is equal to the order of $W$ . These are simplicial cones defined as

$$
C _ {w} = \{\lambda | (w \lambda , \alpha_ {i}) \geq 0, i = 1, \dots , r \}, \quad w \in W \tag {13.68}
$$

These chambers intersect only at their boundaries $(w\lambda, \alpha_i) = 0$ , the reflecting hyperplanes of the $s_i$ 's. The chamber corresponding to the identity element of the Weyl group is called the fundamental chamber, and it will be denoted by $C_0$ . An obvious but fundamental consequence of this splitting is that for any weight $\lambda \notin C_0$ , there exists a $w \in W$ such that $w\lambda \in C_0$ . More precisely, the $W$ orbit of every weight has exactly one point in the fundamental chamber. The $W$ orbit of $\lambda$ is the set of all weights $\{w\lambda | w \in W\}$ . A weight in the fundamental chamber and whose Dynkin labels are all integers, $\lambda_i \in \mathbb{Z}_{+}$ , is said to be dominant. (A dominant weight is thus understood to be integral.) $\theta$ is an example of a dominant weight.

To conclude this section, we present some notation that will be used extensively in the sequel. The modified Weyl reflection

$$
w \cdot \lambda \equiv w (\lambda + \rho) - \rho \tag {13.69}
$$

denoted by a dot, will be referred to as a shifted Weyl reflection. Here $\rho$ is the Weyl vector. It is simple to verify that

$$
w \cdot \left(w ^ {\prime} \cdot \lambda\right) = \left(w w ^ {\prime}\right) \cdot \lambda \tag {13.70}
$$

The length of $w$ , denoted $\ell(w)$ , is the minimum number of $s_i$ among all possible decompositions of $w = \prod_{i} s_i$ . The signature of $w$ is defined as

$$
\epsilon (w) = (- 1) ^ {\ell (w)} \tag {13.71}
$$

In the linear representation of $w$ , this is simply $\det(w)$ (cf. Ex. 13.3). Finally, the longest element of the Weyl group will be denoted by $w_0$ . It is the unique element of $W$ that maps $\Delta_{+}$ to $\Delta_{-}$ .

# 13.1.9. Lattices

In terms of a basis $(\epsilon_1, \dots, \epsilon_d)$ of the $d$ -dimensional Euclidean space $\mathbb{R}^d$ , a lattice is the set of all points whose expansion coefficients, in terms of the specified basis, are all integers:

$$
\mathbb {Z} \epsilon_ {1} + \mathbb {Z} \epsilon_ {2} + \dots + \mathbb {Z} \epsilon_ {d} \tag {13.72}
$$

In other words, it is the $\mathbb{Z}$ span of $\{\epsilon_i\}$ . Three $r$ -dimensional lattices are important for Lie algebras. These are the weight lattice

$$
P = \mathbb {Z} \omega_ {1} + \dots + \mathbb {Z} \omega_ {r} \tag {13.73}
$$

the root lattice

$$
Q = \mathbb {Z} \alpha_ {1} + \dots + \mathbb {Z} \alpha_ {r} \tag {13.74}
$$

and the coroot lattice

$$
Q ^ {\vee} = \mathbb {Z} \alpha_ {1} ^ {\vee} + \dots + \mathbb {Z} \alpha_ {r} ^ {\vee} \tag {13.75}
$$

The relevance of the weight lattice lies in that the weights in finite-dimensional representations have integer Dynkin labels (cf. Eq. (13.27)), hence they belong to $P$ . The connection between $P$ and the generators of $\mathbf{g}$ is twofold. First, the integers specifying the position of a weight in $P$ are the eigenvalues of the Chevalley generators $h^i$ . Second, the effect of the other generators is to shift the eigenvalues by an element of the root lattice $Q$ . Since roots are weights in a particular finite-dimensional representation, $Q \subseteq P$ . Hence, upon the action of $E^\alpha$ , a point of $P$ is translated to another point of $P$ . In the following, we denote by $P_+$ the set of dominant weights

$$
P _ {+} = \mathbb {Z} _ {+} \omega_ {1} + \dots + \mathbb {Z} _ {+} \omega_ {r} \tag {13.76}
$$

For the algebras $G_2, F_4$ , and $E_8$ , it turns out that $Q = P$ . In all other cases, $Q$ is a proper subset of $P$ , and the ratio $P/Q$ is a finite group. Its order, $|P/Q|$ , is equal to the determinant of the Cartan matrix. Actually, it is isomorphic to the center of the group of the algebra under consideration (whose structure will be studied in more detail later). The distinct elements of the coset $P/Q$ define the so-called congruence classes (often called conjugacy classes). A weight $\lambda$ lies in exactly one congruence class. For instance, for $su(2)$ there are two congruence classes given by $\lambda_1 \bmod 2$ (integer or half-integer spins). For $su(3)$ , there are three classes, defined by the triality: $\lambda_1 + 2\lambda_2 \bmod 3$ . The $su(N)$ generalization is

$$
\lambda_ {1} + 2 \lambda_ {2} + \dots + (N - 1) \lambda_ {N - 1} \bmod N \tag {13.77}
$$

For any algebra $\mathbf{g}$ , the congruence classes take the form

$$
\lambda \cdot v = \sum_ {i = 1} ^ {r} \lambda_ {i} v _ {i} \mod | P / Q | \quad (\mathrm {m o d} \mathbb {Z} _ {2} \quad \text {f o r} \quad g = D _ {2 \ell}) \tag {13.78}
$$

where the vector $(\nu_{1},\dots ,\nu_{r})$ , equal to $(1,2,\dots ,N - 1)$ for $su(N)$ , is called the congruence vector. The congruence classes are tabulated in App. 13.A for all simple Lie algebras.

On the other hand, since the bases $\{\omega_i\}$ and $\{\alpha_i^\vee\}$ are dual, $P$ and $Q^\vee$ are dual lattices. A lattice is said to be self-dual if it is equal to its dual. For simple Lie algebras, the weight lattice is self-dual only for $E_8$ .

# 13.1.10. Normalization Convention

Up to now, all the normalizations have been fixed with respect to the root square lengths. In order to fully fix the normalization, it is necessary to give a specific value to these lengths. We follow the standard convention in which the square length of the long roots is set equal to two. Given that $\theta$ is necessarily a long root, we thus fix our normalization by setting

$$
\left| \theta \right| ^ {2} = 2 \tag {13.79}
$$

With $|\alpha_i|^2 \leq 2$ , it follows from Eq. (13.34) that

$$
a _ {i} \geq a _ {i} ^ {\vee} \quad \Rightarrow \quad a _ {i} ^ {\vee} = 1 \quad \text {i f} \quad a _ {i} = 1 \tag {13.80}
$$

and similarly

$$
\alpha_ {i} ^ {\vee} = \alpha_ {i} \frac {a _ {i}}{a _ {i} ^ {\vee}} \quad \Rightarrow \quad \alpha_ {i} ^ {\vee} \geq \alpha_ {i} \tag {13.81}
$$

# 13.1.11. Examples

# EXAMPLE 1: $su(2)$

This is the only simple Lie algebra of rank 1. Its Cartan matrix is $A = (2)$ , meaning that the simple root $\alpha_{1}$ is related to the fundamental weight $\omega_{1}$ by

$$
\alpha_ {1} = 2 \omega_ {1} \tag {13.82}
$$

Since $|\alpha_1|^2 = 2$ , it follows that

$$
\left(\omega_ {1}, \omega_ {1}\right) = \frac {1}{2} \tag {13.83}
$$

The Weyl group is generated by the simple reflection $s_1$ , whose action on a weight $\lambda = \lambda_1\omega_1$ is

$$
s _ {1} \left(\lambda_ {1} \omega_ {1}\right) = \lambda_ {1} \omega_ {1} - \lambda_ {1} \alpha_ {1} = - \lambda_ {1} \omega_ {1} \tag {13.84}
$$

Because $s_1^2 = 1$ , $W$ contains only the two elements $\{1, s_1\}$ . The full system of roots is then seen to be given by $\Delta = \{\alpha_1, -\alpha_1\}$ .

The weight and the root lattices are displayed in Fig. 13.1. The weight lattice is composed of all the nodes, whereas the root lattice contains only those with a cross. The fundamental Weyl chamber is the positive part of the weight lattice (here one-dimensional).

![](images/cc2c013501057fbbda9a4486afb7e88c0175489244c115cc1665682af7c2e095.jpg)  
Figure 13.1. Weight and root lattices for $su(2)$ .

For subsequent reference, we give the explicit form of the commutation relations in different bases. In the Chevalley basis, it reads (dropping the superscript 1):

$$
[ e, f ] = h \quad , \quad [ h, e ] = 2 e \quad , \quad [ h, f ] = - 2 f \tag {13.85}
$$

On a state $|\lambda \rangle$ of weight $\lambda$ , the action of $h$ is:

$$
h | \lambda \rangle = \lambda_ {1} | \lambda \rangle \tag {13.86}
$$

In the Cartan-Weyl basis, the generators are (cf. Eq. (13.36) with $\alpha_{1} = \sqrt{2}$ ):

$$
H = h / \sqrt {2}, \quad E ^ {+} = e, \quad E ^ {-} = f \tag {13.87}
$$

with $E^{\pm} \equiv E^{\pm \alpha_{1}}$ . The commutation relations are thus

$$
\left[ E ^ {+}, E ^ {-} \right] = \sqrt {2} H, \quad \left[ H, E ^ {\pm} \right] = \pm \sqrt {2} E ^ {\pm} \tag {13.88}
$$

and

$$
H | \lambda \rangle = \lambda^ {1} | \lambda \rangle = (\lambda_ {1} / \sqrt {2}) | \lambda \rangle \tag {13.89}
$$

Another frequently used basis in the case of $su(2)$ , which we call the spin basis, is defined by

$$
J ^ {0} = H / \sqrt {2}, \quad J ^ {\pm} = E ^ {\pm} \tag {13.90}
$$

This yields

$$
\left[ J ^ {+}, J ^ {-} \right] = 2 J ^ {0}, \quad \left[ J ^ {0}, J ^ {\pm} \right] = \pm J ^ {\pm} \tag {13.91}
$$

and on the state $|\lambda \rangle = |j,m\rangle$ , the action of the generators is

$$
\begin{array}{l} J ^ {0} | j, m \rangle = m | j, m \rangle \\ J ^ {\pm} | j, m \rangle = \sqrt {(j (j + 1) - m (m \pm 1)} | j, m \pm 1 \rangle \tag {13.92} \\ \end{array}
$$

# EXAMPLE 2: $su(3)$

The Cartan matrix for this rank-2 algebra is

$$
A = \left( \begin{array}{c c} 2 & - 1 \\ - 1 & 2 \end{array} \right) \tag {13.93}
$$

The simple roots $\alpha_{1}$ and $\alpha_{2}$ have the same length (the algebra is simply laced) and they are related to the fundamental weights by

$$
\begin{array}{l} \alpha_ {1} = \alpha_ {1} ^ {\vee} = 2 \omega_ {1} - \omega_ {2} = (2, - 1) \\ \alpha_ {2} ^ {\prime} = \alpha_ {2} ^ {\vee} = - \omega_ {1} + 2 \omega_ {2} = (- 1, 2) \tag {13.94} \\ \end{array}
$$

The scalar products between fundamental weights are

$$
\left(\omega_ {1}, \omega_ {1}\right) = \left(\omega_ {2}, \omega_ {2}\right) = \frac {2}{3}, \quad \left(\omega_ {1}, \omega_ {2}\right) = \frac {1}{3}. \tag {13.95}
$$

The full Weyl group is given by

$$
W = \{1, s _ {1}, s _ {2}, s _ {1} s _ {2}, s _ {2} s _ {1}, s _ {1} s _ {2} s _ {1} \} \tag {13.96}
$$

This follows from the relation

$$
(s _ {1} s _ {2}) ^ {3} = 1 \quad \Longrightarrow \quad s _ {1} s _ {2} s _ {1} = s _ {2} s _ {1} s _ {2} \tag {13.97}
$$

a consequence of Eq. (13.59), which implies that there are no strings of $s_i$ with more than three elements. This identity can also be checked directly by acting on

an arbitrary weight:

$$
s _ {1} \left(\lambda_ {1}, \lambda_ {2}\right) = \left(\lambda_ {1}, \lambda_ {2}\right) - \lambda_ {1} \alpha_ {1} = \left(- \lambda_ {1}, \lambda_ {1} + \lambda_ {2}\right)
$$

$$
s _ {2} \left(\lambda_ {1}, \lambda_ {2}\right) = \left(\lambda_ {1}, \lambda_ {2}\right) - \lambda_ {2} \alpha_ {2} = \left(\lambda_ {1} + \lambda_ {2}, - \lambda_ {2}\right)
$$

$$
s _ {1} s _ {2} \left(\lambda_ {1}, \lambda_ {2}\right) = \left(- \lambda_ {1} - \lambda_ {2}, \lambda_ {1}\right) \tag {13.98}
$$

$$
s _ {2} s _ {1} \left(\lambda_ {1}, \lambda_ {2}\right) = \left(\lambda_ {2}, - \lambda_ {1} - \lambda_ {2}\right)
$$

$$
s _ {1} s _ {2} s _ {1} \left(\lambda_ {1}, \lambda_ {2}\right) = s _ {2} s _ {1} s _ {2} \left(\lambda_ {1}, \lambda_ {2}\right) = \left(- \lambda_ {2}, - \lambda_ {1}\right)
$$

The action of the different elements of the Weyl group on the two simple roots gives all possible roots. For instance, $-\alpha_{1}$ and $\alpha_{1} + \alpha_{2}$ are roots because

$$
s _ {1} \alpha_ {1} = - \alpha_ {1}, \quad s _ {1} \alpha_ {2} = \alpha_ {1} + \alpha_ {2} \tag {13.99}
$$

In this way, $\Delta$ is found to be

$$
\Delta = \left\{\alpha_ {1}, \alpha_ {2}, \alpha_ {1} + \alpha_ {2}, - \alpha_ {1}, - \alpha_ {2}, - \alpha_ {1} - \alpha_ {2} \right\} \tag {13.100}
$$

The highest root is

$$
\theta = \alpha_ {1} + \alpha_ {2} \quad \Longrightarrow \quad a _ {i} = a _ {i} ^ {\vee} = 1, i = 1, 2 \tag {13.101}
$$

The root system and the Weyl chambers are presented in Fig. 13.2. The Weyl chambers are the regions separated by the dashed lines and they are specified here in terms of the elements of the Weyl group.

![](images/f2a6f0bfa19bb973be53a0e66d17d17bf243e535c5388e04485e9095cf9ab914.jpg)  
Figure 13.2. Root system and Weyl chambers $su(3)$ .

# EXAMPLE 3: $sp(4)$

This is again a rank-2 algebra, but it is not simply laced. The Cartan matrix is

$$
A = \left( \begin{array}{c c} 2 & - 1 \\ - 2 & 2 \end{array} \right) \tag {13.102}
$$

so that

$$
\cdot \alpha_ {1} = \frac {1}{2} \alpha_ {1} ^ {\vee} = 2 \omega_ {1} - \omega_ {2} = (2, - 1) \tag {13.103}
$$

$$
\alpha_ {2} = \alpha_ {2} ^ {\vee} = - 2 \omega_ {1} + 2 \omega_ {2} = (- 2, 2)
$$

Because the long root is $\alpha_{2}$

$$
\left| \alpha_ {2} \right| ^ {2} = 2 \quad \Longrightarrow \quad \left| \alpha_ {1} \right| ^ {2} = 1 \tag {13.104}
$$

The components of the quadratic form matrix are

$$
\left(\omega_ {1}, \omega_ {1}\right) = \left(\omega_ {1}, \omega_ {2}\right) = \frac {1}{2}, \quad \left(\omega_ {2}, \omega_ {2}\right) = 1 \tag {13.105}
$$

On the other hand, the complete structure of the Weyl group is easily recovered from the equality

$$
(s _ {1} s _ {2}) ^ {4} = 1 \tag {13.106}
$$

meaning that the longest element is $s_1s_2s_1s_2$ ; hence,

$$
W = \left\{1, s _ {1}, s _ {2}, s _ {1} s _ {2}, s _ {2} s _ {1}, s _ {1} s _ {2} s _ {1}, s _ {2} s _ {1} s _ {2}, s _ {1} s _ {2} s _ {1} s _ {2} \right\} \tag {13.107}
$$

Having determined the Weyl group, the set $\Delta$ can be constructed

$$
\Delta = \left\{\alpha_ {1}, \alpha_ {2}, \alpha_ {1} + \alpha_ {2}, 2 \alpha_ {1} + \alpha_ {2}, - \alpha_ {1}, - \alpha_ {2}, - \alpha_ {1} - \alpha_ {2}, - 2 \alpha_ {1} - \alpha_ {2} \right\} \tag {13.108}
$$

The highest root is thus

$$
\theta = 2 \alpha_ {1} + \alpha_ {2} = 2 \alpha_ {1} ^ {\vee} + \alpha_ {2} ^ {\vee} \quad \Longrightarrow \quad a _ {1} = 2, a _ {2} = a _ {1} ^ {\vee} = a _ {2} ^ {\vee} = 1 \tag {13.109}
$$

In this case, the root vectors separate the Weyl chambers, as can be seen in Fig. 13.3.

![](images/0723a736c8f47a6f33244edf1b47a62fceaa5c897b1778cd204d8cd70d9a0961.jpg)  
Figure 13.3. Root system and Weyl chambers $sp(4)$ .

# §13.2. Highest-Weight Representations

Any finite-dimensional irreducible representation has a unique highest-weight state $|\lambda \rangle$ . Being nondegenerate, $|\lambda \rangle$ is completely specified by its eigenvalues (Dynkin labels) $\lambda(h^i) = \lambda_i$ . Among all the weights in the representation, the highest weight is the one for which the sum of the coefficient expansions in the basis of simple roots is maximal. As a result, for any $\alpha > 0$ , $\lambda + \alpha$ cannot be a weight in the representation, so that

$$
E ^ {\alpha} | \lambda \rangle = 0, \quad \forall \alpha > 0 \tag {13.110}
$$

From Eq. (13.27), it is clear that the highest weight of a finite-dimensional representation is necessarily dominant (i.e., with positive-integer Dynkin labels). Moreover, to each dominant weight $\lambda$ there corresponds a unique irreducible finite-dimensional representation $L_{\lambda}$ whose highest weight is $\lambda$ . By abuse of notation, we will often specify a representation by its highest weight. The highest weight for the adjoint representation is $\theta$ .

# 13.2.1. Weights and Their Multiplicities

Starting from the highest-weight state $|\lambda \rangle$ , all the states in the representation space (or irreducible module) $\mathsf{L}_{\lambda}$ can be obtained by the action of the lowering operators of $\mathbf{g}$ as

$$
E ^ {- \beta} E ^ {- \gamma} \dots E ^ {- \eta} | \lambda \rangle \quad \text {f o r} \quad \beta , \gamma , \eta \in \Delta_ {+} \tag {13.111}
$$

The set of eigenvalues of all the states in $\mathsf{L}_{\lambda}$ is the weight system, written $\Omega_{\lambda}$ . Any weight $\lambda'$ in the set $\Omega_{\lambda}$ is such that $\lambda - \lambda' \in \Delta_{+}$ . An immediate consequence is that all the weights of a given representation lie in exactly one congruence class, that is, one element of the coset $P/Q$ .

In order to find all the weights $\lambda' \in \Omega_{\lambda}$ , the key relation is again Eq. (13.27), which can be rewritten as

$$
\left(\lambda^ {\prime}, \alpha_ {i} ^ {\vee}\right) = \lambda_ {i} ^ {\prime} = - \left(p _ {i} - q _ {i}\right), \quad p _ {i}, q _ {i} \in \mathbb {Z} _ {+} \tag {13.112}
$$

As already mentioned, $\lambda'$ is necessarily of the form $\lambda - \sum n_i \alpha_i$ , with $n_i \in \mathbb{Z}_+$ . If we call $\sum n_i$ the level of the weight $\lambda'$ in the representation $\lambda$ , proceeding level by level, we know at each step the value of $p_i$ . Clearly, $\lambda' - \alpha_i$ is also a weight if $q_i$ is nonzero, that is, if $\lambda_i' - p_i > 0$ .

With this criterion, the systematic construction of all the weights in the representation can be done by means of the following algorithm. We start with the highest weight $\lambda = (\lambda_1, \dots, \lambda_r)$ . For each positive Dynkin label $\lambda_i > 0$ , we construct the sequence of weights $\lambda - \alpha_i, \lambda - 2\alpha_i, \dots, \lambda - \lambda_i\alpha_i$ , which all belong to $\Omega_\lambda$ . The process is then repeated with $\lambda$ replaced by each of the weights just obtained, and iterated until no more weights with positive Dynkin labels are produced. Simple examples will clarify the method. Consider the adjoint representation of $su(3)$ , whose highest weight is (1, 1). The weights obtained at each step can be read from

![](images/9d8117d9ead5bd09579e3ebe20259c2c6b2c62ddefd2f3d2bbf352bce83ea6f6.jpg)  
Figure 13.4. Weights in the adjoint representation of $su(3)$ .

![](images/7f73f7c4fb22f22cc3d896e3167e30438b67e70c3c5e948d7209f1b416591a8f.jpg)  
Figure 13.5. Weights in the adjoint representation of $sp(4)$ .

Fig. 13.4. Similarly, Fig. 13.5 displays the weights in the adjoint representation of $sp(4)$ .

However, this procedure does not keep track of multiplicities. For this, one can use the Freudenthal recursion formula, whose origin will be indicated in Sect. 13.2.3, and which gives the multiplicity of $\lambda'$ in the representation $\lambda$ in

terms of the multiplicity of all the weights above it:

$$
\boxed {[ | \lambda + \rho | ^ {2} - | \lambda^ {\prime} + \rho | ^ {2} ] \operatorname {m u l t} _ {\lambda} (\lambda^ {\prime}) = 2 \sum_ {\alpha > 0} \sum_ {k = 1} ^ {\infty} \left(\lambda^ {\prime} + k \alpha , \alpha\right) \operatorname {m u l t} _ {\lambda} \left(\lambda^ {\prime} + k \alpha\right)} \tag {13.113}
$$

To illustrate the formula, we calculate the multiplicity of the weight $(0,0)$ in the adjoint representation of $su(3)$ . Having proceeded recursively, we know that $k$ can only be 1 and the three weights above $(0,0)$ have multiplicity 1. Furthermore, $(\lambda' + \alpha, \alpha) = 2$ for the three positive roots. Then, using $\lambda = \theta = \rho = \alpha_1 + \alpha_2$ , we easily find that

$$
(8 - 2) \operatorname {m u l t} _ {\theta} (0, 0) = 2 (2 + 2 + 2) \quad \Longrightarrow \quad \operatorname {m u l t} _ {\theta} (0, 0) = 2 \tag {13.114}
$$

Indeed, the zero eigenvalue in the adjoint representation always has multiplicity $r$ , being associated with the generators of the Cartan subalgebra (whereas the nonzero weights (roots) are nondegenerate). Another multiplicity formula is presented in Ex. 13.17.

We note that all the weights in a given $W$ orbit have the same multiplicity:

$$
\operatorname {m u l t} _ {\lambda} \left(w \lambda^ {\prime}\right) = \operatorname {m u l t} _ {\lambda} \left(\lambda^ {\prime}\right) \quad \text {f o r a l l} \quad w \in W \tag {13.115}
$$

This ultimately reflects the arbitrariness of the basis of simple roots, that is, that any set $\{w\alpha_{i}\}$ with $w$ fixed, could serve as a basis.

Finally, we mention that a finite-dimensional irreducible module $\mathsf{L}_{\lambda}$ is always unitary. This means that, with $(H^{i})^{\dagger} = H^{i}$ and $(E^{\alpha})^{\dagger} = E^{-\alpha}$ , the norm of any state $|\lambda^{\prime}\rangle$ in $\mathsf{L}_{\lambda}$ is positive definite:

$$
| \lambda^ {\prime} \rangle = E ^ {- \beta} \dots E ^ {- \gamma} | \lambda \rangle \Longrightarrow \langle \lambda^ {\prime} | \lambda^ {\prime} \rangle = \langle \lambda | E ^ {\gamma} \dots E ^ {\beta} E ^ {- \beta} \dots E ^ {- \gamma} | \lambda \rangle > 0 \tag {13.116}
$$

with $\beta, \gamma \in \Delta_{+}$ . This also holds for linear combinations of such states.

# 13.2.2. Conjugate Representations

In an irreducible finite-dimensional representation, there is obviously a lowest state, also unique. It lies in the $W$ orbit of the highest state, in the chamber exactly opposite to the fundamental one. This chamber is specified by the longest element of the Weyl group $\mathcal{w}_0$ . In terms of the highest state $\lambda$ , the lowest state is thus given by $\mathcal{w}_0\lambda$ . Turning a representation "upside down" produces the conjugate representation, indicated by $\lambda^*$ . Its highest-weight state is the negative of the lowest state of the original representation

$$
\lambda^ {*} = - \left(w _ {0} \lambda\right) = \left(- w _ {0}\right) \cdot \lambda \tag {13.117}
$$

since $\rho$ is the highest weight of a self-conjugate representation: $\rho = -w_0\rho$ . More generally, all the weights in $\Omega_{\lambda^*}$ are the negatives of those in $\Omega_{\lambda}$ . For $su(N), w_0$ is given by

$$
w _ {0} = s _ {1} s _ {2} \dots s _ {N - 1} s _ {1} s _ {2} \dots s _ {N - 2} \dots s _ {1} s _ {2} s _ {1} \tag {13.118}
$$

With $N = 3$ , it yields

$$
\left(- w _ {0}\right) \cdot \left(\lambda_ {1}, \lambda_ {2}\right) = - s _ {1} s _ {2} s _ {1} \left(\lambda_ {1} + 1, \lambda_ {2} + 1\right) - (1, 1) = \left(\lambda_ {2}, \lambda_ {1}\right) \tag {13.119}
$$

The conjugation is related to the reflection symmetry of the Dynkin diagram. This readily shows that for $su(N)$ , the conjugation amounts to reversing the order of the finite Dynkin labels. Because the Dynkin diagram of $so(2r + 1), sp(2r), so(4r), G_2, F_4, E_7$ , and $E_8$ have no symmetry, all representations of these algebras are self-conjugate. For the other algebras, self-conjugate representations are those with highest weight satisfying:

$$
\begin{array}{l} s u (r + 1): \quad \lambda_ {i} = \lambda_ {r - i} \\ \operatorname {s o} (4 r + 2): \quad \lambda_ {r} = \lambda_ {r - 1} \tag {13.120} \\ E _ {6}: \quad \lambda_ {1} = \lambda_ {5}, \quad \lambda_ {2} = \lambda_ {4} \\ \end{array}
$$

# 13.2.3. Quadratic Casimir Operator

A generalization of the $su(2)$ quadratic Casimir operator $\mathcal{Q}$ can be constructed for any semisimple Lie algebra. Up to a scale factor, it is uniquely characterized by its commutativity with all the generators of the algebra. In a generic basis $\{\mathcal{L}^a\}$ , it can be checked to be given by

$$
\boxed {Q = \sum_ {a, b} [ K (\mathcal {L} ^ {a}, \mathcal {L} ^ {b}) ] ^ {- 1} \mathcal {L} ^ {a} \mathcal {L} ^ {b}} \tag {13.121}
$$

where $K$ is the Killing form (which, as already mentioned, is nondegenerate for semisimple Lie algebras). In the orthonormal $\{J^a\}$ basis, it is thus

$$
Q = \sum_ {a} J ^ {a} J ^ {a} \tag {13.122}
$$

On the other hand, in the Cartan-Weyl basis, it reads

$$
Q = \sum_ {i} H ^ {i} H ^ {i} + \sum_ {\alpha > 0} \frac {| \alpha | ^ {2}}{2} \left(E ^ {\alpha} E ^ {- \alpha} + E ^ {- \alpha} E ^ {\alpha}\right) \tag {13.123}
$$

We note that $\mathcal{Q}$ is not an element of $\mathbf{g}$ itself; it lies in its universal enveloping algebra, which is the set of all formal power series in elements of $\mathbf{g}$ .

Since $\mathcal{Q}$ commutes with all the generators of the algebra, its eigenvalue is the same on all the states of an irreducible representation. It is most easily evaluated on the highest-weight state, using the Cartan-Weyl basis. First, we have

$$
\sum_ {i} H ^ {i} H ^ {i} | \lambda \rangle = \sum_ {i} \lambda^ {i} \lambda^ {i} | \lambda \rangle = (\lambda , \lambda) | \lambda \rangle \tag {13.124}
$$

Because $E^{\alpha}|\lambda \rangle = 0$ for $\alpha > 0$ , the term $E^{-\alpha}E^{\alpha}$ does not contribute. For the remaining term, we move $E^{\alpha}$ to the right of $E^{-\alpha}$ using

$$
\left[ E ^ {\alpha}, E ^ {- \alpha} \right] | \lambda \rangle = \frac {2}{| \alpha | ^ {2}} \alpha \cdot H | \lambda \rangle = \frac {2}{| \alpha | ^ {2}} (\alpha , \lambda) | \lambda \rangle \tag {13.125}
$$

The result is

$$
\mathcal {Q} | \lambda \rangle = [ (\lambda , \lambda) + \sum_ {\alpha > 0} (\alpha , \lambda) ] | \lambda \rangle \tag {13.126}
$$

By using the definition (13.47) of the Weyl vector, we can write

$$
\mathcal {Q} | \lambda \rangle = (\lambda , \lambda + 2 \rho) | \lambda \rangle \tag {13.127}
$$

In the adjoint representation, the eigenvalue of the Casimir operator is

$$
\begin{array}{l} (\theta , \theta + 2 \rho) = 2 + 2 (\theta , \rho) = 2 + 2 \sum_ {i, j} a _ {i} ^ {\vee} \left(\alpha_ {i} ^ {\vee}, \omega_ {j}\right) \\ = 2 + 2 \sum_ {i} \alpha_ {i} ^ {\vee} = 2 + 2 (g - 1) = 2 g \tag {13.128} \\ \end{array}
$$

The quadratic Casimir operator does not distinguish a representation from its conjugate

$$
\mathcal {Q} | \lambda^ {*} \rangle = \mathcal {Q} | \lambda \rangle \tag {13.129}
$$

This follows from the equality

$$
\left| \lambda^ {*} + \rho \right| ^ {2} = \left| \lambda + \rho \right| ^ {2} \tag {13.130}
$$

which is itself a simple consequence of Eq. (13.117):

$$
\lambda^ {*} + \rho = (- w _ {0}) \cdot \lambda + \rho = - w _ {0} (\lambda + \rho) \tag {13.131}
$$

and of the invariance of the scalar product with respect to the Weyl group: $|w\mu|^2 = |\mu|^2$ .

The Freudenthal formula (13.113) is obtained by evaluating the trace of $\mathcal{Q}$ in the subspace associated with the weight $\lambda'$ , first using the eigenvalue just obtained and then using the explicit form of $\mathcal{Q}$ in the Cartan-Weyl basis.

For $su(2)$ , the quadratic Casimir operator is the unique operator that commutes with all the generators. However, we mention that for higher-rank algebras there exist Casimir operators of higher degree. Their degrees minus one are called the exponents of the algebras (tabulated in App. 13.A).<sup>8</sup>

# 13.2.4. Index of a Representation

The quadratic Casimir operator enters in the definition of an important quantity, the index of a representation, which gives the relative normalization of invariant bilinear products taken in different representations.

As already stressed, once a normalization is fixed for the length of the long roots, every product is uniquely determined. In particular, the normalization of the invariant bilinear form $\operatorname{Tr}_{\lambda}(\mathcal{R}(J^a)\mathcal{R}(J^b))$ , for $\mathcal{R}(J^a)$ standing for a matrix

representation of the generator $J^a$ , must be fixed. Here the trace is evaluated in $\mathsf{L}_{\lambda}$ . The relative normalization of this product with respect to $|\theta|^2$ defines the Dynkin index $x_{\lambda}$ of the representation $\lambda$

$$
\operatorname {T r} _ {\lambda} \left(\mathcal {R} \left(J ^ {a}\right) \mathcal {R} \left(J ^ {b}\right)\right) = | \theta | ^ {2} x _ {\lambda} \delta_ {a b} = 2 x _ {\lambda} \delta_ {a b} \tag {13.132}
$$

An explicit expression for $x_{\lambda}$ can be easily obtained by setting $a = b$ and summing over all values of $a$ . The l.h.s. becomes equal to the trace of the quadratic Casimir, so that

$$
x _ {\lambda} = \frac {\dim | \lambda | (\lambda , \lambda + 2 \rho)}{2 \dim g} \tag {13.133}
$$

We note that the Dynkin index of the adjoint representation, $\lambda = \theta$ , is simply the dual Coxeter number

$$
x _ {\theta} = g \tag {13.134}
$$

since $\dim |\theta| = \dim g$ and $(\theta, \rho) = g - 1$ (cf. Eq (13.128)).

# §13.3. Tableaux and Patterns (su(N))

In this section, we introduce a useful diagrammatic representation of highest weights, which will also turn out to be a powerful combinatorial tool, particularly efficient in tensor-product calculations. A refinement of this diagrammatic representation leads to the simple construction of a complete basis of states in a finite-dimensional representation. This will be shown to be equivalent to a description of states in terms of triangular arrays of numbers, the so-called Gelfand-Tsetlin patterns. For simplicity, we restrict the whole discussion to $su(N)$ .

# 13.3.1. Young Tableaux

A $su(N)$ integrable highest weight $\lambda$ , with Dynkin labels

$$
\lambda = \left(\lambda_ {1}, \dots , \lambda_ {N - 1}\right) \tag {13.135}
$$

can equally well be specified in terms of its partition

$$
\lambda = \left\{\ell_ {1}; \ell_ {2}; \dots ; \ell_ {N - 1} \right\} \tag {13.136}
$$

where

$$
\ell_ {i} = \lambda_ {i} + \lambda_ {i + 1} + \dots + \lambda_ {N - 1} \tag {13.137}
$$

To a partition, we associate a Young tableau, which is a box array of rows lined up on the left, such that the length of the $i$ -th row is equal to $\ell_i$ . For example, to the $su(5)$ weight

$$
\lambda = (2, 0, 2, 0) = \{4; 2; 2 \} \tag {13.138}
$$

(zero entries in partitions being generally omitted) corresponds the Young tableau

$$
\begin{array}{c c c c} \hline & & & \\ \hline & & \\ \hline & & \end{array} \tag {13.139}
$$

Dynkin labels provide a dual description of the tableau: $\lambda_{i}$ gives the number of columns of $i$ boxes. The fundamental representation $\omega_{\ell}$ is described by a single column of $\ell$ boxes. To the scalar representation corresponds a void tableau or, equivalently, a single column of $N$ boxes. Allowing for columns of $N$ boxes, partitions are fixed by $N$ integers. But clearly, when $\ell_{N} \neq 0$ , we can always subtract $\{\ell_{N}; \ell_{N}; \dots; \ell_{N}\}$ ( $N$ entries) from the partition, which just amounts to eliminating columns of $N$ boxes; for instance,

$$
\{5; 3; 3; 1; 1 \} = \{4; 2; 2 \} \tag {13.140}
$$

Tableaux with $\ell_N = 0$ will be referred to as reduced tableaux, and likewise for partitions.

The transpose of a Young tableau is obtained by interchanging rows and columns. We denote by $\lambda^t$ the corresponding weight. For instance, the transpose of $\lambda = \{4; 2; 2\}$ is

$$
\begin{array}{c c c} \hline & & \\ \hline & & \\ \hline & & \leftrightarrow \lambda^ {t} = \{3; 3; 1; 1 \} \\ \hline & & \end{array} \tag {13.141}
$$

With

$$
\lambda^ {t} = \{\tilde {\ell} _ {1}; \tilde {\ell} _ {2}; \dots \} \tag {13.142}
$$

it is not difficult to see that

$$
\tilde {\ell} _ {i} = \text {n u m b e r o f} \ell_ {j} \text {s u c h t h a t} \ell_ {j} \geq i \tag {13.143}
$$

# 13.3.2. Partitions and Orthonormal Bases

The entries of the partition (13.136) are the expansion coefficients of a dominant weight in a certain basis, which we now describe. We also indicate how partitions can be associated with nondominant integral weights, providing a rationale for the construction of the next section.

Elements of $su(N)$ can be represented by $N \times N$ traceless matrices. In this representation, the Cartan subalgebra is spanned by the set of all diagonal traceless matrices. We let $e_{ij}$ stand for the matrix with 0 everywhere, except for a single 1 at position $(i,j)$ ( $i$ -th row, $j$ -th column). With this notation, the elements of the Cartan subalgebra are of the form $\sum_{i=1}^{N} \epsilon_i e_{ii}$ with $\sum_{i=1}^{N} \epsilon_i = 0$ . The ladder operators are represented by the matrices $e_{ij}, i \neq j$ . The roots are then given by $\epsilon_i - \epsilon_j, i \neq j$ , and a basis of simple roots is

$$
\alpha_ {i} = \epsilon_ {i} - \epsilon_ {i + 1}, \quad i = 1, \dots , N - 1 \tag {13.144}
$$

Generalizing this point of view, we can consider the $\epsilon_{i}$ as orthonormal vectors in an $(r + 1)$ -dimensional space, and in terms of these vectors, the root lattice is simply

$$
Q = \sum_ {i = 1} ^ {N} n _ {i} \epsilon_ {i} \quad \text {w i t h} \quad n _ {i} \in \mathbb {Z} \quad \text {a n d} \quad \sum_ {i = 1} ^ {N} n _ {i} = 0 \tag {13.145}
$$

With $\epsilon_i^2 = 1$ and $\epsilon_i \cdot \epsilon_j = 0$ for $i \neq j$ , we see that $|\alpha|^2 = 2$ for any root. The fundamental weights are related to the simple roots by the quadratic form matrix (since here roots are the same as coroots), which leads to

$$
\omega_ {i} = \epsilon_ {1} + \epsilon_ {2} + \dots + \epsilon_ {i} - \frac {i}{N} \sum_ {i = 1} ^ {N} \epsilon_ {i} \tag {13.146}
$$

Hence, the expansion coefficients of a highest weight in the $\{\epsilon_i\}$ basis are exactly the entries of the partition

$$
\cdot \lambda = \sum_ {i = 1} ^ {N - 1} \lambda_ {i} \omega_ {i} = \sum_ {i = 1} ^ {N} \left(\ell_ {i} - \kappa\right) \epsilon_ {i} \tag {13.147}
$$

where the $\ell_i$ are related to the Dynkin labels by Eq. (13.137) and $\kappa$ is

$$
\kappa = \frac {1}{N} \sum_ {j = 1} ^ {N - 1} j \lambda_ {j} \tag {13.148}
$$

A well-defined partition is thus associated with the highest weight $\lambda$ of each representation. The other weights in the representation are obtained by subtracting from $\lambda$ the positive roots $\epsilon_{i} - \epsilon_{j}, i < j$ . This construction gives directly their expansion coefficients in the $\{\epsilon_{i}\}$ basis. A weight $\lambda'$ can thus be described by a partition $\{\ell_1', \ell_2', \dots, \ell_N'\}$ . We stress that the partition of a weight that is not a highest weight is not related to the shape of a Young tableau. In particular, such a partition is no longer bound to satisfy $\ell_i' \geq \ell_{i+1}'$ .

# 13.3.3. Semistandard Tableaux

We now indicate how tableau techniques can be used to explicitly describe all the states in a representation. This involves filling the boxes of a Young tableau with positive integers, generating the so-called semistandard tableaux. They are defined as follows. We let $c_{i,j}$ be the integer appearing in the box on the $i$ -th row (from top) and the $j$ -th column (from left), and satisfying

$$
1 \leq c _ {i, j} \leq N, \quad c _ {i, j} \leq c _ {i, j + 1}, \quad c _ {i, j} <   c _ {i + 1, j} \tag {13.149}
$$

In other words, the numbers are nondecreasing from left to right and strictly increasing from top to bottom.

Semistandard tableaux of shape $\lambda$ are in one-to-one correspondence with the states in the module $\mathsf{L}_{\lambda}$ . The numbering in the semistandard tableaux encodes the partition of the corresponding weight. We can think of a box with number $i$ as

representing $\epsilon_{i}$ . The number of $i$ 's in the semistandard tableau of weight $\lambda' = \{\ell_1', \dots, \ell_N'\}$ is given by $\ell_i'$ . In the semistandard tableau representing the highest weight $\lambda$ , all boxes of the $i$ -th row have number $i$ . The weight of a semistandard tableau is clearly obtained by adding the weights of all its boxes. The weight of a box marked with a $i$ is

$$
\epsilon_ {i} = \omega_ {i} - \omega_ {i - 1}, \quad i = 1, \dots , N \tag {13.150}
$$

(modulo $\sum \epsilon_{i} / N$ ) with $\omega_0 = \omega_N = 0$

The number of semistandard tableaux of a fixed shape that can be constructed with a given partition gives the multiplicity of the corresponding weight in the representation. In other words, the rules for constructing semistandard tableaux provide a combinatorial realization of the Freudenthal multiplicity formula (13.113).

We consider for example $su(3)$ . The semistandard tableaux of the three states in the representation $\omega_{1}$ are

$$
\boxed {1} \leftrightarrow (1, 0), \quad \boxed {2} \leftrightarrow (- 1, 1), \quad \boxed {3} \leftrightarrow (0, - 1) \tag {13.151}
$$

whereas those in the representation $\omega_{2}$ are

$$
\begin{array}{l} \boxed {1} \\ \hline 2 \end{array} \leftrightarrow (0, 1), \quad \boxed {1} \\ \hline 3 \end{array} \leftrightarrow (1, - 1), \quad \boxed {2} \\ \hline 3 \end{array} \leftrightarrow (- 1, 0) \tag {13.152}
$$

The adjoint representation $(1,1)$ contains the 8 semistandard tableaux:

$$
\begin{array} { c } \framebox { 1 } \framebox { 1 } \\ \framebox { 2 } \end{array}
$$

$$
(1, 1)
$$

$$
\begin{array}{c c} \hline 1 & 2 \\ \hline 2 & \end{array}
$$

$$
(- 1, 2)
$$

$$
\begin{array} { c } \framebox { 1 } \framebox { 3 } \\ \framebox { 2 } \end{array}
$$

$$
(0, 0)
$$

$$
\begin{array} { c } \framebox { 1 } \framebox { 1 } \\ \framebox { 3 } \end{array}
$$

$$
(2, - 1)
$$

$$
\begin{array}{c c} \hline 1 & 2 \\ \hline 3 & \end{array}
$$

$$
(\bar {0}, 0)
$$

$$
\begin{array}{c c} \hline 1 & 3 \\ \hline 3 & \end{array}
$$

$$
(1, - 2)
$$

$$
\begin{array} { c } \framebox { 2 } \framebox { 2 } \\ \framebox { 3 } \end{array}
$$

$$
(- 2, 1)
$$

$$
\begin{array} { c } \framebox { 2 } \framebox { 3 } \\ \framebox { 3 } \end{array}
$$

$$
(- 1, - 1)
$$

(13.153)

Two distinct semistandard tableaux, that is, two distinct states, correspond to the doubly degenerate weight $(0,0)$ . Similarly, to the weight $(0,0,0)$ of the $su(4)$ representation $(2,0,2)$ (with partition $\{2;2;2;2\}$ ) correspond 6 semistandard tableaux:

![](images/be7c6b7d263c377ea84af6ec2bec9554181c6c0d55c5e6eac343069628aff05a.jpg)

![](images/8afd2effc44c92136255da33c5a6834feb65f32e28aa2be685cb4335590b0616.jpg)

![](images/013b289689d3d4c98b23fafeef9c0930ac2b6be3401ec62f1036433ac94e5233.jpg)

![](images/5c8fc4db9038c64ba805d742d129e5df2297f7ef7d0eac48d55e11d5f9b66405.jpg)

![](images/498f9fd23732279f2ea44ae859c62f1fb927e996416dbf3be9d0c21a16755552.jpg)

![](images/8e2260cb71d8f122f51337c189a6294f24b766f630b8db0d763c590a9b35448a.jpg)

# 13.3.4. Gelfand-Tsetlin Patterns

An equivalent representation of the basis of semistandard tableaux is given by the Gelfand-Tsetlin patterns. To a given semistandard tableau we associate the following triangular array of numbers:

$$
\begin{array}{l} \beta_ {1} ^ {(N)} \beta_ {2} ^ {(N)} \dots \dots \beta_ {N} ^ {(N)} \\ \beta_ {1} ^ {(N - 1)} \dots \beta_ {N - 1} ^ {(N - 1)} \\ \begin{array}{c} \dots \dots \\ \beta_ {1} ^ {(2)} \beta_ {2} ^ {(2)} \\ \beta_ {1} ^ {(1)} \end{array} \tag {13.155} \\ \end{array}
$$

such that $\beta_i^{(j)}$ is the number of boxes containing numbers less or equal to $j$ in the $i$ -th row (from top) of the semistandard tableau. For instance, the following semistandard tableau and Gelfand-Tsetlin pattern corresponding to the weight $(-2,1,0)$ in the representation $(1,2,1)$ of $su(4)$ are equivalent:

$$
\begin{array}{c c} \framebox {1} & \framebox {2} & \framebox {2} & \framebox {4} \\ \framebox {2} & \framebox {3} & \framebox {4} \\ \framebox {3} & \end{array} \quad \leftrightarrow \quad \begin{array}{c c} 4 & 3 & 1 & 0 \\ 3 & 2 & 1 \\ 3 & 1 \\ 1 & \end{array} \tag {13.156}
$$

The first line in the Gelfand-Tsetlin pattern is common to all patterns in the representation, being simply the partition of the tableau

$$
\beta_ {i} ^ {(N)} = \ell_ {i} \tag {13.157}
$$

All the states in a representation are then generated by filling a triangular array of $N$ lines with integers $\beta_{i}^{(j)}$ satisfying

$$
\beta_ {i} ^ {(j)} \geq \beta_ {i + 1} ^ {(j)} \quad \beta_ {i} ^ {(j)} \geq \beta_ {i + 1} ^ {(j + 1)} \tag {13.158}
$$

with the first line fixed by the partition. In this way, the 8 patterns of the $(1,1)$ representation of $su(3)$ are found to be

$$
\begin{array}{c c c c c c c c} 2 1 0 & 2 1 0 & 2 1 0 & 2 1 0 & 2 1 0 & 2 1 0 & 2 1 0 & 2 1 0 \\ 2 1 & 2 1 & 1 1 & 2 0 & 2 0 & 1 0 & 2 0 & 1 0 \\ 2 & 1 & 1 & 2 & 1 & 1 & 0 & 0 \end{array} \tag {13.159}
$$

Their ordering corresponds to the semistandard tableaux (13.153). We note that the Gelfand-Tsetlin pattern of the highest-weight state in the representation is completely fixed by the partition, being

$$
\begin{array}{l} \ell_ {1} \ell_ {2} \dots \dots \ell_ {N} \\ \ell_ {1} \dots \ell_ {N - 1} \\ \dots \dots \tag {13.160} \\ \ell_ {1} \ell_ {2} \\ \ell_ {1} \\ \end{array}
$$

# §13.4. Characters

# 13.4.1. Weyl's Character Formula

A character is a useful functional way of coding the whole content of a representation. The character of the representation of highest weight $\lambda$ is formally defined

as

$$
\chi_ {\lambda} = \sum_ {\lambda^ {\prime} \in \Omega_ {\lambda}} \operatorname {m u l t} _ {\lambda} \left(\lambda^ {\prime}\right) e ^ {\lambda^ {\prime}} \tag {13.161}
$$

where the sum is over all the weights of the representation. $e^{\lambda}$ denotes a formal exponential satisfying

$$
e ^ {\lambda} e ^ {\mu} = e ^ {\lambda + \mu} \tag {13.162}
$$

$$
e ^ {\lambda} (\xi) = e ^ {(\lambda , \xi)}
$$

On the r.h.s. of the last expression, $e$ is a genuine exponential function, and $\xi$ is an arbitrary element of the dual Cartan subalgebra (i.e., an arbitrary weight).

This formal character is related to the familiar character in the representation theory of groups as follows. Let $G$ be the Lie group of $\mathfrak{g}$ and $H$ an element of the Cartan subgroup of $G$ . The character of $H$ in some representation is simply its trace evaluated in the corresponding module $V$ :

$$
\operatorname {T r} _ {V} H = \sum_ {\gamma} \operatorname {m u l t} (\gamma) [ \gamma (H) ] \tag {13.163}
$$

where $\gamma(H)$ denotes the eigenvalues of $H$ . Complete information about the representation is obtained by considering the group character as restricted to the full Cartan subgroup. Since $H$ is associated with an element $h$ of the Cartan subalgebra of $\mathfrak{g}$ by $H = \exp(h)$ , spanning the full Cartan subgroup amounts to replacing the single element $h$ by the vector $\tilde{h} = (h^1, h^2, \dots, h^r)$ . Thus $\gamma(H)$ is replaced by $e^{\lambda'(\tilde{h})} = e^{(\lambda_1', \dots, \lambda_r')}$ where $\lambda'$ is a weight. For a vector exponent, $e$ must be regarded as a formal exponential.

The expression (13.161) for the character can be brought into a more manageable form in two steps (we omit the details). At first, the auxiliary quantity

$$
D _ {\rho} = \prod_ {\alpha > 0} \left(e ^ {\alpha / 2} - e ^ {- \alpha / 2}\right) \tag {13.164}
$$

is introduced, and shown to be expressible as a sum over the elements of the Weyl group

$$
D _ {\rho} = \sum_ {w \in W} \epsilon (w) e ^ {w \rho} \tag {13.165}
$$

The second step (more involved) consists of showing, using the Freudenthal multiplicity formula (13.113), that

$$
D _ {\rho} \chi_ {\lambda} = D _ {\lambda + \rho} \tag {13.166}
$$

where $D_{\lambda + \rho}$ is defined from Eq. (13.165) with $\rho$ replaced by $\lambda + \rho$ . This last result is the famous Weyl character formula

$$
\chi_ {\lambda} = \frac {D _ {\lambda + \rho}}{D _ {\rho}} = \frac {\sum_ {w \in W} \epsilon (w) e ^ {w (\lambda + \rho)}}{\sum_ {w \in W} \epsilon (w) e ^ {w \rho}} \tag {13.167}
$$

For $su(2)$ , with $t = e^{\omega_1}$ , this becomes

$$
\chi_ {\lambda} = \frac {t ^ {\lambda_ {1} + 1} - t ^ {- \lambda_ {1} - 1}}{t - t ^ {- 1}} = t ^ {\lambda_ {1}} + t ^ {\lambda_ {1} - 2} + \dots + t ^ {- \lambda_ {1}} \tag {13.168}
$$

For some manipulations, it is more convenient to work with the character evaluated at a special but arbitrary value $\xi$

$$
\chi_ {\lambda} (\xi) = \frac {\sum_ {w \in W} \epsilon (w) e ^ {(w (\lambda + \rho) , \xi)}}{\sum_ {w \in W} \epsilon (w) e ^ {(w \rho , \xi)}} \tag {13.169}
$$

# 13.4.2. The Dimension and the Strange Formulae

As an immediate application, we derive a formula for the dimension of a representation. From Eq. (13.161), it is clear that this amounts to evaluating the character at the special point $\xi = 0$ . But setting $\xi = 0$ in Eq. (13.169) leads to an indeterminate expression since $\sum \epsilon(w) = 0$ ( $W$ has the same number of even and odd elements). Rather, a limiting process must be used. For this we set $\xi = t\rho$ and consider the limit $t \to 0$ . For $\xi$ proportional to $\rho$ , the character takes the simple form

$$
\chi_ {\lambda} (t \rho) = \frac {D _ {\lambda + \rho} (t \rho)}{D _ {\rho} (t \rho)} = \frac {D _ {\rho} (t (\lambda + \rho))}{D _ {\rho} (t \rho)} = \prod_ {\alpha > 0} \frac {\sinh (\alpha , (\lambda + \rho) t / 2)}{\sinh (\alpha , \rho t / 2)} \tag {13.170}
$$

which yields

$$
\dim | \lambda | = \lim  _ {t \rightarrow 0} \chi_ {\lambda} (t \rho) = \prod_ {\alpha > 0} \frac {(\lambda + \rho , \alpha)}{(\rho , \alpha)} \tag {13.171}
$$

For instance, the application of this formula to $su(2), su(3)$ , and $sp(4)$ gives

$$
s u (2): \quad \dim | \lambda | = \lambda_ {1} + 1
$$

$$
s u (3): \quad \dim | \lambda | = \frac {1}{2} \left(\lambda_ {1} + 1\right) \left(\lambda_ {2} + 1\right) \left(\lambda_ {1} + \lambda_ {2} + 2\right) \tag {13.172}
$$

$$
s p (4): \quad \dim | \lambda | = \frac {1}{6} (\lambda_ {1} + 1) (\lambda_ {2} + 1) (\lambda_ {1} + 2 \lambda_ {2} + 3) (\lambda_ {1} + \lambda_ {2} + 2)
$$

Keeping track of the subleading term in Eq. (13.170) leads to another interesting formula. At first, we have

$$
\begin{array}{l} \chi_ {\lambda} (t \rho) = \prod_ {\alpha > 0} \frac {(\lambda + \rho , \alpha)}{(\rho , \alpha)} \left\{1 + \frac {t ^ {2}}{2 4} [ (\alpha , \lambda + \rho) ^ {2} - (\alpha , \rho) ^ {2} ] \right\} \\ = \dim | \lambda | \left\{\dot {1} + \frac {t ^ {2}}{2 4} \sum_ {\alpha > 0} [ (\alpha , \lambda + \rho) ^ {2} - (\alpha , \rho) ^ {2} ] \right\} \tag {13.173} \\ \end{array}
$$

Now, as demonstrated in the following paragraph, we can always write

$$
(\lambda , \mu) = \frac {1}{y} \sum_ {\alpha \in \Delta} (\lambda , \alpha) (\alpha , \mu) = \frac {2}{y} \sum_ {\alpha \in \Delta_ {+}} (\lambda , \alpha) (\alpha , \mu) \tag {13.174}
$$

where the constant $y$ , evaluated below, depends upon the algebra. Thus we have

$$
\chi_ {\lambda} (t \rho) = \dim | \lambda | \left\{1 + \frac {t ^ {2} y}{4 8} [ | \lambda + \rho | ^ {2} - | \rho | ^ {2} ] \right\} \tag {13.175}
$$

The comparison of this expression with the $t$ expansion of

$$
\chi_ {\lambda} (t \rho) = \sum_ {\lambda^ {\prime} \in \Omega_ {\lambda}} \operatorname {m u l t} _ {\lambda} \left(\lambda^ {\prime}\right) e ^ {\left(\lambda^ {\prime}, \rho\right) t} \tag {13.176}
$$

yields

$$
\frac {1}{2} \sum_ {\lambda^ {\prime} \in \Omega_ {\lambda}} \operatorname {m u l t} _ {\lambda} \left(\lambda^ {\prime}\right) \left(\lambda^ {\prime}, \rho\right) ^ {2} = \frac {y}{4 8} \dim | \lambda | [ | \lambda + \rho | ^ {2} - | \rho | ^ {2} ] \tag {13.177}
$$

For the adjoint representation, $\lambda = \theta$ , and the different nonzero weights $\lambda'$ are the roots, which all have multiplicity 1; the l.h.s. then becomes

$$
\frac {1}{2} \sum_ {\lambda^ {\prime} \in \Omega_ {\theta}} \operatorname {m u l t} _ {\theta} \left(\lambda^ {\prime}, \rho\right) ^ {2} = \frac {1}{2} \sum_ {\alpha \in \Delta} (\rho , \alpha) (\alpha , \rho) = \frac {1}{2} y | \rho | ^ {2}
$$

where in the last step we used Eq. (13.174). The r.h.s. is

$$
\frac {y}{4 8} \dim | \theta | (\theta , \theta + 2 \rho) = \frac {y g}{2 4} \dim g \tag {13.178}
$$

(cf. Eq. (13.128)). This yields the Freudenthal-de Vries strange formula:

$$
\boxed {| \rho | ^ {2} = \frac {g}{1 2} \mathrm {d i m} g} \tag {13.179}
$$

We now return to Eq. (13.174). The product $\sum_{\alpha \in \Delta} \alpha \alpha^t$ (where $\alpha^t$ stands for the transpose of $\alpha$ ) is necessarily proportional to the $r \times r$ identity matrix $I_r$ :

$$
\sum_ {\alpha \in \Delta} \alpha \alpha^ {t} = y I _ {r} \tag {13.180}
$$

Indeed, the l.h.s. commutes with any element of the Weyl group since the action of the latter simply amounts to a permutation of the roots. Since the action of the Weyl group on $\sum_{\alpha \in \Delta} \alpha \alpha^t$ is irreducible, this latter quantity must then be proportional to the identity. The proportionality constant $y$ is evaluated by taking the trace of this equation. With $\operatorname{Tr} \alpha \alpha^t = |\alpha|^2$ , this yields

$$
\sum_ {\alpha \in \Delta} | \alpha | ^ {2} = y r \tag {13.181}
$$

The l.h.s. can be evaluated from Eq. (13.132) restricted to the generators of the Cartan subalgebra:

$$
\operatorname {T r} _ {\theta} H ^ {i} H ^ {j} = 2 g \delta^ {i j} \tag {13.182}
$$

Setting $j = i$ and summing over $i$ yields

$$
\sum_ {i} \sum_ {\alpha \in \Delta} \alpha^ {i} \alpha^ {i} = \sum_ {\alpha \in \Delta} | \alpha | ^ {2} = 2 g r \tag {13.183}
$$

This fixes the value of $y$ :

$$
y = 2 g \tag {13.184}
$$

# 13.4.3. Schur Functions

In the $su(N)$ orthogonal basis $\{\epsilon_i\}$ introduced in Sect. 13.3.2, the characters are called Schur functions. In this basis, there is a simple combinatorial formula for the dimension of a representation.

In the orthonormal basis, the Weyl group acts as the permutation group $S_N$ of the $N$ basis vectors. For instance, the action of $s_\alpha$ with $\alpha = \epsilon_i - \epsilon_j$ simply amounts to interchanging $\epsilon_i$ and $\epsilon_j$ . This observation allows us to rewrite the character as a ratio of matrix determinants, that is, as a Schur function. To this end we introduce the variables

$$
q _ {i} = e ^ {\epsilon_ {i}} \tag {13.185}
$$

subject to the constraint

$$
\prod_ {i = 1} ^ {N} q _ {i} = 1 \tag {13.186}
$$

The formal exponential can thus be written as

$$
e ^ {\lambda} = q _ {1} ^ {\ell_ {1}} q _ {2} ^ {\ell_ {2}} \dots q _ {N} ^ {\ell_ {N}} \tag {13.187}
$$

and we have

$$
D _ {\lambda} = \sum_ {w \in W} ^ {\prime} e ^ {w \lambda} = \sum_ {\sigma \in S _ {N}} \epsilon (\sigma) \prod_ {i = 1} ^ {N} q ^ {\ell_ {\sigma (i)}} = \det q _ {j} ^ {\ell_ {i}} \tag {13.188}
$$

Since the $i$ -th entry of the partition of $\rho$ is $N - i$ , the character can be written as

$$
\chi_ {\lambda} = S _ {\lambda} \left(q _ {1}, \dots , q _ {N}\right) = \frac {\det q _ {j} ^ {\ell_ {i} + N - i}}{\det q _ {j} ^ {N - i}} \tag {13.189}
$$

where $S_{\lambda}$ stands for the Schur function. The above denominator is the ubiquitous Vandermonde determinant:

$$
\det q _ {j} ^ {N - i} = \det q _ {j} ^ {i - 1} = \det \left( \begin{array}{c c c c} 1 & 1 & \dots & 1 \\ q _ {1} & q _ {2} & \dots & q _ {N} \\ \vdots & & & \\ q _ {1} ^ {N - 1} & q _ {2} ^ {N - 1} & \dots & q _ {N} ^ {N - 1} \end{array} \right) \tag {13.190}
$$

which can also be written under the more familiar form

$$
\det q _ {j} ^ {N - i} = \prod_ {1 \leq i <   j \leq N} \left(q _ {i} - q _ {j}\right) \tag {13.191}
$$

The dimension of the representation is calculated by letting all $q_{i}$ approach 1 in the expression for the character, with the result

$$
\dim | \lambda | = \prod_ {1 \leq i <   j \leq N} \frac {\left(\ell_ {i} - \ell_ {j} + j - i\right)}{(j - i)} \tag {13.192}
$$

For $su(2)$ and $su(3)$ reduced tableaux, this is easily seen to reproduce Eq. (13.172). (See also Ex. 13.13 for another dimension formula.)

# §13.5. Tensor Products: Computational Tools

In principle, the problem of calculating tensor products is straightforward. In order to calculate the product $\mathsf{L}_{\lambda} \otimes \mathsf{L}_{\mu}$ , usually written as $\lambda \otimes \mu$ , we simply add together all pairs of weights $\lambda', \mu'$ (belonging respectively to the weight systems $\Omega_{\lambda}$ and $\Omega_{\mu}$ ), taking care of their multiplicities, and reorganize the full set of $\dim |\lambda| \times \dim |\mu|$ resulting weights in irreducible representations. We write the result under the form

$$
\lambda \otimes \mu = \bigoplus_ {\nu \in P _ {+}} \mathcal {N} _ {\lambda \mu} ^ {\nu} \nu \tag {13.193}
$$

where the sum is taken over all dominant weights, and $\mathcal{N}_{\lambda \mu}^{\nu}$ , called a tensor-product coefficient, gives the multiplicity of the representation $\nu$ in the decomposition of the tensor product $\lambda \otimes \mu$ .

In practice, this method is obviously too cumbersome, and more efficient techniques are required. Such methods are described in the next subsections. The first one follows directly from manipulations of the Weyl character formula and it is completely general. Although theoretically important, as a computational tool it is not very powerful. For this reason we introduce two other methods, which, however, are presented more like recipes. These are the famous Littlewood-Richardson rule and the more novel Berenstein-Zelevinsky method of triangles. For these, however, the discussion is again restricted to $su(N)$ . Another motivation for introducing these last two methods is that they allow us to determine precisely those very states that contribute to the tensor product, a point on which we will expand in due time. Furthermore, in their framework, a particular tensor-product coefficient can be studied in isolation, that is without necessarily having to compute the full tensor-product decomposition.

But before turning to techniques, we list some general properties of tensor-product coefficients. It is clear that

$$
\mathcal {N} _ {\lambda 0} ^ {v} = \delta_ {\lambda} ^ {v} \tag {13.194}
$$

where 0 stands for the scalar representation (i.e., the representation whose highest weight has all Dynkin labels equal to zero and which thus contains a single state). On the other hand, the tensor product of a representation with its conjugate always contains the scalar representation. It is obtained from the pairing of all the states of $\lambda$ with their negatives, which necessarily lie in $\lambda^{*}$ . In other words,

$$
\mathcal {N} _ {\lambda \lambda^ {*}} ^ {0} = 1 \tag {13.195}
$$

These two relations show that lower and upper indices in $\mathcal{N}$ can be interchanged by means of the conjugate operation:

$$
\mathcal {N} _ {\lambda \mu} ^ {\nu} = \mathcal {N} _ {\lambda \nu^ {*}} ^ {\mu^ {*}} \tag {13.196}
$$

Let the coefficient $\mathcal{N}_{\lambda \mu \sigma}$ , with three lower indices, correspond to the multiplicity of the scalar representation in the triple product $\lambda \otimes \mu \otimes \sigma$ . We thus have

$$
\mathcal {N} _ {\lambda \mu} ^ {v} = \mathcal {N} _ {\lambda \mu v ^ {*}} \tag {13.197}
$$

# 13.5.1. The Character Method

The first method that will be described is based on the specification of a representation by its character. In consequence, Eq. (13.193) must also hold in character form (since the trace of a tensor product is the product of the trace)

$$
\chi_ {\lambda} \chi_ {\mu} = \sum_ {v \in P _ {+}} \mathcal {N} _ {\lambda \mu} ^ {v} \chi_ {v} \tag {13.198}
$$

Using this character equation, we can derive a simple relation between $\mathcal{N}_{\lambda \mu}^{\nu}$ and the multiplicities of the weights $\mu'$ in the representation $\mu$ , which will lead to an efficient way of calculating tensor-product coefficients. We rewrite Eq. (13.198) under the form

$$
\sum_ {w \in W} \epsilon (w) e ^ {w (\lambda + \rho)} \sum_ {\mu^ {\prime} \in \Omega_ {\mu}} \operatorname {m u l t} _ {\mu} \left(\mu^ {\prime}\right) e ^ {\mu^ {\prime}} = \sum_ {v \in P _ {+}} \mathcal {N} _ {\lambda \mu} ^ {v} \sum_ {w \in W} \epsilon (w) e ^ {w (v + \rho)} \tag {13.199}
$$

(using Eq. (13.167) for $\chi_{\lambda}$ and $\chi_{\nu}$ and Eq. (13.161) for $\chi_{\mu}$ ) and compare the contributions of both sides restricted to the fundamental chamber. Since $\nu \in P_{+}$ , only the identity element of the Weyl group contributes on the r.h.s. If $\lambda + \mu' \in P_{+}$ on the l.h.s., then again only $w = 1$ contributes. Otherwise, we first rewrite the second sum on the l.h.s. as

$$
\sum_ {\mu^ {\prime} \in \Omega_ {\mu}} \operatorname {m u l t} _ {\mu} \left(\mu^ {\prime}\right) e ^ {\mu^ {\prime}} = \sum_ {\mu^ {\prime} \in \Omega_ {\mu}} \operatorname {m u l t} _ {\mu} \left(\mu^ {\prime}\right) e ^ {w _ {\mu^ {\prime}} \mu^ {\prime}} \tag {13.200}
$$

for any $w_{\mu'} \in W$ (the multiplicity being constant along a $W$ orbit). The contributing element is the particular element of the Weyl group $w_{\mu'}$ that reflects the weight $\lambda + \mu'$ in the fundamental chamber; it contributes with the sign $\epsilon(w_{\mu'})$ . This proves the relation

$$
\mathcal {N} _ {\lambda \mu} ^ {\nu} = \sum_ {\substack {\mu^ {\prime} \in \Omega_ {\mu} \\ w \in W}} \epsilon (w) \operatorname {mult} _ {\mu} \left(\mu^ {\prime}\right) \tag{13.201}
$$

where we dropped the index $\mu'$ from $w$ for simplicity. There are two summations here: a sum over all the weights in the representation $\mu$ and a sum over those elements of the Weyl group that satisfy the condition $w \cdot (\lambda + \mu') = \nu \in P_{+}$ . The result can be rewritten more simply, with a single summation, as

$$
\mathcal {N} _ {\lambda \mu} ^ {\nu} = \sum_ {w \in W} \epsilon (w) \operatorname {m u l t} _ {\mu} (w \cdot v - \lambda) \tag {13.202}
$$

This method will be referred to as the character method. Its theoretical interest lies in its generality and in that it has a direct extension for affine fusion rules.

# 13.5.2. Algorithm for the Calculation of Tensor Products

Formula (13.202) can be translated into the following algorithm. In order to calculate the product $\lambda \otimes \mu$ , we first write down all the weights $\mu'$ in the representation $\mu$ and add each of them to $\lambda + \rho$ . Degenerate weights are treated separately. The resulting weights $\lambda + \rho + \mu'$ are of two types:

(i) those that can be reflected into dominant weights by an element $w \in W$ of the finite Weyl group;   
(ii) those in the $W$ orbit of a weight with some vanishing Dynkin labels.

Weights of type (i) contribute $\epsilon(w)$ to the tensor-product coefficient $\mathcal{N}_{\lambda \mu}^{\nu}$ , where $\nu$ is the resulting dominant weight. $\mathcal{N}_{\lambda \mu}^{\nu}$ is obtained from the sum of all these contributions.

By definition, a weight $\xi$ of type (ii) is such that there is a $w \in W$ for which $w\xi$ has at least one vanishing Dynkin label. If, for instance, $(w\xi)_i = 0$ , then $s_i(w\xi) = 0$ . Such weights can be ignored since they could be counted with both $\epsilon(w)$ and $\epsilon(s_i w) = -\epsilon(w)$ . They are located at one boundary, or a Weyl reflection thereof, of the fundamental chamber.

It should be stressed that reflecting a weight in the fundamental chamber is a finite process: at most $|\Delta_{+}|$ (the number of positive roots) reflections are needed.

Reformulated in terms of the shifted action of the Weyl group, the procedure is as follows: If $\lambda + \mu'$ can be reflected into a dominant weight by the shifted action of the Weyl group—that is, if there exists a $w \in W$ such that $w \cdot (\lambda + \mu') \in P_+$ —it contributes $\epsilon(w)$ to $\mathcal{N}_{\lambda \mu}^{\nu}$ ; if it cannot, it is ignored.

# $su(2)$ EXAMPLE

As a simple illustration of this procedure, consider the $su(2)$ tensor product $(2) \otimes (7)$ . We display on the $su(2)$ weight lattice all the weights of the representation $(7), (-7\omega_{1}, -5\omega_{1}, \dots, 7\omega_{1})$ , augmented by $2\omega_{1}$ . A shifted Weyl reflection here is a reflection with respect to the weight $-\omega_{1}$ (as $\rho = \omega_{1}$ ). The weight $-\omega_{1}$ is of type (ii) and it is thus ignored. By reflection, the nondominant weights $-5\omega_{1}$ and $-3\omega_{1}$ are sent respectively onto $3\omega_{1}$ and $\omega_{1}$ , and contribute with a minus sign, which cancels the contribution of the representations (1) and (3). This is illustrated in Fig. 13.6, from which the result of the tensor-product decomposition is directly read off:

$$
(2) \otimes (7) = (5) \oplus (7) \oplus (9) \tag {13.203}
$$

This agrees with the familiar rules of angular-momentum addition.

# $su(3)$ EXAMPLE

Consider the $su(3)$ tensor product $(1,0) \otimes (2,0)$ . The six weights in the representation $(2,0)$ are $\{(2,0), (0,1), (1,-1), (-2,2), (-1,0), (0,-2)\}$ . Adding $(1,0)$

![](images/7816ab81741b83d841da49770f772695b38c1478f81fe7647cff09a5a29b66c1.jpg)  
Figure 13.6. The $su(2)$ tensor product (2) $\otimes$ (7). The weights of the representation (7) are centered around $2\omega_{1}$ and the nondominant weights are Weyl reflected back into the dominant sector.

to each of them yields:

$$
(3, 0), (1, 1), (2, - 1), (- 1, 2), (0, 0), (1, - 2) \tag {13.204}
$$

The third and fourth of the weights (13.204) are ignored since they are respectively invariant under the shifted action of $s_2$ and $s_1$ (and are therefore of type (ii)). Acting on the sixth one with $s_2 \cdot y$ yields

$$
s _ {2} \cdot (1, - 2) = s _ {2} (2, - 1) - (1, 1) = (2, - 1) + (- 1, 2) - (1, 1) = (0, 0)
$$

Hence the reflection of the sixth weight into the fundamental chamber contributes to $\epsilon(s_2)(0,0) = -(0,0)$ , and consequently cancels the contribution of the fifth weight in (13.204). The final result is

$$
(1, 0) \otimes (2, 0) = (3, 0) \oplus (1, 1) \tag {13.205}
$$

as illustrated on Fig. 13.7.

![](images/1fb469adb28e4624d8fbd5c8e3ae1e3f082f1a745c32fa3244e462b9811071b9.jpg)  
Figure 13.7. The $su(3)$ tensor product $(1,0) \otimes (2,0)$ by the method of Weyl reflections.

In these two examples, it would have been wiser to interchange the roles of the two representations. For instance, adding $(2,0)$ to the three weights of the representation $(1,0)$ gives directly $(3,0)$ , $(1,1)$ , $(2,-1)$ , and the last one is ignored. Choosing for $\mu$ the highest weight of the smallest of the two representations simplifies the calculation in two respects: fewer states need to be considered and most of the weights in the representation $\mu$ , when added to $\lambda$ , are dominant.

# 13.5.3. The Littlewood-Richardson Rule

The Littlewood-Richardson rule is a simple and powerful algorithm, formulated in terms of the product of Young tableaux. This algorithm proceeds as follows: In the second tableau, we fill the first row with 1's, the second row with 2's, and so on. Then we add all the boxes with a 1 to the first tableau and keep only the resulting tableaux that satisfy the following two conditions:

(i) They must be regular: the number of boxes in a given row must be smaller or equal to the number of boxes in the row just above.   
(ii) They must not contain two boxes marked by 1 in the same column.

Tableaux that do not satisfy these conditions are ignored. To the resulting tableaux, we then add all the boxes marked by a 2 and again we keep only the tableaux that satisfy (i) and (ii), where in (ii), 1 is replaced by 2. We continue until all the boxes of the second tableau in the original product have been used. In this process an additional rule must be respected:

(iii) In counting from right to left and top to bottom, the number of 1's must always be greater or equal to the number of 2's, the number of 2's must always be greater or equal to the number of 3's, and so on.

The resulting Littlewood-Richardson tableaux are the Young tableaux of the irreducible representations occurring in the decomposition.

A warning: In this process, we do not construct semistandard tableaux! However, in Littlewood-Richardson tableaux it is clear that the numbers are strictly increasing in each column and they are nondecreasing in rows.

For example, consider the $su(3)$ tensor product $(2,0) \otimes (1,1)$ :

![](images/1ce9cb2fdf6c740e7fbb3c9a2e2d8d7391b0e19479c39fe8115f6ddd0c2bd312.jpg)

The tableaux obtained after the first step are

![](images/0d0a7d8fdcfe47e5260b5c0c25257f02b11b0e0afd930aced7201adf5cc51164.jpg)

Adding now the box marked by a 2 yields<sup>10</sup>

![](images/8a10153bac25d55dde58c6fafe40b103a653095285cdd2a540e1382b536f3bd5.jpg)

from which we read off

$$
(2, 0) \otimes (1, 1) = (3, 1) \oplus (1, 2) \oplus (2, 0) \oplus (0, 1) \tag {13.206}
$$

(for $su(3)$ , columns of three boxes are ignored).

The multiplicity of a given representation $\nu$ in the tensor product $\lambda \otimes \mu$ can be evaluated directly, without necessarily having to calculate the full decomposition. For this we simply add to the Young tableau representing $\lambda$ all boxes of the tableau $\mu$ such that the resulting tableau has weight $\nu$ . The added boxes are then filled

with the following set of numbers: $1 (\mu_1 + \dots + \mu_{N-1} \text{ times}), 2 (\mu_2 + \dots + \mu_{N-1} \text{ times})$ , up to $N - 1 (\mu_{N-1} \text{ times})$ , in a way that respects the Littlewood-Richardson rule. $\mathcal{N}_{\lambda \mu}^{\nu}$ is the number of distinct Littlewood-Richardson tableaux that can be produced in this way.

For instance, to the $su(4)$ tensor product $(1,2,1) \otimes (1,2,1) \supset (1,2,1)$ , there correspond 5 Littlewood-Richardson tableaux:

![](images/ea6e98c7af9d37277ef44f058cff59c16cea7bfc125622f3cc36ada485755153.jpg)

which means that the tensor-product coefficient $\mathcal{N}_{(121)(121)}^{(121)}$ is 5.

In some applications, it is necessary to know which states contribute to the tensor product. It turns out that this information is coded in the Littlewood-Richardson tableaux. More precisely, there is a one-to-one correspondence between a Littlewood-Richardson tableau associated with the product $\lambda \otimes \mu \supset \nu$ and a Gelfand-Tsetlin pattern $\{\beta_j^{(i)}\}$ of weight $\mu' = \nu - \lambda$ in the representation $\mu$ . The entries $\beta_j^{(i)}$ of the Gelfand-Tsetlin pattern can be read off the Littlewood-Richardson tableau as follows:

$$
\beta_ {j} ^ {(i)} = \text {n u m b e r o f j ' s i n t h e f i r s t i r o w s o f} \tag {13.207}
$$

the Littlewood-Richardson tableau

The states associated with each Littlewood-Richardson tableau in the previous example are

![](images/034a3ceab1fcc006328abbc9eeddd64cd93e17f75704b652bd2b40d4bfae4d19.jpg)

![](images/01ee9d9accd96346bdca566091bc4c7f929ae099b2cafddf47ca63226af5b27f.jpg)

![](images/81adeed47323292e4f9fcd1a3de544bfe80bdadc62f60a0e8fd8069b0ef595bb.jpg)

![](images/abc04498892eb3407d31ce61d3dca1e687de926351b72b0996ebd13c85ae8bd1.jpg)

$$
\begin{array}{c c c c c c}\hline&&&&1&1\\\hline&&&2&2\\\hline&1&3\\\hline 1&2\\\hline\end{array}\rightarrow\begin{array}{c c c c c c}4&3&1&0\\3&2&1\\2&2\\2&\end{array}\leftrightarrow\begin{array}{c c c c c c}\hline 1&1&3&4\\\hline 2&2&4\\\hline 3\\\hline\end{array}\tag {13.208}
$$

The weight $\mu' = \nu - \lambda = (0, 0, 0)$ in the representation $(1, 2, 1)$ has multiplicity 7. The two states that do not contribute to the tensor product are

$$
\begin{array}{c c} 4 3 1 0 & 4 3 1 0 \\ 3 2 1 & \boxed {1} \\ 3 1 & \boxed {2} \\ 2 & \boxed {3} \end{array} \quad \begin{array}{c c} 4 3 1 0 & \boxed {1} \\ 4 2 0 & \boxed {3} \\ 4 0 & \boxed {4} \end{array} \leftrightarrow
$$

For completeness, we mention that this relationship between states and Littlewood-Richardson tableaux can be used to obtain an algebraic description of the tensor-product coefficients:

$$
\mathcal {N} _ {\lambda \mu} ^ {\nu} = \text {n u m b e r o f G e l f a n d - T s e t l i n p a t t e r n s} \left\{\beta_ {j} ^ {(i)} \right\} \tag {13.209}
$$

of weight $\mu^{\prime} = \nu -\lambda$ in the representation $\pmb{\mu}$ that satisfy the conditions $d_j^{(i)}\leq \lambda_i$ for all values of $j,1\leq j\leq i\leq N - 1$

where

$$
d _ {j} ^ {(i)} = \sum_ {1 \leq n <   j} \left(\beta_ {n} ^ {(i + 1)} - 2 \beta_ {n} ^ {(i)} + \beta_ {n} ^ {(i - 1)}\right) + \left(\beta_ {j} ^ {(i + 1)} - \beta_ {j} ^ {(i)}\right) \tag {13.210}
$$

For the first noncontributing pattern of the previous example: $d_2^{(3)} = 2 > \lambda_3 = 1$ , and for the other one: $d_1^{(1)} = 2 > \lambda_1 = 1$ .

# 13.5.4. Berenstein-Zelevinsky Triangles

Berenstein-Zelevinsky triangles (BZ) provide a powerful way to calculate the multiplicity of a triple product, that is, the multiplicity of the scalar representation in $\lambda \otimes \mu \otimes \nu$ . (We point out the slight change in the notation for the third weight: we take it to be $\nu$ instead of $\nu^{*}$ .) They also contain information on the states contributing to the product. We first describe the construction for $su(3)$ .

We consider the set of three $su(3)$ highest weights $(\lambda_1, \lambda_2)$ , $(\mu_1, \mu_2)$ , and $(\nu_1, \nu_2)$ . We construct triangles according to the following rules:

$$
m _ {1 3}
$$

$$
\begin{array}{c c} n _ {1 2} & l _ {2 3} \\ 3 & m _ {1 2} \end{array} \tag {13.211}
$$

$$
\begin{array}{c c c c} n _ {1 3} & l _ {1 2} & n _ {2 3} & l _ {1 3} \end{array}
$$

where the nine nonnegative integers $l_{ij}, m_{ij}, n_{ij}$ are related to the Dynkin labels of the three integrable weights by

$$
m _ {1 3} + n _ {1 2} = \lambda_ {1} \quad n _ {1 3} + l _ {1 2} = \mu_ {1} \quad l _ {1 3} + m _ {1 2} = v _ {1} \tag {13.212}
$$

$$
m _ {2 3} + n _ {1 3} = \lambda_ {2} \quad n _ {2 3} + l _ {1 3} = \mu_ {2} \quad l _ {2 3} + m _ {1 3} = \nu_ {2}
$$

They must further satisfy the so-called hexagon conditions

$$
\begin{array}{l} n _ {1 2} + m _ {2 3} = n _ {2 3} + m _ {1 2} \\ l _ {1 2} + m _ {2 3} = l _ {2 3} + m _ {1 2} \tag {13.213} \\ l _ {1 2} + n _ {2 3} = l _ {2 3} + n _ {1 2} \\ \end{array}
$$

This means that the length of opposite sides in the hexagon formed by $n_{12}, l_{23}, m_{12}, n_{23}, l_{12}$ and $m_{23}$ in (13.211) are equal, the length of a segment being defined as the sum of its two vertices.

The number of such triangles gives the value of $\mathcal{N}_{\lambda \mu \nu}$ . If it is not possible to construct such a triangle, it means that $\nu^{*}$ does not occur in the tensor product $\lambda \otimes \mu$ .

The integers in the BZ triangles have the following origin. Each pair of indices $ij$ , $i < j$ , on the labels of the triangle is related to a positive root of $su(3)$ . We recall that the positive roots of $su(N)$ can be written as $\epsilon_i - \epsilon_j$ , $1 \leq i < j \leq N$ in terms of orthonormal vectors $\epsilon_i$ in $\mathbb{R}^N$ . The triangle encodes three sums of positive roots:

$$
\begin{array}{l} \mu^ {\prime} + v - \lambda^ {*} = \sum_ {i <   j} l _ {i j} \left(\epsilon_ {i} - \epsilon_ {j}\right) \\ \nu + \lambda - \mu^ {*} = \sum_ {i <   j} m _ {i j} \left(\epsilon_ {i} - \epsilon_ {j}\right) \tag {13.214} \\ \lambda + \mu - v ^ {*} = \sum_ {i <   j} n _ {i j} \left(\epsilon_ {i} - \epsilon_ {j}\right) \\ \end{array}
$$

The hexagon relations (13.213) can be seen as consistency conditions for these three expansions.

The four triangles for the example (13.206) are

$$
\begin{array}{c c c c c c c c c} 2 & & 1 & & 1 & & 0 \\ 0 & 1 & & 1 & 0 & & 2 & 0 \\ 0 & 0 & & 0 & 1 & & 0 & 1 \\ 0 & 1 & 0 & 1 & 0 & 1 & 0 & 0 \end{array} \tag {13.215}
$$

On the other hand, corresponding to the coupling $(2,2) \otimes (2,2) \otimes (2,2)$ , three triangles can be constructed:

$$
\begin{array}{r r r r r r r r} & 0 & & & 1 & & & 2 \\ & 2 & 2 & & & 1 & 1 & \\ 2 & & 2 & & & 1 & & 0 \\ 0 & 2 & 2 & 0 & & 1 & 1 & 1 \end{array} \quad \begin{array}{r r r r r r r r} & 0 & 0 \\ & 0 & 0 \\ 2 & 0 & 0 & 2 \end{array} \tag {13.216}
$$

and accordingly the multiplicity of the scalar representation in this triple product is 3.

The states involved in a specific coupling can be read off a triangle as follows. Consider the product $\lambda \otimes \mu \supset \nu^{*}$ associated with the triangle (13.211). The state

of weight $\mu' = \nu^* - \lambda$ in this coupling is described by the Gelfand-Tsetlin pattern

$$
\mu_ {1} + \mu_ {2} \quad \mu_ {2} \quad 0
$$

$$
\mu_ {1} + \mu_ {2} - n _ {1 3} \quad \mu_ {2} - n _ {2 3} \tag {13.217}
$$

$$
\mu_ {1} + \mu_ {2} - n _ {1 3} - n _ {1 2}
$$

For example, the Gelfand-Tsetlin patterns and corresponding semistandard tableaux in the representation $\mu = (2,2)$ associated with the three triangles of the last example (13.216) are (in the same order)

$$
\begin{array}{c c} 4 2 0 & \\ 4 0 & \leftrightarrow \begin{array}{c c c c} \hline 1 & 1 & 2 & 2 \\ \hline 3 & 3 \end{array} \\ 2 & \end{array} \qquad \begin{array}{c c c c} 4 2 0 & \\ 3 & 1 & \leftrightarrow \begin{array}{c c c c} \hline 1 & 1 & 2 & 3 \\ \hline 2 & 3 \end{array} \\ 2 & \end{array} \qquad \begin{array}{c c c c} 4 2 0 & \\ 2 & \leftrightarrow \begin{array}{c c c c} \hline 1 & 1 & 3 & 3 \\ \hline 2 & 2 \end{array} \\ 2 & \end{array}
$$

For $su(4)$ , the BZ triangles are defined in a similar way, in terms of eighteen nonnegative integers:

$$
\begin{array}{c c c c c} & m _ {1 4} \\ & n _ {1 2} & l _ {3 4} \\ & m _ {2 4} & m _ {1 3} \\ & n _ {1 3} & l _ {2 3} & n _ {2 3} & l _ {2 4} \\ m _ {3 4} & m _ {2 3} & m _ {1 2} \\ n _ {1 4} & l _ {1 2} & n _ {2 4} & l _ {1 3} & n _ {3 4} \end{array} \tag {13.218}
$$

related to the Dynkin labels by

$$
m _ {1 4} + n _ {1 2} = \lambda_ {1} \quad n _ {1 4} + l _ {1 2} = \mu_ {1} \quad l _ {1 4} + m _ {1 2} = \nu_ {1}
$$

$$
m _ {2 4} + n _ {1 3} = \lambda_ {2} \quad n _ {2 4} + l _ {1 3} = \mu_ {2} \quad l _ {2 4} + m _ {1 3} = v _ {2} \tag {13.219}
$$

$$
m _ {3 4} + n _ {1 4} = \lambda_ {3} \quad n _ {3 4} + l _ {1 4} = \mu_ {3} \quad l _ {3 4} + m _ {1 4} = \nu_ {3}
$$

Furthermore, a $su(4)$ BZ triangle has 3 hexagons:

$$
n _ {1 2} + m _ {2 4} = m _ {1 3} + n _ {2 3} \quad n _ {1 3} + l _ {2 3} = l _ {1 2} + n _ {2 4} \quad l _ {2 4} + n _ {2 3} = l _ {1 3} + n _ {3 4}
$$

$$
n _ {1 2} + l _ {3 4} = l _ {2 3} + n _ {2 3} \quad n _ {1 3} + m _ {3 4} = n _ {2 4} + m _ {2 3} \quad n _ {2 3} + m _ {2 3} = m _ {1 2} + n _ {3 4}
$$

$$
m _ {2 4} + l _ {2 3} = l _ {3 4} + m _ {1 3} \quad m _ {3 4} + l _ {1 2} = l _ {2 3} + m _ {2 3} \quad l _ {1 3} + m _ {2 3} = l _ {2 4} + m _ {1 2} \tag {13.220}
$$

The $su(N)$ generalization is straightforward; the triangles are built out of $(N - 1)(N - 2)/2$ hexagons and three corner points. On the other hand, for $su(2)$ there are no hexagons: the tensor products are described by the simple triangles

$$
\begin{array}{c} m _ {1 2} \\ n _ {1 2} \end{array} l _ {1 2} \tag {13.221}
$$

written in terms of three nonnegative integers constrained by

$$
m _ {1 2} + n _ {1 2} = \lambda_ {1}
$$

$$
n _ {1 2} + l _ {1 2} = \mu_ {1} \tag {13.222}
$$

$$
l _ {1 2} + m _ {1 2} = v _ {1}
$$

With $\lambda_{1}$ and $\mu_{1}$ fixed, $\nu_{1}$ satisfies

$$
v _ {1} = \lambda_ {1} + \mu_ {1} - 2 n _ {1 2} \tag {13.223}
$$

which reproduces the rule for $su(2)$ tensor products in a very simple way.

With the last two methods described, it is possible to study a particular triple product in isolation, that is, without necessarily computing the full product $\lambda \otimes \mu$ . This is a clear advantage when reasonably large representations are involved. The BZ triangles have the further advantage of preserving most of the symmetries of the tensor-product coefficients. In fact, the only symmetry that is not manifest is $\mathcal{N}_{\lambda \mu \nu} = \mathcal{N}_{\mu \lambda \nu}$ .

We note finally that, in contradistinction with the Littlewood-Richardson rule, the generalization of the BZ triangles to $so(N)$ and $sp(N)$ is unknown at this time.

# §13.6. Tensor Products: A Fusion-Rule Point of View

In this section we discuss tensor products from a point of view close in spirit to the approach used in fusion-rule calculations. At first, we indicate how generic tensor-product coefficients are fixed by associativity in terms of tensor-product coefficients involving the fundamental representations. We recall that for minimal models, the fusion ring was found to be generated by $\phi_{(1,2)}$ and $\phi_{(2,1)}$ . In that context, Chebyshev polynomials appeared naturally. These polynomials and their generalizations resurface here. In a second step, we derive the Lie algebra version of the Verlinde formula (10.201).

The associativity of tensor products translates into the following condition

$$
\sum_ {\sigma} \mathcal {N} _ {\lambda \mu} ^ {\sigma} \mathcal {N} _ {\sigma v \xi} = \sum_ {\zeta} \mathcal {N} _ {\mu v} ^ {\zeta} \mathcal {N} _ {\zeta \lambda \xi} \tag {13.224}
$$

It is clear that if $\lambda$ and $\mu$ are fundamental representations, any general coefficient $\mathcal{N}_{\sigma \nu \xi}$ can be deduced from this condition whenever all tensor-product coefficients involving at least one fundamental representation are known. Again, by introducing a matrix $N_{\lambda}$ with entries

$$
\left(N _ {\lambda}\right) _ {\mu} ^ {\sigma} = \mathcal {N} _ {\lambda \mu} ^ {\sigma} \tag {13.225}
$$

we see that the associativity requirement boils down to the commutativity of the matrices $N$ :

$$
\left(N _ {\lambda} N _ {v}\right) _ {\mu \xi} = \left(N _ {v} N _ {\lambda}\right) _ {\mu \xi} \tag {13.226}
$$

These matrices provide a representation of the tensor-product algebra:

$$
N _ {\lambda} N _ {v} = \sum_ {\sigma \in P _ {+}} \mathcal {N} _ {\lambda v} ^ {\sigma} N _ {\sigma} \tag {13.227}
$$

Here we did nothing but rewrite (13.224) in matrix form. We note that these matrices are infinite.

We will now see how Chebyshev-like polynomials arise in this picture. We consider first the $su(2)$ case. From the Littlewood-Richardson rule (or the angular-momentum addition theory), we easily see that

$$
(1) \otimes (n) = (n + 1) \oplus (n - 1) \tag {13.228}
$$

where it is understood that if the Dynkin label $n - 1$ is negative, the second representation on the r.h.s. is omitted. The comparison of this product rule with Eq. (13.227) shows that the matrix $N_{1}$ is simply:

$$
\left(N _ {1}\right) _ {j} ^ {k} = \delta_ {j, k + 1} + \delta_ {j, k - 1} \tag {13.229}
$$

Mutually commuting matrices associated with other representations can be constructed as follows. One first observes that Eq. (13.228) translates into the following relation

$$
N _ {1} N _ {n} = N _ {n + 1} + N _ {n - 1} \tag {13.230}
$$

which can be regarded as a recurrence relation to be solved for $N_{n}$ in terms of $N_{1}$ . This becomes clearer if we replace $N_{1}$ by $x$ , $N_{n}$ by $U_{n}(x)$ and rewrite the above equation in the form

$$
x U _ {n} = U _ {n + 1} + U _ {n - 1} \tag {13.231}
$$

With $U_0 = 1$ , $U_1 = x$ , this is the defining relation for Chebyshev polynomials of the second kind, which already arose in the context of minimal-model fusion rules (cf. Eq. (8.101)). The desired expression for $N_n$ is thus

$$
N _ {n} = U _ {n} \left(N _ {1}\right) \tag {13.232}
$$

For instance, given the matrix $N_{1}$ ——which fixes all tensor products with the fundamental representation—the matrix $N_{2}$ describing the products with the adjoint representation is

$$
N _ {2} = N _ {1} ^ {2} - 1 \tag {13.233}
$$

That is,

$$
N _ {1} = \left( \begin{array}{c c c c c c} 0 & 1 & 0 & 0 & 0 & \dots \\ 1 & 0 & 1 & 0 & 0 & \dots \\ 0 & 1 & 0 & 1 & 0 & \dots \\ 0 & 0 & 1 & 0 & 1 & \dots \\ \dots & \dots & \dots \end{array} \right) \Longrightarrow N _ {2} = \left( \begin{array}{c c c c c c} 0 & 0 & 1 & 0 & 0 & \dots \\ 0 & 1 & 0 & 1 & 0 & \dots \\ 1 & 0 & 1 & 0 & 1 & \dots \\ 0 & 1 & 0 & 1 & 0 & \dots \\ \dots & \dots & \dots \end{array} \right) \tag {13.234}
$$

from which we read off directly that

$$
\begin{array}{l} (2) \otimes (0) = (2) \\ (2) \otimes (1) = (1) \oplus (3) \tag {13.235} \\ (2) \otimes (2) = (0) \oplus (2) \oplus (4) \\ \end{array}
$$

and so on. Because they are all constructed out of polynomials in $N_{1}$ , the matrices $N_{n}$ necessarily commute among themselves.

It is interesting to construct the generating function of the Chebyshev polynomials. This is done in the standard way: one multiplies Eq. (13.231) by $t^n$ and sums the result from $n = 0$ to $n = \infty$ ; by simple manipulations, each term can be reexpressed in terms of

$$
F (x; t) = \sum_ {n = 0} ^ {\infty} U _ {n} t ^ {n} \tag {13.236}
$$

with the result

$$
x F = (F - 1) / t + t F \tag {13.237}
$$

that is,

$$
F (x; t) = \frac {1}{1 - x t + t ^ {2}} \tag {13.238}
$$

A similar analysis can be done for any Lie algebra. For instance, for $su(3)$ , the Littlewood-Richardson rule immediately tells us that

$$
\begin{array}{l} (1, 0) \otimes (\lambda_ {1}, \lambda_ {2}) = (\lambda_ {1} + 1, \lambda_ {2}) \oplus (\lambda_ {1}, \lambda_ {2} - 1) \oplus (\lambda_ {1} - 1, \lambda_ {2} + 1) \\ (0, 1) \otimes (\lambda_ {1}, \lambda_ {2}) = (\lambda_ {1}, \lambda_ {2} + 1) \oplus (\lambda_ {1} - 1, \lambda_ {2}) \oplus (\lambda_ {1} + 1, \lambda_ {2} - 1) \tag {13.239} \\ \end{array}
$$

Again, we can replace the representations in these expressions by their corresponding tensor-product matrices $N_{(\lambda_1,\lambda_2)}$ . These matrices turn out to be expressible in terms of some generalized Chebyshev polynomials $U_{(\lambda_1,\lambda_2)}$ , a function of two variables $x_1, x_2$ associated respectively with $N_{(1,0)}$ and $N_{(0,1)}$ , as follows:

$$
N _ {(\lambda_ {1}, \lambda_ {2})} = U _ {(\lambda_ {1}, \lambda_ {2})} \left(N _ {(1, 0)}, N _ {(0, 1)}\right) \equiv U _ {(\lambda_ {1}, \lambda_ {2})} \left(x _ {1}, x _ {2}\right) \tag {13.240}
$$

These polynomials are defined in terms of the generating function

$$
\begin{array}{l} F \left(x _ {1}, x _ {2}; t, s\right) = \sum_ {\lambda_ {1}, \lambda_ {2} = 0} ^ {\infty} U _ {\left(\lambda_ {1}, \lambda_ {2}\right)} \left(x _ {1}, x _ {2}\right) t ^ {\lambda_ {1}} s ^ {\lambda_ {2}} \tag {13.241} \\ = \frac {1 - t s}{\left(1 - t x _ {1} + t ^ {2} x _ {2} - t ^ {3}\right) \left(1 - s x _ {2} + s ^ {2} x _ {1} - s ^ {3}\right)} \\ \end{array}
$$

The details of this analysis are left to the reader (cf. Ex. 13.20).

We now turn to a Lie algebra version of the Verlinde formula (10.201). The starting point is the character product form (13.198), in which all the characters are supposed to be evaluated at the particular point

$$
X = - 2 \pi i \sum_ {i = 1} ^ {r} t _ {i} \alpha_ {i} ^ {\vee} \tag {13.242}
$$

where the $t_i$ 's are real numbers valued in the range [0, 1]. With $\chi_{\lambda} = D_{\lambda + \rho} / D_{\rho}$ , this becomes

$$
\frac {D _ {\lambda + \rho} (X) D _ {\mu + \rho} (X)}{D _ {\rho} (X)} = \sum_ {v} \mathcal {N} _ {\lambda \mu} ^ {v} D _ {v + \rho} (X) \tag {13.243}
$$

The $D_{\lambda + \rho}(X)$ 's satisfy the following orthogonality relation

$$
\int_ {0} ^ {1} \left(\prod_ {i} d t _ {i}\right) D _ {\lambda + \rho} (X) D _ {\mu + \rho} (X) = | W | \delta_ {\mu , \lambda}. \tag {13.244}
$$

where $|W|$ is the order of the Weyl group. This follows directly from the definition (13.165) and the expression (13.117) for the conjugate of a representation, which implies that $D_{\lambda^{\bullet} + \rho}(X)$ is the complex conjugate of $D_{\lambda + \rho}(X)$ . To proceed, we multiply Eq. (13.243) by $D_{\sigma^{\bullet} + \rho}(X)$ and integrate the result over $X$ (i.e., integrate over all $t_i$ 's from 0 to 1). This gives the desired Verlinde-type formula

$$
\mathcal {N} _ {\lambda \mu} ^ {\sigma} = \int_ {0} ^ {1} \left(\prod_ {i} d t _ {i}\right) \frac {S _ {\lambda} (X) S _ {\mu} (X) \bar {S} _ {\sigma} (X)}{S _ {0} (X)} \tag {13.245}
$$

where

$$
\mathcal {S} _ {\lambda} (X) = \frac {1}{\sqrt {| W |}} D _ {\lambda + \rho} (X) \tag {13.246}
$$

Such an $S$ matrix is analogous to the one obtained for a finite group in Ex. 10.18: it is indexed by $r$ discrete numbers, the Dynkin labels $\lambda_i$ , and $r$ continuous ones, the $t_i$ . This immediately tells us that such an $S$ matrix cannot be a transformation matrix of characters into themselves, like the modular transformation matrix. For $su(2)$ , it takes the simple form

$$
S _ {\lambda_ {1}} \left(t _ {1}\right) = - i \sqrt {2} \sin \left[ 2 \pi t _ {1} \left(\lambda_ {1} + 1\right) \right] \tag {13.247}
$$

# §13.7. Algebra Embeddings and Branching Rules

As mentioned in the introduction, we will often encounter "affine" generalizations of simple Lie algebra embeddings. This fact motivates the general remarks of this section.

# 13.7.1. Embedding Index

We first present different ways of characterizing an embedding $\mathbf{p} \subset \mathbf{g}$ , deferring classification issues to the next subsection.

# i) Branching rules:

Viewed from the standpoint of the smaller algebra $\mathfrak{p}$ , an irreducible representation of $\mathfrak{g}$ usually breaks down into many irreducible representations of $\mathfrak{p}$ . Such decompositions are called branching rules and are noted as

$$
\mathsf {L} _ {\lambda} \mapsto \bigoplus_ {\mu \in P _ {+}} b _ {\lambda \mu} \mathsf {L} _ {\mu} \tag {13.248}
$$

or simply as

$$
\lambda \mapsto \bigoplus_ {\mu \in P _ {+}} b _ {\lambda \mu} \mu \tag {13.249}
$$

The branching coefficient $b_{\lambda \mu}$ gives the multiplicity of the irreducible representation $\mu$ of $\mathfrak{p}$ in the decomposition of the irreducible representation $\lambda$ of $\mathfrak{g}$ . The decomposition of the lowest-dimensional nontrivial representation is sufficient to characterize an embedding. To each of its inequivalent branching rules corresponds a distinct embedding.

ii) Projection matrix:

A projection matrix $\mathcal{P}$ gives the explicit projection of every weight of $\mathbf{g}$ onto a weight of $\mathbf{p}$ . Hence, to calculate the branching rules one first projects all the weights of a given irreducible representation of $\mathbf{g}$ into $\mathbf{p}$ -weights and reorganizes them into irreducible representations. Projection matrices are not unique: a Weyl reflection of the root diagram modifies them without affecting the embedding.

iii) Embedding index:

The embedding index $x_{e}$ is defined as the ratio of the square length of the projection of $\theta$ , the highest root of $\mathbf{g}$ , to the square length of the highest root of $\mathbf{p}$ , which is denoted by $\vartheta$ :

$$
\boxed {x _ {e} = \frac {| \mathcal {P} \theta | ^ {2}}{| \vartheta | ^ {2}}} \tag {13.250}
$$

Given a branching rule, the embedding index can also be calculated from

$$
\boxed {x _ {e} = \sum_ {\mu \in P _ {+}} b _ {\lambda \mu} \frac {x _ {\mu}}{x _ {\lambda}}} \tag {13.251}
$$

where $x_{\lambda}$ is the index of the representation $\lambda$ of $\mathbf{g}$ defined in Eq. (13.133). The proof of this relation is left as an exercise (cf. Ex. 13.22).

![](images/82a595a5b85d1cacb250ee2e15bf07925030626cd32667bbea7398dd929d33ff.jpg)  
Figure 13.8. Projection of the $su(3)$ adjoint representation onto $su(2)$ .

As an example, we show how $su(2)$ can be embedded into $su(3)$ . Fig. 13.8 shows how the $su(3)$ root system is projected along the highest root vector and gives a possible assignment of the $su(2)$ weights. The representation $(1,1)$ of $su(3)$ decomposes into the $su(2)$ representations $(2) \oplus (4)$ (of respective dimension 3 and 5); it is thus characterized by the branching rule:

$$
(1, 1) \mapsto (4) \oplus (2) \tag {13.252}
$$

The embedding index is easily found by noticing that the highest weight of $su(3)$ , $\alpha_{1} + \alpha_{2}$ , is projected onto $2\alpha_{1}$ . The ratio of highest roots is thus 4, and $x_{e} = 4$ . This can also be seen from Eq. (13.251) using the branching rule (13.252). The required representation indices are:

$$
s u (3): x _ {(1, 1)} = 3 \tag {13.253}
$$

$$
s u (2): x _ {(4)} = 1 0, \quad x _ {(2)} = 2
$$

Their substitution in Eq. (13.251) reproduces the value $x_{e} = 4$ . The projection matrix for this embedding can be chosen as

$$
\mathcal {P} _ {(4)} = (2, 2) \tag {13.254}
$$

if the $su(3)$ weight is written in a column matrix whose entries are its Dynkin labels. Hence, the $su(3)$ weight $(\lambda_1,\lambda_2)$ is projected into the $su(2)$ weight of Dynkin label $2\lambda_{1} + 2\lambda_{2}$ (which is thus always even):

$$
(2, 2) \binom {\lambda_ {1}} {\lambda_ {2}} = (2 \lambda_ {1} + 2 \lambda_ {2}) \tag {13.255}
$$

Using this matrix, the branching rules $(1,0) \mapsto (2)$ , $(0,1) \mapsto (2)$ are easily derived.

Dividing all $su(2)$ Dynkin labels by 2 in Fig 13.8 leads to another possible assignment for the $su(2)$ weights. The branching rule specifying this embedding is

$$
(1, 1) \mapsto (2) \oplus 2 (1) \oplus (0) \tag {13.256}
$$

Because $\alpha_{1} + \alpha_{2}$ is projected onto $\alpha_{1}$ , the embedding index is equal to 1. A candidate projection matrix for this embedding is

$$
\mathcal {P} _ {(1)} = (1, 1) \tag {13.257}
$$

The basic branching rule is $(1,0)\mapsto (1) + (0)$ .

These two embeddings are most conveniently described by means of the following generating functions:

$$
F _ {(1)} = \frac {1}{(1 - L _ {1} M) (1 - L _ {2} M) (1 - L _ {1}) (1 - L _ {2})} \tag {13.258}
$$

$$
F _ {(4)} = \frac {\left(1 + L _ {1} L _ {2} M ^ {2}\right)}{\left(1 - L _ {1} M ^ {2}\right) \left(1 - L _ {2} M ^ {2}\right) \left(1 - L _ {1} ^ {2}\right) \left(1 - L _ {2} ^ {2}\right)}
$$

where the subscript indicates the embedding index. To obtain the decomposition of the $su(3)$ weight $(\lambda_1, \lambda_2)$ , we expand $F$ and collect all the terms multiplying $L_1^{\lambda_1} L_2^{\lambda_2}$ ; its coefficient, of the form $aM^m + bM^n + \dots$ , codes the decomposition of $(\lambda_1, \lambda_2)$ :

$$
\left(\lambda_ {1}, \lambda_ {2}\right) \mapsto a (m) \oplus b (n) \oplus \dots \tag {13.259}
$$

We take for instance the embedding with $x_{e} = 4$ . In the power expansion of $F_{(4)}$ , the term $L_1^3 L_2^0$ is multiplied by $M^6 + M^2$ , so that

$$
(3, 0) \mapsto (6) \oplus (2) \tag {13.260}
$$

A few remarks complete this discussion. An obvious necessary condition for the branching coefficient $b_{\lambda \mu}$ to be nonzero is

$$
\mathcal {P} \lambda - \mu \in \mathcal {P} Q \tag {13.261}
$$

where $Q$ is the root lattice of $g$ . This simply means that the integrable weight $\mu$ must lie somewhere in the integrable representation $\lambda$ , after projection; since any weight in $\Omega_{\lambda}$ can be obtained from $\lambda$ by subtracting an appropriate number of positive roots, the condition follows. In the examples above, the root lattice projects as follows. For the first embedding,

$$
x _ {e} = 4: \quad \mathcal {P} Q _ {s u (3)} = Q _ {s u (2)} \tag {13.262}
$$

since the $su(3)$ roots $\alpha_{1}^{\prime}$ and $\alpha_{2}$ are projected onto the $su(2)$ weight $2\omega_{1}$ , that is, onto the $su(2)$ simple root $\alpha_{1}$ . The condition (13.261) forces then

$$
\left(\mathcal {P} \lambda\right) _ {1} = \mu_ {1} \bmod 2 \tag {13.263}
$$

where $(\mathcal{P}\lambda)_1$ is the Dynkin label of the projected $su(3)$ weight. For the other embedding, we find that both $\alpha_{1}$ and $\alpha_{2}$ are mapped onto $\omega_{1}$ of $su(2)$ , so that

$$
x _ {e} = 1: \quad \mathcal {P} Q _ {s u (3)} = P _ {s u (2)} \tag {13.264}
$$

where, as usual, (noncalligraphic) $P$ stands for the weight lattice. As a result Eq. (13.261) gives no constraint: both $\lambda$ and $\mu$ are integrable weights.

We note finally that a useful tool for the computation of branching rules uses tensor products. If

$$
\lambda \mapsto \bigoplus_ {\mu} b _ {\lambda \mu} \mu \quad \text {a n d} \quad \xi \mapsto \bigoplus_ {\nu} b _ {\xi \nu} \nu \tag {13.265}
$$

then

$$
\lambda \otimes \xi \mapsto \bigoplus_ {\mu , v} b _ {\lambda \mu} b _ {\xi v} \mu \otimes v \tag {13.266}
$$

For instance, given the branching rules $(1,0) \mapsto (1) \oplus (0)$ and $(0,1) \mapsto (1) \oplus (0)$ for $su(2) \subset su(3)$ with $x_{e} = 1$ , we find the branching rule for $(2,0)$ from

$$
(1, 0) \otimes (1, 0) = (2, 0) \oplus (0, 1) \mapsto [ (1) \oplus (0) ] \otimes [ (1) \oplus (0) ] = (2) \oplus 2 (1) \oplus 2 (0) \tag {13.267}
$$

that is, $(2,0)\mapsto (2)\oplus (1)\oplus (0)$

# 13.7.2. Classification of Embeddings

We now briefly address the question of classifying the possible embeddings. The following discussion is restricted to maximal embeddings; these are embeddings

$\mathfrak{p} \subset \mathfrak{g}$ for which there is no $\mathfrak{p}'$ such that $\mathfrak{p} \subset \mathfrak{p}' \subset \mathfrak{g}$ . All nonmaximal embeddings can be obtained from a chain of maximal ones. We also suppose that $\mathfrak{g}$ is semisimple; that also makes $\mathfrak{p}$ semisimple up to a possible $u(1)$ factor.

The simplest embeddings are those for which there exists a basis of $\mathbf{g}$ in which a subset of generators form the generators of $\mathfrak{p}$ . In other words, if the $\mathfrak{p}$ generators are denoted by a tilde, we have $\{\tilde{E}^{\alpha}\} \subset \{E^{\alpha}\}$ and $\{\tilde{H}^i\} \subset \{H^i\}$ . These are called the regular subalgebras. The maximal regular subalgebras have the same rank as the algebra $\mathbf{g}$ and they are easily described in terms of the root system of $\mathbf{g}$ .

We first construct the extended Dynkin diagram of $\mathfrak{g}$ by adding an extra node, associated with $-\theta$ . The extended Dynkin diagrams are displayed in Fig. 14.1 of Chap. 14. Promoting $-\theta$ to a "simple root" preserves the characteristic property that the difference between two simple roots is not a root (i.e., $\alpha_{i} + \theta$ cannot be a root since $\theta$ is the highest root). However in order to restore the linear independence of the simple roots, at least one $\alpha_{i}$ has to be removed from this augmented set of simple roots. All semisimple maximal regular subalgebras are obtained by removing from the extended Dynkin diagram of $\mathfrak{g}$ any node whose mark is a prime number.[12] Maximal regular algebras that are not semisimple are constructed from the removal of two nodes with mark $a_{i} = 1$ and the addition of a $u(1)$ factor. The embedding $su(2) \subset su(3)$ with $x_{e} = 1$ is a regular embedding because it can be obtained from $su(3)$ by dropping one simple root (hence two simple roots with unit mark from the extended $su(3)$ diagram). However, it is not a maximal embedding, being associated with the regular chain $su(2) \subset su(2) \oplus u(1) \subset su(3)$ . As another example, consider the extended Dynkin diagram of $E_{8}$ out of which one of the simple roots $\{\alpha_{1}, \alpha_{2}, \alpha_{4}, \alpha_{7}, \alpha_{8}\}$ is removed. The resulting algebras in each case are, respectively, $su(2) \oplus E_{7}$ , $su(3) \oplus E_{6}$ , $su(4) \oplus so(11)$ , $so(16)$ , and $su(9)$ . Since $E_{8}$ has no simple root with unit mark, its maximal regular subalgebras are all semisimple.

The calculation of branching rules in a regular embedding proceeds as follows. We first add to all the weights in the representation $\mathsf{L}_{\lambda}$ an extra Dynkin label, associated with the extra simple root $-\theta$ . Since the decomposition of $\theta$ in terms of the simple coroots is known (the expansion coefficients being the comarks), this extra Dynkin label is simply

$$
\lambda_ {- \theta} = - \sum_ {i} a _ {i} ^ {\vee} \lambda_ {i} \tag {13.268}
$$

$$
\widehat {F} _ {4} / \alpha_ {3} \subset \widehat {F} _ {4} / \alpha_ {4}, \quad \widehat {E} _ {7} / \alpha_ {3} \subset \widehat {E} _ {7} / \alpha_ {1}, \quad \widehat {E} _ {8} / \alpha_ {6} \subset \widehat {E} _ {8} / \alpha_ {1}
$$

$$
\widehat {E} _ {8} / \alpha_ {5} \subset \widehat {E} _ {8} / \alpha_ {2}, \quad \widehat {E} _ {8} / \alpha_ {3} \subset \widehat {E} _ {8} / \alpha_ {7}
$$

If the regular subalgebra $\mathfrak{p}$ is obtained by deleting the simple root $\alpha_{j}$ , we simply delete the Dynkin label $\lambda_{j}$ from all the weights. The resulting weights are exactly the projected weights, and they can be reorganized into irreducible representations of $\mathfrak{p}$ . This is illustrated in Ex. 13.25. The same procedure works for the semisimple algebra obtained from the removal of two nodes.

Table 13.1. Maximal semisimple special algebras of exceptional Lie algebras. The upper index gives the value of the embedding index.   

<table><tr><td>Exceptional g</td><td>Maximal special p</td></tr><tr><td>G2</td><td>A1(28)</td></tr><tr><td>F4</td><td>A1(156), G2(1) ⊕ A1(8)</td></tr><tr><td>E6</td><td>A1(9), G2(3), C4(1), G2(1) ⊕ A2(2), F4(1)</td></tr><tr><td>E7</td><td>A1(399), A1(231), A2(21), G2(1) ⊕ C3(1), F4(1) ⊕ A2(3), G2(2) ⊕ A1(7), A1(24) ⊕ A1(15)</td></tr><tr><td>E8</td><td>A1(1240), A1(760), A1(520), G2(1) ⊕ F4(1), A2(6) ⊕ A1(16), B2(12)</td></tr></table>

Nonregular subalgebras are called special subalgebras. There is still a general method for obtaining the special embeddings of the classical algebras, but the exceptional ones require a case-by-case analysis, whose result is given in Table 13.1. For the classical algebras, we use the realization of the corresponding compact groups as a group of matrices. This makes the following embeddings almost immediate:

$$
s u (p) \oplus s u (q) \subset s u (p q)
$$

$$
s o (p) \oplus s o (q) \subset s o (p q)
$$

$$
s p (2 p) \oplus s p (2 q) \subset s o (4 p q) \tag {13.269}
$$

$$
s p (2 p) \oplus s o (q) \subset s p (2 p q)
$$

$$
s o (p) \oplus s o (q) \subset s o (p + q) \quad \text {f o r} p \text {a n d} q \text {o d d}
$$

On the other hand, if the algebra $\mathfrak{p}$ has an $N$ -dimensional representation with an invariant bilinear form, $\mathfrak{p}$ can be embedded in $so(N)$ (resp. $sp(N)$ ) if this bilinear form is symmetric (resp. antisymmetric). If the representation has no invariant bilinear form, it realizes an embedding into $su(N)$ .<sup>13</sup> A necessary condition for $\mathsf{L}_{\lambda}$ to have an invariant bilinear form is that $-\lambda \in \Omega_{\lambda}$ , which means that $\mathsf{L}_{\lambda}$ must be self-conjugate. These representations have already been identified in Sect. 13.2.2. The symmetry of the bilinear form is determined by the height of the representation, defined in terms of a height vector $u$ , tabulated in App. 13.A. The form is symmetric

(resp. antisymmetric) if $\lambda \cdot u = \sum_{i} \lambda_{i} u_{i} = 0$ (resp. 1) mod 2. For instance, all representations of $su(2)$ are self-conjugate ( $u_{1} = 1$ ), so the representation $\mathsf{L}_{\lambda}$ is symmetric when $\lambda_{1}$ is even, in which case it leads to the special embeddings $su(2) \subset so(\lambda_{1} + 1)$ . An interesting generic example is $\mathfrak{p} \subset so(\dim \mathfrak{p})$ . Indeed, the adjoint representation is always self-conjugate (i.e., the highest weight is $\theta$ and $-\theta$ is also a root) and it has a symmetric bilinear form (the Killing form).

In the following chapters we often encounter a special embedding, which, although not maximal, deserves a particular mention. This is the diagonal embedding $\mathbf{g} \subset \mathbf{g} \oplus \mathbf{g}$ , in which the two weights $(\lambda, \mu)$ of $\mathbf{g} \oplus \mathbf{g}$ are projected onto the weight $\lambda + \mu$ . Because the highest root is of the same length for the algebra and its subalgebra, the embedding index of a diagonal embedding is always equal to 1.

# Appendix 13.A. Properties of Simple Lie Algebras

The following summaries present the essential information needed for all simple Lie algebras. The Cartan notation is used, and, for the classical algebras, the compact real form is also given in parentheses. For each algebra, we present the Dynkin diagram, a short list of basic properties, the Cartan matrix and the quadratic form matrix. Black nodes in the Dynkin diagrams refer to short roots. The numbers appearing beside the nodes of the Dynkin diagrams give (in this order) the numbering of the corresponding simple root, its mark, and its comark. For simply laced algebras, the third entry is omitted (marks and comarks are identical). The numbering of the simple roots also gives the numbering of the fundamental weights; this is the numbering used when a weight is specified in terms of a sequence of Dynkin labels as $\lambda = (\lambda_1,\dots ,\lambda_r)$ . Marks and comarks are defined in Eq. (13.33). The list of properties includes the dimension of the algebra, the dual Coxeter number $g$ , the order of the Weyl group $|W|$ , the highest root $\theta$ (in Dynkin label notation), the finite group that corresponds to the ratio of the weight lattice $P$ to the root lattice $Q$ , the associated congruence vector $(\nu$ , defined in Sect. 13.1.9), the height vector $(u$ , defined in Sect. 13.7.2) and the exponents (defined at the end of Sect. 13.2.3). For some algebras, the entry " $P/Q"$ and "congruence vector" do not appear; for those cases, $P/Q = I$ .

$$
\mathbf {A} _ {\mathrm {r} \geq 2} (s u (r + 1))
$$

$$
\begin{array}{l l} \bigcirc & (1; 1) \\ \bigcirc & (2; 1) \\ \bigcirc & (3; 1) \\ \bigcirc & (r; 1) \end{array} \qquad \begin{array}{l l} \dim g = r ^ {2} + 2 r \\ g = r + 1 \\ | W | = (r + 1)! \\ \theta = (1, 0, \dots , 1) \\ P / Q = Z _ {r + 1} \\ v = (1, 2, \dots , r) \\ u = (r, 2 (r - 1), \dots , r) \\ \text {e x p o n e n t s} = 1, 2, \dots , r \end{array}
$$

Cartan matrix:

$$
\left( \begin{array}{r r r r r r} 2 & - 1 & 0 & \dots & 0 & 0 \\ - 1 & 2 & - 1 & \dots & 0 & 0 \\ 0 & - 1 & 2 & \dots & 0 & 0 \\ \cdot & \cdot & \cdot & \dots & \cdot & \cdot \\ 0 & 0 & 0 & \dots & 2 & - 1 \\ 0 & 0 & 0 & \dots & - 1 & 2 \end{array} \right)
$$

Quadratic form matrix:

$$
\frac {1}{r + 1} \left( \begin{array}{c c c c c c} r & r - 1 & r - 2 & \dots & 2 & 1 \\ r ^ {\prime} - 1 & 2 (r - 1) & 2 (r - 2) & \dots & 4 & 2 \\ r - 2 & 2 (r - 2) & 3 (r - 2) & \dots & 6 & 3 \\ . & . & . & \dots & . & . \\ 2 & 4 & 2 (r - 1) & \dots & 2 (r - 1) & r - 1 \\ 1 ^ {\prime} & 2 & r - 1 & \dots & r - 1 & r \end{array} \right)
$$

$$
\mathbf {B} _ {\mathrm {r} \geq 3} (s o (2 r + 1))
$$

![](images/34abb1b68da4e9f93b627940cdf3d5da62c03692e14f8852641a9680acdfb655.jpg)

(1;1;1) (2;2;2)

$$
\dim \mathbf {g} = 2 r ^ {2} + r
$$

$$
g = 2 r - 1
$$

$$
| W | = 2 ^ {r} r!
$$

$$
\theta = (0, 1, \dots , 0)
$$

$$
P / Q = \mathbb {Z} _ {2}
$$

![](images/f2ff519776a58b9929f4897919b137e28630cea2511e5f1fbf13d7690c6e1418.jpg)

$(r - 1;2;2)$

$(r;2;1)$

$$
v _ {i} = (0, \dots , 0, 1)
$$

$$
u = (2 r, 2 (2 r - 1), \dots , (r - 1) (r + 2), \frac {1}{2} r (r + 1))
$$

$$
\text {e x p o n e n t s} = 1, 3, \dots , 2 r - 1
$$

Cartan matrix:

$$
\left( \begin{array}{c c c c c c} 2 & - 1 & 0 & \dots & 0 & 0 \\ - 1 & 2 & - 1 & \dots & 0 & 0 \\ 0 & - 1 & 2 & \dots & 0 & 0 \\ \cdot & \cdot & \cdot & \dots & \cdot & \cdot \\ 0 & 0 & 0 & \dots & 2 & - 2 \\ 0 & 0 & 0 & \dots & - 1 & 2 \end{array} \right)
$$

Quadratic form matrix:

$$
\frac {1}{2} \left( \begin{array}{c c c c c c} 2 & 2 & 2 & \dots & 2 & 1 \\ 2 & 4 & 4 & \dots & 4 & 2 \\ 2 & 4 & 6 & \dots & 6 & 3 \\ . & . & . & \dots & . & . \\ 2 & 4 & 6 & \dots & 2 (r - 1) & r - 1 \\ 1 & 2 & 3 & \dots & r - 1 & r / 2 \end{array} \right)
$$

$$
\mathbf {C} _ {r \geq 2} (s p (2 r))
$$

![](images/531aa6ef194df05d379e14d3dcdafd2b3a91f868dbaa0a20664fe12387266aad.jpg)

$$
\dim \mathbf {g} = 2 r ^ {2} + r
$$

$$
g = r + 1
$$

$$
| W | = 2 ^ {r} r!
$$

$$
\theta = (2, 0, \dots , 0)
$$

$$
P / Q = \mathbb {Z} _ {2}
$$

$$
v = (1, 2, \dots , r - 1, r)
$$

$$
u = ((2 r - 2), 2 (2 r - 2), \dots , (r - 1) (r + 1), r ^ {2})
$$

$$
\text {e x p o n e n t s} = 1, 3, \dots , 2 r - 1
$$

Cartan matrix:

$$
\left( \begin{array}{c c c c c c} 2 & - 1 & 0 & \dots & 0 & 0 \\ - 1 & 2 & - 1 & \dots & 0 & 0 \\ 0 & - 1 & 2 & \dots & 0 & 0 \\ \cdot & \cdot & \cdot & \dots & \cdot & \cdot \\ 0 & 0 & 0 & \dots & 2 & - 1 \\ 0 & 0 & 0 & \dots & - 2 & 2 \end{array} \right)
$$

Quadratic form matrix:

$$
\frac {1}{2} \left( \begin{array}{c c c c c c} 1 & 1 & 1 & \dots & 1 & 1 \\ 1 & 2 & 2 & \dots & 2 & 2 \\ 1 & 2 & 3 & \dots & 3 & 3 \\ . & . & . & \dots & . & . \\ 1 & 2 & 3 & \dots & r - 1 & r - 1 \\ 1 & 2 & 3 & \dots & r - 1 & r \end{array} \right)
$$

$\mathbf{D}_{r \geq 4}$ (so(2r))

![](images/4e67000fddcef9ad8dabe3cc8bdcf6ede74693a718aeedf5dd1d4217456ef6b1.jpg)

Cartan matrix:

$$
\left( \begin{array}{c c c c c c c c} 2 & - 1 & 0 & \dots & 0 & 0 & 0 \\ - 1 & 2 & - 1 & \dots & 0 & 0 & 0 \\ 0 & - 1 & 2 & \dots & 0 & 0 & 0 \\ \cdot & \cdot & \cdot & \dots & \cdot & \cdot & \cdot \\ 0 & 0 & 0 & \dots & 2 & - 1 & - 1 \\ 0 & 0 & 0 & \dots & - 1 & 2 & 0 \\ 0 & 0 & 0 & \dots & - 1 & 0 & 2 \end{array} \right)
$$

# $\S 13.\mathbf{A}$ .Properties of Simple Lie Algebras

Quadratic form matrix:

$$
\frac {1}{2} \left( \begin{array}{c c c c c c c} 2 & 2 & 2 & \dots & 2 & 1 & 1 \\ 2 & 4 & 4 & \dots & 4 & 2 & 2 \\ 2 & 4 & 6 & \dots & 6 & 3 & 3 \\ \cdot & \cdot & \cdot & \dots & \cdot & \cdot & \cdot \\ 2 & 4 & 6 & \dots & 2 (r - 2) & r - 2 & r - 2 \\ 1 & 2 & 3 & \dots & r - 2 & r / 2 & (r - 2) / 2 \\ 1 & 2 & 3 & \dots & r - 2 & (r - 2) / 2 & r / 2 \end{array} \right)
$$

E

![](images/288cfda2fd84b32deef5fec670ebd66fe16012a73f37c42ee1208c1ec520a878.jpg)

$$
\begin{array}{l} \dim \mathbf {g} = 2 4 8 \\ g = 3 0 \\ | W | = 6 9 6 7 2 9 6 0 0 \\ \theta = (1, 0, \dots , 0) \\ u = (5 8, 1 1 4, 1 6 8, 2 2 0, 2 7 0, 1 8 2, 9 2, 1 3 6) \\ \text {e x p o n e n t s} = 1, 7, 1 1, 1 3, 1 7, 1 9, 2 3, 2 9 \\ \end{array}
$$

Cartan matrix:

$$
\left( \begin{array}{r r r r r r r r} 2 & - 1 & 0 & 0 & 0 & 0 & 0 & 0 \\ - 1 & 2 & - 1 & 0 & 0 & 0 & 0 & 0 \\ 0 & - 1 & 2 & - 1 & 0 & 0 & 0 & 0 \\ 0 & 0 & - 1 & 2 & - 1 & 0 & 0 & 0 \\ 0 & 0 & 0 & - 1 & 2 & - 1 & 0 & - 1 \\ 0 & 0 \cdot & 0 & 0 & - 1 & 2 & - 1 & 0 \\ 0 & 0 & 0 & 0 & 0 & - 1 & 2 & 0 \\ 0 & 0 & 0 & 0 & - 1 & 0 & 0 & 2 \end{array} \right)
$$

Quadratic form matrix:

$$
\left( \begin{array}{c c c c c c c c} 2 & 3 & 4 & 5 & 6 & 4 & 2 & 3 \\ 3 & 6 & 8 & 1 0 & 1 2 & 8 & 4 & 6 \\ 4 & 8 & 1 2 & 1 5 & 1 8 & 1 2 & 6 & 9 \\ 5 & 1 0 & 1 5 & 2 0 & 2 4 & 1 6 & 8 & 1 2 \\ 6 & 1 2 & 1 8 & 2 4 & 3 0 & 2 0 & 1 0 & 1 5 \\ 4 & 8 & 1 2 & 1 6 & 2 0 & 1 4 & 7 & 1 0 \\ 2 & 4 & 6 & 8 & 1 0 & 7 & 4 & 5 \\ 3 & 6 & . 9 & 1 2 & 1 5 & 1 0 & 5 & 8 \end{array} \right)
$$

E7

![](images/a5e714b97781447c9fac0f968f0bca36792d2db38e16c54ea55be3c9c9e2265e.jpg)

$$
\dim \mathbf {g} = 1 3 3
$$

$$
g = 1 8
$$

$$
| W | = 2 9 0 3 0 4 0
$$

$$
\theta = (1, 0, \dots , 0)
$$

$$
P / Q = \mathbb {Z} _ {2}
$$

$$
\boldsymbol {v} = (0, 0, 0, 1, 0, 1, 1)
$$

$$
u = (3 4, 6 6, 9 6, 7 5, 5 2, 2 7, 4 9)
$$

$$
\text {e x p o n e n t s} = 1, 5, 7, 9, 1 1, 1 3, 1 7
$$

Cartan matrix:

$$
\left( \begin{array}{r r r r r r r} 2 & - 1 & 0 & 0 & 0 & 0 & 0 \\ - 1 & 2 & - 1 & 0 & 0 & 0 & 0 \\ 0 & - 1 & 2 & - 1 & 0 & 0 & - 1 \\ 0 & 0 & - 1 & 2 & - 1 & 0 & 0 \\ 0 & 0 & 0 & - 1 & 2 & - 1 & 0 \\ 0 & 0 & 0 & 0 & - 1 & 2 & 0 \\ 0 & 0 & - 1 & 0 & 0 & 0 & 2 \end{array} \right)
$$

Quadratic form matrix:

$$
\frac {1}{2} \left( \begin{array}{c c c c c c c} 4 & 6 & 8 & 6 & 4 & 2 & 4 \\ 6 & 1 2 & 1 6 & 1 2 & 8 & 4 & 8 \\ 8 & 1 6 & 2 4 & 1 8 & 1 2 & 6 & 1 2 \\ 6 & 1 2 & 1 8 & 1 5 & 1 0 & 5 & 9 \\ 4 & 8 & 1 2 & 1 0 & 8 & 4 & 6 \\ 2 & 4 & 6 & 5 & 4 & 3 & 3 \\ 4 & 8 & 1 2 & 9 & 6 & 3 & 7 \end{array} \right)
$$

E6

![](images/3f921e8ff0486d267e7822417d07a92a0e0309ac2b2843070f070434f8b2007d.jpg)

$$
\dim \mathbf {g} = 7 8
$$

$$
g = 1 2
$$

$$
| W | = 5 1 8 4 0
$$

$$
\theta = (0, 0, \dots , 1)
$$

$$
P / Q = \mathbb {Z} _ {3}
$$

$$
\nu = (1, 2, 0, 1, 2, 0)
$$

$$
u = (1 6, 3 0, 4 2, 3 0, 1 6, 2 2)
$$

$$
\text {e x p o n e n t s} = 1, 4, 5, 7, 8, 1 1
$$

Cartan matrix:

$$
\left( \begin{array}{r r r r r r} 2 & - 1 & 0 & 0 & 0 & 0 \\ - 1 & 2 & - 1 & 0 & 0 & 0 \\ 0 & - 1 & 2 & - 1 & 0 & - 1 \\ 0 & 0 & - 1 & 2 & - 1 & 0 \\ 0 & 0 & 0 & - 1 & 2 & 0 \\ 0 & 0 & - 1 & 0 & 0 & 2 \end{array} \right)
$$

# $\S 13.\mathbf{A}$ .Properties of Simple Lie Algebras

Quadratic form matrix:

$$
\frac {1}{3} \left( \begin{array}{c c c c c c} 4 & 5 & 6 & 4 & 2 & 3 \\ 5 & 1 0 & 1 2 & 8 & 4 & 6 \\ 6 & 1 2 & 1 8 & 1 2 & 6 & 9 \\ 4 & 8 & 1 2 & 1 0 & 5 & 6 \\ 2 & 4 & 6 & 5 & 4 & 3 \\ 3 & 6 & 9 & 6 & 3 & 6 \end{array} \right)
$$

F4

![](images/03011d3772518cd200641b3f44054ad1ccc75a78edabf5134a52df982f593e0b.jpg)

$$
\dim \mathbf {g} = 5 2
$$

$$
g = 9
$$

$$
| W | = 1 1 5 2
$$

$$
\theta = (1, 0, 0, 0)
$$

$$
u = (2 2, 4 2, 3 0, 1 6)
$$

$$
\text {e x p o n e n t s} = 1, 5, 7, 1 1
$$

Cartan matrix:

$$
\left( \begin{array}{r r r r} 2 & - 1 & 0 & 0 \\ - 1 & 2 & - 2 & 0 \\ 0 & - 1 & 2 & - 1 \\ 0 & 0 & - 1 & 2 \end{array} \right)
$$

Quadratic form matrix:

$$
\left( \begin{array}{c c c c} 2 & 3 & 2 & 1 \\ 3 & 6 & 4 & 2 \\ 2 & 4 & 3 & \frac {3}{2} \\ 1 & 2 & \frac {3}{2} & 1 \end{array} \right)
$$

G2

![](images/73c6f0b35f31a9a9ac2e5090ebd8b83f5c2be0d3ee43021208b3c99634cf1a96.jpg)

(1;2;2)

(2;3;1)

$$
\dim \mathbf {g} = 1 4
$$

$$
g = 4
$$

$$
| W | = 1 2
$$

$$
\theta = (1, 0)
$$

$$
u = (1 0, 6)
$$

$$
e x p o n e n t s = 1, 5
$$

Cartan matrix:

$$
\left( \begin{array}{c c} 2 & - 3 \\ - 1 & 2 \end{array} \right)
$$

Quadratic form matrix:

$$
\frac {1}{2} \left( \begin{array}{l l} 6 & 3 \\ 2 & 0 \end{array} \right)
$$

# Appendix 13.B. Notation for Simple Lie Algebras

g, h: finite Lie algebras

$G, H$ : corresponding Lie groups

dim $\mathbf{g}$ : dimension of the algebra $\mathbf{g}$

$\pmb{r}$ :rank

$J^{a}$ $(a = 1,\dots ,\dim \mathbf{g})$ : generators of $\mathbf{g}$

$H^{i}$ $(i = 1,\dots ,r)$ : generators of the Cartan subalgebra in the Cartan-Weyl basis

$E^{\alpha}$ : ladder generators in the Cartan-Weyl basis

$h^i (i = 1,\dots ,r)$ : generators of the Cartan subalgebra in the Chevalley basis

$e^{i},f^{i}$ $(i = 1,\dots ,r)$ : raising and lowering operators associated with the simple roots in the Chevalley basis

$\alpha ,\beta :$ roots

$\alpha^i: i$ -th component of $\alpha$ in the Cartan-Weyl basis

$\alpha_{i}$ : simple roots

$\alpha_{i}^{\vee}$ : simple coroots $= 2\alpha_{i} / \alpha_{i}^{2}$

$A_{ij}$ : Cartan matrix element $= 2(\alpha_i, \alpha_j^{\vee})$

$\Delta$ : set of roots

$\Delta_{+}, \Delta_{-}$ : set of positive, negative roots

$|\Delta |,|\Delta_{+}|:$ number of roots, number of positive roots

$\theta$ : highest root

$a_{i}, a_{i}^{\vee}$ : marks and comarks; $\theta = \sum_{i=1}^{r} a_{i} \alpha_{i} = \sum_{i=1}^{r} a_{i}^{\vee} \alpha_{i}^{\vee}$

$g$ : dual Coxeter number $= \sum_{i = 1}^{r}a_{i}^{\vee} + 1$

$\lambda, \mu, \nu$ : finite weights (usually highest weights)

$\dim |\lambda |:$ dimension of the representation of highest weight $\pmb{\lambda}$

$\Omega_{\lambda}$ : weight system of the representation of highest weight $\pmb{\lambda}$

$\mathsf{L}_{\lambda}$ : irreducible module of highest weight $\pmb{\lambda}$

$\lambda^{\prime}, \lambda^{\prime \prime}$ : particular weights in $\Omega_{\lambda}$

$\left|\lambda^{\prime}\right\rangle$ : particular state, of weight $\lambda^\prime$ , in the module $\mathsf{L}_{\lambda}$

$\mathrm{mult}_{\lambda}(\lambda')$ : multiplicity of $\lambda'$ in the highest-weight representation $\lambda$

$\rho$ : Weyl vector (half-sum of positive roots)

$\omega_{i}$ : fundamental weights

$F_{ij}$ : quadratic form matrix $= (\omega_{i},\omega_{j})$

$\lambda_{i}$ : Dynkin labels: $\lambda = \sum_{i=1}^{r} \lambda_{i} \omega_{i} = (\lambda_{1}, \ldots, \lambda_{r})$ ; $h^{i}$ eigenvalues of $|\lambda\rangle$

$\lambda^i:H^i$ eigenvalues of $|\lambda \rangle$

$\epsilon_{i}$ : orthonormal vectors

$\{\ell_1; \ell_2; \ldots; \ell_r\} : \text{partition of the Young tableau associated with the } su(N) \text{ weight }$ $\lambda = \sum_{i} \ell_i \epsilon_i = \sum_{i} \lambda_i \omega_i$ so that $\ell_i = \lambda_i + \lambda_{i+1} + \ldots + \lambda_r = \text{length of the } i$ -th row (from top)

$\{\tilde{\ell}_1;\tilde{\ell}_2;\ldots ;\tilde{\ell}_s\} :$ transposed partition (change rows and columns)

$\lambda^t$ : transposed weight

$\lambda^{*}$ : conjugate of $\pmb{\lambda}$

$\chi_{\lambda}$ : character of the representation $\lambda$

W: Weyl group

$|W|$ : order of the Weyl group

$s_{\alpha}$ : reflection with respect to the root $\alpha$

$s_i$ : reflection with respect to the simple root $\alpha_i$ (a simple Weyl reflection)

$\pmb{w}$ : element of the Weyl group (a Weyl reflection)

$w_0$ : longest element of the Weyl group

$\epsilon (\pmb {w})$ : signature of $\pmb{w}$

$\ell (w)$ : length of $\pmb{w}$

$\pmb{w}:$ shifted Weyl reflection: $\pmb {w}\cdot \lambda = \pmb {w}(\lambda +\rho) - \rho$

$C_w$ : Weyl chamber associated with the element $w$

$Q$ : root lattice

$Q^{\vee}$ : coroot lattice

$\pmb{P}$ : weight lattice

$P_{+}$ : set of dominant weights ( $\equiv$ set of highest weights for irreducible representations)

$\mathcal{N}_{\lambda \mu \nu} = \mathcal{N}_{\lambda \mu}{}^{\nu^{*}}$ : tensor-product coefficients

$\mathcal{Q}$ :quadraticCasimir operator

ad : adjoint operator; $\operatorname {ad}(X)Y = [X,Y]$

$K(\mathbf{\theta},\mathbf{\theta}):$ (normalized) Killing form; $K(X,Y) = \mathrm{Tr}(\mathrm{ad}X,\mathrm{ad}Y) / 2g$

$\kappa$ : Kostant partition function

$\pmb{x}_{\lambda}$ : Dynkin index of the representation $\pmb{\lambda}$

$\pmb{x}_{e}$ : embedding index

$\pmb{b}_{\lambda \mu}$ : branching coefficient (multiplicity of $\mathsf{L}_{\mu}$ in $\mathsf{L}_{\lambda}$ )

$\pmb{\nu}$ : congruence vector

$\pmb{u}$ : height vector

$B(G)$ : center of the group $G^{14}$

# Exercises

# 13.1 The Killing form

a) Verify Eq. (13.18) and check that the only nonzero Killing norms are $K(H^i, H^i)$ and $K(E^\alpha, E^{-\alpha})$ .

b) Calculate the $su(2)$ Killing form $\tilde{K}$ in the Chevalley basis (13.85). Result: With the ordering $e, h, f$ , it reads

$$
\tilde {K} = \left( \begin{array}{c c c} 0 & 0 & 4 \\ 0 & 8 & 0 \\ 4 & 0 & 0 \end{array} \right)
$$

A rescaling by a factor $\frac{1}{4} = 1 / (2g)$ yields the standard normalization:

$$
K (e, f) = K (f, e) = \frac {1}{2} K (h, h) = 1
$$

# 13.2 Weyl group for $G_{2}$ and $su(4)$

Starting from the corresponding Cartan matrix given in App. 13.A, find the Weyl group and the set of all roots for:

a) $G_{2}$   
b) $su(4)$

# 13.3 Linear representation of the Weyl group

The linear representation of the simple Weyl reflection $s_j$ is the $r \times r$ matrix that maps the column vector with components $\lambda_i$ to that with components $(s_j\lambda)_i$ .

a) Show that $\det s_j = -1$ . Deduce that for a general Weyl reflection $w$ ,

$$
\det w = (- 1) ^ {\ell}
$$

where $\ell$ is the number of simple reflections in the decomposition of $w$ .

b) Find the matrix representation of the simple reflections of $G_{2}$ and verify the relations (13.59).   
c) Same as (b) for the algebra $F_4$ .

# 13.4 Order of the Weyl group

Verify the following formula for the order of the Weyl group of a simple Lie algebra of rank $r$ with marks $\{a_i\}$ :

$$
| W | = | P / Q | r! \prod_ {i = 1} ^ {r} a _ {i}
$$

Proceed case by case, using the data of App. 13.A.

# 13.5 Weight systems

Write all weights in the representation of highest weight:

a) $(1,0)$ of $G_{2}$   
b) $(0,0,1)$ of $so(7)$ ,   
c) $(0,0,0,0,1)$ of so(10).

# 13.6 Weight multiplicities

Find the multiplicity of the $su(4)$ weight $(-2,3,0)$ in the representation $(3,1,1)$ using:

a) the Freudenthal formula (13.113);   
b) semistandard tableaux (cf. Sect. 13.3.3).

Hint: The calculation in (a) is greatly simplified if the weight is first transformed into a dominant one.

# 13.7 su(3) Gelfand-Tsetlin patterns

a) Write all the Gelfand-Tsetlin patterns for the $su(3)$ representation of highest weight (2,2).   
b) For a $su(3)$ weight $\lambda' \in \Omega_{\lambda}$ , there corresponds $\mathrm{mult}_{\lambda}(\lambda')$ Gelfand-Tsetlin patterns of the form

$$
\begin{array}{c c c} \lambda_ {1} + \lambda_ {2} & \lambda_ {2} & 0 \\ a & b \\ c & \end{array}
$$

Relate the parameters $a, b, c$ to the two Dynkin labels $\lambda_1', \lambda_2'$ . Find inequalities satisfied by the free parameter of the Gelfand-Tsetlin pattern, and deduce a simple formula for $\mathrm{mult}_{\lambda}(\lambda')$ . Compare with the example of part (a).

# 13.8 The Demazure character formula

An expression equivalent to Eq. (13.167) is given by

$$
\chi_ {\lambda} = M _ {w _ {0}} \left(e ^ {\lambda}\right)
$$

where $w_0$ is the longest element of the Weyl group, and for $w_0 = s_i \cdots s_j$ , $M_{w_0}(e^\lambda)$ is defined by

$$
M _ {w _ {0}} = M _ {i} \dots M _ {j}
$$

with

$$
M _ {i} (e ^ {\lambda}) = \frac {e ^ {\lambda} - e ^ {s _ {i} \cdot \lambda}}{1 - e ^ {- \alpha_ {i}}}
$$

(notice that the Weyl reflection is shifted), where, as usual, $\alpha_{i}$ stands for a simple root and

$$
M _ {i} M _ {j} \left(e ^ {\lambda}\right) \equiv M _ {i} \left(M _ {j} \left(e ^ {\lambda}\right)\right)
$$

This is called the Demazure character formula.

a) Verify the following properties of $M_{i}$ :

$$
\begin{array}{l} M _ {i} \left(e ^ {\lambda}\right) = e ^ {\lambda} + e ^ {\lambda + 1} + \dots + e ^ {\lambda - \lambda_ {i} \alpha_ {i}} \quad \text {i f} \quad \lambda_ {i} \geq 0 \\ = 0 \quad \text {i f} \quad \lambda_ {i} = - 1 \\ = - e ^ {\lambda + \alpha_ {i}} - e ^ {\lambda + \alpha_ {i} + 1} - \dots - e ^ {\lambda - (\lambda_ {i} + 1) \alpha_ {i}} \quad \text {i f} \quad \lambda_ {i} \leq - 2 \\ \end{array}
$$

and

$$
(M _ {i}) ^ {2} = M _ {i}
$$

b) For $su(2)$ , show that the Demazure formula is equivalent to the Weyl character formula.  
c) Check the formula for the $su(3)$ representation (1,2) (compare the result with Eq. (13.161)). For this representation, verify also that

$$
M _ {s _ {1} s _ {2} s _ {1}} \left(e ^ {\lambda}\right) = M _ {s _ {2} s _ {1} s _ {2}} \left(e ^ {\lambda}\right)
$$

d) Another version of the Demazure formula is

$$
\chi_ {\lambda} = \sum_ {w \in W} N _ {w} (e ^ {\lambda})
$$

where, in terms of a (minimal) decomposition of $w$ in simple Weyl reflections, e.g., if $w = s_{l} \cdots s_{k}, N_{w}$ is given by

$$
N _ {w} (e ^ {\lambda}) = N _ {l} \dots N _ {k} (e ^ {\lambda})
$$

and

$$
N _ {i} (e ^ {\lambda}) = \frac {e ^ {s _ {i} \lambda} - e ^ {\lambda}}{1 - e ^ {\alpha_ {i}}}
$$

Express $N_{i}$ as a sum, as done in part (a) for $M_{i}$ .

e) Evaluate the different $N_{w}(e^{\lambda})$ 's for the $su(3)$ highest weight $\lambda = (1,2)$ . Observe that each $N_{w}(e^{\lambda})$ is a positive sum.

f) Prove the relation:

$$
(1 + N _ {i}) (e ^ {\lambda}) = M _ {i} (e ^ {\lambda})
$$

# 13.9 Dimension of $G_{2}$ representations

Derive the dimension formula for the irreducible representations of $G_2$ and check that $L_{(0,1)}$ and $L_{(1,0)}$ have respective dimensions 7 and 14.

# 13.10 Another expression for the dual Coxeter number

Equations (13.181) and (13.184) lead to the following expression for the dual Coxeter number:

$$
\begin{array}{l} g = (2 n _ {L} + n _ {S}) / 2 r \quad \text {f o r} \quad g \neq G _ {2} \\ = \left(3 n _ {L} + n _ {S}\right) / 3 r \quad \text {f o r} \quad G _ {2} \\ \end{array}
$$

where $n_{L,S}$ denotes the number of long and short roots, respectively. Verify this result for $sp(4)$ and $G_2$ .

Remark: For simply laced algebras, this reduces to the relation: $|\Delta| = gr$ .

# 13.11 Cauchy determinant and Schur functions

a) Show that

$$
\phi (\{x \}, \{y \}) = \frac {\Delta (y)}{\prod_ {1 \leq i , j \leq N} (1 - x _ {i} y _ {j})}
$$

where $\Delta (x) = \prod_{1\leq i < j\leq N}(x_i - x_j)$ , is a generating function for the Schur functions (13.189), namely

$$
\phi (\{x \}, \{y \}) = \sum_ {m _ {1}, m _ {2}, \dots , m _ {N} \geq 0} y _ {1} ^ {m _ {1}} \dots y _ {N} ^ {m _ {N}} S _ {\lambda} (x _ {1}, \dots , x _ {N})
$$

where $\lambda = \{\ell_i\}$ , and $\ell_{i} = m_{i} + i - N$

b) By means of the Cauchy determinant formula (see Ex. 12.12 for a proof; take $z_{i} = 1 / x_{i}$ and $w_{j} = y_{j}$ in the formula (12.195))

$$
\det  \left[ \frac {1}{1 - x _ {i} y _ {j}} \right] _ {1 \leq i, j \leq N} = \frac {\Delta (x) \Delta (y)}{\prod_ {1 \leq i , j \leq N} (1 - x _ {i} y _ {j})}
$$

rewrite the generating function $\phi (\{x\} ,\{y\})$ as the single determinant

$$
\begin{array}{l} \phi (\{x \}, \{y \}) = \frac {\Delta (y)}{\prod_ {1 \leq i , j \leq N} (1 - x _ {i} y _ {j})} \\ = \det  \left[ \frac {y _ {i} ^ {N - j}}{\prod_ {k = 1} ^ {N} \left(1 - y _ {i} x _ {k}\right)} \right] _ {1 \leq i, j \leq N} \\ \end{array}
$$

Hint: Represent the quantity $\Delta(y)$ as a determinant (13.191).

c) The Schur polynomials of the variables $t_1, t_2, \ldots$ are defined through the generating function

$$
F (y) = \sum_ {m \geq 0} y ^ {m} P _ {m} (t.) = e ^ {\sum_ {k = 1} ^ {\infty} y ^ {k} \frac {t _ {k}}{k}}
$$

This definition is supplemented by the convention that $P_{m}(t_{\cdot}) = 0$ for $m \leq -1$ . Show that

$$
F (y) = \prod_ {k = 1} ^ {N} \frac {1}{(1 - y x _ {k})}
$$

iff the $t_k$ are expressed as

$$
t _ {k} = \sum_ {i = 1} ^ {N} x _ {i} ^ {k}
$$

for some integer $N$ .

d) Prove the following properties of the Schur polynomials

$$
\frac {\partial}{\partial t _ {k}} P _ {m} (t.) = P _ {m - k} (t.)
$$

$$
P _ {m} (1) = \frac {1}{m !}
$$

where 1 stands for $t_k = 1$ for all $k \geq 1$ .

e) Express the generating function $\phi(\{x\}, \{y\})$ in terms of Schur polynomials. Deduce the following formula expressing the Schur functions as determinants of Schur polynomials of the variable $t_k = \sum_{i=1}^{N} x_i^k$ .

$$
S _ {\lambda} \left(x _ {1}, \dots , x _ {N}\right) = \det  \left[ P _ {\ell_ {i} + j - i} (t.) \right] _ {1 \leq i, j \leq N}
$$

# 13.12 Partitions and Schur functions

a) Work out the details of the derivation of Eqs. (13.189) and (13.192).   
b) Prove directly the equivalence of Eqs. (13.192) and (13.172) by evaluating the scalar products in Eq. (13.172) in the orthogonal basis.   
c) Find the action of the $s_i$ 's on the partitions.

# 13.13 Dimension of $su(N)$ representations and hooks

The dimension of a representation can be read off a Young tableau in a rather simple way using hooks. The hook associated with the box at position $(i,j)$ ( $i$ -th row, $j$ -th column) is composed of two lines joined at right angle in the box $(i,j)$ and leaving the tableau downward and toward the right. Its length, denoted by $h_{i,j}$ , is the number of boxes it crosses. The following tableau is filled with the numbers $h_{i,j}$

$$
h _ {i, j}: \begin{array}{c c c c} \hline 6 & 5 & 2 & 1 \\ \hline 3 & 2 \\ \hline 2 & 1 \\ \hline \end{array}
$$

In terms of hooks, the dimension of a $su(N)$ representation reads

$$
\dim | \lambda | = \prod_ {i, j} \frac {(N - i + j)}{h _ {i , j}}
$$

where the product is taken over all the boxes of the tableau.

a) Verify the equivalence of this formula with Eq. (13.192) for the above $su(4)$ tableau.   
b) Using this expression, reproduce the $su(2)$ and $su(3)$ dimension formulae (13.172).

# 13.14 sp(4) tensor product: character method

Calculate the $sp(4)$ tensor product $(1,1) \otimes (2,0)$ using the character method and check the result by calculating the total dimension of each sides.

# 13.15 Weyl-group folding in the character method

Extending the validity of Eq. (13.171) to nondominant weights, prove that

$$
\dim | w \cdot \lambda | = \epsilon (w) \dim | \lambda |
$$

In the character method for tensor-product calculations, this shows that weights that are ignored have zero dimension, and two weights cancel each other if their dimensions add up to zero. Check this explicitly for the $su(3)$ example $(3,2) \otimes (2,4)$ , to be worked out graphically using the algorithm underlying the character method.

# 13.16 Littlewood-Richardson and Berenstein-Zelevinsky methods

a) Using the Littlewood-Richardson method once and then the BZ triangles, calculate the following tensor products:

$$
\begin{array}{l} s u (3): (3, 2) \otimes (0, 3) \\ s u (4): (1, 0, 1) \otimes (1, 0, 1) \\ \end{array}
$$

b) Using Littlewood-Richardson tableaux once and then the BZ triangles, find the multiplicity of the scalar representation in the following triple tensor products

$$
\begin{array}{l} s u (3): (4, 4) \otimes (4, 4) \otimes (4, 4) \\ s u (4): (2, 1, 1) \otimes (1, 2, 1) \otimes (1, 1, 2) \\ \end{array}
$$

c) Observe that all the $su(3)$ triangles in (b) are related to each other by addition or subtraction of the "basic" triangle

$$
\Omega = \begin{array}{c c c c} & 1 \\ & - 1 & - 1 \\ & - 1 & - 1 \\ 1 & - 1 & - 1 & 1 \end{array}
$$

Hence, once a triangle is found, all the others are readily generated. Relate this to a one-parameter indeterminacy in (13.212). Find the analogous result for $su(4)$ and compare with the example worked out in (b).

d) Prove, using either Littlewood-Richardson tableaux or BZ triangles, that the $su(3)$ tensor-product coefficient $\mathcal{N}_{\lambda \mu \nu}$ is at most 1 if one of the three weights has at least one vanishing Dynkin label.

# 13.17 Kostant's multiplicity formula

The Weyl character formula leads directly to a new expression for weight multiplicities, Kostant's formula. For this, we introduce the partition function $\kappa(\mu)$ defined to be the number of distinct decompositions of $\mu$ in terms of positive roots. In other words, $\kappa(\mu)$ is

the number of solutions $\{k_{\alpha}\}$ , $\alpha \in \Delta_{+}$ of the equation $\sum_{\alpha > 0} k_{\alpha} \alpha = \mu$ , with all $k_{\alpha} \geq 0$ . Of course, if there is no such decomposition, $\mathcal{K}(\mu) = 0$ . Setting $\mathcal{K}(0) = 1$ , we have

$$
\prod_ {\alpha > 0} \frac {1}{1 - e ^ {\alpha}} = \sum_ {\mu} \mathcal {K} (\mu) e ^ {\mu}
$$

In terms of this partition function, show that the multiplicity of the weight $\lambda'$ in the representation $\lambda$ is given by

$$
\operatorname {m u l t} _ {\lambda} \left(\lambda^ {\prime}\right) = \sum_ {w \in W} \epsilon (w) \mathcal {K} \left(w (\lambda + \rho) - \left(\lambda^ {\prime} + \rho\right)\right)
$$

Hint: Use the product form of $D_{\rho}^{-1}$ to relate it to the partition function $\kappa$ .

The advantage of Kostant's formula over Freudenthal's is that a given weight can be treated in isolation. The price that has to be paid is a sum over the whole Weyl group. Nevertheless, in favorable circumstances only a few terms contribute. Illustrate this by calculating the multiplicity of the weight $(0,0)$ in the adjoint representation of $su(3)$ .

# 13.18 Steinberg formula for tensor products

Use the Kostant multiplicity formula to obtain the Steinberg formula for tensor-product coefficients:

$$
\mathcal {N} _ {\lambda \mu} ^ {\nu} = \sum_ {w, w ^ {\prime} \in W} \epsilon (w w ^ {\prime}) \kappa (w \cdot \lambda + w ^ {\prime} \cdot \mu - \nu)
$$

# 13.19 Associativity in tensor products

Tensor product coefficients can be calculated from the fusion coefficients involving fundamental weights, that is, $\{\mathcal{N}_{\lambda \mu}^{\omega_i}\}$ for $i = 1,\dots ,r$ and any $\lambda ,\mu$ , and the associativity condition (13.224). Illustrate this by calculating, from these data, the $su(3)$ coefficient $\mathcal{N}_{(1,1)(1,1)}^{(1,1)}$

# 13.20 Generalized Chebyshev polynomials and tensor products

a) Verify the relations (13.239), regarded as the defining recursion relations for the generalized Chebyshev polynomials $U_{(\lambda_1,\lambda_2)}$ , associated with the tensor-product matrix $N_{(\lambda_1,\lambda_2)}$ . Check further that

$$
U _ {(\lambda_ {1}, \lambda_ {2})} = U _ {(\lambda_ {1}, 0)} U _ {(0, \lambda_ {2})} - U _ {(\lambda_ {1} - 1, 0)} U _ {(0, \lambda_ {2} - 1)}
$$

for $\lambda_1, \lambda_2 > 1$ . Argue that the matrices $N_{(1,0)}$ and $N_{(0,1)}$ must commute. Use these relations to obtain the generating function (13.241).

b) Derive analogous results for $sp(4)$ . With $N_{(1,0)} = x_1$ and $N_{(0,1)} = x_2$ , the generating function $F(x_1, x_2; t, s)$ is

$$
\frac {1 + s \left(t ^ {2} + 1\right) + s ^ {2} t ^ {2} - t s x _ {1}}{\left(1 + t ^ {2} + t ^ {4} - x _ {1} \left(t ^ {3} + t\right) + t ^ {2} x _ {2}\right) \left(1 + s + s ^ {3} + s ^ {4} - x _ {2} \left(s + 2 s ^ {2} + s ^ {3}\right) + s ^ {2} x _ {1} ^ {2}\right)}
$$

# 13.21 Verlinde formula for a Lie algebra

Check carefully the derivation of the orthogonality relation (13.244). Use the Verlinde formula (13.245) to recover the $su(2)$ tensor-product matrices $N_{1}$ and $N_{2}$ .

# 13.22 Embedding index

a) Prove the relation (13.251).

b) For the embedding $E_{8} \supset su(2) \oplus su(3)$ , calculate the embedding index, using the branching rule:

$$
(1, 0, 0, 0, 0, 0, 0, 0) \mapsto \{(6) \otimes (1, 1) \} \oplus \{(4) \otimes (3, 0) \}
$$

c) For the embedding $so(7) \supset su(4)$ , calculate the embedding index, using the projection matrix:

$$
\mathcal {P} = \left( \begin{array}{c c c} 0 & 1 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{array} \right)
$$

# 13.23 Embeddings of $su(2)$

a) Describe all possible embeddings of $su(2)$ in $sp(4)$ . In each case, find the branching rule for $(1,0)$ , the projection matrix and the embedding index.

b) Same as (a) for the embeddings $su(2)\subset G_2$ , using the representation $(0,1)$ .

# 13.24 Regular maximal subalgebras

Find all regular maximal subalgebras of $F_4, E_6$ , and $E_7$ .

# 13.25 Branching rules in regular embeddings

a) Consider the regular embedding $su(3) \subset G_2$ . Draw the extended Dynkin diagram of $G_2$ (i.e., calculate the number of links between the new root $-\theta$ and $\alpha_1, \alpha_2$ ). Identify the node that must be deleted to recover the $su(3)$ Dynkin diagram. Write all the weights in the $(0, 1)$ representation of $G_2$ and their extended Dynkin labels $[\lambda_{-\theta}, \lambda_1, \lambda_2]$ , where

$$
\lambda_ {- \theta} = - 2 \lambda_ {1} - \lambda_ {2}
$$

(cf. Eq. (13.268)). Delete the Dynkin label appropriate for the $su(3)$ embedding and reorganize the resulting $su(3)$ weights in irreducible representations. This gives the branching of the $(0,1)G_{2}$ representation into $su(3)$ ones.

b) By proceeding similarly for the regular embedding $su(4) \subset so(7)$ , find the branching of the $so(7)$ representation $(1,0,0)$ .

# Notes

Except for some aspects of tensor-product calculations and tableaux techniques, the content of this chapter is rather standard. It is covered, for instance, in Cahn [61], Wybourne [361], Fulton and Harris [155], Jacobson [209], Humphreys [196], Bourbaki [56], and Zelobenko [368]. The book of Cahn provides a clear and concise first introduction to the subject, and that of Fulton and Harris is a particularly readable mathematical textbook; tableaux techniques are well covered there. A sharp focus on the material presented in Sects. 13.1 and 13.2 can be found in those sections of Kass et al. [228] related to finite Lie algebras. The theory of semisimple Lie algebras is also well summarized in the first chapter of Fuchs [148]. The proof of the strange formula follows Freudenthal and de Vries [138]. The relation between semistandard tableaux and Gelfand-Tsetlin patterns can be found in Ref. [193].

The character method for tensor products is presented in Racah [301], Speiser [329], and Klimyk [239]. The relation between Littlewood-Richardson tableaux and Gelfand-Tsetlin patterns can be found in Gelfand and Zelevinsky [164]. It is equivalent to the method for calculating tensor-product coefficients by means of semistandard tableaux, which is

presented in [257, 354, 278]. Berenstein-Zelevinsky triangles were introduced in Ref. [38] and further developed in Refs. [74, 39].

The basics of algebra embeddings are explained in Cahn [61]. For a more detailed discussion, the reader is referred to the original articles of Dynkin [117, 118]. The generating functions for the embeddings of $su(2)$ into $su(3)$ (and many others) can be found in Patera and Sharp [291].

The Demazure formula of Ex. 13.8 is proved in Ref. [90] (see also Ref. [163]).

Our conventions and most of our notations follow mainly that of Patera and collaborators [268, 59], which makes easier the consultation of these extensive and very useful tables of weight multiplicities, dimensions of representations, branching rules, and so forth: