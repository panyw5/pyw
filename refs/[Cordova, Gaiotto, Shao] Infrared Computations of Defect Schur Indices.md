# Infrared Computations of Defect Schur Indices

Clay Córdova<sup>1*</sup>, Davide Gaiotto<sup>2†</sup>, and Shu-Heng Shao<sup>3‡</sup>

$^{1}$ School of Natural Sciences, Institute for Advanced Study, Princeton, NJ, USA

$^{2}$ Perimeter Institute for Theoretical Physics, Waterloo, Ontario, Canada N2L 2Y5

$^{3}$ Jefferson Physical Laboratory, Harvard University, Cambridge, MA, USA

# Abstract

We conjecture a formula for the Schur index of four-dimensional  $\mathcal{N} = 2$  theories in the presence of boundary conditions and/or line defects, in terms of the low-energy effective Seiberg-Witten description of the system together with massive BPS excitations. We test our proposal in a variety of examples for  $SU(2)$  gauge theories, either conformal or asymptotically free. We use the conjecture to compute these defect-enriched Schur indices for theories which lack a Lagrangian description, such as Argyres-Douglas theories. We demonstrate in various examples that line defect indices can be expressed as sums of characters of the associated two-dimensional chiral algebra and that for Argyres-Douglas theories the line defect OPE reduces in the index to the Verlinde algebra.

# June, 2016

# Contents

# 1 Introduction 2

1.1 The IR Formula for the Schur Index 5  
1.2 The IR Formula for the Schur Index with Line Defects 7  
1.3 The IR Formula for the Schur Half-Index 9

# 2 A Review of the Schur Index and its IR Formulation 10

2.1 The Schur Index 10  
2.2 An IR Formula for the Schur Index 12

2.2.1 Convergence of the IR Formula 15

2.3  $SU(2)$  Gauge Theory with  $N_{f} = 4$  Flavors 16

# 3 Line Defects and Their Schur Indices 19

3.1 Defect Junctions and the Schur Index 20  
3.2 Examples of the Line Defect Schur Index 24

3.2.1 Wilson Lines in  $SU(2)$  Gauge Theory 24  
3.2.2 't Hooft Lines in  $SU(2)$  Gauge Theory 25  
3.2.3 Wilson Lines in  $SU(2)$  Gauge Theory with  $N_{f} = 4$  Flavors 26

3.3 An IR Formula for the Line Defect Schur Index 27  
3.4 Examples of the IR Formula 29

3.4.1 Wilson Lines in  $SU(2)$  Gauge Theory 29  
3.4.2 't Hooft Lines in  $SU(2)$  Gauge Theory 31  
3.4.3 Wilson Lines in  $SU(2)$  Gauge Theory with  $N_{f} = 4$  Flavors 31

# 4 Half-Indices and Boundary Conditions 34

4.1 Indices and RG Interfaces 37  
4.2 A Free Hypermultiplet 38  
4.3 Pure  $SU(2)$  Gauge Theory 39

# 5 Defect Indices in Argyres-Douglas Theories 42

5.1 Line Defect OPEs and Schur Indices 43  
5.2  $A_{2}$  Argyres-Douglas Theory 44  
5.3  $A_{3}$  Argyres-Douglas Theory 46  
5.4  $A_{4}$  Argyres-Douglas Theory 49

# 6 Chiral Algebra and Line Defects 50

6.1  $SU(2)$  Gauge Theory with  $N_{f} = 4$  Flavors 51  
6.2 Verlinde Algebra from Line Defects 55

6.2.1  $A_{2}$  Argyres-Douglas Theory 57  
6.2.2  $A_{3}$  Argyres-Douglas Theory 59  
6.2.3  $A_{4}$  Argyres-Douglas Theory 61

# A Supercharges of Line Defects and Chiral Algebras 63

A.1 Supercharges Preserved by Full Lines 64  
A.2 Supercharges Preserved by Half Lines 66  
A.3 Supercharges Shared by the Chiral Algebra and Line Defects 66

# B Framed Quivers for the Argyres-Douglas Theories 68

B.1  $A_{3}$  Argyres-Douglas Theory 69  
B.2  $A_{4}$  Argyres-Douglas Theory 72

# C Affine Characters of Kac-Moody Algebra at Negative Level 75

C.1 Generalities on Affine Lie Algebra 75  
C.2 Affine Characters and the Kazhdan-Lusztig Polynomials 79

C.2.1 The Kazhdan-Lusztig Conjecture 81

C.3 Affine Characters of  $\widehat{su(2)}_{-\frac{4}{3}}$  81  
C.4 Affine Characters of  $\widehat{so(8)}_{-2}$  82

# 1 Introduction

The Schur index, introduced in [1-3], is a specialization of the superconformal index of four-dimensional  $\mathcal{N} = 2$  theories. It depends on a single fugacity  $q$ , and counts only superconformal multiplets which are quarter BPS. As a trace over the Hilbert space on  $S^3$  it is defined as

$$
\mathcal {I} (q) = \mathrm {T r} \left[ (- 1) ^ {F} q ^ {\Delta - R} \right], \tag {1.1}
$$

where  $\Delta$  is the scaling dimension and  $R$  is the Cartan of the  $SU(2)_R$  symmetry. Because of the enhanced supersymmetry, the Schur index is highly computable even for non-Lagrangian theories [4-8]. For instance, in the context of class  $\mathcal{S}$  theories, the Schur index equals the  $q$ -deformed topological Yang-Mills partition function [2,3,9-12]. More recently, it was

demonstrated that the local operators counted by this index secretly form a two-dimensional chiral algebra and correspondingly,  $\mathcal{I}(q)$  is the vacuum character of that algebra [13] (see also [14-21] for further developments).

One consequence of the extra supersymmetry of the Schur index is that it can be modified by adding to the system a variety of half-BPS defects, including in particular line defects and interfaces or boundary conditions. Supersymmetric line defects have been previously studied in [22-31], while BPS boundary conditions in  $\mathcal{N} = 2$  theories were investigated in [32-38]. Indices enriched by these defects have been computed in [36,39-41] and are the main objects of study in this work.

The structure that results strongly resembles that of the four-dimensional ellipsoid  $(S_b^4)$  partition function decorated by defects [42-47]. Moreover, the relationship between the two indices is analogous to the relationship between the three-dimensional ellipsoid  $(S_b^3)$  partition function and the superconformal index of three-dimensional  $\mathcal{N} = 2$  theories [36]. A striking feature of the ellipsoid partition function in both three and four dimensions is that although the original partition function can be defined for superconformal field theories by a conformal transformation from flat space, these partition functions can also be defined for non-conformal theories as well.

The situation with the Schur index appears to be similar, and recent work by some of the authors [8] strongly indicates that:

- The Schur index is a meaningful quantity for non-conformal  $\mathcal{N} = 2$  field theories.  
- The Schur index can be computed in the IR using the Seiberg-Witten description [48, 49] as an Abelian gauge theory enriched by BPS particles. The details of the calculation depend on the choice of a chamber in the Coulomb branch but the result does not.

Although the two statements are logically independent, they are closely related. Once we accept that the index can be computed on the Coulomb branch, which spontaneously breaks both conformal symmetry and  $U(1)_r$ , it is natural to expect that these symmetries are not needed to define the index. $^2$

The general idea that wall-crossing invariant generating functions of BPS states in four-dimensional field theories are related to the local operators at the UV superconformal fixed point originated in [51, 52] following related ideas in [53]. For instance, [51] found a relationship between the BPS spectrum and the  $U(1)_r$  charges of chiral operators. Our IR formulation of the Schur index and its generalization to line defects draws heavily from these works and their prescriptions for constructing such generating functions.

For non-conformal field theories there are two possible interpretations of the Schur index, which are not obviously equivalent. It may be a counting function of certain supersymmetric local operators, or an index for an  $S^3$  compactification of the theory. We cannot use a conformal transformation to directly relate these two perspectives. However, it is likely that although the physical theory depends on the conformal factor, the Witten index does not and thus both perspectives are valid.

The operator counting perspective leads to an immediate UV definition for Lagrangian theories. It should be straightforward to reproduce the UV calculation by localization in a judicious supersymmetric compactification of the UV Lagrangian. On the other hand, as we discuss, the sphere compactification perspective gives an intuitive motivation for the IR formulation of the index using BPS states. It would be interesting to justify the IR formula directly at the level of operator counting, perhaps by employing some map from local operators to fans of BPS particles as in [54, 55]. In any case, in [8] the agreement between the UV and IR calculations it was checked for various asymptotically free examples.

There is a strong analogy between many computations in this paper and calculations in two-dimensional (2,2) supersymmetric quantum field theories. Indeed, this analogy was a central motivation of [51] to define chamber independent combinations of BPS states. The dictionary proceeds as follows:

- The elliptic genus  $\chi(y, p; \alpha) = \operatorname{Tr}(-1)^F y^{F_L} p^{L_0} \bar{p}^{\bar{L}_0} \alpha^{J_0}$  is the  $2d$  version of the general superconformal index. It does not depend on  $\bar{p}$  and counts holomorphic operators.  
- The specialization to  $y = 1$  given as  $\chi(\alpha) = \operatorname{Tr}(-1)^F p^{L_0} \bar{p}^{\bar{L}_0} \alpha^{J_0}$  only receives contributions from chiral operators. It is analogous to the Schur index.  
- The specialization  $\chi(\alpha)$  can be computed in mass-deformed or even asymptotically free theories. The Cecotti-Vafa formula [53] expresses  $\chi(\alpha)$  in terms of the spectrum of IR BPS solitons.  
- Although  $\chi(\alpha)$  naturally counts chiral operators, in a non-conformal theory it can also be computed as the Witten index of a special twisted  $S^1$  compactification, which gives an intuitive understanding of the Cecotti-Vafa formula [53-55].

In a separate work we will discuss hybrid  $2d/4d$  systems at some length [56, 57] which involve coupling the  $4d$  theory to surface defects.

# 1.1 The IR Formula for the Schur Index

In the absence of defects, the formal IR expression for the Schur index is [8]

$$
\mathcal {I} (q) = (q) _ {\infty} ^ {2 r} \operatorname {T r} \left[ \mathcal {O} (q) \right], \tag {1.2}
$$

where  $r$  is the rank of the Coulomb branch and the various symbols are defined as:

-  $X$  stands for a quantum torus algebra

$$
X _ {\gamma} X _ {\gamma^ {\prime}} = q ^ {\frac {1}{2} \langle \gamma , \gamma^ {\prime} \rangle} X _ {\gamma + \gamma^ {\prime}}, \tag {1.3}
$$

of operators labelled by the IR charge lattice  $\Gamma$  of the  $\mathcal{N} = 2$  gauge theory.

- The charge lattice  $\Gamma$  has a flavor sub-lattice  $\Gamma_f$ . The trace  $\mathrm{Tr}$  sets to 0 all the  $X_{\gamma}$  variables which carry non-zero gauge charge, i.e., such that  $\gamma \notin \Gamma_f$ . The surviving  $X_{\gamma_f}$  variables commute and are identified with flavor fugacities.  
-  $\mathcal{O}(q)$  is defined as a product of quantum Kontsevich-Soibelman (KS) factors over the whole BPS spectrum, ordered along the phases of their central charges:

$$
\mathcal {O} (q) = \prod_ {\gamma \in \Gamma} ^ {\hat {\cap}} K \left(q; X _ {\gamma}; \Omega_ {j} (\gamma)\right). \tag {1.4}
$$

We will refer to the operator  $\mathcal{O}(q)$  as the quantum monodromy operator. An important property of this formal IR expression for the Schur index is that it is invariant under wall-crossing: the quantum KS wall-crossing formula [59, 60] guarantees that the ordered product of quantum KS factors remains unchanged across walls of marginal stability. This moduli invariance was also noted in [51] where traces of this type were first considered. Wall-crossing invariant is clearly necessary for (1.2) to make sense. The generalization to traces of higher powers of  $\mathcal{O}(q)$  and their  $4d$  interpretation have been considered in [18], and an extension to five dimensions was explored in [61].

The infrared formulation (1.2) of the index has a simple heuristic interpretation. It is the Schur index as naively computed from the infrared QED description of the Coulomb branch physics. Indeed, the factors  $(q)_{\infty}^{2}$  are Schur indices of Abelian vector multiplets while KS factor contributions of BPS particles are the Schur indices for massive hypermultiplets with appropriate spin. Finally, the trace selects the gauge invariant states which carry vanishing electric and magnetic charge.

Although this is a suggestive caricature of the physics encapsulated by the IR formula (1.2), it ignores a crucial feature: the BPS particles carry both electric and magnetic charges so there is no local field description which includes these objects as fields. This observation is intimately related to the appearance of non-commutative variables  $X_{\gamma}$  and the need for a prescribed ordering of the KS factors to resolve ambiguities.

To understand these issues it is helpful to unpack (1.4) and interpret it on the sphere. The quantum KS factors  $K(q;X_{\gamma};\Omega_{j}(\gamma))$  are graded Witten indices of Fock spaces of BPS particles of charge  $\gamma$ . We can decompose them by particle number

$$
K (q; X _ {\gamma}; \Omega_ {j} (\gamma)) = \sum_ {n \geq 0} N (q; n, \Omega_ {j} (\gamma)) X _ {n \gamma}. \tag {1.5}
$$

Thus the Schur index can be written, at least formally, as a sum over Fock spaces of BPS particles

$$
\mathcal {I} (q) = \sum_ {n _ {\gamma} \geq 0} \left[ \prod_ {\gamma \in \Gamma} N (q; n, \Omega_ {j} (\gamma)) \right] (q) _ {\infty} ^ {2 r} \operatorname {T r} \left[ \prod_ {\gamma} ^ {\cap} X _ {n _ {\gamma} \gamma} \right], \tag {1.6}
$$

The expression multiplying the Fock space degeneracies is just the Schur index of the free IR Abelian gauge theory, decorated by a collection of Abelian 't Hooft-Wilson line defects (the  $X_{\gamma}$ ) which describe the coupling of the corresponding massive BPS particles to the low-energy Abelian gauge theory. The interpretation as line defects also helps to clarify their non-commutative nature. A mutually non-local pair of line defects sources an electromagnetic field which contains angular momentum accounted for by  $q$  in (1.3).

Another advantage of viewing the BPS particles as line defects on the sphere is that we can now understand the ordering prescription in (1.4). Like BPS particles, half-BPS line defects  $L$  preserve a combination of supercharges controlled by a phase  $\vartheta$  (See Appendix A for details). When we conformally map to the sphere, this phase determines where they sit along a great circle. See Figure 1. Therefore, the ordered product simply takes into account the position of the effective line defects along this circle.

We refer to these expressions as formal because the Schur index is expected to be a power series in  $q$ . On the other hand, the crucial factor  $\mathrm{Tr}\left[\prod_{\gamma \in \Gamma}^{\widehat{\gamma}}X_{n_{\gamma}\gamma}\right]$  may produce powers of  $q$  which can be arbitrarily negative as the gauge charges increase. These negative powers are supposed to cancel out in the final answer, but in general the summation in (1.6) is conditionally convergent and hence the expression, and claimed cancellations, are ill-defined without a precise prescription for the order in which the  $n_{\gamma}$  are summed.

This important technical point limited the checks of the IR calculation of the Schur index in [8] to UV non-conformal gauge theories and non-Lagrangian examples, in special chambers only. Here we resolve this issue and in Section 2.3 explicitly evaluate the IR formula for the Schur index of  $SU(2)$  superconformal QCD. We find a precise match with

![](images/852342627301d8517e5b5830385a041a1bb1ab2aa56c44f58ca0ebc2cc2ef4f7.jpg)  
Figure 1: The geometry of line defects. When we conformally map to  $S^3 \times S^1$ , each (half) line defect  $L_i$  wraps around the  $S^1$  and sits at a point on a great circle (blue line) on the  $S^3$  according to their phases  $\vartheta_i$ . Here  $\vartheta_{ij} = \vartheta_i - \vartheta_j$ . The worldline of the line defect is shown in red.

![](images/c2f280056331fc4a78cf29c0ff09db79ea12e6f62659dce93ca509c8835625eb.jpg)

the index obtained by localization.

Specifically, we first by introduce the quantum spectrum generator  $S(q)$  defined as the product in the  $[0,\pi)$  sector of the quantum KS operator [59]

$$
\mathcal {S} (q) = \prod_ {\arg \left(Z _ {\gamma}\right) \in [ 0, \pi)} ^ {\hat {\cap}} K \left(q; X _ {\gamma}; \Omega_ {j} (\gamma)\right). \tag {1.7}
$$

Similarly the conjugate  $\bar{S}(q)$  is defined as the product in the  $[\pi, 2\pi)$  sector.

We then refine (1.4), and propose that for all well-defined  $\mathcal{N} = 2$  theories the coefficients of  $X_{\gamma}$  in the quantum spectrum generator  $\mathcal{S}(q)$  can be expanded as a power series in  $q$ , starting from a non-negative power which grows fast enough as  $\gamma$  grows so that the trace

$$
\mathcal {I} (q) = (q) _ {\infty} ^ {2 r} \operatorname {T r} \left[ \mathcal {S} (q) \bar {\mathcal {S}} (q) \right], \tag {1.8}
$$

is defined as a power series in  $q$ . We discuss further aspects of this proposal in Section 2.

# 1.2 The IR Formula for the Schur Index with Line Defects

In Section 3 we introduce line defects and study aspects of their indices from both UV and IR points of view. We concentrate on line defects inserted at a single point on  $S^3$  which

wrap the  $S^1$  (see Figure 3a). After a decompactification of the  $S^1$  and a conformal map to  $\mathbb{R}^4$ , these defects are rays which end at the origin. Indices in the presence of such a defect may be thought of heuristically as counting gauge non-invariant local operators which can absorb the charge carried by the defect.

To extend our conjecture to include the insertions of line defects is straightforward. As we have already discussed, Abelian line defects are represented in our setup by the quantum torus variables  $X_{\gamma}$ . Thus, from the IR point of view it is trivial to generalize to include additional insertions of Abelian line defects. Indeed, as described in [51] these are obtained by simply inserting a quantum torus variable  $X_{\gamma}$  into the trace. These insertions of IR line defects and their traces were computed [51] for a variety of examples.

From this point of view, the main physical input required to compute the Schur index in the presence of a line defect is the IR description of a UV line defect. In general, a given UV line defect  $L$  will map in the IR to a superposition of Abelian 't Hooft-Wilson line defects  $X_{\gamma}$ , each associated to some local ground state of the system in the presence of the UV defect. The Witten indices of these spaces of ground states are dubbed framed BPS degeneracies  $\overline{\Omega}(L, \gamma, q)$  defined in [27], and further studied in [62-68]. These degeneracies are naturally collected into a generating function  $F(L, \vartheta)$ ,

$$
F (L, \vartheta) = \sum_ {\gamma \in \Gamma} \underline {{\Omega}} (L, \vartheta , \gamma , q) X _ {\gamma}. \tag {1.9}
$$

We conjecture that the Schur index decorated by the defect  $L$  with central charge phase  $\vartheta$  can be computed by introducing  $S_{\vartheta}(q)$  defined as a the quantum spectrum generator for the phase range  $[\vartheta, \vartheta + \pi)$  and computing the trace with the generating function  $F(L, \vartheta)$  inserted

$$
\mathcal {I} _ {L} (q) = (q) _ {\infty} ^ {2 r} \operatorname {T r} \left[ F (L, \vartheta) \mathcal {S} _ {\vartheta} (q) \mathcal {S} _ {\vartheta + \pi} (q) \right]. \tag {1.10}
$$

This is well-defined if the coefficients of the spectrum generators have the expected  $q$ -expansion. As a test of this conjecture, in Section 3.4 we check directly in  $SU(2)$  superconformal QCD that this IR formula reproduces UV localization results for non-Abelian Wilson lines inserted at a point on  $S^3$  and wrapping the  $S^1$  circle.

A crucial feature of our conjecture (1.2) is that it is wall-crossing invariant. Indeed, the framed BPS states which govern the decomposition of a UV defect  $L$  jump as moduli are varied. However the framed wall-crossing formula [27] ensures that these line defect indices are invariant.

In Section 5 we use (1.10) as a tool to compute line defect indices in the Argyres-Douglas theories [69-71]. These are non-Lagrangian theories arising from special loci on the moduli space of more familiar  $\mathcal{N} = 2$  gauge theories. They can also be engineered from M5-brane compactifications [72-75]. Their BPS spectra on the moduli space have been studied

extensively in [27, 51, 76-79]. The superconformal indices of these theories were considered in [7,8,12,80,81]. The application to these examples illustrates the power of the conjecture: although these models are strongly interacting, in the IR it is possible to reconstruct the properties of their UV line operators using their known framed BPS spectra.

We conclude our discussion of line defects in Section 6 by discussing a variety of experimental connections between line defect indices and chiral algebras. This relationship was first found in [51], where it was pointed out that the insertion of IR line defects into the trace can yield characters of chiral algebras, and moreover that there is a connection between IR line defect OPEs and the  $2d$  Verlinde algebra. Indeed, this connection between chiral algebras and traces of the quantum monodromy was an important clue toward formulating the conjecture of [8] relating the Schur index and the BPS spectrum.

We generalize these observations using our calculations of UV line defect indices. In all models where we have obtained explicit expressions, the Schur index in the presence of a UV line defect produces a sum of characters of the chiral algebra associated to the  $4d$  theory. Moreover, we also find that for Argyres-Douglas theories the UV line defect operator product expansions, when inserted into the Schur index, reduce to the associated Verlinde algebra in the  $q \to 1$  limit as anticipated in [51].

As a specific example of these results, consider  $SU(2)$  superconformal QCD and let  $L$  be a half Wilson line in the doublet. According to [13], the chiral algebra is the affine Kac-Moody algebra  $\widehat{so(8)}_{-2}$ , and we find that line defect index  $\mathcal{I}_L(q)$  can be written as the following linear combination of characters of  $\widehat{so(8)}_{-2}$ ,

$$
\mathcal {I} _ {L} (q) = \sum_ {k = 1} ^ {\infty} (- 1) ^ {k} q ^ {\frac {k ^ {2} + k - 1}{2}} \left(1 - q ^ {k}\right) \chi_ {[ - 2 k - 1, 2 k - 1, 0, 0, 0 ]} (q), \tag {1.11}
$$

where  $\chi_{[a_0,a_1,a_2,a_3,a_4]}(q)$  is the affine character of  $\widehat{so(8)}_{-2}$  with affine Dynkin labels  $[a_0,a_1,a_2,a_3,a_4]$ , and we have normalized the  $\widehat{so(8)}_{-2}$  affine characters to start from order  $q^0$ .

We are able to explain aspects of these results in theories, like  $SU(2)$  superconformal QCD, which are continuously connected to free theories, but leave a complete explanation of these phenomena as an open problem.

# 1.3 The IR Formula for the Schur Half-Index

Finally, in Section 4 we propose an IR expression for the hemisphere index in the presence of some UV boundary condition. As before the key idea is to describe the boundary condition in the IR. Typically the effective description involves a  $3d\mathcal{N} = 2$  theory with a  $U(1)^r$  global symmetry which is coupled at the boundary to the bulk IR Abelian gauge fields. Note that the choice of  $3d\mathcal{N} = 2$  subalgebra of the  $4d$  theory also selects a phase  $\vartheta$ .

The  $3d$  index for the IR boundary degrees of freedom can be expanded in a charge basis including both electric charges and magnetically charged monopole operators. We thus obtain a collection of formal  $q$  power series  $Z_{\gamma}(q)$  labelled by a bulk charge  $\gamma$  and we collect them into a generating function

$$
Z ^ {I R} (q) [ X ] = \sum_ {\gamma \in \Gamma} Z _ {\gamma} ^ {I R} (q) X _ {\gamma}. \tag {1.12}
$$

Given this input, we conjecture that the Schur half-index in the presence of the given boundary condition is

$$
\mathcal {I I} (q) = (q) _ {\infty} ^ {r} \operatorname {T r} \left[ Z ^ {I R} (q) [ X ] \mathcal {S} _ {\vartheta + \pi} (q) \right]. \tag {1.13}
$$

Again, we demonstrate that this formula is wall-crossing invariant: as  $\vartheta$  crosses a BPS ray, the IR boundary condition changes in a known way [38] and the  $3d$  index varies in the opposite way as  $S_{\vartheta}$ . Our formula can be decorated further by line defects in an obvious manner and extended to the case of interfaces. We verify this formula for the examples of the Dirichlet and the RG boundary conditions [38] for the pure  $SU(2)$  gauge theory.

The RG boundary conditions of [38] play a key role in our formula. For any theory, the RG boundary condition has the property that it flows in the IR to simple Dirichlet boundary conditions. In particular, this means that the  $3d$  index of the RG interface theory can be interpreted as an invertible kernel which directly relates IR and UV Schur index calculations. Moreover, this kernel satisfies functional relations which imply the equality of IR and UV formulas for the bulk theory, either bare or decorated by any set of defects.

As a consequence of this we obtain a novel interpretation of the quantum spectrum generator  $S(q)$ . The Schur half-index in the presence of the RG boundary condition can be expanded in a charge basis as in (1.12), and the resulting generating function is simply  $(q)_{\infty}^{r}S(q)$ . Thus the expansion of  $S(q)$  into quantum torus variables can be interpreted as the contribution of the bulk BPS hypermultiplets to the half-index with this boundary condition. In this way the IR formula for the Schur half-index gives the abstract quantum wall-crossing formalism a direct physical meaning.

# 2 A Review of the Schur Index and its IR Formulation

# 2.1 The Schur Index

Our discussion of the Schur index follows [2,3]. In general for a superconformal field theory with flavor symmetry of rank  $n_f$  the Schur index is defined as a trace over the Hilbert space on  $S^3$ . It depends on a single universal fugacity  $q$  and may be refined to include flavor

fugacities  $z_{i}$

$$
\mathcal {I} \left(q, z _ {1}, \dots , z _ {n _ {f}}\right) = \operatorname {T r} \left[ e ^ {2 \pi i R} q ^ {\Delta - R} \prod_ {i = 1} ^ {n _ {f}} z _ {i} ^ {f _ {i}} \right]. \tag {2.1}
$$

The state operator correspondence implies that the same quantity may be computed by counting local operators. These operators are quarter-BPS (annihilated by two  $Q$ 's and two  $S$ 's) and obey the following restrictions on their quantum numbers

$$
\frac {1}{2} \left(\Delta - j _ {1} - j _ {2}\right) - R = 0, \quad r + j _ {1} - j _ {2} = 0, \tag {2.2}
$$

where  $j_{i}$  are spins for the Lorentz group,  $\Delta$  is the scaling dimension, and  $R,r$  are Cartans of the  $SU(2)_R\times U(1)_r$  symmetry.5

Note that in our definition, we have chosen a slightly unconventional fermion number

$$
(- 1) ^ {F} = e ^ {2 \pi i R}, (2. 3)
$$

compared to [2,3] in which  $(-1)^{F} = e^{2\pi i(j_{1} + j_{2})}$ . The two conventions of the Schur index are related by an  $\mathbb{Z}_2$  flavor charge insertion  $e^{2\pi i(j_1 + j_2 + R)}$ . As a consequence of the shortening conditions (2.2) the two conventions are related by shifting  $q^{\frac{1}{2}} \rightarrow -q^{\frac{1}{2}}$ .

For theories with a Lagrangian description, the Schur index may be computed by simply counting gauge invariant local operators built out of the free fields. The fact that it is an index then ensures that the result is correct even for an interacting theory. This yields a simple matrix integral expression for the index. The objects entering the expression are single letter partition functions for vector multiplets and hypermultiplets

$$
f ^ {V} (q) = - \frac {2 q}{1 - q}, \quad f ^ {\frac {1}{2} H} = - \frac {q ^ {1 / 2}}{1 - q}, \tag {2.4}
$$

as well as the plethyestic exponential

$$
P. E. [ f (q, u, z) ] = \exp \left[ \sum_ {n = 1} ^ {\infty} \frac {1}{n} f \left(q ^ {n}, u ^ {n}, z ^ {n}\right) \right]. \tag {2.5}
$$

Note that the sign in  $f^{\frac{1}{2} H}$ , compared to [2,3], comes from our choice of the fermion number  $(-1)^{F} = e^{2\pi iR}$ . We can also write

$$
P. E. \left[ f ^ {V} (q) u \right] = \left(q u; q\right) _ {\infty} ^ {2}, \quad P. E. \left[ f ^ {\frac {1}{2} H} (q) u \right] = \left(- q ^ {1 / 2} u; q\right) _ {\infty} ^ {- 1}, \tag {2.6}
$$

where the Pochhammer symbol is defined as

$$
(a; q) _ {n} = \left\{ \begin{array}{l l} 1 & n = 0, \\ \prod_ {j = 0} ^ {n - 1} \left(1 - a q ^ {j}\right) & n > 0. \end{array} \right. \tag {2.7}
$$

We also define  $(q)_{n}\equiv (q;q)_{n}$

For the a Lagrangian theory with gauge group  $G$  and matter in a representation  $\mathbf{R}$  of  $G$  and in a representation  $\mathbf{F}$  of the flavor symmetry, the Schur index is

$$
\mathcal {I} (q, z) = \int [ d u ] P. E. \left[ f ^ {V} (q) \chi_ {G} (u) + f ^ {\frac {1}{2} H} (q) \chi_ {\mathbf {R}} (u) \chi_ {\mathbf {F}} (z) \right], \tag {2.8}
$$

where  $\chi_{\alpha}$  are characters of the gauge and flavor group and  $[du]$  is the Haar measure on the maximal torus of  $G$ . Here  $z$  collectively denotes the flavor fugacities.

Strictly speaking, our discussion so far involves superconformal field theories. However, as elaborated on in the introduction the consistency of the IR formulation of the Schur index reviewed below strongly suggests that the Schur index may be defined for non-conformal  $\mathcal{N} = 2$  theories as well. When we discuss such examples in the following, we take the operator counting formula (3.5) as a working definition which applies to models with Lagrangians.

# 2.2 An IR Formula for the Schur Index

We now turn to the IR formula for the Schur index conjectured in [8]. This formulation can be made in any generic vacuum on the Coulomb branch where the theory is IR free. At such a point the theory is described by a  $U(1)^r$  gauge theory ( $r$  is called the rank). There is an integral charge lattice  $\Gamma$  which is equipped with three structures:

- A Dirac pairing  $\langle \cdot, \cdot \rangle$  which is bilinear, antisymmetric, integer-valued.  
- A linear central charge function  $\mathcal{Z}:\Gamma \to \mathbb{C}$ . The central charge function is the main output of the Seiberg-Witten solution [48, 49] of the low-energy dynamics.  
- A sublattice  $\Gamma_f$  of "flavor charges" which has zero Dirac pairing with other charges.

$$
\left(q ^ {\frac {1}{2}} z; q\right) _ {\infty} ^ {- 1} = \sum_ {n = 0} ^ {\infty} \frac {\left(q ^ {\frac {1}{2}} z\right) ^ {n}}{\left(q\right) _ {n}}. \tag {2.9}
$$

The Dirac pairing is non-degenerate on the quotient lattice  $\Gamma_g = \Gamma / \Gamma_f$  of gauge charges.

Associated to the lattice is a quantum torus algebra. For each charge vector  $\gamma \in \Gamma$  we introduce a variable  $X_{\gamma}$  which obey

$$
X _ {\gamma} X _ {\gamma^ {\prime}} = q ^ {\frac {1}{2} \langle \gamma , \gamma^ {\prime} \rangle} X _ {\gamma + \gamma^ {\prime}}. \tag {2.10}
$$

The torus algebra variables have a simple physical interpretation: they are line defects in the IR abelian gauge theory modeling infinitely massive source dyons with charge  $\gamma$ . Note that this explains the algebra of these variables as well, since a pair of dyons (the left-hand side above) sources a electromagnetic fields carrying angular momentum, while a single dyon (the right-hand side) does not. The variable  $q$  is thus a fugacity for rotations and keeps track of this difference.

In order to compute the Schur index we require knowledge of the spectrum of supersymmetric massive excitations of the low-energy theory described by the BPS states. Each massive BPS particle is a representation of the super little group which is  $SU(2)_J \times SU(2)_R$ . After factoring out the center of mass degrees of freedom the one-particle Hilbert space for the charge sector  $\gamma$  may be written as

$$
H _ {\gamma} = \left[ \left(\mathbf {2}, \mathbf {1}\right) \oplus \left(\mathbf {1}, \mathbf {2}\right) \right] \otimes h _ {\gamma}. \tag {2.11}
$$

The degeneracies we require are integers  $\Omega_{n}(\gamma)$  that are encoded in  $h_\gamma$  as

$$
\operatorname {T r} _ {h _ {\gamma}} \left[ y ^ {J} (- y) ^ {R} \right] = \sum_ {n \in \mathbb {Z}} \Omega_ {n} (\gamma) y ^ {n}. \tag {2.12}
$$

From the above physical data we can now formulate the index. We introduce the  $q$ -exponential, sometimes also called the quantum dilogarithm

$$
E _ {q} (z) = \left(- q ^ {\frac {1}{2}} z; q\right) _ {\infty} ^ {- 1} = \prod_ {i = 0} ^ {\infty} \left(1 + q ^ {i + \frac {1}{2}} z\right) ^ {- 1} = \sum_ {n = 0} ^ {\infty} \frac {\left(- q ^ {\frac {1}{2}} z\right) ^ {n}}{\left(q\right) _ {n}}. \tag {2.13}
$$

For each charge vector  $\gamma$  we then define a KS factor as

$$
K (q; X _ {\gamma}; \Omega_ {j} (\gamma)) = \prod_ {n \in \mathbb {Z}} E _ {q} ((- 1) ^ {n} q ^ {n / 2} X _ {\gamma}) ^ {(- 1) ^ {n} \Omega_ {n} (\gamma)}. \tag {2.14}
$$

Naively, we define the quantum KS operator as a product of these factors

$$
\mathcal {O} (q) = \prod_ {\gamma \in \Gamma} ^ {\hat {\cap}} K (q; X _ {\gamma}; \Omega_ {j} (\gamma)). \tag {2.15}
$$

Here the ordering in the product of non-commutative KS factors is defined using the central charge  $\mathcal{Z}$ . If  $\arg(\mathcal{Z}(\gamma_1)) < \arg(\mathcal{Z}(\gamma_2))$  then  $K(X_{\gamma_1})$  appears to the left of  $K(X_{\gamma_2})$  in the product. This definition has severe convergence issues, which we will address momentarily. Up to these issues, according to the wall-crossing formula [59], as moduli are varied the individual factors  $K(q; X_{\gamma}; \Omega_j(\gamma))$  may jump, but their formal product defining  $\mathcal{O}(q)$  is invariant.

Finally, we can extract the Schur index from these ingredients. Observe that flavor charges are those elements of  $\Gamma$  with trivial Dirac pairings. It then follows from the relations (2.10) that the associated variables  $X_{\gamma}$  are central elements of the torus algebra. We define a trace operation which projects the torus algebra onto these central flavor elements:

$$
\operatorname {T r} \left[ X _ {\gamma} \right] = \left\{ \begin{array}{l l} X _ {\gamma_ {f}} & \gamma = \gamma_ {f} \in \Gamma_ {f}, \\ 0 & \text {e l s e ,} \end{array} \right. \tag {2.16}
$$

The trace operation is then extended linearly to sums of the  $X_{\gamma}$ .

The flavor torus algebra generators in turn are identified with flavor fugacities appearing in the index according with the map between UV and IR global symmetries. If we pick some basis of flavor fugacities  $z_{a}$  and denote the corresponding components of the flavor charge as  $\gamma_{f}^{a}$ , we can write:

$$
\operatorname {T r} \left(X _ {\gamma_ {f}}\right) = z _ {\gamma_ {f}} \equiv \prod_ {a} z _ {a} ^ {\gamma_ {f} ^ {a}}. \tag {2.17}
$$

At last we can precisely state the conjecture of [8]. It states that the Schur index may be calculated from the infrared as

$$
\mathcal {I} (q, z) = (q) _ {\infty} ^ {2 r} \operatorname {T r} [ \mathcal {O} (q) ]. \tag {2.18}
$$

By construction this formula is wall-crossing invariant. In [8] it was tested against nonconformal  $SU(2)$  gauge theories (using the operator counting formula (3.5) as a UV definition) as well as non-Lagrangian Argyres-Douglas models using comparisons with chiral algebra techniques.

For instance, the simplest example of (2.18) is the case of a free hypermultiplet. The IR formula for the Schur index is (taking into account our choice of fermion number (2.3))

$$
\mathcal {I} (q, X _ {\gamma}) = E _ {q} \left(X _ {\gamma}\right) E _ {q} \left(X _ {- \gamma}\right) = \left(- q ^ {\frac {1}{2}} X _ {\gamma}; q\right) _ {\infty} ^ {- 1} \left(- q ^ {\frac {1}{2}} X _ {- \gamma}; q\right) _ {\infty} ^ {- 1}, \tag {2.19}
$$

which is equal to the UV answer (2.6) upon the identification  $\mathrm{Tr}(X_{\gamma}) = z$

# 2.2.1 Convergence of the IR Formula

Although the examples considered in [8] provide significant evidence towards the validity of this conjecture, the expression (2.18) suffers from an important problem: in general, as discussed in Section 1.1, the trace of the operator  $\mathcal{O}(q)$  is not well-defined.

This difficulty is closely related to the problem of giving a meaning to partial products of KS factors

$$
\mathcal {S} _ {V} = \prod_ {\gamma \in \Gamma_ {V}} ^ {\hat {\cap}} K (q; X _ {\gamma}; \Omega_ {j} (\gamma)) \tag {2.20}
$$

restricted to charges such that the central charge  $\mathcal{Z}$  lies in a radial sector  $V$  of the complex plane. If the sector  $V$  has width less than  $\pi$  (or has width  $\pi$  but is open on one side), the product can be interpreted as an element in a group of formal power series in  $X_{\gamma}$ , as  $\gamma$  and  $-\gamma$  never belong simultaneously to  $V$  and  $V$  is closed under addition. The coefficients in the sum are rational functions of  $q$ . This does not work if the width of  $V$  is greater than  $\pi$ .

In particular,  $\mathcal{O}(q)$  itself cannot be understood as a formal power series in  $X_{\gamma}$ , but the quantum spectrum generator  $\mathcal{S}(q)$  defined as the product in the  $[0,\pi)$  sector is well-defined, and so is the conjugate  $\bar{\mathcal{S}}(q)$  defined as the product in the  $[\pi,2\pi)$  sector. To give meaning to the conjecture (2.18) in general, we therefore write instead

$$
\mathcal {I} (q) = (q) _ {\infty} ^ {2 r} \operatorname {T r} \left[ \bar {\mathcal {S}} (q) \mathcal {S} (q) \right]. \tag {2.21}
$$

To be even more explicit about how this regulates the trace of  $\mathcal{O}(q)$ , let us fix a basis of charges  $\gamma_{i}$  for  $\Gamma$  such that  $\mathcal{Z}(\gamma_i)$  lies in the upper half-plane. We can then define a truncated version of the quantum spectrum generator  $S_N(q)$  by setting to zero all torus algebra variables  $X_{\gamma}$  such that their coefficients in this basis expansion are larger than  $N$

$$
X _ {a _ {i} \gamma_ {i}} \mapsto 0 \quad \text {i f} \quad \sum_ {i} a _ {i} > N. \tag {2.22}
$$

Plugging into (2.21) we find a truncated version of the Schur index  $\mathcal{I}_N(q)$ . Then, we conjecture that the coefficient of  $q^k$  in the Schur index  $\mathcal{I}(q)$  can be obtained from that  $\mathcal{I}_N(q)$  provided that  $N$  is sufficiently large compared to  $k$ . In particular, the limit as  $N$  tends to infinity of  $\mathcal{I}_N(q)$  is the Schur index  $\mathcal{I}(q)$ . We illustrate this procedure in Section 2.3.

Notice that we could have split the spectrum in two halves in other ways, in terms of quantum spectrum generators  $S_{\vartheta}(q)$  associated to other half planes  $[\vartheta, \vartheta + \pi)$ :

$$
\mathcal {I} (q) = \left(q\right) _ {\infty} ^ {2 r} \operatorname {T r} \left[ \mathcal {S} _ {\vartheta + \pi} (q) \mathcal {S} _ {\vartheta} (q) \right], \tag {2.23}
$$

We expect that for a physical theory the coefficient of  $X_{\gamma}$  in all possible spectrum generators will involve only non-negative, growing powers of  $q$  and the above formula to be correct for all choices of  $\vartheta$ . This appears to be a somewhat non-trivial statement, especially because general  $S_V$  for sectors of width smaller than  $\pi$  definitely involve negative powers of  $q$ .

Given this assumption, standard wall-crossing invariance will be automatic, as  $S_V(q)$  is invariant unless the central charge of BPS particles enters/exists  $V$ .

# 2.3  $SU(2)$  Gauge Theory with  $N_{f} = 4$  Flavors

A simple example to which our formalism applies is  $SU(2)$  gauge theory with  $N_{f} = 4$  hypermultiplets in the fundamental representation. This is a superconformal field theory where both the UV and IR formulas for the Schur index can be computed and compared. As we shall see these two calculations yield perfect agreement. This example is also significant because to properly evaluate the IR formula for the Schur index we require the regularization of the trace discussed in the previous Section.

The Schur index is readily computed from the formula (3.5) and the UV Lagrangian. This results in the following integral

$$
\mathcal {I} (q, z _ {i}) = \frac {1}{\pi} \int_ {0} ^ {2 \pi} d \theta \sin^ {2} \theta P. E. \left[ - 2 \frac {q}{1 - q} (e ^ {2 i \theta} + e ^ {- 2 i \theta} + 1) - \frac {q ^ {\frac {1}{2}}}{1 - q} (e ^ {i \theta} + e ^ {- i \theta}) \sum_ {i = 1} ^ {4} (z _ {i} + z _ {i} ^ {- 1}) \right].
$$

or

$$
\begin{array}{l} \mathcal {I} (q, z _ {i}) \\ = - \frac {1}{4 \pi i} \oint \frac {d u}{u} (u - u ^ {- 1}) ^ {2} \frac {(q) _ {\infty} ^ {2} (q u ^ {2} ; q) _ {\infty} ^ {2} (q u ^ {- 2} ; q) _ {\infty} ^ {2}}{\prod_ {i = 1} ^ {4} (- q ^ {\frac {1}{2}} u z _ {i} ; q) _ {\infty} (- q ^ {\frac {1}{2}} u ^ {- 1} z _ {i} ; q) _ {\infty} (- q ^ {\frac {1}{2}} u z _ {i} ^ {- 1} ; q) _ {\infty} (- q ^ {\frac {1}{2}} u ^ {- 1} z _ {i} ^ {- 1} ; q) _ {\infty}}. \tag {2.25} \\ \end{array}
$$

Here,  $z_{i}$  ( $i = 1, \dots, 4$ ) are the flavor fugacities for the  $SO(2)^{4}$  Cartan subgroup of the  $SO(8)$  flavor symmetry. Note again that the sign in front of each factor of  $q^{\frac{1}{2}}$  comes from our choice of the fermion number  $(-1)^{F} = e^{2\pi i R}$  (2.3).

To make the full  $SO(8)$  flavor symmetry manifest, we perform the following change of

$$
(q) _ {\infty} \left(q u ^ {2}; q\right) _ {\infty} \left(q u ^ {- 2}; q\right) _ {\infty} = \sum_ {j = 0} ^ {\infty} (- 1) ^ {n} \chi_ {j} (u) q ^ {j (j + 1) / 2} \tag {2.26}
$$

to decompose the vector multiplet contribution into  $SU(2)$  characters and expand the hypermultiplet contribution explicitly into powers of gauge and flavor fugacities.

variables,

$$
\eta_ {1} = z _ {1}, \quad \eta_ {2} = z _ {1} z _ {2}, \quad \eta_ {3} = \sqrt {z _ {1} z _ {2} z _ {3} z _ {4}}, \quad \eta_ {4} = \sqrt {\frac {z _ {1} z _ {2} z _ {3}}{z _ {4}}}, \qquad (2. 2 7)
$$

where now  $\eta_{i}$  ( $i = 1, \dots, 4$ ) are the flavor fugacities of  $SO(8)$  in the convention that the power of  $\eta_{i}$  is the Dynkin label of the  $i$ -node. We choose the second node to be the central one in the  $SO(8)$  Dynkin diagram. For example, the character for the  $\mathbf{8}_{v}$  is  $\eta_{1} + \frac{\eta_{2}}{\eta_{1}} + \frac{\eta_{3}\eta_{4}}{\eta_{2}} + \frac{\eta_{4}}{\eta_{3}} + \frac{\eta_{3}}{\eta_{4}} + \frac{\eta_{2}}{\eta_{3}\eta_{4}} + \frac{\eta_{1}}{\eta_{2}} + \frac{1}{\eta_{1}}$ . Incidentally, the Schur index of the  $SU(2)$  with  $N_{f} = 4$  flavors theory equals to the vacuum character of the affine Lie algebra  $\widehat{so(8)}_{-2}$  as established in [13].

The BPS spectrum of this theory has been investigated in [78, 79, 82]. The  $SU(2)$ $N_{f} = 4$

![](images/7160dd35bf93539fc4b693cd8ec8ea15e7495800ec5089fc023e74751d7a6505.jpg)  
Figure 2: The BPS quiver for the  $\mathcal{N} = 2$ $SU(2)$  gauge theory with  $N_{f} = 4$  hypermultiplets in the fundamental representation. The Dirac pairings  $\langle \gamma_i,\gamma_j\rangle$  are given by the arrows between the nodes.

theory has a nice finite chamber where the BPS spectrum consists of 12 hypermultiplets with various gauge and flavor charges, encoded in the BPS quiver in Figure 2. It is worth pointing out that such a convenient the chamber only exists upon mass deformation, which breaks the  $SO(8)$  global symmetry to an  $SU(2) \times SU(2) \times U(1) \times U(1)$  subgroup.

The 12 BPS hypermultiplets in increasing phase order are

$$
\gamma_ {1}, \gamma_ {2}, \gamma_ {1} + \gamma_ {4}, \gamma_ {1} + \gamma_ {6}, \gamma_ {2} + \gamma_ {3}, \gamma_ {2} + \gamma_ {5}, \gamma_ {1} + \gamma_ {4} + \gamma_ {6}, \gamma_ {2} + \gamma_ {3} + \gamma_ {5}, \gamma_ {3}, \gamma_ {4}, \gamma_ {5}, \gamma_ {6}. \tag {2.28}
$$

In addition there are also the antiparticles to these BPS states. The BPS spectrum is

organized into multiplets of that global symmetry, with  $(\gamma_4,\gamma_6)$  and  $(\gamma_{3},\gamma_{5})$  being doublets of the two unbroken  $SU(2)$  global symmetries.

We also have the following identification between the flavor quantum torus generators and  $\eta_{i}$ ,

$$
\eta_ {1} = X _ {\frac {1}{2} (\gamma_ {1} + \gamma_ {2} + \gamma_ {3} + 2 \gamma_ {4} + \gamma_ {6})},
$$

$$
\eta_ {2} = X _ {\frac {1}{2} \left(\gamma_ {1} + \gamma_ {2} + \gamma_ {3} + \gamma_ {4}\right)}, \tag {2.29}
$$

$$
\eta_ {3} = X _ {\frac {1}{2} (\gamma_ {1} + \gamma_ {2} + 2 \gamma_ {3} + \gamma_ {4} + \gamma_ {6})},
$$

$$
\eta_ {4} = X _ {\frac {1}{2} (2 \gamma_ {1} + 2 \gamma_ {2} + \gamma_ {3} + \gamma_ {4} + \gamma_ {5} + \gamma_ {6})}.
$$

The quantum spectrum generator is determined from the spectrum to be

$$
\begin{array}{l} \mathcal {S} (q) = E _ {q} \left(X _ {\gamma_ {1}}\right) E _ {q} \left(X _ {\gamma_ {2}}\right) E _ {q} \left(X _ {\gamma_ {1} + \gamma_ {4}}\right) E _ {q} \left(X _ {\gamma_ {1} + \gamma_ {6}}\right) E _ {q} \left(X _ {\gamma_ {2} + \gamma_ {3}}\right) E _ {q} \left(X _ {\gamma_ {2} + \gamma_ {5}}\right) \tag {2.30} \\ \times E _ {q} \left(X _ {\gamma_ {1} + \gamma_ {4} + \gamma_ {6}}\right) E _ {q} \left(X _ {\gamma_ {2} + \gamma_ {3} + \gamma_ {5}}\right) E _ {q} \left(X _ {\gamma_ {3}}\right) E _ {q} \left(X _ {\gamma_ {4}}\right) E _ {q} \left(X _ {\gamma_ {5}}\right) E _ {q} \left(X _ {\gamma_ {6}}\right). \\ \end{array}
$$

After a somewhat lengthy rearrangement, we can write

$$
\mathcal {S} (q) = \sum_ {\substack {\ell_ {1}, \dots , \ell_ {6}, \\ p _ {1}, \dots , p _ {6} = 0}} ^ {\infty} \frac {(- 1) ^ {\sum_ {i = 1} ^ {6} (\ell_ {i} + p _ {i})} q ^ {\frac {1}{2} A}}{(q) _ {\ell_ {1}} \cdots (q) _ {\ell_ {6}} (q) _ {p _ {1}} \cdots (q) _ {p _ {6}}} X _ {\sum_ {i = 1} ^ {6} a _ {i} \gamma_ {i}}, \tag{2.31}
$$

where

$$
\begin{array}{l} A \equiv (p _ {3} - p _ {4} + p _ {5} - p _ {6}) (\ell_ {1} - \ell_ {2}) + (\ell_ {3} - \ell_ {4} + \ell_ {5} - \ell_ {6}) (p _ {1} - p _ {2}) \\ + \left(\ell_ {3} - \ell_ {4} + \ell_ {5} - \ell_ {6} - 2 p _ {1} + 2 p _ {2}\right) \left(\ell_ {1} - \ell_ {2} - p _ {3} + p _ {4} - p _ {5} + p _ {6}\right) \tag {2.32} \\ + (- p _ {3} + p _ {4} - p _ {5} + p _ {6}) (p _ {1} - p _ {2}) + \sum_ {i = 1} ^ {6} (\ell_ {i} + p _ {i}), \\ \end{array}
$$

and

$$
\begin{array}{l} a _ {1} = \ell_ {1} + p _ {1} + p _ {4} + p _ {6}, \quad a _ {2} = \ell_ {2} + p _ {2} + p _ {3} + p _ {5}, \quad a _ {3} = \ell_ {3} + p _ {2} + p _ {3}, \tag {2.33} \\ a _ {4} = \ell_ {4} + p _ {1} + p _ {4}, \quad a _ {5} = \ell_ {5} + p _ {2} + p _ {5}, \quad a _ {6} = \ell_ {6} + p _ {1} + p _ {6}. \\ \end{array}
$$

To determine the Schur index from these expressions we now use the regularization procedure discussed in Section 2.2. We find that compute terms of order  $q^k$  in the Schur

index we must compute  $S_N(q)$  where  $N \geq 6k$ . For instance  $S_6(q)$  expanded to order  $q$  is

$$
\begin{array}{l} \mathcal{S}_{6}(q) = 1 - q^{\frac{1}{2}}\sum_{i = 1}^{6}X_{\gamma_{i}} + q\Big(X_{\gamma_{1} + \gamma_{2}} + \sum_{\substack{i,j = 3\\ i <   j}}^{6}X_{\gamma_{i} + \gamma_{j}} + \sum_{i = 1}^{6}X_{2\gamma_{i}} + X_{\gamma_{1} + \gamma_{2} + \gamma_{3} + \gamma_{4}} + X_{\gamma_{1} + \gamma_{2} + \gamma_{4} + \gamma_{5}} \\ + \left. X _ {\gamma_ {1} + \gamma_ {2} + \gamma_ {3} + \gamma_ {6}} + X _ {\gamma_ {1} + \gamma_ {2} + \gamma_ {5} + \gamma_ {6}} + X _ {\gamma_ {1} + \gamma_ {2} + \gamma_ {3} + \gamma_ {4} + \gamma_ {5} + \gamma_ {6}}\right) + \mathcal {O} \left(q ^ {\frac {3}{2}}\right), \tag {2.34} \\ \end{array}
$$

and  $\overline{S}_6(q)$  is simply given by replacing every  $\gamma_i$  by  $-\gamma_i$  in  $S_6(q)$ . The IR formula for the Schur index is then, to order  $q$

$$
\begin{array}{l} \mathcal {I} (q, \eta_ {i}) = (q) _ {\infty} ^ {2} \operatorname {T r} [ \overline {{\mathcal {S}}} _ {6} (q) \mathcal {S} _ {6} (q) ] \\ = 1 + q \left(4 + X _ {\gamma_ {1} + \gamma_ {2}} + X _ {- \gamma_ {1} - \gamma_ {2}} + X _ {\gamma_ {3} + \gamma_ {4}} + X _ {- \gamma_ {3} - \gamma_ {4}} + X _ {\gamma_ {3} - \gamma_ {5}} + X _ {- \gamma_ {3} + \gamma_ {5}} \right. \\ + X _ {\gamma_ {4} + \gamma_ {5}} + X _ {- \gamma_ {4} - \gamma_ {5}} + X _ {\gamma_ {3} + \gamma_ {6}} + X _ {- \gamma_ {3} - \gamma_ {6}} + X _ {\gamma_ {4} - \gamma_ {6}} + X _ {- \gamma_ {4} + \gamma_ {6}} + X _ {\gamma_ {5} + \gamma_ {6}} + X _ {- \gamma_ {5} - \gamma_ {6}} \\ + X _ {\gamma_ {1} + \gamma_ {2} + \gamma_ {3} + \gamma_ {4}} + X _ {- \gamma_ {1} - \gamma_ {2} - \gamma_ {3} - \gamma_ {4}} + X _ {\gamma_ {1} + \gamma_ {2} + \gamma_ {4} + \gamma_ {5}} + X _ {- \gamma_ {1} - \gamma_ {2} - \gamma_ {4} - \gamma_ {5}} \tag {2.35} \\ + X _ {\gamma_ {1} + \gamma_ {2} + \gamma_ {3} + \gamma_ {6}} + X _ {- \gamma_ {1} - \gamma_ {2} - \gamma_ {3} - \gamma_ {6}} + X _ {\gamma_ {1} + \gamma_ {2} + \gamma_ {5} + \gamma_ {6}} + X _ {- \gamma_ {1} - \gamma_ {2} - \gamma_ {5} - \gamma_ {6}} \\ \left. + X _ {\gamma_ {1} + \gamma_ {2} + \gamma_ {3} + \gamma_ {4} + \gamma_ {5} + \gamma_ {6}} + X _ {- \gamma_ {1} - \gamma_ {2} - \gamma_ {3} - \gamma_ {4} - \gamma_ {5} - \gamma_ {6}}\right) + \mathcal {O} (q ^ {2}), \\ = 1 + \chi_ {\mathbf {2 8}} (\eta_ {i}) q + \mathcal {O} (q ^ {2}), \\ \end{array}
$$

where we have used the relations between the flavor  $X_{\gamma}$  with the  $SO(8)$  fugacities (2.29). Here  $\chi_{\mathbf{28}}(\eta_i)$  is the character of the 28 of  $SO(8)$ . Including the higher order terms in  $\mathcal{S}(q)$ , we have computed the trace of the quantum monodromy operator to  $q^4$  order,

$$
\begin{array}{l} (q) _ {\infty} ^ {2} \operatorname {T r} [ \overline {{\mathcal {S}}} _ {2 4} (q) \mathcal {S} _ {2 4} (q) ] = 1 + \chi_ {\mathbf {2 8}} q + (\chi_ {\mathbf {1}} + \chi_ {\mathbf {2 8}} + \chi_ {\mathbf {3 0 0}}) q ^ {2} + (\chi_ {\mathbf {1}} + 2 \chi_ {\mathbf {2 8}} + \chi_ {\mathbf {3 0 0}} + \chi_ {\mathbf {3 5 0}} + \chi_ {\mathbf {1 9 2 5}}) q ^ {3} \\ + \left(2 \chi_ {1} + 3 \chi_ {2 8} + \chi_ {3 5 _ {v}} + \chi_ {3 5 _ {s}} + \chi_ {3 5 _ {c}} + 3 \chi_ {3 0 0} + \chi_ {3 5 0} + \chi_ {1 9 2 5} + \chi_ {4 0 9 6} + \chi_ {8 9 1 8}\right) q ^ {4} + \mathcal {O} \left(q ^ {5}\right) \tag {2.36} \\ \end{array}
$$

This agrees perfectly with the UV integral expression for the Schur index (2.24).

# 3 Line Defects and Their Schur Indices

In this section we study supersymmetric line defects and their indices in  $\mathcal{N} = 2$  field theories. These include 't Hooft-Wilson lines in gauge theories as well as their generalizations to non-Lagrangian field theories. See [22-30] for further background.

# 3.1 Defect Junctions and the Schur Index

The class of defects of interest can be characterized by the symmetries that they preserve. The most symmetric situation occurs when the defect is point-like in space and extended along time. It is then stabilized by the following odd and even generators:

- Four supercharges (thus the defect is half-BPS).  
- The group  $SU(2)_J$  of spatial rotations about the defect, the R-symmetry  $SU(2)_R$ , and time translations (but no other translations).

We refer to any object preserving these symmetries as a full line defect  $L$ .

Implicit in this definition is a parameter  $\zeta \in \mathbb{C}^*$  which characterizes which four supercharges are preserved by the defect. The  $U(1)_r$  symmetry (which is broken by the defect) rotates  $\zeta$  by a phase. In the special case where  $|\zeta| = 1$  we express it in terms of a phase as  $\zeta = \exp(-i\vartheta)$ . In this case the symmetry algebra of the line defect has a simple physical interpretation: it is the symmetries of a massive BPS particle at rest, where  $\vartheta$  is the phase of the central charge. This interpretation is implicit in the following.

In the special case of a conformal field theory, we can strengthen the requirements on line defects to promote the translations along the defect to a full  $SL(2,\mathbb{R})$  symmetry. Alternatively we may characterize the same objects as supersymmetric boundary condition on  $AdS_2 \times S^2$  [22, 83].

In addition to these line defects extended along time, there are other configurations of defects which will be significant to us. Specifically, it is useful to also consider line defects which extend along a ray in  $\mathbb{R}^4$  and terminate at the origin. We sometimes refer to these configurations as half line defects to distinguish them from the full line defects defined above. These half defects are also supersymmetric and the preserved supercharges are given in Appendix A. The origin where the half line defect terminates can support a variety of endpoint operators and we seek to count these in an index.

It is instructive to consider both the full and half defects on  $S^3 \times \mathbb{R}$  via conformal mapping. The latter is simplest, it marks the sphere at a single point associated to the defect  $L$ , and thereby modifies the radially quantized Hilbert space. Similarly, in the case of a full defect one modifies the  $S^3$  at two antipodal points by insertion of  $L$  and its CPT conjugate defect  $\bar{L}$ . See Figure 3 for the distinction between a full and a half line defect.

Because of this geometry, a full line defect can be thought simply as a junction of two half line defects. More generally, we can consider a junction of an arbitrary number of radial half line defects. There is however an important constraint on such junctions. As we demonstrate in Appendix A, in order to preserve supersymmetry all of the half lines must lie in a common two-plane in Euclidean space. Moreover, the angle of a given half line defect in the plane is exactly the same as the central charge phase  $\vartheta$  of the defect.

![](images/9efdec71ff38c34da619ae400bd29e9799235bd41a68c0e6fc402310159052b0.jpg)

![](images/b61b01061d0f9e46c639656ff1d6f2b1c9dd929259b5b49ff4a8e2f5c07285b6.jpg)  
(a)

![](images/0f2935eb33bfaef20a304d9d60fb949d4578fb14540b9cbb4f828ca578c0fbad.jpg)

![](images/30f9db497a96ba8322dafcfaef0f3b0344bc78434dd48552151053b42ba844d6.jpg)  
Figure 3: The conformal map (together with the compactification of  $\mathbb{R}$  to  $S^1$ ) of (a) a half line defect and (b) a full line defect from  $\mathbb{R}^4$  to  $S^3 \times S^1$ . The worldline of the line defect is shown in red.

![](images/92c421e43e0ccbb062d93b83a436e89a55559354731d3ee47a7e3fe21f17a171.jpg)  
(b)

![](images/4c9ed0be23aef723a6dca06be9ea2508cf146870cb161572da610869ce1c9077.jpg)

Conformally mapping to the sphere the defect insertions then lie along a fixed great circle. See Figure 1.

The symmetry preserved by these junctions consists of  $SU(2)_R$ , as well as a  $U(1)$  rotation in the plane transverse to the rays defining the junction. The supercharges preserved by this configuration are compatible with those use to define the Schur index (see Appendix A) and allow us to extend the definition Schur index to include these insertions [36]:

$$
\mathcal {I} _ {L _ {1} \left(\vartheta_ {1}\right) L _ {2} \left(\vartheta_ {2}\right) \dots L _ {n} \left(\vartheta_ {n}\right)} (q) = \operatorname {T r} \left[ e ^ {2 \pi i R} q ^ {\Delta - R} \right]. \tag {3.1}
$$

Here the trace is over the Hilbert space on  $S^3$  with defects  $L_{i}$  inserted at angle  $\vartheta_{i}$  along a great circle. Note that this index does not depend continuously on the parameters  $\vartheta_{i}$ , but does depend on the relative ordering of the points along the circle. In practice we will mostly focus on the case of a single half line defect insertion in the index in which case the  $\vartheta$  dependence can be suppressed.

In theories with a UV Lagrangian formulation, the localization formula (3.5) can be simply extended to include line defects [40]. Consider first the case of a Wilson line in a representation  $\mathbf{R}$  of the gauge group, and let  $\chi_{\mathbf{R}}(u)$  denote the character of this representation. The index in the presence of the half Wilson line  $L_{\mathbf{R}}$  is then

$$
\mathcal {I} _ {L _ {\mathbf {R}}} (q, z) = \int [ d u ] \chi_ {\mathbf {R}} (u) Z (q, u, z), \tag {3.2}
$$

where  $Z(q,u,z)$  is the integrand in the absence of the line defect.

The above may be readily generalized to include multiple half Wilson lines in general representations. In this case the various half lines are all mutually local and thus their relative ordering on the great circle in  $S^3$  does not effect the index. To add these to the index we simply include a separate character factor  $\chi_{\mathbf{R}}$  for each of the half Wilson lines. In particular, for the specific case of two half lines, which is equivalent to the insertion of a full unbroken line defect in a representation  $\mathbf{R}$ , we add the character of the representation  $\mathbf{R}$ , associated to a defect at the north pole of  $S^3$ , and the character of  $\overline{\mathbf{R}}$ , associated to the defect insertion at the south pole.

One interesting aspect of the localization formula (3.2) is that it gives a more intuitive description of what exactly is being counted in the line defect index. In the absence of the character, the integral (3.2) counts gauge invariant local operators satisfying the Schur shortening conditions (2.2). With the character  $\chi_{\mathbf{R}}(u)$  it counts "gauge non-invariant local operators" (i.e. words in the free field variables) which satisfy the same shortening conditions and transform in the representation  $\mathbf{R}$ . Indeed, these are exactly the objects that may end on the defect and absorb its charge.

We can generalize from junctions of Wilson lines to a localization formula for the Schur index in the presence of a 't Hooft-Wilson line half-defect  $L$  [36]. This requires introducing some new notations.

It is useful to interpret the integral expression for the Schur index in terms of an inner product in a space of functions of gauge fugacities and magnetic charges, denoted as

$$
(\mathcal {A}, \mathcal {B}) \equiv \sum_ {\vec {m}} \int [ d u ] _ {\vec {m}} \mathcal {A} _ {- \vec {m}} (u) \mathcal {B} _ {\vec {m}} (u), \tag {3.3}
$$

where  $[du]_{\vec{m}}$  is a certain shifted Haar measure with magnetic charge  $\vec{m}$ . The usual Schur index is written as an inner product

$$
\mathcal {I} (q, z) = \left(\Pi^ {N}, \Pi^ {S}\right) \tag {3.4}
$$

of two half-indices  $\Pi_{\vec{m}}^{N,S}(q,u,z)$  associated with the two hemispheres.

Concretely, the gauge theory integrand is evenly distributed between the two half-indices

$$
\Pi_ {\vec {m}} ^ {N, S} (q, u, z) = \delta_ {\vec {m}, 0} P. E. \left[ \frac {1}{2} f ^ {V} (q) \chi_ {G} (u) + f ^ {\frac {1}{2} H} (q) \chi_ {\mathbf {R} \times \mathbf {F}} ^ {N, S} (u, z) \right], \tag {3.5}
$$

where we pick any Lagrangian splitting of the hypermultiplets into two sets of half-hypermultiplets with characters  $\chi_{\mathbf{R}\times \mathbf{F}}^{N,S}(u,z)$ .

The Schur index decorated by a line defect  $L$  can be written as follows:

$$
\mathcal {I} _ {L} (q, z) = \left(\Pi^ {N}, \hat {O} _ {L} \Pi^ {S}\right), \tag {3.6}
$$

where  $\hat{O}_L$  is a certain difference operator acting on functions of  $q, u, \vec{m}$ . The specific form of  $\hat{O}_L$  follows from localization computations as in [46]. It can also be obtained with the help of the AGT correspondence [44, 84] (see also the relation to quantization of the Coulomb branch of the circle-compactified theory [85, 86]). To include more half line defects, we include more difference operators  $\hat{O}_{L_i}$ ,

$$
\mathcal {I} _ {(L _ {i})} (q, z) = \left(\Pi^ {N}, \hat {O} _ {L _ {1}} \dots \hat {O} _ {L _ {n}} \Pi^ {S}\right), \tag {3.7}
$$

Now the order of insertion on the circle matters and is captured by the order in which the operators act.

In this formula, the line defects are inserted along a southern quarter of the great circle which goes from the equator to the south pole of the three-sphere. There is a second set of difference operators,  $\hat{O}_L^{\prime}$ , which represents an insertion along the other southern quarter of the great circle and commute with the first set. It is convenient to represent the action of the second set of operators from the right, writing

$$
\mathcal {I} _ {(L _ {i}), (L _ {i} ^ {\prime})} (q, z) = \left(\Pi^ {N}, \hat {O} _ {L _ {1}} \dots \hat {O} _ {L _ {n}} \Pi^ {S} \hat {O} _ {L _ {1}} ^ {\prime} \dots \hat {O} _ {L _ {n ^ {\prime}}} ^ {\prime}\right), \tag {3.8}
$$

The inner product is defined in such a way that

$$
(\mathcal {A}, \hat {\mathcal {O}} _ {L} \mathcal {B}) = (\mathcal {A} \hat {\mathcal {O}} _ {L} ^ {\prime}, \mathcal {B}), \quad (\mathcal {A}, \mathcal {B} \hat {\mathcal {O}} _ {L} ^ {\prime}) = (\hat {\mathcal {O}} _ {L} \mathcal {A}, \mathcal {B}), \tag {3.9}
$$

while the half-indices satisfy

$$
\hat {O} _ {L} \Pi^ {S} = \Pi^ {S} \hat {O} _ {L} ^ {\prime}, \quad \hat {O} _ {L} \Pi^ {N} = \Pi^ {N} \hat {O} _ {L} ^ {\prime}. \tag {3.10}
$$

These relations encode the fact that the location of a line defect can be moved freely along the great circle.

In an Abelian theory, writing the gauge fugacity as  $u$ , the 't Hooft-Wilson lines are

represented by monomials in the operators

$$
x = q ^ {\frac {m}{2}} u, \quad p = (m \rightarrow m + 1, u \rightarrow q ^ {\frac {1}{2}} u), \quad x ^ {\prime} = q ^ {- \frac {m}{2}} u, \quad p ^ {\prime} = (m \rightarrow m + 1, u \rightarrow q ^ {- \frac {1}{2}} u). \tag {3.11}
$$

We can map a function  $\mathcal{A}_m(u)$  to a generating function

$$
\mathcal {A} [ X ] = \sum_ {m \in \mathbb {Z}}: \mathcal {A} _ {m} \left(X _ {\gamma}\right) X _ {- m \gamma^ {\prime}}: \tag {3.12}
$$

where:  $X_{a}X_{b}\coloneqq X_{a + b}$  and  $\langle \gamma^{\prime},\gamma \rangle = 1$  . Then

$$
(x \mathcal {A}) [ X ] = X _ {\gamma} \mathcal {A} [ X ], \qquad (p \mathcal {A}) [ X ] = X _ {\gamma^ {\prime}} \mathcal {A} [ X ],
$$

$$
\left(x ^ {\prime} \mathcal {A}\right) [ X ] = \mathcal {A} [ X ] X _ {\gamma}, \quad \left(p ^ {\prime} \mathcal {A}\right) [ X ] = \mathcal {A} [ X ] X _ {\gamma^ {\prime}}. \tag {3.13}
$$

Furthermore,  $\Pi^{N,S} = \delta_{m,0}$  for a pure Abelian gauge theory. This explains why Abelian Schur index calculations can be expressed as traces over the quantum torus algebra.

As with the Schur index, the localization formulas, (3.2)-(3.6) and their interpretation as traces strictly speaking apply only to conformal field theories. Consistency of the conjectures to follow strongly suggests that these concepts have a more universal definition with the localization formulas describing Lagrangian non-conformal theories as well. We thus continue to apply these formulas to non-conformal systems.

# 3.2 Examples of the Line Defect Schur Index

In this subsection we compute various line defect indices using the UV localization formula in the pure  $SU(2)$  gauge theory and in the  $SU(2)$  superconformal QCD. Our methods follow directly from [36, 40]. We will later reproduce these line defect indices from an IR calculation in Section 3.3 and Section 3.4.

# 3.2.1 Wilson Lines in  $SU(2)$  Gauge Theory

Let us consider various explicit calculations of the index in the presence of half line defects in  $SU(2)$  gauge theory. The Schur index of the pure  $SU(2)$  gauge theory with a half Wilson line defect  $L_{0,n}$  (in the representation of dimension  $n + 1$ ) is

$$
\mathcal {L} _ {L _ {0, n}} (q) = \frac {1}{\pi} \int_ {0} ^ {2 \pi} d \theta \sin^ {2} \theta \left(\frac {e ^ {i (n + 1) \theta} - e ^ {- i (n + 1) \theta}}{e ^ {i \theta} - e ^ {- i \theta}}\right) P. E. \left[ - \frac {2 q}{1 - q} (e ^ {2 i \theta} + e ^ {- 2 i \theta} + 1) \right]. (3. 1 4)
$$

If  $n$  is odd,  $\mathcal{I}_{L_{0,n}}(q)$  vanishes, while if  $n$  is even we obtain non-trivial results. For example,

$$
\mathcal {I} _ {L _ {0, 2}} = - 2 q + q ^ {2} - 2 q ^ {4} + q ^ {6} - 2 q ^ {9} + q ^ {1 2} + \mathcal {O} \left(q ^ {1 3}\right),
$$

$$
\mathcal {I} _ {L _ {0, 4}} = q ^ {2} + 2 q ^ {3} - 2 q ^ {4} + q ^ {6} + 2 q ^ {7} - 2 q ^ {9} + q ^ {1 2} + \mathcal {O} \left(q ^ {1 3}\right), \tag {3.15}
$$

$$
\mathcal {I} _ {L _ {0, 6}} = - 2 q ^ {4} - q ^ {6} + 2 q ^ {7} - 2 q ^ {9} - 2 q ^ {1 1} + q ^ {1 2} + \mathcal {O} \left(q ^ {1 3}\right).
$$

# 3.2.2 't Hooft Lines in  $SU(2)$  Gauge Theory

Let us consider the Schur index with the presence of a full 't Hooft line with minimal magnetic charge in the pure  $SU(2)$  gauge theory. The magnetic charge  $m$  rangess over non-negative half integers,  $0, \frac{1}{2}, 1, \ldots$ . The shifted Haar measure  $[du]_m$  is

$$
[ d u ] _ {m} = \frac {1}{2 \pi} \left(1 - \frac {1}{2} \delta_ {m, 0}\right) q ^ {- m} (1 - q ^ {m} e ^ {2 i \theta}) (1 - q ^ {m} e ^ {- 2 i \theta}), \tag {3.16}
$$

where we have written  $u = e^{i\theta}$ . The half-index  $\Pi_{m}$  is

$$
\Pi_ {m} ^ {N, S} (q, \theta) = \delta_ {m, 0} P. E. \left[ - \frac {q}{1 - q} (e ^ {2 i \theta} + e ^ {- 2 i \theta} + 1) \right]. \tag {3.17}
$$

The difference operator  $\hat{O}_{1,0}$  for a 't Hooft line with minimal magnetic charge can be read off from that of the  $\mathcal{N} = 4$ $SU(2)$  gauge theory by decoupling the adjoint hypermultiplet [36]. In terms of the Abelian difference operators  $x,p$  defined above it is

$$
\hat {O} _ {1, 0} = \frac {i}{x - x ^ {- 1}} \left(p ^ {\frac {1}{2}} - p ^ {- \frac {1}{2}}\right). \tag {3.18}
$$

The Schur index with a full 't Hooft line is

$$
\mathcal {I} _ {L _ {1, 0}} (q) = \sum_ {m} \int [ d u ] _ {m} \Pi_ {- m} ^ {N} (q, \theta) (\hat {O} _ {1, 0}) ^ {2} \Pi_ {m} ^ {S} (q, \theta). \tag {3.19}
$$

Since  $\Pi_{-m}^{N}$  is nonzero only when  $m = 0$ , it suffices to compute the term in  $\left(\hat{O}_{1,0}\right)^2\Pi_m(q,\theta)$  that comes with  $\delta_{m,0}$ ,

$$
\left(\hat {O} _ {1, 0}\right) ^ {2} \Pi_ {m} (q, \theta) = - \delta_ {m, 0} \frac {q ^ {\frac {1}{2}} (1 + q)}{(1 - q e ^ {2 i \theta}) (1 - q e ^ {- 2 i \theta})} \Pi_ {0} (q, \theta) + \dots . \tag {3.20}
$$

It follows that

$$
\begin{array}{l} \mathcal {I} _ {L _ {1, 0}} (q) = - \frac {1}{\pi} \int_ {0} ^ {2 \pi} d \theta \sin^ {2} \theta \frac {q ^ {\frac {1}{2}} (1 + q)}{(1 - q e ^ {2 i \theta}) (1 - q e ^ {- 2 i \theta})} P. E. \left[ - \frac {2 q}{1 - q} \left(e ^ {2 i \theta} + e ^ {- 2 i \theta} + 1\right) \right] \tag {3.21} \\ = - q ^ {\frac {1}{2}} + q ^ {\frac {5}{2}} - q ^ {\frac {7}{2}} - q ^ {\frac {9}{2}} + q ^ {\frac {1 3}{2}} + q ^ {\frac {1 5}{2}} - q ^ {\frac {1 7}{2}} - q ^ {\frac {1 9}{2}} + \dots . \\ \end{array}
$$

# 3.2.3 Wilson Lines in  $SU(2)$  Gauge Theory with  $N_{f} = 4$  Flavors

In the conformal case  $N_{f} = 4$ , the half Wilson line defect  $L_{0,n}$  in the  $(n + 1)$ -dimensional irreducible representation of  $SU(2)$  is [40]

$$
\begin{array}{l} \mathcal {I} _ {L _ {0, n}} (q, z _ {i}) = \frac {1}{\pi} \int_ {0} ^ {2 \pi} d \theta \sin^ {2} \theta \left(\frac {e ^ {i (n + 1) \theta} - e ^ {- i (n + 1) \theta}}{e ^ {i \theta} - e ^ {- i \theta}}\right) (3.22) \\ \times P. E. \left[ - 2 \frac {q}{1 - q} \left(e ^ {2 i \theta} + e ^ {- 2 i \theta} + 1\right) - \frac {q ^ {\frac {1}{2}}}{1 - q} \left(e ^ {i \theta} + e ^ {- i \theta}\right) \sum_ {i = 1} ^ {4} \left(z _ {i} + z _ {i} ^ {- 1}\right) \right], (3.22) \\ \end{array}
$$

where  $z_{i}$  ( $i = 1, \dots, 4$ ) are the flavor fugacities for  $SO(2)^{4}$ . Cartan subgroup of the  $SO(8)$  flavor symmetry. They are related to the  $SO(8)$  fugacities  $\eta_{i}$  ( $i = 1, \dots, 4$ ) by (2.27). The first few terms of the doublet half Wilson line defect index are

$$
\mathcal {I} _ {L _ {0, 1}} (q, \eta_ {i}) = - \chi_ {[ 1, 0, 0, 0 ]} q ^ {\frac {1}{2}} - \chi_ {[ 1, 1, 0, 0 ]} q ^ {\frac {3}{2}} - \left(\chi_ {[ 1, 0, 0, 0 ]} + \chi_ {[ 0, 0, 1, 1 ]} + \chi_ {[ 1, 1, 0, 0 ]} + \chi_ {[ 1, 2, 0, 0 ]}\right) q ^ {\frac {5}{2}} \dots , \tag {3.23}
$$

where  $\chi_{[a_1,a_2,a_3,a_4]}(\eta_i)$  is the (finite)  $SO(8)$  character for the representation  $[a_1,a_2,a_3,a_4]$ .

It is instructive to enumerate the operators that are counted by the doublet half Wilson line defect index (3.23) in the free limit of zero gauge coupling. The defect operators living at the end of a half Wilson line transform as doublets under the gauge  $SU(2)$  and satisfy the Schur operator conditions. For the  $SU(2)$ $N_{f} = 4$  theory at the free coupling point, the single-letter operators that contribute to the Schur index are the complex scalars  $H^{ia}$  of the 8 half-hypermultiplets, 2 components of the gauginos  $\rho^{\pm A}$  (the  $\pm$  denotes their  $U(1)_r$  charges) in the vector multiplet, and 1 derivative  $\partial \equiv \partial_{-\dot{+}}$  [3]. Here  $i = 1,\dots ,8$  is the index for the  $\mathbf{8}_v$  of the flavor  $SO(8)$ , while  $a = 1,2$  and  $A = 1,2,3$  are the indices for the 2 and 3 of the gauge  $SU(2)$ , respectively. The representations of the single-letter Schur operators under the flavor  $SO(8)$  and gauge  $SU(2)$  and their contributions to the Schur index are summarized below

<table><tr><td></td><td>SO(8)</td><td>SU(2)</td><td>Schur index</td></tr><tr><td>Hia</td><td>8v</td><td>2</td><td>q1/2</td></tr><tr><td>ρ±A</td><td>1</td><td>3</td><td>-q</td></tr><tr><td>δ</td><td>1</td><td>1</td><td>q</td></tr></table>

At the zero coupling point, the line defect index (3.23) is counting all the composite operators of the single-letter Schur operators that are in the doublet of the gauge  $SU(2)$ .

At  $q^{\frac{1}{2}}$  order in the line defect index (3.23), the only contributing operator is  $H^{ia}$ , which transforms in the  $\mathbf{8}_v$  of the flavor  $SO(8)$ . Hence the coefficient of  $q^{\frac{1}{2}}$  in the defect index is minus the (finite)  $SO(8)$  character of of  $\mathbf{8}_v$ . The sign comes from our our choice of the fermion number  $(-1)^F = e^{2\pi iR}$  and the fact that  $H^{ia}$  has  $R = 1/2$ . At  $q^{\frac{3}{2}}$  order, the doublet Schur operators are

<table><tr><td></td><td>SO(8)</td><td>SU(2)</td></tr><tr><td>∂Hia</td><td>8v</td><td>2</td></tr><tr><td>(ρ±A Hib)(TA)a b</td><td>8v</td><td>2</td></tr><tr><td>HiaHjbHkc</td><td>8v ⊕ 160v</td><td>2</td></tr></table>

(3.25)

where  $(T^A)_b$  are the generators of  $SU(2)$  in the doublet representation, i.e., the Pauli matrices. Adding their contributions together, we indeed find the  $q^{\frac{3}{2}}$  coefficient to be minus the character of  $\mathbf{160}_v$  of  $SO(8)$ .

We also record the defect indices  $\mathcal{I}_{L_{0,n}}$  of half Wilson lines in the  $(n + 1)$ -dimensional representation with flavor fugacities set to be 1,

$$
\begin{array}{l} \mathcal {I} _ {L _ {0, 1}} (q, \eta_ {i} = 1) = - \left(8 q ^ {\frac {1}{2}} + 1 6 0 q ^ {\frac {3}{2}} + 1 6 2 4 q ^ {\frac {5}{2}} + 1 1 7 6 8 q ^ {\frac {7}{2}} + 6 8 3 7 6 q ^ {\frac {9}{2}} + 3 3 9 4 0 8 q ^ {\frac {1 1}{2}} + 1 4 9 3 0 6 4 q ^ {\frac {1 3}{2}} \right. \\ \left. + 5 9 6 5 1 9 2 q ^ {\frac {1 5}{2}} + 2 2 0 1 5 9 3 6 q ^ {\frac {1 7}{2}} + 7 6 0 0 7 9 0 4 q ^ {\frac {1 9}{2}}\right) + \mathcal {O} \left(q ^ {\frac {2 1}{2}}\right), \\ \end{array}
$$

$$
\begin{array}{l} \mathcal {I} _ {L _ {0, 2}} (q, \eta_ {i} = 1) = 3 4 q + 5 6 7 q ^ {2} + 5 2 3 6 q ^ {3} + 3 5 4 7 6 q ^ {4} + 1 9 6 0 7 2 q ^ {5} + 9 3 5 3 3 4 q ^ {6} + 3 9 8 2 5 9 8 q ^ {7} \\ + 1 5 4 8 0 6 1 8 q ^ {8} + 5 5 8 0 4 0 0 8 q ^ {9} + 1 8 8 7 3 8 9 7 8 q ^ {1 0} + 6 0 4 2 6 9 7 5 8 q ^ {1 1} + 1 8 4 4 1 8 2 0 6 3 q ^ {1 2} \\ + 5 3 9 5 2 1 2 2 7 2 q ^ {1 3} + 1 5 1 9 9 5 8 2 9 3 9 q ^ {1 4} + \mathcal {O} (q ^ {1 5}), \\ \end{array}
$$

$$
\begin{array}{l} \mathcal {I} _ {L _ {0, 3}} (q, \eta_ {i} = 1) = - \left(1 0 4 q ^ {\frac {3}{2}} + 1 5 6 0 q ^ {\frac {5}{2}} + 1 3 5 1 2 q ^ {\frac {7}{2}} + 8 7 3 1 2 q ^ {\frac {9}{2}} + 4 6 5 0 7 2 q ^ {\frac {1 1}{2}} + 2 1 5 2 5 8 4 q ^ {\frac {1 3}{2}} \right. \\ + 8 9 3 6 5 3 6 q ^ {\frac {1 5}{2}} + 3 3 9 9 0 7 0 4 q ^ {\frac {1 7}{2}} + 1 2 0 2 3 2 2 1 6 q ^ {\frac {1 9}{2}} + 3 9 9 9 0 8 1 2 0 q ^ {\frac {2 1}{2}} \\ + 1 2 6 1 4 0 1 3 6 0 q ^ {\frac {2 3}{2}}) + \mathcal {O} \left(q ^ {\frac {2 5}{2}}\right), \tag {3.26} \\ \end{array}
$$

$$
\begin{array}{l} \mathcal {I} _ {L _ {0, 4}} (q, \eta_ {i} = 1) = 2 5 9 q ^ {2} + 3 6 3 4 q ^ {3} + 3 0 1 1 2 q ^ {4} + 1 8 7 9 9 4 q ^ {5} + 9 7 4 0 5 0 q ^ {6} + 4 4 0 5 0 1 4 q ^ {7} + 1 7 9 2 8 4 9 8 q ^ {8} \\ + 6 7 0 2 3 5 0 2 q ^ {9} + 2 3 3 4 8 4 0 5 0 q ^ {1 0} + 7 6 6 0 7 8 4 8 6 q ^ {1 1} + 2 3 8 6 8 6 0 9 5 5 q ^ {1 2} + \mathcal {O} (q ^ {1 3}). \\ \end{array}
$$

# 3.3 An IR Formula for the Line Defect Schur Index

In this section, we generalize the IR formula for the Schur index to include line defects. The basic intuition is easy to explain. The IR formula for the Schur index (2.21) can be interpreted as the index of an Abelian gauge theory with independent fields for each BPS hypermultiplet. It is straightforward to generalize such a formula to include infrared line

defects: we simply include the expression  $X_{\gamma}$  in the trace. Thus to extract a correct UV line defect index all that we require is the data of how a UV line defect  $L$  is decomposed in the IR into a sum of Abelian defects  $X_{\gamma}$ . As we now review, this is precisely the data captured by framed BPS states.

Consider a full line defect  $L$  at a point in space and extended along time. On the Coulomb branch of the theory,  $L$  modifies the Hilbert space and there is a new class of BPS states, so-called framed BPS states, which may be viewed as ordinary particles bound to the defect. The framed BPS Hilbert space  $\mathcal{H}_L$  is graded by electromagnetic charges

$$
\mathcal {H} _ {L} = \bigoplus_ {\gamma \in \Gamma} \mathcal {H} _ {L, \gamma}. \tag {3.27}
$$

We encode the framed BPS states in a framed protected spin character

$$
\underline {{\Omega}} (L, \gamma , q) = \mathrm {T r} _ {\mathcal {H} _ {L, \gamma}} q ^ {J} (- q) ^ {R}, \tag {3.28}
$$

The framed BPS Hilbert spaces, as well as the framed protected spin characters, jump at walls of marginal stability.

The framed BPS states also characterize the defect renormalization group flow. In the infrared, the defect  $L$  is described by a collection of defects in the IR abelian gauge theory. These are simply dyons characterized by their electromagnetic charges and are represented exactly by the quantum torus variables  $X_{\gamma}$ . The  $X_{\gamma}$  and the OPE (1.3) provide a convenient way of describing the ultraviolet line defects in terms of the infrared data. For each  $L$  we introduce the generating function [27, 65]

$$
F (L, \vartheta) = \sum_ {\gamma \in \Gamma} \underline {{\Omega}} (L, \gamma , q) X _ {\gamma}. \tag {3.29}
$$

The above expression specifies how the ultraviolet defect  $L$  is decomposed into infrared pieces. Again as a consequence of wall crossing, this decomposition will jump.

With these ingredients, we now formulate our conjecture for the Schur index in the presence of line defect  $L$  with central charge phase  $\vartheta$ . It is simply a modified trace with the defect generating function inserted at the appropriate phase, $^{10}$

$$
\mathcal {I} _ {L} (q) = (q) _ {\infty} ^ {2 r} \operatorname {T r} \left[ F (L, \vartheta) \mathcal {S} _ {\vartheta} (q) \mathcal {S} _ {\vartheta + \pi} (q) \right], \tag {3.30}
$$

where  $S_{\vartheta}(q)$  is the quantum spectrum generator associated to the half plane  $[\vartheta, \vartheta + \pi)$ .

As advocated before, since the line defect index is originally defined in the UV, its IR

formula (3.30) should be framed wall-crossing invariant. The framed wall-crossing phenomenon occurs when the central charge  $Z_{\gamma} = |Z_{\gamma}|e^{i\vartheta_{\gamma}}$  of an ordinary BPS state crosses the central charge of the line defect  $\zeta$  as we vary the moduli. The framed wall-crossing formula for the generating function  $F(L,\vartheta)$  is [27]

$$
F (L, \vartheta_ {\gamma} - \epsilon) = E _ {q} \left(X _ {\gamma}\right) F (L, \vartheta_ {\gamma} + \epsilon) E _ {q} \left(X _ {\gamma}\right) ^ {- 1}, \tag {3.31}
$$

where  $\epsilon > 0$ . As the central charge  $Z_{\gamma}$  of a BPS state crosses the line defect central charge  $\zeta$  from above, the quantum spectrum generator  $S_{\vartheta}(q)$  also jumps

$$
\mathcal {S} _ {\vartheta_ {\gamma} - \epsilon} (q) = E _ {q} \left(X _ {\gamma}\right) S _ {\vartheta_ {\gamma} + \epsilon} (q) E _ {q} \left(X _ {- \gamma}\right) ^ {- 1}. \tag {3.32}
$$

Similarly, the quantum spectrum generator for the other half plane  $S_{\vartheta + \pi}(q)$  also jumps discontinuously,  $S_{\vartheta_{\gamma} + \pi - \epsilon}(q) = E_q(X_{-\gamma}) S_{\vartheta_{\gamma} + \pi + \epsilon}(q) E_q(X_{\gamma})^{-1}$ . Combining the above transformations, we conclude that our IR formula (3.30) for the Schur index is framed wall-crossing invariant,

$$
(q) _ {\infty} ^ {2 r} \operatorname {T r} \left[ F (L, \vartheta) \mathcal {S} _ {\vartheta} (q) \mathcal {S} _ {\vartheta + \pi} (q) \right] \Big | _ {\vartheta = \vartheta_ {\gamma} - \epsilon} = (q) _ {\infty} ^ {2 r} \operatorname {T r} \left[ F (L, \vartheta) \mathcal {S} _ {\vartheta} (q) \mathcal {S} _ {\vartheta + \pi} (q) \right] \Big | _ {\vartheta = \vartheta_ {\gamma} + \epsilon}. \quad (3. 3 3)
$$

More generally, for an arbitrary junction of half line defects  $L_{i}$  with central charge phases  $\vartheta_{1} < \vartheta_{2} < \dots < \vartheta_{n}$  all chosen to be on the right of the ordinary BPS particles and to the left of the anti-particles, we propose the following IR formula for the line defect index (3.1),

$$
\mathcal {I} _ {L _ {1} \left(\vartheta_ {1}\right) L _ {2} \left(\vartheta_ {2}\right) \dots L _ {n} \left(\vartheta_ {n}\right)} (q) = (q) _ {\infty} ^ {2 r} \operatorname {T r} \left[ F \left(L _ {1}, \vartheta_ {1}\right) F \left(L _ {2}, \vartheta_ {2}\right) \dots F \left(L _ {n}, \vartheta_ {n}\right) \mathcal {S} _ {\vartheta_ {n}} (q) \mathcal {S} _ {\vartheta_ {n} + \pi} (q) \right]. \tag {3.34}
$$

The IR formula for the more general central charge phases  $\vartheta_{i}$  assignment can be obtained by applying the framed wall-crossing formula.

# 3.4 Examples of the IR Formula

In this subsection we will apply the IR formula introduced in Section 3.3 to various line defects in the pure  $SU(2)$  gauge theory and in the  $SU(2)$  gauge theory with  $N_{f} = 4$  flavors. In particular, we will reproduce the line defect indices obtained in Section 3.2 from the trace of the quantum KS operator and the generating function  $F(L,\vartheta)$  of the line defect.

# 3.4.1 Wilson Lines in  $SU(2)$  Gauge Theory

We begin by testing our IR formula in the case of Wilson lines in the pure  $SU(2)$  gauge theory. Let  $L_{0,n}$  denote the Wilson line in the irreducible  $n + 1$ -dimensional representation

![](images/fa0a97cfeead7f89aa93910b437e7102ad09780d05eca148e25ff0bc0ebeb8df.jpg)  
Figure 4: The BPS quiver for the  $\mathcal{N} = 2$  pure  $SU(2)$  gauge theory.

of  $SU(2)$ . To describe the resulting generating functions  $F(L_{0,n},\vartheta)$ , it is helpful to introduce centrally symmetric  $q$ -binomial coefficients as [65]

$$
\binom {m} {r} _ {q} \equiv q ^ {- \frac {1}{2} r (m - r)} \frac {(1 - q ^ {m}) (1 - q ^ {m - 1}) \cdots (1 - q ^ {m - r + 1})}{(1 - q ^ {r}) (1 - q ^ {r - 1}) \cdots (1 - q)}. \qquad \qquad (3. 3 5)
$$

Then (assuming the defect phase  $\vartheta$  is chosen to lie to the right of all ordinary BPS states) the generating functions are

$$
F \left(L _ {0, n}, \vartheta\right) = \sum_ {s = 0} ^ {n} \sum_ {r = 0} ^ {s} \binom {n - r} {n - s} _ {q} \binom {s} {r} _ {q} X _ {\left(r - \frac {n}{2}\right) \gamma_ {1} + \left(s - \frac {n}{2}\right) \gamma_ {2}}. \tag {3.36}
$$

Here, charges  $\gamma_{1}$  and  $\gamma_{2}$  are chosen such that  $\langle \gamma_1,\gamma_2\rangle = 2$ , and the standard electric charge in the lattice is  $\frac{1}{2} (\gamma_1 + \gamma_2)$ . Thus in the sum above, the terms with  $r = s$  describe the expected decomposition of the representation of  $SU(2)$  into electrically charged states on the Coulomb branch. The terms with  $r\neq s$  carry non-vanishing magnetic charge.

The general formula (1.2) for the defect index reduces to

$$
\mathcal {I} _ {L _ {0, n}} (q) = (q) _ {\infty} ^ {2} \operatorname {T r} \left[ F (L _ {0, n}, \vartheta) E _ {q} (X _ {\gamma_ {1}}) E _ {q} (X _ {\gamma_ {2}}) E _ {q} (X _ {\gamma_ {1}} ^ {- 1}) E _ {q} (X _ {\gamma_ {2}} ^ {- 1}) \right]. \tag {3.37}
$$

The trace involves a linear combination of the following quantity

$$
\begin{array}{l} \mathrm {T r} \left[ X _ {\gamma} E _ {q} (X _ {\gamma_ {1}}) E _ {q} (X _ {\gamma_ {2}}) E _ {q} (X _ {\gamma_ {1}} ^ {- 1}) E _ {q} (X _ {\gamma_ {2}} ^ {- 1}) \right] \\ = \sum_ {k _ {1}, k _ {2}, \ell_ {1}, \ell_ {2} = 0} ^ {\infty} \frac {(- 1) ^ {k _ {1} + k _ {2} + \ell_ {1} + \ell_ {2}} q ^ {\frac {1}{2} (k _ {1} + k _ {2} + \ell_ {1} + \ell_ {2})}}{(q) _ {k _ {1}} (q) _ {k _ {2}} (q) _ {\ell_ {1}} (q) _ {\ell_ {2}}} \operatorname {T r} \left[ X _ {\gamma} X _ {\gamma_ {1}} ^ {\ell_ {1}} X _ {\gamma_ {2}} ^ {\ell_ {2}} X _ {\gamma_ {1}} ^ {- k _ {1}} X _ {\gamma_ {2}} ^ {- k _ {2}} \right] \\ = \sum_ {k _ {1}, k _ {2}, \ell_ {1}, \ell_ {2} = 0} ^ {\infty} \frac {(- 1) ^ {k _ {1} + k _ {2} + \ell_ {1} + \ell_ {2}} q ^ {\frac {1}{2} (k _ {1} + k _ {2} + \ell_ {1} + \ell_ {2})}}{(q) _ {k _ {1}} (q) _ {k _ {2}} (q) _ {\ell_ {1}} (q) _ {\ell_ {2}}} q ^ {a _ {1} a _ {2} + 2 (\ell_ {1} + a _ {1}) (k _ {2} - a _ {2})} \delta_ {\ell_ {1} + a _ {1}, k _ {1}} \delta_ {\ell_ {2}, k _ {2} - a _ {2}} \tag {3.38} \\ = (-1)^{a_{1} + a_{2}}q^{\frac{1}{2} (a_{1} + a_{2}) + a_{1}a_{2}}\sum_{\substack{\ell_{1} = \max (\lceil -a_{1}\rceil ,0)\\ \ell_{2} = \max (\lceil -a_{2}\rceil ,0)}}^{\infty}\frac{q^{\ell_{1} + \ell_{2} + 2\ell_{1}\ell_{2} + 2a_{1}\ell_{2}}}{(q)_{\ell_{1}}(q)_{\ell_{2}}(q)_{\ell_{1} + a_{1}}(q)_{\ell_{2} + a_{2}}}, \\ \end{array}
$$

where  $\gamma = a_{1}\gamma_{1} + a_{2}\gamma_{2}$ . Thus, all that remains is to sum these quantities as dictated by generating functions  $F(L_{0,n})$  which may be found in (3.36). Carrying out this straightforward

calculation we find for instance

$$
\mathcal {I} _ {L _ {0, 2}} = - 2 q + q ^ {2} - 2 q ^ {4} + q ^ {6} - 2 q ^ {9} + q ^ {1 2} + \mathcal {O} \left(q ^ {1 3}\right), \tag {3.39}
$$

$$
\mathcal {I} _ {L _ {0, 4}} = q ^ {2} + 2 q ^ {3} - 2 q ^ {4} + q ^ {6} + 2 q ^ {7} - 2 q ^ {9} + q ^ {1 2} + \mathcal {O} \left(q ^ {1 3}\right),
$$

in exact agreement with the UV defect indices (3.15) computed from localization.

# 3.4.2 't Hooft Lines in  $SU(2)$  Gauge Theory

For the 't Hooft line  $L_{1,0}$  with minimal magnetic charge in the pure  $SU(2)$  gauge theory, the generating function of the frame BPS states is very simple

$$
F \left(L _ {1, 0}\right) = X _ {- \frac {1}{2} \gamma_ {2}}. \tag {3.40}
$$

The general formula (1.2) for a full 't Hooft line  $L_{1,0}$  can be computed in a similar way as in the Wilson line defect cases,

$$
\begin{array}{l} \mathcal {I} _ {L _ {1, 0}} (q) = \left(q\right) _ {\infty} ^ {2} \operatorname {T r} \left[ F \left(L _ {1, 0}\right) ^ {2} E _ {q} \left(X _ {\gamma_ {1}}\right) E _ {q} \left(X _ {\gamma_ {2}}\right) E _ {q} \left(X _ {\gamma_ {1}} ^ {- 1}\right) E _ {q} \left(X _ {\gamma_ {2}} ^ {- 1}\right) \right] \tag {3.41} \\ = - q ^ {\frac {1}{2}} + q ^ {\frac {5}{2}} - q ^ {\frac {7}{2}} - q ^ {\frac {9}{2}} + q ^ {\frac {1 3}{2}} + q ^ {\frac {1 5}{2}} - q ^ {\frac {1 7}{2}} - q ^ {\frac {1 9}{2}} + \dots , \\ \end{array}
$$

which equals the 't Hooft line defect index (3.21) computed from UV localization.

# 3.4.3 Wilson Lines in  $SU(2)$  Gauge Theory with  $N_{f} = 4$  Flavors

Let us move on to the half Wilson line index in the  $SU(2)$  superconformal QCD. The framed BPS degeneracies for the line defect can be read off, say, from the class  $\mathcal{S}$  description, given i.e. in [27] in a slightly different chamber. For illustrative purposes, we will reproduce the answer for the Wilson line in the doublet using the representation theory of the framed BPS quiver, i.e., an extended BPS quiver with one extra node representing the defect [65].<sup>11</sup>

In the chamber shown in Figure 2, the core charge for a Wilson line in the 2 of  $SU(2)$  is

$$
\gamma_ {c} = - \frac {1}{2} \left(\gamma_ {1} + \gamma_ {2} + \gamma_ {3} + \gamma_ {5}\right). \tag {3.42}
$$

The framed BPS quiver for this Wilson line is shown in Figure 5, with the square node

representing the doublet Wilson line.

![](images/b317f574bc07b745aad5af8f0c6cc58ed288485957b2af557fca8d65b513eafe.jpg)  
Figure 5: The framed quiver for a Wilson line in the  $\mathbf{2}$  in the  $SU(2)$  gauge theory with  $N_{f} = 4$  flavors. The core charge of the doublet Wilson line is given by  $\gamma_{c} = -\frac{1}{2} (\gamma_{1} + \gamma_{2} + \gamma_{3} + \gamma_{5})$ .

The generating function  $F(L_{0,1})$  for the Wilson line in the 2 of  $SU(2)$  is determined from the framed BPS states. The framed BPS states of the framed quiver can in turn be obtained by mutations as in Figure 6,

$$
F \left(L _ {0, 1}\right) = X _ {\gamma_ {c}} + X _ {\gamma_ {c} + \gamma_ {2}} + X _ {\gamma_ {c} + \gamma_ {2} + \gamma_ {3}} + X _ {\gamma_ {c} + \gamma_ {2} + \gamma_ {5}} + X _ {\gamma_ {c} + \gamma_ {2} + \gamma_ {3} + \gamma_ {5}} + X _ {\gamma_ {c} + \gamma_ {1} + \gamma_ {2} + \gamma_ {3} + \gamma_ {5}}. \tag {3.43}
$$

Where the central charge of the framed node is to the right of the vanilla BPS states.

The line defect index (3.23) can be reproduced from the trace of  $F(L_{0,1})$ . For example, the leading term  $-\chi_{[1,0,0,0]}(\eta_i)q^{\frac{1}{2}}$  in (3.23) can be reproduced from the  $q^{\frac{1}{2}}$  term in  $S(q)$  (2.34) with the help of (2.29),

$$
\begin{array}{l} (q) _ {\infty} ^ {2} \operatorname {T r} [ F (L _ {0, 1}) \mathcal {S} (q) \overline {{\mathcal {S}}} (q) ] = (q) _ {\infty} ^ {2} \operatorname {T r} \left[ F (L _ {0, 1}) \left(1 - q ^ {\frac {1}{2}} \sum_ {i = 1} ^ {6} X _ {\gamma_ {i}}\right) \left(1 - q ^ {\frac {1}{2}} \sum_ {j = 1} ^ {6} X _ {- \gamma_ {j}}\right) \right] + \mathcal {O} (q) \\ = - q ^ {\frac {1}{2}} \Big (X _ {\frac {1}{2} (\gamma_ {1} + \gamma_ {2} + \gamma_ {3} - \gamma_ {5})} + X _ {\frac {1}{2} (- \gamma_ {1} - \gamma_ {2} + \gamma_ {3} - \gamma_ {5})} + X _ {\frac {1}{2} (\gamma_ {1} + \gamma_ {2} - \gamma_ {3} + \gamma_ {5})} + X _ {\frac {1}{2} (- \gamma_ {1} - \gamma_ {2} - \gamma_ {3} + \gamma_ {5})} \\ \left. + X _ {\frac {1}{2} \left(\gamma_ {1} + \gamma_ {2} + \gamma_ {3} + 2 \gamma_ {4} + \gamma_ {5}\right)} + X _ {\frac {1}{2} \left(- \gamma_ {1} - \gamma_ {2} - \gamma_ {3} - 2 \gamma_ {4} - \gamma_ {5}\right)} + X _ {\frac {1}{2} \left(\gamma_ {1} + \gamma_ {2} + \gamma_ {3} + \gamma_ {5} + 2 \gamma_ {6} \right.} + X _ {\frac {1}{2} \left(- \gamma_ {1} - \gamma_ {2} - \gamma_ {3} - \gamma_ {5} - 2 \gamma_ {6}\right)}\right) \\ + \mathcal {O} \left(q ^ {\frac {3}{2}}\right). \tag {3.44} \\ \end{array}
$$

![](images/af8355d8e05735bced64c6cfac0c35937ff6bd7582092e11fe05591e6c452a88.jpg)

![](images/5cd4b99cde76e9b57aeeaddd9be811a77ae1ec9ecbf552dd58867be59a07a578.jpg)

![](images/df5b4d1bf16ee69adcbfb25796a8a768636e087eaf8e5dbef9d09693a429b5cc.jpg)

![](images/4611f2cbcd20add8aaece66c9907ad5713d264e9b5ad769a18ced6fa12bc6eaf.jpg)

![](images/5c12775a607b02f6c4793b23b7a6b3e2dd306f27af5e17a8ca93c52468808928.jpg)

![](images/ff7a7cac5c0a44a5512e497ef30e575b000d1e22b2a675772dc9e9000bb1d99c.jpg)

![](images/59154cae9470f0aac53d4a60a6fd19b208431b5f585b8bf4fab6aa7c5688fdd7.jpg)  
Figure 6: Mutations of the framed quiver for the doublet Wilson line defect in the  $SU(2)$  gauge theory with  $N_{f} = 4$  flavors. The core charge for the doublet Wilson line is  $\gamma_{c} = -\frac{1}{2} (\gamma_{1} + \gamma_{2} + \gamma_{3} + \gamma_{5})$ . The crossed denotes the node that is about to be right-mutated.

We have further computed the trace of  $L_{0,1}$  to  $q^{\frac{3}{2}}$  order,

$$
(q) _ {\infty} ^ {2} \operatorname {T r} [ F (L _ {0, 1}) \mathcal {S} (q) \overline {{\mathcal {S}}} (q) ] = - \chi_ {[ 1, 0, 0, 0 ]} (\eta_ {i}) q ^ {\frac {1}{2}} - \chi_ {[ 1, 1, 0, 0 ]} (\eta_ {i}) q ^ {\frac {3}{2}} + \mathcal {O} (q ^ {\frac {5}{2}}), \tag {3.45}
$$

where we have used (2.29) to express the flavor  $X_{\gamma}$  in terms of the  $SO(8)$  fugacities  $\eta_{i}$ . Indeed, we see that the trace of  $L_{0,1}$  nicely agrees with the line defect index (3.23). In the case when the flavor fugacities are off,  $\eta_{i} = 1$ , we have computed the trace of  $L_{0,1}$  up to order  $q^{\frac{7}{2}}$ ,

$$
(q) _ {\infty} ^ {2} \mathrm {T r} [ F (L _ {0, 1}) \mathcal {S} (q) \overline {{\mathcal {S}}} (q) ] \Big | _ {\eta_ {i} = 1} = - 8 q ^ {\frac {1}{2}} - 1 6 0 q ^ {\frac {3}{2}} - 1 6 2 4 q ^ {\frac {5}{2}} - 1 1 7 6 8 q ^ {\frac {7}{2}} + \mathcal {O} (q ^ {\frac {9}{2}}), \quad (3. 4 6)
$$

which is equal to the doublet Wilson line defect index  $\mathcal{I}_{L_{0,1}}(q,\eta_i = 1)$  in (3.26) obtained from localization.

# 4 Half-Indices and Boundary Conditions

The Schur index can be generalized further by the insertion of half-BPS boundaries or interfaces along the equator of the sphere. These boundary conditions will preserve  $3d$ $\mathcal{N} = 2$  supersymmetry. $^{12}$

We have already encountered the simplest example in the form of the half-index

$$
\Pi_ {\vec {m}} ^ {S} (q, u, z), \tag {4.1}
$$

which corresponds to a choice of Dirichlet boundary conditions for the UV gauge fields. Remember that we made a choice of Lagrangian splitting of the hypermultiplet scalar fields, which determines which half of the scalars has Dirichlet boundary condition and which half has Neumann boundary condition [37]. The bulk gauge symmetry becomes a  $3d$  global symmetry at a Dirichlet boundary and thus the fugacity  $u$  and magnetic charge  $\vec{m}$  should be interpreted as associated to that  $3d$  global symmetry.

Remember that the general index is written as

$$
\mathcal {I} (q, z) = \left(\Pi^ {N}, \Pi^ {S}\right) = \sum_ {\vec {m}} \int [ d u ] _ {\vec {m}} \Pi_ {- \vec {m}} ^ {N} (q, u, z) \Pi_ {\vec {m}} ^ {S} (q, u, z). \tag {4.2}
$$

We can interpret the a sum over magnetic fluxes and the integral over gauge fugacities as following from the fact that the index is glued from two hemisphere indices with Dirichlet boundary condition by restoring the gauge fields at the equator. Of course, only the term

$\vec{m} = 0$  contributes unless we add line defects in the two hemispheres.

The half index for a more general boundary condition, with Neumann boundary condition for the gauge group and general boundary matter fields is written analogously as

$$
\mathcal {I I} _ {\vec {n}} (q, z, \xi) = (\Pi^ {N}, Z) = \sum_ {\vec {m}} \int [ d u ] _ {\vec {m}} \Pi_ {- \vec {m}} ^ {N} (q, u, z) Z _ {\vec {m}, \vec {n}} (q, u, z, \xi), \tag {4.3}
$$

where  $Z_{\vec{m},\vec{n}}(q,u,z,\xi)$  is the 3d index of the boundary matter fields. We have included fugacities  $\xi$  and magnetic flux  $\vec{n}$  for possible other global symmetries of the 3d matter fields. It is also possible to gauge at the boundary only a subgroup of the original gauge group, by restricting appropriately the  $u$  integral and  $\vec{m}$  sum.

Again, only the term  $\vec{m} = 0$  contributes unless we add bulk line defects in the hemisphere, as in<sup>13</sup>

$$
\mathcal {I} \mathcal {I} _ {\vec {n}} ^ {L} (q, z, \xi) = (\Pi^ {N}, \hat {O} _ {L} Z). \tag {4.4}
$$

The line defect can be brought to any position along the great circle of the sphere, which intersects the equator at two points  $N$  and  $S$ , the poles of the boundary  $S^2$ . In a purely 3d context, the index  $Z_{\vec{m},\vec{n}}(q,u,z,\xi)$  often satisfies difference equations which arise from the insertion of line defects at  $N$  or  $S$ . More precisely,  $Z$  satisfies two sets of difference equations built from the two commuting sets of difference operators  $x,p$  and  $x',p'$ .

We propose the following IR description of the Schur half-index,

$$
\mathcal {I} \mathcal {I} _ {\vec {n}} (q, z, \xi) = (q) _ {\infty} ^ {r} \mathrm {T r} \left[ Z _ {\vec {n}} ^ {I R} (q, \xi) [ X ] \overline {{\mathcal {S}}} (q) \right], \tag {4.5}
$$

and for the half-index with line defect insertion,

$$
\mathcal {I} \mathcal {I} _ {\vec {n}} ^ {L} (q, z, \xi) = (q) _ {\infty} ^ {r} \operatorname {T r} \left[ F (L) Z _ {\vec {n}} ^ {I R} (q, \xi) [ X ] \overline {{\mathcal {S}}} (q) \right], \tag {4.6}
$$

where

$$
Z _ {\vec {n}} ^ {I R} (q, \xi) [ X ] \equiv \sum_ {\gamma} Z _ {\gamma , \vec {n}} ^ {I R} (q, \xi) X _ {\gamma} \tag {4.7}
$$

is a formal generating function for the  $3d$  indices of the IR degrees of freedom living on the domain wall, expressed in a charge basis for the Abelian symmetries which are coupled to the bulk Abelian gauge fields. And as usual,  $r$  is the rank of the Coulomb branch.

We can give an intuitive interpretation of formula (4.5) by interpreting it in the IR effective QED description of the Coulomb branch. Indeed, If we pretend the BPS particles are free and mutually local, this would be a Lagrangian splitting for the bulk hypermultiplets,

which is expected to be such that the bulk fields which survive at the boundary are those whose charges appear in the  $S(q)$  product.

The full story is likely more complex and requires some way to define some kind of effective action, both in the bulk and boundary, analogous to what is done in  $2d$  in [54, 55]. At the level of the index, though, this approximate perspective is expected to be sufficient.

With this caveat, the wall-crossing behavior of the IR boundary conditions is well-understood [38]. Across a wall of marginal stability for a BPS particle, where some BPS ray exits the half plane associated to  $S(q)$  and the opposite ray enters it, the canonical choice of boundary condition also flips. The boundary degrees of freedom change in such a way as to compensate the change in boundary condition.

For BPS hypermultiplets, the flip of boundary condition adds an extra chiral field to the boundary degrees of freedom. This multiplies  $Z_{\gamma}^{IR}$  by the Fourier modes of

$$
Z _ {m} ^ {\text {c h i r a l}} (z) = \frac {\prod_ {n = 0} ^ {\infty} \left(1 + z ^ {- 1} q ^ {- m / 2 + n + \frac {1}{2}}\right)}{\prod_ {n = 0} ^ {\infty} \left(1 + z q ^ {- m / 2 + n + \frac {1}{2}}\right)}, \tag {4.8}
$$

i.e. changes  $Z^{IR}[X]$  to

$$
\tilde {Z} ^ {I R} [ X ] = \prod_ {n = 0} ^ {\infty} (1 + X _ {- \gamma} q ^ {n + \frac {1}{2}}) Z ^ {I R} [ X ] \frac {1}{\prod_ {n = 0} ^ {\infty} (1 + X _ {\gamma} q ^ {n + \frac {1}{2}})} = E _ {q} ^ {- 1} (X _ {- \gamma}) Z ^ {I R} [ X ] E _ {q} (X _ {\gamma}).
$$

Then the candidate Schur index is invariant

$$
\mathcal {I} \mathcal {I} _ {\vec {n}} (q, z, \xi) = (q) _ {\infty} ^ {r} \mathrm {T r} \left[ \mathcal {S} (q) Z _ {\vec {n}} ^ {I R} (q, \xi) [ X ] \right] = (q) _ {\infty} ^ {r} \mathrm {T r} \left[ \tilde {\mathcal {S}} (q) \tilde {Z} _ {\vec {n}} ^ {I R} (q, \xi) [ X ] \right], \tag {4.10}
$$

as  $\mathcal{S}(q) = E_q(X_\gamma)\tilde{\mathcal{S}}(q)E_q^{-1}(X_{-\gamma})$ . Vice versa, adding a chiral of the opposite charge to cross the wall backwards gives

$$
Z ^ {I R} [ X ] = E _ {q} (X _ {- \gamma}) \tilde {Z} ^ {I R} [ X ] E _ {q} ^ {- 1} (X _ {\gamma}). \tag {4.11}
$$

We expect these relations to hold for BPS particles of every spin. It would be interesting to understand which boundary degrees of freedom are added or removed in that case, but as higher spin BPS particles usually come together with infinite cohorts of hypermultiplet particles (see e.g. [87]), individual wall-crossing events are perhaps less physically meaningful.

# 4.1 Indices and RG Interfaces

There is a special class of boundary conditions/ interfaces which is very useful in relating BPS quantities in the UV and IR description of  $\mathcal{N} = 2$  gauge theories:  $RG$  interfaces. These are special interfaces between the UV theory and its IR effective description, obtained by applying the IR effective description to one side only of the identity interface in the UV [38,88].

A very useful properties of RG interfaces is that they intertwine between the IR and UV description of several BPS objects, included line defects and boundary conditions. Thus the IR description of a UV boundary condition is obtained by acting on it with the RG interface, and vice versa. This implies a precise relation between the indices  $Z_{\vec{m}}^{UV}(q,u)$  and  $Z^{IR}[X]$  of the UV and IR boundary degrees of freedom. Similar considerations apply to line defects. We will write down these relations momentarily.

We will call the interface degrees of freedom defining an RG interface the  $RG$  theory. Formally, the interface degrees of freedom can be obtained by starting from UV Dirichlet boundary conditions and flowing to the IR: the result should be the RG theory coupled to the IR Abelian gauge theory. Vice versa, the UV boundary condition defined by the RG theory will flow to Dirichlet boundary conditions for the IR theory. Of course, the RG theory depends on a choice of hypermultiplet splitting in the UV and a chamber as well as an electromagnetic duality frame for the IR theory. The RG theory transforms appropriately as these choices are modified. Explicit examples of conjectural RG theories were described in [38].

At the level of the indices, the UV Dirichlet boundary condition gives us the half-index and thus we expect

$$
\Pi_ {\vec {m}} ^ {S} (q, u, z) = (q) _ {\infty} ^ {r} \operatorname {T r} \left[ K _ {\vec {m}} (q, u, z) [ X ] \mathcal {S} _ {\vartheta + \pi} (q) \right], \tag {4.12}
$$

where  $K_{\gamma,\vec{m}}(q,u,z)$  is the 3d index of the RG theory. Building the UV boundary condition which flows to a simple Dirichlet IR boundary condition gives us an inverse relation:

$$
(q) _ {\infty} ^ {r} \mathcal {S} _ {\vartheta} (q) = (\Pi^ {N}, K [ X ]) = \sum_ {\vec {m}} \int [ d u ] _ {\vec {m}} \Pi_ {- \vec {m}} ^ {N} (q, u, z) K _ {\vec {m}} (q, u, z) [ X ]. \tag {4.13}
$$

In particular, this gives a direct physical meaning of the quantum spectrum generator  $S(q)$  as a generating function for Schur indices in the presence of the RG interface boundary condition.

Notice also that these relations are compatible and imply the identity between the UV and IR bulk Schur indices:

$$
\left(\Pi^ {N}, \Pi^ {S}\right) = (q) _ {\infty} ^ {r} \operatorname {T r} \left[ \left(\Pi^ {N}, K [ X ]\right) \mathcal {S} _ {\vartheta + \pi} (q) \right] = (q) _ {\infty} ^ {2 r} \operatorname {T r} \left[ \mathcal {S} _ {\vartheta} (q) \mathcal {S} _ {\vartheta + \pi} (q) \right]. \tag {4.14}
$$

Furthermore, the existence of the Kernel  $K_{\vec{m}}(q,u,z)[X]$  also implies our formulae for line defects: the index of the RG interface satisfies difference equations of the form

$$
\hat {O} _ {L} K [ X ] = F (L, \vartheta) K [ X ], \tag {4.15}
$$

and

$$
K [ X ] \hat {O} _ {L} ^ {\prime} = K [ X ] F (L, \vartheta + \pi), \tag {4.16}
$$

The above two relations imply our proposal for the line defect indices (3.30),

$$
\begin{array}{l} (\Pi^ {N}, \hat {O} _ {L} \Pi^ {S}) = (q) _ {\infty} ^ {r} \mathrm {T r} \left[ (\Pi^ {N}, \hat {O} _ {L} K [ X ]) \mathcal {S} _ {\vartheta + \pi} (q) \right] = (q) _ {\infty} ^ {r} \mathrm {T r} \left[ F (L, \vartheta) (\Pi^ {N}, K [ X ]) \mathcal {S} _ {\vartheta + \pi} (q) \right] \\ = (q) _ {\infty} ^ {2 r} \operatorname {T r} \left[ F (L, \vartheta) \mathcal {S} _ {\vartheta} (q) \mathcal {S} _ {\vartheta + \pi} (q) \right]. \tag {4.17} \\ \end{array}
$$

Indeed, one can argue that equations (4.15) and (4.16) are the truly crucial relationships. For example, they imply a recursion relation

$$
F (L, \vartheta) (\Pi^ {N}, K [ X ]) = (\Pi^ {N}, K [ X ]) F (L, \vartheta + \pi), \tag {4.18}
$$

which in turn implies its identification with  $(q)_{\infty}^{r}\mathcal{S}_{\vartheta}(q)$  up to an overall function of  $q$  and similarly for  $(q)_{\infty}^{r}\mathrm{Tr}\left[K_{\vec{m}}(q,u,z)[X]\mathcal{S}_{\vartheta +\pi}(q)\right]$ .

Similar considerations apply in the presence of boundary conditions. If  $Z^{UV}$  is the partition function in the presence of a UV boundary condition and  $Z^{IR}$  the partition function of the IR boundary condition to which it flows, we expect the relations

$$
(Z ^ {U V}, K [ X ]) = Z ^ {I R} [ X ], \qquad \qquad Z ^ {U V} = \mathrm {T r} \left[ Z ^ {I R} [ X ] \bar {K} [ X ] \right]. \qquad (4. 1 9)
$$

These then imply

$$
(\Pi^ {N}, Z ^ {U V}) = \mathrm {T r} \left[ Z ^ {I R} [ X ] (\Pi^ {N}, \bar {K} [ X ]) \right] = (q) _ {\infty} ^ {r} \mathrm {T r} \left[ Z ^ {I R} [ X ] \mathcal {S} _ {\vartheta + \pi} (q) \right]. \qquad (4. 2 0)
$$

# 4.2 A Free Hypermultiplet

The free hypermultiplet index in IR conventions is

$$
\mathcal {I} _ {\text {h y p e r}} = \frac {1}{\prod_ {n = 0} ^ {\infty} \left(1 + z q ^ {n + \frac {1}{2}}\right) \left(1 + z ^ {- 1} q ^ {n + \frac {1}{2}}\right)} = E _ {q} (z) E _ {q} \left(z ^ {- 1}\right). \tag {4.21}
$$

This is an obvious example of the BPS formula.

The half-indices for the two possible half-BPS boundary conditions are

$$
\mathcal {I I} _ {\text {h y p e r}, \pm} = \frac {1}{\prod_ {n = 0} ^ {\infty} \left(1 + z ^ {\pm} q ^ {n + \frac {1}{2}}\right)} = E _ {q} \left(z ^ {\pm}\right). \tag {4.22}
$$

This is a neat example of the BPS formula for half-indices, which is already somewhat non-trivial.

If we pick the phase of the hypermultiplet mass in such a way that  $\mathcal{S}(q) = E_q(z)$ , we have that  $\mathcal{I}\mathcal{I}_{\mathrm{hyper},+}$  simply equals  $\mathcal{S}(q)$ . This makes sense: we are already using the "canonical" choice of IR boundary condition for the bulk hypermultiplet.

On the other hand, we have

$$
\mathcal {I I} _ {\text {h y p e r , -}} = E _ {q} \left(z ^ {- 1}\right) = \mathcal {S} (q) E _ {q} ^ {- 1} (z) E _ {q} \left(z ^ {- 1}\right). \tag {4.23}
$$

We recognize the expression  $E_{q}^{-1}(z)E_{q}(z^{-1})$  as the 3d index of a 3d chiral field with no magnetic flux on the two-sphere. Indeed, the two boundary conditions can be related by adding a boundary chiral field with appropriate boundary superpotential coupling to the hypermultiplet [36-38].

# 4.3 Pure  $SU(2)$  Gauge Theory

The RG theory for pure  $SU(2)$  gauge theory consists of a  $SU(2)$  doublet of chiral fields, transforming with charge  $-1$  under a  $U(1)$  global symmetry [38].

We can immediately compute

$$
\begin{array}{l} \mathcal {I I} _ {n} ^ {R G} (q, \xi) = - \frac {1}{4 \pi i} \int \frac {d u}{u} (u - u ^ {- 1}) ^ {2} \left[ (q; q) _ {\infty} (q u ^ {2}; q) _ {\infty} (q u ^ {- 2}; q) _ {\infty} \right] \\ \times \left[ \frac {\left(- q ^ {\frac {1}{2} + \frac {n}{2}} \xi u ^ {- 1} ; q\right) _ {\infty}}{\left(- q ^ {\frac {1}{2} + \frac {n}{2}} \xi^ {- 1} u ; q\right) _ {\infty}} \frac {\left(- q ^ {\frac {1}{2} + \frac {n}{2}} \xi u ; q\right) _ {\infty}}{\left(- q ^ {\frac {1}{2} + \frac {n}{2}} \xi^ {- 1} u ^ {- 1} ; q\right) _ {\infty}} \right] , \tag {4.24} \\ \end{array}
$$

where  $\xi$  and  $n$  are the fugacity and the magnetic flux, respectively, for the 3d  $U(1)$  global symmetry of the RG theory. We have separated in the integrand the contributions from the hemisphere and from the boundary doublet. The superscript  $RG$  indicates that we are computing the half index in the RG boundary condition.

Explicit calculation at finite order in  $q$  suggests that this complicated contour integral has a dramatically simple answer:

$$
\mathcal {I I} _ {n} ^ {R G} (q, \xi) = (q) _ {\infty} \sum_ {e = \max  (0, - n)} ^ {\infty} \xi^ {2 e} \frac {q ^ {e (e + n)}}{(q) _ {e} (q) _ {e + n}}, \tag {4.25}
$$

The corresponding generating function is

$$
\begin{array}{l} \mathcal {I I} ^ {R G} (q; X) \equiv : \sum_ {n} (- 1) ^ {n} q ^ {\frac {n}{2}} \mathcal {I I} _ {n} ^ {R G} (q, \xi = q ^ {\frac {1}{2}} X _ {\gamma}) X _ {- n \gamma^ {\prime}}: \\ = (q) _ {\infty} \sum_ {n} \sum_ {e = \max  (0, n)} ^ {\infty} (- 1) ^ {n} q ^ {e - \frac {n}{2}} \frac {q ^ {e (e - n)}}{(q) _ {e} (q) _ {e - n}} X _ {n \gamma^ {\prime} + 2 e \gamma}, \tag {4.26} \\ \end{array}
$$

where we denote the electric charge and the magnetic charge by  $\gamma$  and  $\gamma'$ , respectively. The insertion of  $(-1)^n$  is interpreted as a convenient shift of the fermion number of monopole operators, while the insertion of  $q^{\frac{n}{2}}$  and the factor of  $q^{\frac{1}{2}}$  in  $\xi$  are useful re-definitions of the boundary R-charge. The latter essentially assigns trivial R-charge and fermion number to the bosonic components of the boundary chiral multiplet. $^{14}$  Importantly, the R-charge shift  $\xi \rightarrow q^{1/2}\xi$  has to be performed after the contour integral in  $u$ , not before.

The generating function for the half index can be written as

$$
\mathcal {I I} ^ {R G} (q; X) = (q) _ {\infty} \sum_ {e \geq 0} \sum_ {e ^ {\prime} \geq 0} (- q ^ {\frac {1}{2}}) ^ {e + e ^ {\prime}} \frac {q ^ {e e ^ {\prime}}}{(q) _ {e} (q) _ {e ^ {\prime}}} X _ {- e ^ {\prime} \gamma^ {\prime} + e \left(\gamma^ {\prime} + 2 \gamma\right)} = (q) _ {\infty} E _ {q} \left(X _ {2 \gamma + \gamma^ {\prime}}\right) E _ {q} \left(X _ {- \gamma^ {\prime}}\right), \tag {4.27}
$$

where we set the Dirac pairing  $\langle \gamma', \gamma \rangle = 1$ . This is the same as  $(q)_{\infty} S(q)$  if we choose the following electromagnetic duality frame  $\gamma_1 = 2\gamma + \gamma'$  and  $\gamma_2 = -\gamma'$  for the two nodes in the BPS quiver of the pure  $SU(2)$  theory! We can thus identify

$$
K [ X ] \equiv : \sum_ {n} (- 1) ^ {n} q ^ {\frac {n}{2}} \left[ \frac {\left(- q ^ {1 + \frac {n}{2}} X _ {\gamma} u ^ {- 1} ; q\right) _ {\infty}}{\left(- q ^ {\frac {n}{2}} X _ {- \gamma} u ; q\right) _ {\infty}} \frac {\left(- q ^ {1 + \frac {n}{2}} X _ {\gamma} u ; q\right) _ {\infty}}{\left(- q ^ {\frac {n}{2}} X _ {- \gamma} u ^ {- 1} ; q\right) _ {\infty}} \right] X _ {- n \gamma^ {\prime}}: . \tag {4.28}
$$

Similarly, if we insert a Wilson loop operators

$$
\begin{array}{l} \mathcal {I I} _ {n} ^ {R G; W} (q, \xi) \equiv - \frac {1}{4 \pi i} \int \frac {d u}{u} (u - u ^ {- 1}) ^ {2} (u + u ^ {- 1}) \left[ (q) _ {\infty} (q u ^ {2}; q) _ {\infty} (q u ^ {- 2}; q) _ {\infty} \right] \\ \times \left[ \frac {\left(- q ^ {\frac {1}{2} + \frac {n}{2}} \xi u ^ {- 1} ; q\right) _ {\infty}}{\left(- q ^ {\frac {1}{2} + \frac {n}{2}} \xi^ {- 1} u ; q\right) _ {\infty}} \frac {\left(- q ^ {\frac {1}{2} + \frac {n}{2}} \xi u ; q\right) _ {\infty}}{\left(- q ^ {\frac {1}{2} + \frac {n}{2}} \xi^ {- 1} u ^ {- 1} ; q\right) _ {\infty}} \right] , \tag {4.29} \\ \end{array}
$$

we get a neat expression

$$
\mathcal {I I} _ {- n} ^ {R G; W} (q, \xi) = (q) _ {\infty} \sum_ {e = \max  (0, n)} ^ {\infty} \xi^ {1 - 2 e} \frac {q ^ {e ^ {2} - e n - e + \frac {n + 1}{2}} - q ^ {e ^ {2} - e n + \frac {n + 1}{2}} - q ^ {e ^ {2} - e n - \frac {n - 1}{2}}}{(q) _ {e} (q) _ {e - n}}, \tag {4.30}
$$

The corresponding generating function is

$$
\mathcal {I I} ^ {R G; W} (q; X) = (q) _ {\infty} \sum_ {e \geq 0} \sum_ {e ^ {\prime} \geq 0} (- q ^ {\frac {1}{2}}) ^ {e + e ^ {\prime} - 1} \frac {q ^ {(e - \frac {1}{2}) (e ^ {\prime} - \frac {1}{2}) + \frac {1}{4}} \left(1 - q ^ {e} - q ^ {e ^ {\prime}}\right)}{(q) _ {e} (q) _ {e ^ {\prime}}} X _ {e ^ {\prime} \gamma^ {\prime} - e (\gamma^ {\prime} + 2 \gamma) + \gamma}, \tag {4.31}
$$

which can be manipulated to

$$
\mathcal {I I} ^ {R G; W} (q; X) = (q) _ {\infty} E _ {q} (X _ {\gamma^ {\prime}}) E _ {q} (X _ {- 2 \gamma - \gamma^ {\prime}}) X _ {\gamma} - (q) _ {\infty} E _ {q} (X _ {\gamma^ {\prime}}) X _ {- \gamma} E _ {q} (X _ {- 2 \gamma - \gamma^ {\prime}}), \quad (4. 3 2)
$$

which yields a reasonable framed BPS degeneracy:

$$
\mathcal {I I} ^ {R G; W} (q; X) = (q) _ {\infty} \mathcal {S} (q) \left[ X _ {\gamma} - X _ {- \gamma} - X _ {- 3 \gamma - \gamma^ {\prime}} \right]. \tag {4.33}
$$

We can deal in a similar manner with 't Hooft lines. Consider for example  $\hat{O}_{1,0}$

$$
\frac {i}{q ^ {\frac {m}{2}} u - q ^ {- \frac {m}{2}} u ^ {- 1}} \left(p ^ {\frac {1}{2}} - p ^ {- \frac {1}{2}}\right) \tag {4.34}
$$

so that (  $n$  has to be half-integral here)

$$
\begin{array}{l} \mathcal {I I} _ {n} ^ {R G; \hat {O} _ {1, 0}} (q, \xi) \equiv - \frac {1}{4 \pi} \int \frac {d u}{u} (u - u ^ {- 1}) \left[ (q) _ {\infty} (q u ^ {2}; q) _ {\infty} (q u ^ {- 2}; q) _ {\infty} \right] \\ \left[ \frac {\left(- q ^ {\frac {n}{2}} \xi u ^ {- 1} ; q\right) _ {\infty}}{\left(- q ^ {\frac {1}{2} + \frac {n}{2}} \xi^ {- 1} u ; q\right) _ {\infty}} \frac {\left(- q ^ {1 + \frac {n}{2}} \xi u ; q\right) _ {\infty}}{\left(- q ^ {\frac {1}{2} + \frac {n}{2}} \xi^ {- 1} u ^ {- 1} ; q\right) _ {\infty}} - \frac {\left(- q ^ {1 + \frac {n}{2}} \xi u ^ {- 1} ; q\right) _ {\infty}}{\left(- q ^ {\frac {1}{2} + \frac {n}{2}} \xi^ {- 1} u ; q\right) _ {\infty}} \frac {\left(- q ^ {\frac {n}{2}} \xi u ; q\right) _ {\infty}}{\left(- q ^ {\frac {1}{2} + \frac {n}{2}} \xi^ {- 1} u ^ {- 1} ; q\right) _ {\infty}} \right] (4. 3 5) \\ \end{array}
$$

which becomes

$$
\mathcal {I I} _ {n} ^ {R G; \hat {O} _ {1, 0}} (q, \xi) \equiv \xi q ^ {\frac {n}{2}} \mathcal {I I} _ {n + \frac {1}{2}} ^ {R G} (q, q ^ {\frac {1}{4}} \xi), \tag {4.36}
$$

i.e.

$$
\mathcal {I I} ^ {R G; \hat {O} _ {1, 0}} (q; X) \equiv (q) _ {\infty} \mathcal {S} (q) X _ {\frac {\gamma^ {\prime}}{2} + \gamma}. \tag {4.37}
$$

Similarly, for the 't Hooft-Wilson line, the difference operator  $\hat{O}_{1,1}$  is

$$
\frac {i}{q ^ {\frac {m}{2}} u - q ^ {- \frac {m}{2}} u ^ {- 1}} \left(q ^ {\frac {m}{2}} u p ^ {\frac {1}{2}} - q ^ {- \frac {m}{2}} u ^ {- 1} p ^ {- \frac {1}{2}}\right) \tag {4.38}
$$

so that (  $n$  has to be half-integral here)

$$
\begin{array}{l} \mathcal {I I} _ {n} ^ {R G; \hat {O} _ {1, 1}} (q, \xi) \equiv - \frac {1}{4 \pi i} \int \frac {d u}{u} (u - u ^ {- 1}) \left[ (q) _ {\infty} (q u ^ {2}; q) _ {\infty} (q u ^ {- 2}; q) _ {\infty} \right] \\ \left\{u \left[ \frac {(- q ^ {\frac {n}{2}} \xi u ^ {- 1} ; q) _ {\infty}}{(- q ^ {\frac {1}{2} + \frac {n}{2}} \xi^ {- 1} u ; q) _ {\infty}} \frac {(- q ^ {1 + \frac {n}{2}} \xi u ; q) _ {\infty}}{(- q ^ {\frac {1}{2} + \frac {n}{2}} \xi^ {- 1} u ^ {- 1} ; q) _ {\infty}} \right] - u ^ {- 1} \left[ \frac {(- q ^ {1 + \frac {n}{2}} \xi u ^ {- 1} ; q) _ {\infty}}{(- q ^ {\frac {1}{2} + \frac {n}{2}} \xi^ {- 1} u ; q) _ {\infty}} \frac {(- q ^ {\frac {n}{2}} \xi u ; q) _ {\infty}}{(- q ^ {\frac {1}{2} + \frac {n}{2}} \xi^ {- 1} u ^ {- 1} ; q) _ {\infty}} \right] \right\}, \\ \end{array}
$$

which becomes

$$
\mathcal {I I} _ {n} ^ {R G; \hat {O} _ {1, 1}} (q, \xi) \equiv \mathcal {I I} _ {n + \frac {1}{2}} ^ {R G} (q, q ^ {\frac {1}{4}} \xi), \tag {4.40}
$$

i.e.

$$
\mathcal {I I} ^ {R G; \hat {O} _ {1, 1}} (q, X) \equiv q ^ {- \frac {1}{4}} (q) _ {\infty} \mathcal {S} (q) X _ {\frac {\gamma^ {\prime}}{2}}. \tag {4.41}
$$

Next, we should compute the IR index for Dirichlet boundary condition

$$
\begin{array}{l} (q) _ {\infty} \mathrm {T r} \left[ \overline {{\mathcal {S}}} (q) K _ {m} (q, u) [ X ] \right] = (q) _ {\infty} \sum_ {n = - \infty} ^ {\infty} \sum_ {e = \max (0, n)} ^ {\infty} q ^ {2 e - n} \frac {q ^ {e (e - n)}}{(q) _ {e} (q) _ {e - n}} \\ \times \frac {1}{2 \pi i} \oint \frac {d \xi}{\xi} \xi^ {2 e} \left[ \frac {\left(- q ^ {\frac {1}{2} - \frac {n}{2} + \frac {m}{2}} \xi^ {- 1} u ^ {- 1} ; q\right) _ {\infty}}{\left(- q ^ {\frac {1}{2} - \frac {n}{2} + \frac {m}{2}} \xi u ; q\right) _ {\infty}} \frac {\left(- q ^ {\frac {1}{2} - \frac {n}{2} - \frac {m}{2}} \xi^ {- 1} u ; q\right) _ {\infty}}{\left(- q ^ {\frac {1}{2} - \frac {n}{2} - \frac {m}{2}} \xi u ^ {- 1} ; q\right) _ {\infty}} \right] \tag {4.42} \\ \end{array}
$$

The  $q^{2e - n}$  factor is due again to the choice of quantum numbers for the boundary theory.

Amazingly, the sum vanishes unless  $m = 0$  and at  $m = 0$  it reproduces the expected answer:

$$
(q) _ {\infty} \mathrm {T r} \left[ \overline {{\mathcal {S}}} (q) K _ {m} (q, u) [ X ] \right] = \delta_ {m, 0} (q) _ {\infty} (q u ^ {2}; q) _ {\infty} (q u ^ {- 2}; q) _ {\infty} = \Pi_ {m} ^ {S} (q, u) \tag {4.43}
$$

Of course, all these miraculous-looking relations are somewhat demystified by the recursion relations satisfied by  $K[X]$ , which can be easily seen to intertwine between the difference operators associated to UV line defects and the corresponding generating functions of IR framed BPS degeneracies.

# 5 Defect Indices in Argyres-Douglas Theories

One important application of our IR formula in Section 3.3 is a prediction for the line defect Schur indices in the strongly-coupled Argyres-Douglas theories, where a UV localization calculation is not available. In this section, we will compute the Schur indices of the  $A_{2}, A_{3}, A_{4}$  Argyres-Douglas theories with the presence of line defects using the conjectural

formula (3.30). We will also discuss how the OPEs between the defects are respected by the line defect indices  $\mathcal{I}_L(q)$ . We will defer the discussion on the relation between these defect OPEs with the Verlinde algebra of the associated chiral algebra to Section 6.

# 5.1 Line Defect OPEs and Schur Indices

The UV line defects satisfy a non-commutative defect OPE that takes the form [27, 65]

$$
L _ {\alpha} L _ {\beta} \equiv \lim  _ {\epsilon \rightarrow 0} L _ {\alpha} (\vartheta) L _ {\beta} (\vartheta + \epsilon) = \sum_ {\gamma} c _ {\alpha \beta} ^ {\gamma} (q) L _ {\gamma} (\vartheta), \tag {5.1}
$$

where the coefficients  $c_{\alpha \beta}^{\gamma}(q)$  are valued in  $\mathbb{Z}_{\geq 0}[q^{\frac{1}{2}}, q^{-\frac{1}{2}}]$ . We have restored the central charge phase  $\vartheta$  dependence to indicate their positions on  $S^3 \times S^1$ . The defect OPE can be intuitively understood as bringing two line defects, which are points on a great circle in  $S^3 \times S^1$ , close to each other to form a composite line defect, which then admits the above expansion in terms of simple defects. This configuration preserves the  $U(1)$  rotation transverse to the great circle and the  $SU(2)_R$  symmetry, and we can turn on the variable  $q$  to keep track of these quantum numbers. The resulting OPE is non-commutative because the first line defect can approach the second one either from above or from below on the great circle. See Figure 7.

![](images/2852be3f883eb44e4961473ba964a629cca108bf699d73fea49f5282dd4502cf.jpg)  
Figure 7: The OPE for line defects. As we bring two line defects near each other on the  $S^3$  while they both wrap around the  $S^1$ , the composite defect can be expanded into a sum of simple defects.

The simplest example of the line defect OPE is that between the IR Abelian defects  $X_{\gamma}$  given in (1.3),

$$
X _ {\gamma} X _ {\gamma^ {\prime}} = q ^ {\frac {1}{2} \langle \gamma , \gamma^ {\prime} \rangle} X _ {\gamma + \gamma^ {\prime}},
$$

where the Dirac pairing  $\langle \gamma, \gamma' \rangle$  captures the angular momentum of the composite defect.

Because of supersymmetry, the defect OPE is independent of the distance separated between the two defects, and thus the OPE computed in the UV should agree with that computed in the IR. This implies that the IR description of the UV line defect, i.e., the generating function (3.29)

$$
F (L) = \sum_ {\gamma \in \Gamma} \underline {{\Omega}} (L, \gamma , q) X _ {\gamma},
$$

has to obey the same OPE

$$
\lim  _ {\epsilon \rightarrow 0} F (L _ {\alpha}, \vartheta) F (L _ {\beta}, \vartheta + \epsilon) = \sum_ {\gamma} c _ {\alpha \beta} ^ {\gamma} (q) F (L _ {\gamma}, \vartheta), \tag {5.2}
$$

where the product on the lefthand side is given by the non-commutative product of IR Abelian line defects in (1.3). This provides a strong consistency check on the framed BPS state degeneracies  $\overline{\Omega} (L,\gamma ,q)$ . Since the defect OPE can be computed in the UV, the coefficients  $c_{\alpha \beta}^{\gamma}(q)$  must be wall-crossing invariant, while the individual generating functions  $F(L_{\alpha},\vartheta)$  are not.

One consequence of our IR formula for the line defect Schur index (3.30),

$$
\mathcal {I} _ {L} (q) = \left(q\right) _ {\infty} ^ {2 r} \operatorname {T r} \left[ F (L, \vartheta) \mathcal {S} _ {\vartheta} (q) \mathcal {S} _ {\vartheta + \pi} (q) \right], \tag {5.3}
$$

is that it manifestly respect the defect OPE because  $F(L, \vartheta)$  does. That is,

$$
\mathcal {I} _ {L _ {\alpha} L _ {\beta}} (q) = \sum_ {\gamma} c _ {\alpha \beta} ^ {\gamma} (q) \mathcal {I} _ {L _ {\gamma}} (q). \tag {5.4}
$$

In addition, we will obtain more relations between the defect indices than those descended from the defect OPEs in the case of the Argyres-Douglas theories.

# 5.2  $A_{2}$  Argyres-Douglas Theory

![](images/1f4a7d02fe305ae83b53a8c3e4d66ffb4575595bd05219ea27b028a3dd7394c1.jpg)  
Figure 8: The BPS quiver for the  $A_{2}$  Argyres-Douglas theory.

The  $A_{2}$  Argyres-Douglas theory arises from special points on the moduli space of the pure  $SU(3)$  gauge theory or from the  $SU(2)$  SQCD with  $N_{f} = 1$  flavor [69,70]. The UV line defects in the  $A_{2}$  Argyres-Douglas theory are generated by five defects  $L_{i}$ 's together with the unit operator. Assuming the defect phase  $\vartheta$  is chosen to lie to the right of all

ordinary BPS states, the generating functions for  $L_{i}$ 's are [27,65],

$$
F (L _ {1}) = X _ {\gamma_ {1}},
$$

$$
F (L _ {2}) = X _ {\gamma_ {2}} + X _ {\gamma_ {1} + \gamma_ {2}},
$$

$$
F \left(L _ {3}\right) = X _ {- \gamma_ {1}} + X _ {- \gamma_ {1} + \gamma_ {2}} + X _ {\gamma_ {2}}, \tag {5.5}
$$

$$
F (L _ {4}) = X _ {- \gamma_ {1} - \gamma_ {2}} + X _ {- \gamma_ {1}},
$$

$$
F (L _ {5}) = X _ {- \gamma_ {2}}.
$$

It follows that the five  $L_{i}$ 's satisfy the OPE algebra,

$$
L _ {i} L _ {i + 2} = 1 + q ^ {\frac {1}{2}} L _ {i + 1}, \tag {5.6}
$$

which also implies  $L_{i}L_{i - 2} = 1 + q^{-\frac{1}{2}}L_{i - 1}$ . We have taken the index  $i$  to be periodic mod 5.

Following our proposal (1.10), the Schur index with the insertion of  $L_{i}$  is computed by the trace of  $F(L_{i})$ , i.e.,  $\mathcal{I}_{L_i} = (q)_{\infty}^2\mathrm{Tr}[F(L_i)\mathcal{S}(q)\overline{\mathcal{S}} (q)]$ . After a similar calculation as in the pure  $SU(2)$  gauge theory,  $\mathcal{I}_L$  can be computed to be

$$
\begin{array}{l} \mathcal {I} _ {L _ {i}} (q) = (q) _ {\infty} ^ {2} \operatorname {T r} \left[ F (L _ {i}) E _ {q} \left(X _ {\gamma_ {1}}\right) E _ {q} \left(X _ {\gamma_ {2}}\right) E _ {q} \left(X _ {\gamma_ {1}} ^ {- 1}\right) E _ {q} \left(X _ {\gamma_ {2}} ^ {- 1}\right) \right] \tag {5.7} \\ = - q ^ {\frac {1}{2}} \left(1 + q ^ {3} + q ^ {4} + q ^ {5} + q ^ {6} + q ^ {7} + 2 q ^ {8} + 2 q ^ {9} + 3 q ^ {1 0} + \dots\right), \quad \text {f o r a l l} i = 1, \dots , 5. \\ \end{array}
$$

Notice that the dependence on the index  $i$  is washed out inside the trace, reflecting the  $\mathbb{Z}_5$  symmetry of the  $A_{2}$  Argyres-Douglas theory. We can therefore define  $\mathcal{I}_L$  unambiguously as,

$$
\mathcal {I} _ {L} \equiv \mathcal {I} _ {L _ {1}} = \mathcal {I} _ {L _ {2}} = \dots = \mathcal {I} _ {L _ {5}}, \tag {5.8}
$$

More explicitly, for example when  $i = 1$ , the trace of  $L_{i}$  can written as

$$
\begin{array}{l} \mathcal {I} _ {L} = (q) _ {\infty} ^ {2} \sum_ {\ell_ {1}, \ell_ {2} = 0} ^ {\infty} \frac {(- 1) q ^ {\ell_ {1} + 2 \ell_ {2} + \ell_ {1} \ell_ {2} + \frac {1}{2}}}{(q) _ {\ell_ {1}} [ (q) _ {\ell_ {2}} ] ^ {2} (q) _ {\ell + 1}} \tag {5.9} \\ = - q ^ {- \frac {1}{2}} \left(1 + q ^ {3} + q ^ {4} + q ^ {5} + q ^ {6} + q ^ {7} + 2 q ^ {8} + 2 q ^ {9} + 3 q ^ {1 0} + \dots\right), \quad \text {f o r a l l} i = 1, \dots , 5. \\ \end{array}
$$

We further observe the following relations among the line defect Schur indices (no sum

in the indices),

$$
\mathcal {I} _ {L _ {i} L _ {i - 2}} = \mathcal {I} + q ^ {- \frac {1}{2}} \mathcal {I} _ {L},
$$

$$
\mathcal {I} _ {L _ {i} L _ {i - 1}} = q ^ {- 1} \mathcal {I} + q ^ {- \frac {3}{2}} \mathcal {I} _ {L},
$$

$$
\mathcal {I} _ {L _ {i} L _ {i}} = q ^ {- 1} \mathcal {I} + q ^ {- \frac {3}{2}} \mathcal {I} _ {L}, \tag {5.10}
$$

$$
\mathcal {I} _ {L _ {i} L _ {i + 1}} = \mathcal {I} + q ^ {- \frac {1}{2}} \mathcal {I} _ {L},
$$

$$
\mathcal {I} _ {L _ {i} L _ {i + 2}} = \mathcal {I} + q ^ {\frac {1}{2}} \mathcal {I} _ {L},
$$

which hold true for any  $i = 1, \dots, 5$ . Here  $\mathcal{I}$  is the Schur index without the insertion of line defects. Note that the first and the last relations simply follow from the UV line defect OPE (5.6). The third relation was already noticed in [51] in the case of the inverse of the quantum KS operator, and a connection with the Verlinde algebra was observed. We will give a similar proposal in Section 6.

# 5.3  $A_{3}$  Argyres-Douglas Theory

![](images/508676d9c7f068c90708044d11183f23a69d1070af6de6b5922fd969e8c5dd55.jpg)  
Figure 9: The BPS quiver for the  $A_{3}$  Argyres-Douglas theory.

The  $A_{3}$  Argyres-Douglas theory arises from special points on the moduli space of the  $SU(2)$  SQCD with  $N_{f} = 2$  flavors [70]. The  $A_{3}$  Argyres-Douglas theory has an  $SU(2)$  flavor symmetry which corresponds to the direction  $\gamma_{1} - \gamma_{3}$  in the charge lattice.

The core charges of the line defects in the  $A_{3}$  Argyres-Douglas theory can be derived from the seeds of the BPS quiver and their dual cones [27, 65]. We present the details of the derivation of the framed BPS quivers in Appendix B.1 following the logic of [65]. The result is that the line defects in the  $A_{3}$  Argyres-Douglas theory are generated by six defects,  $A_{i}$ ,  $B_{i}$ ,  $i = 1, 2, 3$ , one flavor defect  $C$ , and the unit operator.

The generating functions for these line defects can then be obtained straightforwardly from the associated framed BPS quivers, say, by the mutation method [79]. Assuming the defect phase  $\vartheta$  is to the right of all the ordinary BPS state phases, the generating functions

for the above six defects are

$$
F (A _ {1}) = X _ {\gamma^ {\prime}},
$$

$$
F (A _ {2}) = X _ {- \gamma^ {\prime} - \gamma_ {2}} + X _ {- \gamma^ {\prime}},
$$

$$
F (A _ {3}) = X _ {- \gamma^ {\prime}} + (z + z ^ {- 1}) X _ {\gamma_ {2}} + X _ {\gamma^ {\prime} + \gamma_ {2}} + X _ {- \gamma^ {\prime} + \gamma_ {2}},
$$

$$
F (B _ {1}) = X _ {- \gamma_ {2}},
$$

$$
F (B _ {2}) = X _ {\gamma_ {2}} + (z + z ^ {- 1}) X _ {- \gamma^ {\prime}} + (z + z ^ {- 1}) X _ {- \gamma^ {\prime} + \gamma_ {2}} + (q ^ {\frac {1}{2}} + q ^ {- \frac {1}{2}}) X _ {- 2 \gamma^ {\prime}} + X _ {- 2 \gamma^ {\prime} + \gamma_ {2}} + X _ {- 2 \gamma^ {\prime} - \gamma_ {2}},
$$

$$
F (B _ {3}) = X _ {\gamma_ {2}} + (z + z ^ {- 1}) X _ {\gamma^ {\prime} + \gamma_ {2}} + X _ {2 \gamma^ {\prime} + \gamma_ {2}},
$$

$$
F (C) = z + z ^ {- 1}. \tag {5.11}
$$

Here  $z$  is the  $SU(2)$  flavor fugacity that is related to the flavor generator as

$$
\operatorname {T r} \left[ X _ {\frac {\gamma_ {1} - \gamma_ {3}}{2}} \right] = z. \tag {5.12}
$$

$\gamma^\prime$  is defined as

$$
\gamma^ {\prime} \equiv \frac {\gamma_ {1} + \gamma_ {3}}{2}. \tag {5.13}
$$

Note that  $\gamma^\prime$  has the same Dirac pairings with every charge vector as those of  $\gamma_{1}$  and  $\gamma_{3}$ . Note that the UV flavor defect  $C$  is a Wilson line in the 2 of the flavor  $SU(2)$  symmetry, whose insertion into the path integral is just an overall multiplication of  $z + z^{-1}$ .

The defect OPE can be readily derived from the above generating functions,

$$
A _ {i} A _ {i + 1} = 1 + q ^ {- \frac {1}{2}} B _ {i},
$$

$$
B _ {i} B _ {i + 1} = 1 + q ^ {- \frac {1}{2}} C A _ {i + 1} + q ^ {- 1} A _ {i + 1} ^ {2}, \tag {5.14}
$$

$$
A _ {i} B _ {i + 1} = C + q ^ {- \frac {1}{2}} A _ {i + 1} + q ^ {\frac {1}{2}} A _ {i + 2},
$$

with the flavor defect  $C$  commuting with everything. Here we view the index  $i$  as periodic mod 3.

To compute the line defect indices of  $A_{i},B_{i}$  , we will repeatedly encounter the insertion

of an IR line defect  $X_{a\gamma' + b\gamma_2}$

$$
\begin{array}{l} (q) _ {\infty} ^ {2} \operatorname {T r} \left[ X _ {a \gamma^ {\prime} + b \gamma_ {2}} \mathcal {S} (q) \overline {{\mathcal {S}}} (q) \right] \\ = (q) _ {\infty} ^ {2} \operatorname {T r} \left[ X _ {a \gamma^ {\prime} + b \gamma_ {2}} E _ {q} \left(X _ {\gamma^ {\prime}}\right) E _ {q} \left(X _ {\gamma_ {3}}\right) E _ {q} \left(X _ {\gamma_ {2}}\right) E _ {q} \left(X _ {\gamma^ {\prime}} ^ {- 1}\right) E _ {q} \left(X _ {\gamma_ {3}} ^ {- 1}\right) E _ {q} \left(X _ {\gamma_ {2}} ^ {- 1}\right) \right] \\ = (q) _ {\infty} ^ {2} \sum_ {\substack {\ell_ {1}, \ell_ {2}, \ell_ {3}, \\ k _ {1}, k _ {2}, k _ {3} = 0}} ^ {\infty} \frac {(- 1) ^ {a + b} q ^ {\ell_ {1} + \ell_ {2} + \ell_ {3} + \ell_ {2} (\ell_ {1} + \ell_ {3}) + \frac {1}{2} (a + b + a b) + a \ell_ {2}}}{(q) _ {\ell_ {1}} (q) _ {\ell_ {2}} (q) _ {\ell_ {3}} (q) _ {k _ {1}} (q) _ {k _ {2}} (q) _ {k _ {3}}} z ^ {2 (l _ {1} - k _ {1}) + a} \delta_ {k _ {2}, \ell_ {2} + b} \delta_ {k _ {1} + k _ {3}, \ell_ {1} + \ell_ {3} + a}. \tag{5.15} \\ \end{array}
$$

Following the conjecture (5.3), together with (5.11) and (5.15), we obtain the line defect Schur indices for  $A_{i}$  and  $B_{i}$ ,

$$
\begin{array}{l} \mathcal {I} _ {A _ {i}} (q, z) = (q) _ {\infty} ^ {2} \operatorname {T r} [ F (A _ {i}) \mathcal {S} (q) \overline {{\mathcal {S}}} (q) ] \\ = - q ^ {\frac {1}{2}} \left[ \chi_ {2} + \chi_ {4} q + \left(\chi_ {2} + \chi_ {4} + \chi_ {6}\right) q ^ {2} + \left(2 \chi_ {2} + 2 \chi_ {4} + \chi_ {6} + \chi_ {8}\right) q ^ {3} \right. \\ \left. + \left(3 \chi_ {\mathbf {2}} + 3 \chi_ {\mathbf {4}} + 3 \chi_ {\mathbf {6}} + \chi_ {\mathbf {8}} + \chi_ {\mathbf {1 0}}\right) q ^ {4} + \dots \right], \quad \text {f o r a l l} \quad i = 1, 2, 3, \tag {5.16} \\ \end{array}
$$

$$
\begin{array}{l} \mathcal {I} _ {B _ {i}} (q, z) = (q) _ {\infty} ^ {2} \operatorname {T r} [ F (B _ {i}) \mathcal {S} (q) \overline {{\mathcal {S}}} (q) ] \\ = - q ^ {\frac {1}{2}} \left[ \chi_ {1} + \chi_ {3} q ^ {2} + \left(\chi_ {1} + \chi_ {3}\right) q ^ {3} + \left(\chi_ {1} + \chi_ {3} + \chi_ {5}\right) q ^ {4} + \left(\chi_ {1} + 2 \chi_ {3} + \chi_ {5}\right) q ^ {5} + \dots \right], \\ \text {f o r a l l} i = 1, 2, 3, \tag {5.17} \\ \end{array}
$$

where  $\chi_{\mathbf{n}}$  is the character for the  $n$ -dimensional irreducible representation of  $SU(2)$ , normalized such that  $\chi_{\mathbf{2}} = z + z^{-1}$ . As in the  $A_{2}$  Argyres-Douglas theory, the dependence on the index  $i$  is washed out inside the trace. We can therefore define  $\mathcal{I}_A$  and  $\mathcal{I}_B$  unambiguously as

$$
\mathcal {I} _ {A} \equiv \mathcal {I} _ {A _ {1}} = \mathcal {I} _ {A _ {2}} = \mathcal {I} _ {A _ {3}},
$$

$$
\mathcal {I} _ {B} \equiv \mathcal {I} _ {B _ {1}} = \mathcal {I} _ {B _ {2}} = \mathcal {I} _ {B _ {3}}.
$$

We observe the following relations between the line defect indices (no sum in the indices),

$$
\mathcal {I} _ {A _ {i} A _ {i}} = \mathcal {I} _ {A _ {j} A _ {j + 1}} = \mathcal {I} + q ^ {- \frac {1}{2}} \mathcal {I} _ {B},
$$

$$
\mathcal {I} _ {B _ {i} B _ {i}} = \mathcal {I} _ {B _ {j} B _ {j + 1}} = \mathcal {I} + q ^ {- \frac {1}{2}} (z + z ^ {- 1}) \mathcal {I} _ {A} + q ^ {- 1} \mathcal {I} _ {A _ {k} A _ {k}}, \tag {5.18}
$$

$$
q \mathcal {I} _ {A _ {i} B _ {i}} = \mathcal {I} _ {A _ {j} B _ {j + 2}} = \mathcal {I} _ {A _ {k} B _ {k + 1}} = (z + z ^ {- 1}) \mathcal {I} + (q ^ {\frac {1}{2}} + q ^ {- \frac {1}{2}}) \mathcal {I} _ {A},
$$

which hold true for all values of  $i, j, k = 1, 2, 3$ . Note that the rightmost equalities in the above relations are implied by the UV line defect OPE (5.14).

# 5.4  $A_{4}$  Argyres-Douglas Theory

![](images/04804d96b05d0cf68dacaeec50e613e7cc18f21f6270616ac61e5f7cbc1244d2.jpg)  
Figure 10: The BPS quiver for the  $A_4$  Argyres-Douglas theory.

For the  $A_4$  Argyres-Douglas theory, we present a similar (but much more tedious) derivation of the core charges of the line defects in Appendix B.2 again following the method of [65]. The result is that there are fourteen generators  $A_i$ ,  $B_i$  ( $i = 1, \dots, 7$ ) together with the unit operator for the defects. Assuming the defect phase  $\vartheta$  is chosen to lie to the right of all ordinary BPS states, the generating functions for the line defects are,

$$
F (A _ {1}) = X _ {- \gamma_ {2} + \gamma_ {4}},
$$

$$
F (A _ {2}) = X _ {\gamma_ {2} - \gamma_ {4}} + X _ {\gamma_ {1} + \gamma_ {2} - \gamma_ {4}},
$$

$$
F (A _ {3}) = X _ {- \gamma_ {1} - \gamma_ {2} + \gamma_ {4}} + X _ {- \gamma_ {1} + \gamma_ {4}} + X _ {- \gamma_ {1} + \gamma_ {3} + \gamma_ {4}},
$$

$$
F (A _ {4}) = X _ {\gamma_ {1} - \gamma_ {3} - \gamma_ {4}} + X _ {\gamma_ {1} - \gamma_ {3}},
$$

$$
F (A _ {5}) = X _ {- \gamma_ {1} + \gamma_ {3}},
$$

$$
F (A _ {6}) = X _ {\gamma_ {1} - \gamma_ {3}} + X _ {\gamma_ {1} - \gamma_ {3} + \gamma_ {4}} + X _ {\gamma_ {1} + \gamma_ {4}},
$$

$$
F (A _ {7}) = X _ {- \gamma_ {1} - \gamma_ {4}} + X _ {- \gamma_ {1} + \gamma_ {2} - \gamma_ {4}} + X _ {\gamma_ {2} - \gamma_ {4}},
$$

$$
F \left(B _ {1}\right) = X _ {- \gamma_ {1}} + X _ {\gamma_ {2}} + X _ {\gamma_ {2} - \gamma_ {1}} + X _ {\gamma_ {2} + \gamma_ {3}} + X _ {- \gamma_ {1} + \gamma_ {2} + \gamma_ {3}}, \tag {5.19}
$$

$$
F (B _ {2}) = X _ {- \gamma_ {3} - \gamma_ {4}} + X _ {\gamma_ {2}} + X _ {\gamma_ {1} + \gamma_ {2}} + X _ {\gamma_ {2} - \gamma_ {3}} + X _ {\gamma_ {1} + \gamma_ {2} - \gamma_ {3}} + X _ {- \gamma_ {3}} + X _ {\gamma_ {2} - \gamma_ {3} - \gamma_ {4}} + X _ {\gamma_ {1} + \gamma_ {2} - \gamma_ {3} - \gamma_ {4}},
$$

$$
F \left(B _ {3}\right) = X _ {- \gamma_ {2} - \gamma_ {3}} + X _ {- \gamma_ {3}} + X _ {\gamma_ {4}} + X _ {\gamma_ {4} - \gamma_ {3}} + X _ {- \gamma_ {2} - \gamma_ {3} + \gamma_ {4}},
$$

$$
F (B _ {4}) = X _ {- \gamma_ {1} - \gamma_ {2}} X _ {- \gamma_ {1}},
$$

$$
F (B _ {5}) = X _ {- \gamma_ {4}},
$$

$$
F (B _ {6}) = X _ {\gamma_ {1}},
$$

$$
F (B _ {7}) = X _ {\gamma_ {4}} + X _ {\gamma_ {3} + \gamma_ {4}}.
$$

The fourteen generators for the line defects satisfy the following defect OPE,

$$
A _ {i} A _ {i + 1} = 1 + q ^ {\frac {1}{2}} B _ {2 i + 4},
$$

$$
B _ {i} B _ {i + 2} = 1 + q ^ {\frac {1}{2}} A _ {4 i - 1} B _ {i + 1},
$$

$$
B _ {i} B _ {i + 3} = A _ {4 i - 4} A _ {4 i - 1} + A _ {4 i + 1}, \tag {5.20}
$$

$$
A _ {i} B _ {2 i} = q ^ {\frac {1}{2}} A _ {i - 2} + B _ {2 i + 1},
$$

$$
A _ {i} B _ {2 i - 1} = q ^ {- \frac {1}{2}} A _ {i + 2} + B _ {2 i - 2}.
$$

We have taken the index  $i$  to be periodic mod 7.

The Schur indices with insertions of  $A_{i}$  and  $B_{i}$  computed from our IR formula (5.3) are,

$$
\begin{array}{l} \mathcal {I} _ {A _ {i}} (q) = (q) _ {\infty} ^ {4} \operatorname {T r} \left[ F (A _ {i}) \mathcal {S} (q) \overline {{\mathcal {S}}} (q) \right] \\ = q + q ^ {4} + q ^ {5} + q ^ {6} + 2 q ^ {7} + 2 q ^ {8} + 3 q ^ {9} + 4 q ^ {1 1} + 5 q ^ {1 2} + 7 q ^ {1 3} + 8 q ^ {1 4} + \mathcal {O} (q ^ {1 5}), \\ \end{array}
$$

for all  $i = 1,\dots ,7$  (5.21)

$$
\begin{array}{l} \mathcal {I} _ {B _ {i}} (q) = (q) _ {\infty} ^ {4} \operatorname {T r} \left[ F (B _ {i}) \mathcal {S} (q) \overline {{\mathcal {S}}} (q) \right] \\ = - q ^ {\frac {1}{2}} \left(1 + q ^ {2} + q ^ {3} + q ^ {4} + 2 q ^ {5} + 3 q ^ {6} + 3 q ^ {7} + 4 q ^ {8} + 5 q ^ {9} + 7 q ^ {1 0} + 8 q ^ {1 1} + 1 1 q ^ {1 2} + \mathcal {O} \left(q ^ {1 3}\right)\right), \\ \end{array}
$$

for all  $i = 1,\dots ,7$  (5.22)

where for the  $A_4$  Argyres-Douglas theory

$$
\mathcal {S} (q) \overline {{\mathcal {S}}} (q) = E _ {q} \left(X _ {\gamma_ {1}}\right) E _ {q} \left(X _ {\gamma_ {2}}\right) E _ {q} \left(X _ {\gamma_ {3}}\right) E _ {q} \left(X _ {\gamma_ {4}}\right) E _ {q} \left(X _ {\gamma_ {1}} ^ {- 1}\right) E _ {q} \left(X _ {\gamma_ {2}} ^ {- 1}\right) E _ {q} \left(X _ {\gamma_ {3}} ^ {- 1}\right) E _ {q} \left(X _ {\gamma_ {4}} ^ {- 1}\right). \tag {5.23}
$$

Notice that the dependence on the index  $i = 1,\dots ,7$  is washed out inside the trace, reflecting the  $\mathbb{Z}_7$  symmetry of the  $A_4$  Argyres-Douglas theory. We can therefore define  $\mathcal{I}_A(q)$  and  $\mathcal{I}_B(q)$  unambiguously as  $(q)_{\infty}^{6}\mathrm{Tr}\left[F(A_{i})\mathcal{S}(q)\overline{\mathcal{S}} (q)\right]$  and  $(q)_{\infty}^{6}\mathrm{Tr}\left[F(B_{i})\mathcal{S}(q)\overline{\mathcal{S}} (q)\right]$  for any choice of  $i$ , respectively.

# 6 Chiral Algebra and Line Defects

To every  $4d\mathcal{N} = 2$  superconformal field theory, we can associate a  $2d$  chiral algebra a la the work of [13]. The states in the vacuum module of the  $2d$  chiral algebra are in one-to-one correspondence with the protected operators in the  $4d$  theory that contribute to the Schur index. It is then natural to ask whether one has access to the states in the other modules of the chiral algebra from the  $4d$  physics.

In Appendix A, we show that when the line defects are extended on a plane, say the 12-plane, transverse to the chiral algebra plane, say, the 34-plane, the combined system preserves two supercharges (A.17). The incidence geometry of the line defects and the chiral algebra plane is shown in Figure 11. Given that the line defects share some common supercharges with the chiral algebra plane, it is tempting to speculate that the defect operators are related to the states in the other modules of the chiral algebra.

In this section we demonstrate in several examples, including the Argyres-Douglas theories and  $SU(2)$  gauge theory with  $N_{f} = 4$  flavors, that the line defect indices can indeed be written as linear combinations of the characters for the other modules in the chiral algebra. This should not come as a surprise in the case of Lagrangian theories. Indeed, in Section 6.1 we show explicitly how the defect Schur operators at the zero coupling point of

![](images/8b46a9a26af19dee8f9c83f643076a9fd71562af839d0bfe830ee56e94480864.jpg)  
Figure 11: The incidence geometry of line defects and the chiral algebra plane. We suppress one dimension of the chiral algebra plane (the 34-plane) and represent it as the red line above. The black lines are the line defects lying on the 12-plane. The line defects are oriented on the 12-plane by their central charge phases  $\vartheta_{i}$ . Here  $\vartheta_{ij} = \vartheta_{i} - \vartheta_{j}$ .

$SU(2)$  SQCD organize themselves into modules of the associated chiral algebra. However a general abstract derivation which holds also for non-Lagrangian theories is still lacking and is an open problem for future research. It would also be interesting to explore the modular properties of the resulting sums of characters such as (6.1) and (6.20) as in [19-21,89].

In the case of Argyres-Douglas theories, we go further and describe how the fusion rule, or the Verlinde algebra, of the  $2d$  chiral algebra can be realized from line defect indices in  $4d$  in the  $q \to 1$  limit. We will start with a general proposal in Section 6.2 and demonstrate it with the examples of the  $A_2$ ,  $A_3$ , and  $A_4$  Argyres-Douglas theories. This connection between the  $2d$  Verlinde algebra and the  $4d$  defect indices was first observed in [51].

# 6.1  $SU(2)$  Gauge Theory with  $N_{f} = 4$  Flavors

The chiral algebra associated to the  $4d$ $SU(2)$  gauge theory with  $N_{f} = 4$  flavors is  $\widehat{so(8)}_{-2}$  [13]. The Schur index without any insertion of defect is reproduced by the vacuum character of  $\widehat{so(8)}_{-2}$ . It is natural to speculate that the defect indices discussed in Section 3.2.3 can be related to the other characters of  $\widehat{so(8)}_{-2}$ . Indeed, we observe that the line defect index for a half Wilson in the 2 can be written as the following linear combinations of the characters

for  $\widehat{so(8)}_{-2}$ ,

$$
\mathcal {I} _ {L _ {0, 1}} (q, \eta_ {i} = 1) = \sum_ {k = 1} ^ {\infty} (- 1) ^ {k} q ^ {\frac {k ^ {2} + k - 1}{2}} (1 - q ^ {k}) \chi_ {[ - 2 k - 1, 2 k - 1, 0, 0, 0 ]} (q, \eta_ {i} = 1) \tag {6.1}
$$

where  $\chi_{[a_0,a_1,a_2,a_3,a_4]}(q,\eta_i)$  is the affine character $^{15}$  of  $\widehat{so(8)}_{-2}$  with affine Dynkin labels  $[a_0,a_1,a_2,a_3,a_4]$ . We have normalized the  $\widehat{so(8)}_{-2}$  affine characters to start from 1. We present the details of the calculation for the affine characters of  $\widehat{so(8)}_{-2}$  in Appendix C. The  $SO(8)$  flavor fugacities have been set to 1 for simplicity and we have checked the above relation to  $\mathcal{O}(q^{\frac{19}{2}})$ .

The vacuum representation of  $\widehat{so(8)}_{-2}$  has certain null states that can be translated into null relations among the Kac-Moody currents  $J^{[ij]}(z)$ . It is perhaps more natural to consider representations that respect these null relations. Among all the highest weight irreducible representations of  $\widehat{so(8)}_{-2}$ , there are four of them, in addition to the vacuum representation, that respect these null relations [19]. Their highest weights are  $[0, -2, 0, 0, 0]$ ,  $[0, 0, -1, 0, 0]$ ,  $[0, 0, 0, -2, 0]$ , and  $[0, 0, 0, 0, -2]$ . We do not know whether the line defect indices can be decomposed into sum of these four characters together with the vacuum character. It would be very interesting to pursue this direction further.

The appearance of non-vacuum characters in the Wilson line indices of Lagrangian theories is not a surprise. Indeed, by going to the free theory point on the moduli space, it is easy to see why the line defect Schur operators organize themselves into modules of  $\widehat{so(8)}_{-2}$ . Let us denote the holomorphic operators in the chiral algebra by the same symbols as their  $4d$  Schur operators. In the case of  $SU(2)$  gauge theory with  $N_{f} = 4$  flavors, we have a fermionic current  $\rho_{\pm}^{A}(z)$  with dimension 1 from the vector multiplet, and a bosonic current  $H^{ia}(z)$  with dimension  $\frac{1}{2}$  from the hypermultiplets. Here  $a$ ,  $A$ , and  $i$  are the indices for the  $SU(2)$  doublet,  $SU(2)$  triplet, and the  $\mathbf{8}_v$  of  $SO(8)$ , respectively. The dimensions of  $H^{ia}$  and  $\rho^{\alpha A}$  are determined by their  $4d$  quantum numbers  $\Delta - R$ , the exponent of  $q$  in the Schur index. Their  $2d$  OPEs can be obtained from the  $4d$  OPEs following the construction of the chiral algebra in [13],

$$
\rho^ {\alpha A} (z) \rho^ {\beta B} (0) \sim \frac {\delta^ {A B} \epsilon^ {\alpha \beta}}{z ^ {2}},
$$

$$
H ^ {i a} (z) H ^ {j b} (0) \sim \frac {\delta^ {i j} \epsilon^ {a b}}{z}, \tag {6.2}
$$

$$
H ^ {i a} (z) \rho^ {\beta B} (0) \sim 0,
$$

where  $\epsilon^{\alpha \beta}$  is the two-by-two antisymmetric tensor normalized such that  $\epsilon^{+ - } = +1$ . We will

also use  $\epsilon_{\alpha \beta}$  which is the inverse of  $\epsilon^{\alpha \beta}$ , i.e.  $\epsilon_{+-} = -1$ .

We can then construct the dimension 1  $\widehat{so(8)}_{-2}$  currents  $J^{[ij]}(z)$ ,

$$
J ^ {[ i j ]} (z) = \mathcal {N} \epsilon_ {a b}: H ^ {i a} H ^ {j b}: (z), \tag {6.3}
$$

where  $\mathcal{N}$  is a normalization constant. Using (6.2), one can show that  $J^{[ij]}(z)$  satisfies the so(8) current algebra at level  $-2$ .

The half Wilson line index  $\mathcal{I}_{L_{0,1}}$  counts "gauge non-invariant" Schur operators in the doublets of  $SU(2)$ . As discussed in Section 3.2.3, we can enumerate these operators at the free theory point, which are normal-ordered operators in the doublets made out of  $H^{ia}$  and  $\rho^{\alpha A}$ . Now since the  $\widehat{so(8)}_{-2}$  current  $J^{[ij]}(z)$  is a singlet under  $SU(2)$ , any descendant  $J_{-n_1}^{[i_1j_1]}\cdots J_{-n_k}^{[i_kj_k]}\cdot \mathcal{O}^a (z)$  of a doublet  $\mathcal{O}^a (z)$  is still a doublet made out of  $H^{ia}$  and  $\rho^{\alpha A}$ . Hence if  $\mathcal{O}^a (z)$  is counted by the line defect index, so are all its current algebra descendants. It follows that the line defect Schur operators fall into modules of the chiral algebra  $\widehat{so(8)}_{-2}$ .

We can see explicitly how this works for the low dimension operators. At dimension  $1/2$ , the only line defect Schur operator is  $H^{ia}(z)$ . It is a current algebra primary since it is annihilated by all the  $J_{+n}^{[ij]}$  with  $n > 0$ .

At dimension  $3/2$ , the Schur operators are  $\partial H^{ia}$ ,  $(\rho^{\pm A}H^{ib})(T^{A})_{b}^{a}$ ,  $(H^{3})_{8v}$ , and  $(H^{3})_{160v}$  listed in (3.25). Here  $(H^{3})_{\mathbf{R}}$  denotes the part of the normal-ordered operator<sup>16</sup>:  $H^{ia}H^{jb}H^{kc}$ : that transforms in the representation  $\mathbf{R}$  of  $SO(8)$  and in the doublet of  $SU(2)$ . Let us organize these dimension  $3/2$  operators into primaries and descendants of the current algebra. The level 1 current algebra descendant of the dimension  $1/2$  primary  $H^{ia}(z)$  is

$$
\begin{array}{l} J _ {- 1} ^ {[ i j ]} \cdot H ^ {k a} (0) = \oint \frac {d z}{2 \pi i} \frac {1}{z} J ^ {[ i j ]} (z) H ^ {k a} (0) \\ = \mathcal {N} \left[ - \delta^ {i k} \partial H ^ {j a} (0) + \delta^ {j k} \partial H ^ {i a} (0) + \epsilon_ {b c}: H ^ {i b} H ^ {j c} H ^ {k a}: (0) \right]. \\ \end{array}
$$

This descendant, when decomposed into irreducible representations, contains  $(H^3)_{\mathbf{160}_v}$  and a particular linear combination of  $(H^3)_{\mathbf{8}_v}$  and  $\partial H^{ia}$ :

$$
- \epsilon_ {b c}: H ^ {j a} H ^ {j b} H ^ {i c} + 7 \partial H ^ {i a}. \tag {6.5}
$$

After having determined the descendants, we need to identify the primaries at dimension

3/2. To begin with, it is easy to see that  $(\rho^{\pm A}H^{ib})(T^{A})_{b}^{a}$  is annihilated by all the  $J_{+n}^{[ij]}$  with  $n > 0$ , and is hence a primary of the current algebra. Next, one can straightforwardly check that the following combination

$$
\epsilon_ {b c}: H ^ {j a} H ^ {j b} H ^ {i c}: (0) - 3 \partial H ^ {i a} (0) \tag {6.6}
$$

is a current algebra primary.

In summary, we have organized the dimension  $1/2$  and  $3/2$  Wilson line defect Schur operators into current algebra descendants and primaries,

<table><tr><td>Dimension</td><td>Chiral Algebra Operator</td><td>so(8)_2</td><td>SO(8) Rep</td></tr><tr><td>1/2</td><td>Hia</td><td>Primary</td><td>8v</td></tr><tr><td>3/2</td><td>εbcHjaHjbHic-3δHia</td><td>Primary</td><td>8v</td></tr><tr><td>3/2</td><td>(ρ±A Hib)(TA)b</td><td>Primary</td><td>8v</td></tr><tr><td>3/2</td><td>-εbcHjaHjbHic+7δHia</td><td>Descendant</td><td>8v</td></tr><tr><td>3/2</td><td>(H3)160v</td><td>Descendant</td><td>160v</td></tr></table>

Returning to the doublet half Wilson line index (6.1), we have, up to  $q^{3/2}$

$$
\mathcal {I} _ {L _ {0, 1}} = - q ^ {\frac {1}{2}} \chi_ {\mathbf {8} _ {v}} + q ^ {\frac {3}{2}} \chi_ {\mathbf {8} _ {v}} + \mathcal {O} \left(q ^ {\frac {5}{2}}\right). \tag {6.8}
$$

Recall that  $\chi_{\mathbf{8}_v} = 8 + 168q + \dots$ . It is now clear that the first term in  $\mathcal{I}_{L_{0,1}}$  is the contribution from the current algebra primary  $H^{ia}$  and its descendants, while the second term is the contribution from the bosonic primary  $\epsilon_{bc}H^{ja}H^{jb}H^{ic} - 3\partial H^{ia}$  and the two fermionic primaries  $(\rho^{\pm A}H^{ib})(T^{A})_{b}^{a}$  as well as their descendants. This analysis gives a direct understanding of why the non-vacuum characters of the chiral algebra appears in the line defect indices in the case of Lagrangian theories.

Similarly, the other half Wilson line defect indices are also related to the  $\widehat{so(8)}_{-2}$  affine characters,

$$
\mathcal {L} _ {L _ {0, 2}} (q, \eta_ {i} = 1) = - q \chi_ {[ - 2, 0, 0, 0, 0 ]} + (1 - q) \sum_ {k = 1} ^ {\infty} (- 1) ^ {k + 1} q ^ {\frac {k (k + 1)}{2}} (1 - q ^ {2 k}) \chi_ {[ - 2 k - 2, 2 k, 0, 0, 0 ]},
$$

$$
\begin{array}{l} \mathcal {I} _ {L _ {0, 3}} (q, \eta_ {i} = 1) = q ^ {\frac {3}{2}} (1 - q ^ {2}) \chi_ {[ - 3, 1, 0, 0, 0 ]} - q ^ {\frac {3}{2}} (1 - q + q ^ {3} - q ^ {6}) \chi_ {[ - 5, 3, 0, 0, 0 ]} \\ + q ^ {\frac {7}{2}} (1 - q ^ {2} + q ^ {5}) \chi_ {[ - 7, 5, 0, 0, 0 ]} - q ^ {\frac {1 3}{2}} (1 - q ^ {3}) \chi_ {[ - 9, 7, 0, 0, 0 ]} + q ^ {\frac {2 1}{2}} \chi_ {[ - 1 1, 9, 0, 0, 0 ]} + \mathcal {O} (q ^ {\frac {2 5}{2}}), \\ \end{array}
$$

$$
\begin{array}{l} \mathcal {L} _ {L _ {0, 4}} (q, \eta_ {i} = 1) = q ^ {3} \chi_ {[ - 2, 0, 0, 0, 0 ]} - q ^ {2} (1 - q ^ {2} - q ^ {3} + q ^ {5}) \chi_ {[ - 4, 2, 0, 0, 0 ]} \\ + q ^ {2} (1 - q + q ^ {5} - q ^ {6} - q ^ {8} + q ^ {1 0}) \chi_ {[ - 6, 4, 0, 0, 0 ]} - q ^ {4} (1 - 2 q ^ {2} + q ^ {3} + q ^ {8}) \chi_ {[ - 6, 4, 0, 0, 0 ]} \\ + q ^ {7} \left(1 - q ^ {2} - q ^ {3} + q ^ {4}\right) \chi_ {[ - 1 0, 8, 0, 0, 0 ]} - q ^ {1 1} \chi_ {[ - 1 2, 1 0, 0, 0, 0 ]} + \mathcal {O} \left(q ^ {1 3}\right), \tag {6.9} \\ \end{array}
$$

where  $\mathcal{I}_{L_{0,n}}$  is the line defect index for a half Wilson line in the  $(n + 1)$ -dimensional representation of  $SU(2)$ . For the half Wilson line defect indices  $\mathcal{I}_{L_{0,2}}$  in the  $\mathbf{3}$ , we conjectured the above exact relation and checked it to  $\mathcal{O}(q^{14})$ . For the other two cases we recorded our observation above without closed form formulas. Here we have normalized the  $\widehat{so(8)}_{-2}$  characters to start from order  $q^0$ .

# 6.2 Verlinde Algebra from Line Defects

In this subsection we give a precise proposal on how the fusion rules in the two-dimensional chiral algebra can be realized from the four-dimensional line defect indices in the case of Argyres-Douglas theories, generalizing the results of [51].

In all known chiral algebras for the Argyres-Douglas theories, $^{17}$  there is always a distinguished finite set of primaries whose characters form a modular vector. For example, the chiral algebra of the  $A_{2n}$  Argyres-Douglas theory is the  $(2,2n + 3)$  Virasoro minimal model. As another example that we will encounter in this section, the chiral algebra for the  $A_{3}$  Argyres-Douglas theory is the affine Lie algebra  $\widehat{su(2)}_{-\frac{4}{3}}$ . There are three distinguished representations of  $\widehat{su(2)}_{-\frac{4}{3}}$ , called the admissible representation, whose characters are linearly mapped to each other under the modular transformation. We will focus on these primaries  $\Phi_{\alpha}$  with this nice modular property.  $\Phi_0$  is chosen to be the identity.

For the Argyres-Douglas theories, we claim that the generators for the UV electromagnetic line defects can be labeled by the non-identity primaries  $\Phi_{\alpha \neq 0}$  with a (finite) degeneracy labeled by  $i$ ,

$$
L _ {\alpha i}. \tag {6.10}
$$

In particular, this implies that the number of the generators for the electromagnetic line defects is an integer multiple of the number of non-identity primaries in the chiral algebra. We will test this claim explicitly in the  $A_{2}, A_{3}, A_{4}$  Argyres-Douglas theories.

In both the Argyres-Douglas theories and in the  $SU(2)$  superconformal QCD, we observe that the line defect indices  $\mathcal{I}_{L_{\alpha i}}$  are related to the characters  $\chi_{\alpha}(q,z)$  of the chiral algebra by

$$
\mathcal {I} _ {L _ {\alpha i}} (q, z) = \sum_ {\beta \in \text {m o d u l e s}} v _ {\alpha i} ^ {\beta} (q, z) \chi_ {\beta} (q, z), \tag {6.11}
$$

where  $z$  denotes collectively all the flavor fugacities and  $v_{\alpha i}^{\beta}(q,z)$  are some polynomials in  $q$  and  $z$ . The sum in the modules above is finite for the Argyres-Douglas theory, but infinite

for  $SU(2)$  superconformal QCD. Similarly we observe that the line defect index with the insertion of two half lines  $L_{\alpha i}, L_{\beta i}$  (in the phase order, say,  $\vartheta_{\alpha i} < \vartheta_{\beta j}$ ) can be decomposed into characters of the chiral algebra,

$$
\mathcal {I} _ {L _ {\alpha i} L _ {\beta j}} (q, z) = \sum_ {\gamma \in \text {m o d u l e s}} v _ {\alpha i, \beta j} ^ {\gamma} (q, z) \chi_ {\gamma} (q, z), \tag {6.12}
$$

with some polynomials  $v_{\alpha i,\beta j}^{\gamma}(q,z)$ . Note that  $v_{\alpha i,0}^{\gamma}(q,z) = v_{0,\alpha i}^{\gamma}(q,z) = v_{\alpha i}^{\gamma}(q,z)$  by definition.

In all the Argyres-Douglas theories we have investigated, we find that when setting  $q = z = 1$ , the coefficients  $v_{\alpha i}^{\gamma}(1,1)$  and  $v_{\alpha i,\beta j}^{\gamma}(1,1)$  are independent of the degeneracy index  $i$ , while they still depend on the index  $\alpha$ , which labels the primaries in the  $2d$  chiral algebra. We can therefore define

$$
V _ {\alpha} ^ {\gamma} \equiv v _ {\alpha i} ^ {\gamma} (q = 1, z = 1), \tag {6.13}
$$

$$
V _ {\alpha \beta} ^ {\gamma} \equiv v _ {\alpha i, \beta j} ^ {\gamma} (q = 1, z = 1).
$$

Note that  $V_{\alpha 0}^{\gamma} = V_{0\alpha}^{\delta} = V_{\alpha}^{\gamma}$  by definition. We also observe that  $V_{\alpha \beta}^{\gamma}$  is symmetric in  $\alpha$  and  $\beta$ , whereas  $v_{\alpha i,\beta j}^{\gamma}(q,z)$  is not in general. In the following subsections we will determine the coefficients  $V_{\alpha}^{\beta}(q,z)$  and  $V_{\alpha \beta}^{\gamma}(q,z)$  for the  $A_2, A_3, A_4$  Argyres-Douglas theories.

Our main observation is that in the Argyres-Douglas theories, where the index  $\alpha$  runs over finitely many modules, the coefficients  $V_{\alpha}^{\beta}$  and  $V_{\alpha \beta}^{\gamma}$  obey the Verlinde algebra of the associated chiral algebra,

$$
\boxed {V _ {\alpha \beta} ^ {\gamma} = \sum_ {\alpha^ {\prime}, \beta^ {\prime} \in \text {m o d u l e s}} \mathcal {N} _ {\alpha^ {\prime} \beta^ {\prime}} ^ {\gamma} V _ {\alpha} ^ {\alpha^ {\prime}} V _ {\beta} ^ {\beta^ {\prime}},} \tag {6.14}
$$

where  $\mathcal{N}_{\alpha' \beta'}^{\gamma}$  are the fusion coefficients of the  $2d$  chiral algebra. Importantly, we will define the fusion coefficients  $\mathcal{N}_{\alpha \beta}^{\gamma}$  by the Verlinde formula [90]

$$
\mathcal {N} _ {\alpha \beta} ^ {\gamma} = \sum_ {\delta \in \mathrm {m o d u l e s}} \frac {\mathcal {S} _ {\alpha} ^ {\delta} \mathcal {S} _ {\beta} ^ {\delta} \bar {\mathcal {S}} ^ {\delta \gamma}}{\mathcal {S} _ {0} ^ {\delta}}, \tag {6.15}
$$

where  $\mathcal{S}$  is the modular transformation matrix and  $\bar{S}^{\delta \gamma} = S_{\beta}^{\delta}\mathcal{C}^{\beta \gamma}$  with  $\mathcal{C} = \mathcal{S}^2$ . This definition circumvents certain subtleties in the fusion rules in non-rational CFTs. For chiral algebra whose  $\mathcal{S}$  matrix is not well-defined (for example, the  $\widehat{so(8)}_{-2}$  in the  $SU(2)$  with  $N_f = 4$  flavors theory), we will not attempt to define the fusion coefficients and we do

not have a proposal for the relation between the  $4d$  line defect product and the  $2d$  fusion rule.

We can phrase the above relation in a more illuminating way. For each class of  $4d$  line defects labeled by  $\alpha$ ,  $\{L_{\alpha i}\}_{i}$ , we associate it to an element  $[L_{\alpha}]$  of the Verlinde algebra using the coefficients  $V_{\alpha}^{\beta}$ ,

$$
\left[ L _ {\alpha} \right] \equiv \sum_ {\beta \in \text {m o d u l e s}} V _ {\alpha} ^ {\beta} \left[ \Phi_ {\beta} \right]. \tag {6.16}
$$

Similarly we define  $[L_{\alpha}L_{\beta}]$  as

$$
\left[ L _ {\alpha} L _ {\beta} \right] \equiv \sum_ {\gamma \in \text {m o d u l e s}} V _ {\alpha \beta} ^ {\gamma} \left[ \Phi_ {\gamma} \right]. \tag {6.17}
$$

Now our main observation (6.14) can be written as

$$
\boxed {[ L _ {\alpha} L _ {\beta} ] = [ L _ {\alpha} ] \times [ L _ {\beta} ],} \tag {6.18}
$$

where the product  $\times$  is the fusion product in the  $2d$  Verlinde algebra. A similar proposal was made and verified in various examples in [51].

# 6.2.1  $A_{2}$  Argyres-Douglas Theory

The chiral algebra associated to the  $4d$ $A_{2}$  Argyres-Douglas theory is the  $(2,5)$  Virasoro minimal model [8,13,91]. The primaries of the  $(2,5)$  minimal model are the identity  $\Phi_0 = 1$  and a non-identity primary  $\Phi_1\equiv \Phi_{1,2}$  with weight  $-1 / 5$ , hence the index  $\alpha$  in the previous section runs over 0,1. On the other hand there are five non-unity generators  $L_{i}$  for the UV line defects with the same line defect indices, hence the degeneracy labeled by  $i$  in the previous section is five here,  $i = 1,\dots ,5$ .

Recall that the Schur index without any insertion of line defects equals to the vacuum character of the (2,5) Virasoro minimal model [8],

$$
\mathcal {I} (q) = \chi_ {0} (q). \tag {6.19}
$$

We find that the line defect index $^{19}$ $\mathcal{I}_L(q)$  is related to the character for the non-identity primary  $\Phi_1$  with weight  $h_{1,2} = -1/5$  in the following way

$$
\mathcal {I} _ {L} (q) = q ^ {- \frac {1}{2}} \chi_ {0} (q) - q ^ {- \frac {1}{2}} \chi_ {1} (q). \tag {6.20}
$$

We have normalized the characters such that they start from 1. For completeness, we record the two characters in the (2,5) minimal model character (see (6.41))

$$
\begin{array}{l} \chi_ {0} (q) = 1 + q ^ {2} + q ^ {3} + q ^ {4} + q ^ {5} + 2 q ^ {6} + 2 q ^ {7} + 3 q ^ {8} + 3 q ^ {9} + 4 q ^ {1 0} + 4 q ^ {1 1} + \dots , \tag {6.21} \\ \chi_ {1} (q) = 1 + q + q ^ {2} + q ^ {3} + 2 q ^ {4} + 2 q ^ {5} + 3 q ^ {6} + 3 q ^ {7} + 4 q ^ {8} + 5 q ^ {9} + 6 q ^ {1 0} + \dots . \\ \end{array}
$$

A similar observation was made in [51] (see, in particular, (9.51)) in the case of the inverse of the quantum KS operator, but the precise linear combination of line defects are different. Incidentally, the characters for the vacuum  $\Phi_0$  and the non-identity primary  $\Phi_1$  in the (2, 5) minimal model are known to be the two Rogers-Ramanujan functions  $H(q)$  and  $G(q)$ , respectively.

The coefficients  $V_{\beta}^{\alpha}$  can be read off to be

$$
V _ {0} ^ {\alpha} = (1, 0), \quad V _ {1} ^ {\alpha} = (1, - 1). \tag {6.22}
$$

The trace of two line defects are given in (5.10). For example,

$$
\mathcal {I} _ {L _ {i} L _ {i + 2}} = \mathcal {I} + q ^ {\frac {1}{2}} \mathcal {I} _ {L} = 2 \chi_ {0} (q) - \chi_ {1} (q). \tag {6.23}
$$

One can easily verify that the coefficients  $V_{\alpha \beta}^{\gamma}$  are indeed independent of the degeneracy index  $i$  and are given by

$$
V _ {1 1} ^ {\alpha} = (2, - 1), \tag {6.24}
$$

together with  $V_{\alpha 0}^{\beta}$  given by  $V_{\alpha}^{\beta}$ .

On the other hand, the Verlinde algebra for the  $(2,5)$  minimal model is

$$
\left[ \Phi_ {1} \right] \times \left[ \Phi_ {1} \right] = [ 1 ] + \left[ \Phi_ {1} \right]. \tag {6.25}
$$

One can easily check that (6.14) is satisfied. Indeed, from the coefficients  $V_{\alpha}^{\beta}$  and  $V_{\alpha \beta}^{\gamma}$  we have

$$
[ L ] = [ 1 ] - \left[ \Phi_ {1} \right], \tag {6.26}
$$

and

$$
[ L L ] = 2 [ 1 ] - \left[ \Phi_ {1} \right], \tag {6.27}
$$

and the equivalent statement (6.18) is satisfied

$$
[ L L ] = [ L ] \times [ L ]. \tag {6.28}
$$

# 6.2.2  $A_{3}$  Argyres-Douglas Theory

The chiral algebra associated to the  $A_{3}$  Argyres-Douglas theory is the affine Lie algebra  $\widehat{su(2)}_{-\frac{4}{3}}$  [7,8,13,91]. The weights of  $\widehat{su(2)}_k$  are labeled by its Dynkin labels  $[\lambda_0,\lambda_1]$  with  $\lambda_0 + \lambda_1 = k$ . There are three representations in  $\widehat{su(2)}_{-\frac{4}{3}}$ , whose highest weights are

$$
\Phi_ {0} = \left[ - \frac {4}{3}, 0 \right], \quad \Phi_ {1} = \left[ - \frac {2}{3}, - \frac {2}{3} \right], \quad \Phi_ {2} = \left[ 0, - \frac {4}{3} \right]. \tag {6.29}
$$

that are known to be admissible [92] (see also [93]). Admissible representations have the nice property that their characters transform linearly into each other under modular transformation, so the  $S$  matrix of modular transformation is well-defined. Note that the first highest weight above  $\left[-\frac{4}{3},0\right]$  is that for the vacuum module, whose Dynkin label of the finite  $SU(2)$  Lie algebra is zero,  $\lambda_{1} = 0$ . The latter two representations are conjugate to each other. The index  $\alpha$  in (6.10) runs over  $0,1,2$ , which labels the primaries. As we saw in Section 5.3, there are six non-flavor line defects  $L_{1i}\equiv A_i$  and  $L_{2i}\equiv B_i$ , with  $i = 1,2,3$ .

The characters for the three admissible representations can be computed using the Kazhdan-Lusztig formula as reviewed in Appendix C (they can also be found in Chapter 18 of [93]),

$$
\chi_ {0} (q, z) = \frac {\sum_ {m = 0} ^ {\infty} (- 1) ^ {m} \frac {z ^ {2 m + 1} - z ^ {- (2 m + 1)}}{z - z ^ {- 1}} q ^ {\frac {3 m (m + 1)}{2}}}{\prod_ {n = 1} ^ {\infty} (1 - q ^ {n}) (1 - z ^ {2} q ^ {n}) (1 - z ^ {- 2} q ^ {n})},
$$

$$
\chi_ {1} (q, z) = \frac {1 + \sum_ {n = 1} ^ {\infty} (- 1) ^ {n} \left(z ^ {- 2 n} q ^ {\frac {n}{2} (3 n - 1)} + z ^ {2 n} q ^ {\frac {n}{2} (3 n + 1)}\right)}{\left(1 - z ^ {- 2}\right) \prod_ {n = 1} ^ {\infty} \left(1 - q ^ {n}\right) \left(1 - z ^ {2} q ^ {n}\right) \left(1 - z ^ {- 2} q ^ {n}\right)}, \tag {6.30}
$$

$$
\chi_ {2} (q, z) = \frac {1 + \sum_ {n = 1} ^ {\infty} (- 1) ^ {n} \left(z ^ {2 n} q ^ {\frac {n}{2} (3 n - 1)} + z ^ {- 2 n} q ^ {\frac {n}{2} (3 n + 1)}\right)}{(1 - z ^ {- 2}) \prod_ {n = 1} ^ {\infty} (1 - q ^ {n}) (1 - z ^ {2} q ^ {n}) (1 - z ^ {- 2} q ^ {n})}.
$$

Note that for the latter two modules, there are infinitely many states at each grade created by the zero modes  $J_0^A$  of the Kac-Moody algebra, due to the fact that finite  $SU(2)$  Dynkin labels  $\lambda_{1}$  are negative fractional. Hence the two characters diverge as  $1 / (1 - z^{-2})$  as  $z\to 1$ . The vacuum character  $\chi_0(q,z)$  has been computed previously in [7,8].

We find that the line defect indices are related to the characters of  $\widehat{su(2)}_{-\frac{4}{3}}$  as follows,

$$
\mathcal {I} (q, z) = \chi_ {0} (q, z),
$$

$$
\mathcal {I} _ {A} (q, z) = q ^ {- \frac {1}{2}} z ^ {- 1} \left[ - \chi_ {1} (q, z) + \chi_ {2} (q, z) \right], \tag {6.31}
$$

$$
\mathcal {I} _ {B} (q, z) = q ^ {- \frac {1}{2}} \left[ \chi_ {0} (q, z) - \chi_ {1} (q, z) + z ^ {- 2} \chi_ {2} (q, z) \right],
$$

where  $\chi_{\alpha}(q,z)$  is the character for the primary  $\Phi_{\alpha}$ . The first line is the Schur index without any insertion of line defects  $\mathcal{I}(q,z)$ , which equals to the vacuum character  $\chi_0(q,z)$  of  $\widehat{su(2)}_{-\frac{4}{3}}$

[7,8,13,91]. Hence the coefficients  $V_{\alpha}^{\beta}$  are

$$
V _ {0} ^ {\alpha} = (1, 0, 0), \qquad V _ {1} ^ {\alpha} = (0, - 1, 1), \qquad V _ {2} ^ {\alpha} = (1, - 1, 1). \qquad \qquad (6. 3 2)
$$

The Schur indices for two (half) line defects are given in (5.18). We have (no sum in the indices)

$$
\mathcal {I} _ {A _ {i} A _ {i}} (q, z) = (1 + q ^ {- 1}) \chi_ {0} (q, z) - q ^ {- 1} \chi_ {1} (q, z) + q ^ {- 1} z ^ {- 2} \chi_ {2} (q, z),
$$

$$
\begin{array}{l} \mathcal {I} _ {B _ {i} B _ {i}} (q, z) = \left(1 + q ^ {- 1} + q ^ {- 2}\right) \chi_ {0} (q, z) - \left[ q ^ {- 1} \left(1 + z ^ {- 2}\right) + q ^ {- 2} \right] \chi_ {1} (q, z) \tag {6.33} \\ + \left[ q ^ {- 1} (1 + z ^ {- 2}) + q ^ {- 2} z ^ {- 2} \right] \chi_ {2} (q, z), \\ \end{array}
$$

$$
\mathcal {I} _ {A _ {i} B _ {i}} (q, z) = (z + z ^ {- 1}) \chi_ {0} (q, z) - (1 + q ^ {- 1}) z ^ {- 1} \chi_ {1} (q, z) + (1 + q ^ {- 1}) z ^ {- 1} \chi_ {2} (q, z).
$$

Hence the coefficients  $V_{\alpha \beta}^{\gamma}$  are

$$
V _ {1 1} ^ {\alpha} = (2, - 1, 1),
$$

$$
V _ {2 2} ^ {\alpha} = (3, - 3, 3), \tag {6.34}
$$

$$
V _ {1 2} ^ {\alpha} = V _ {2 1} ^ {\alpha} = (2, - 2, 2),
$$

together with  $V_{\alpha 0}^{\beta} = V_{\alpha}^{\beta}$

On the other hand, the  $S$  matrix for these three admissible representations is [93]

$$
\mathcal {S} _ {\beta} ^ {\alpha} = - \frac {1}{\sqrt {3}} \left( \begin{array}{c c c} 1 & - 1 & 1 \\ - 1 & e ^ {\frac {4 \pi i}{3}} & - e ^ {\frac {2 \pi i}{3}} \\ 1 & - e ^ {\frac {2 \pi i}{3}} & e ^ {\frac {4 \pi i}{3}} \end{array} \right). \tag {6.35}
$$

The conjugation matrix  $\mathcal{C} = \mathcal{S}^2$  is given by

$$
\mathcal {C} _ {\beta} ^ {\alpha} = \left( \begin{array}{c c c} 1 & 0 & 0 \\ 0 & 0 & - 1 \\ 0 & - 1 & 0 \end{array} \right). \tag {6.36}
$$

Note that  $\Phi_1 = \left[-\frac{2}{3}, -\frac{2}{3}\right]$  and  $\Phi_2 = \left[0, -\frac{4}{3}\right]$  are conjugate to each other. The fusion rules obtained from the Verlinde formula (6.15) are

$$
\left[ \Phi_ {1} \right] \times \left[ \Phi_ {1} \right] = \left[ \Phi_ {2} \right]
$$

$$
\left[ \Phi_ {2} \right] \times \left[ \Phi_ {2} \right] = - \left[ \Phi_ {1} \right] \tag {6.37}
$$

$$
\left[ \Phi_ {1} \right] \times \left[ \Phi_ {2} \right] = - \left[ \Phi_ {0} \right].
$$

Note that the minus sign in the fusion rule signals the negative central charge of the affine Lie algebra  $\widehat{su(2)}_{-\frac{4}{3}}$ .

From the coefficients  $V_{\alpha}^{\beta}$  and  $V_{\alpha \beta}^{\gamma}$ , we have

$$
[ A ] = - [ \Phi_ {1} ] + [ \Phi_ {2} ] \tag {6.38}
$$

$$
[ B ] = [ \Phi_ {0} ] - [ \Phi_ {1} ] + [ \Phi_ {2} ],
$$

and

$$
[ A A ] = 2 [ \Phi_ {0} ] - [ \Phi_ {1} ] + [ \Phi_ {2} ],
$$

$$
[ B B ] = 3 [ \Phi_ {0} ] - 3 [ \Phi_ {1} ] + 3 [ \Phi_ {2} ], \tag {6.39}
$$

$$
[ A B ] = 2 [ \Phi_ {0} ] - 2 [ \Phi_ {1} ] + 2 [ \Phi_ {2} ].
$$

One can check straightforwardly that our proposal (6.18) is satisfied, i.e.  $[AA] = [A]\times [A]$ ,  $[BB] = [B]\times [B]$ , and  $[AB] = [A]\times [B]$  using the fusion rules (6.37).

# 6.2.3  $A_{4}$  Argyres-Douglas Theory

The chiral algebra associated to the  $4d$ $A_4$  Argyres-Douglas theory is the  $(2,7)$  Virasoro minimal model [8,13,91]. The primaries of the  $(2,7)$  minimal model are the identity  $\Phi_{1,1} = 1$  and two non-identity primaries  $\Phi_{1,2}$  and  $\Phi_{1,3}$  with weight  $-2/7$  and  $-3/7$ , respectively. On the other hand there are 14 non-unity generators for the UV line defects grouped into  $A_i$  and  $B_i$  with  $i = 1,\dots,7$ . The line defect indices are related to the characters of the  $(2,7)$  Virasoro minimal model by

$$
\mathcal {I} (q) = \chi_ {(1, 1)} (q),
$$

$$
\mathcal {I} _ {A} (q) = - q ^ {- 1} \chi_ {(1, 2)} (q) + q ^ {- 1} \chi_ {(1, 3)} (q), \tag {6.40}
$$

$$
\mathcal {I} _ {B} (q) = q ^ {- \frac {1}{2}} \chi_ {(1, 1)} (q) - q ^ {- \frac {1}{2}} \chi_ {(1, 2)} (q).
$$

The characters of the  $\Phi_{s,r}$  primary with  $1\leq s\leq p - 1$  and  $1\leq r\leq p^{\prime} - 1$  in the  $(p,p^{\prime})$  Virasoro minimal model is given by (see, for example, [93])

$$
\chi_ {(s, r)} (q) = q ^ {- \frac {\left(r p - s p ^ {\prime}\right) ^ {2} - \left(p - p ^ {\prime}\right) ^ {2}}{4 p p ^ {\prime}} + \frac {1}{2 4} \left(1 - \frac {6 \left(p - p ^ {\prime}\right) ^ {2}}{p p ^ {\prime}}\right)} \left(K _ {s, r} ^ {(p, p ^ {\prime})} (q) - K _ {- s, r} ^ {(p, p ^ {\prime})} (q)\right) \tag {6.41}
$$

where

$$
K _ {s, r} ^ {(p, p ^ {\prime})} (q) = \frac {q ^ {- \frac {1}{2 4}}}{(q) _ {\infty}} \sum_ {n \in \mathbb {Z}} q ^ {\frac {\left(2 p p ^ {\prime} n + p r - p ^ {\prime} s\right) ^ {2}}{4 p p ^ {\prime}}}. \tag {6.42}
$$

Again we have normalized the character to start from 1. Note that we have the following identification between primaries,  $\Phi_{s,r} = \Phi_{p - s,p' - r}$ .

For the purpose of demonstrating the Verlinde algebra, we only need the following Schur

indices of two line defects (no sum in the indices)

$$
\mathcal {I} _ {A _ {i} A _ {i}} (q) = q ^ {- 3} (1 + q) \chi_ {(1, 1)} (q) - q ^ {- 3} \chi_ {(1, 2)} (q),
$$

$$
\mathcal {I} _ {B _ {i} B _ {i}} (q) = q ^ {- 2} (1 + q) \chi_ {(1, 1)} (q) - q ^ {- 2} (1 + q) \chi_ {(1, 2)} (q) + q ^ {- 1} \chi_ {(1, 3)} (q), \tag {6.43}
$$

$$
\mathcal {I} _ {A _ {5} B _ {6}} (q) = q ^ {- 1 / 2} \chi_ {(1, 1)} (q) - 2 q ^ {- 1 / 2} \chi_ {(1, 2)} (q) + q ^ {- 1 / 2} \chi_ {(1, 3)} (q).
$$

From the above relations between the line defect indices and the characters, we define  $\Phi_{1,3}$

$$
\left[ \begin{array}{l} A \end{array} \right] = - \left[ \Phi_ {1, 2} \right] + \left[ \Phi_ {1, 3} \right], \tag {6.44}
$$

$$
[ B ] = [ \Phi_ {1, 1} ] - [ \Phi_ {1, 2} ].
$$

and

$$
[ A A ] = 2 [ \Phi_ {1, 1} ] - [ \Phi_ {1, 2} ],
$$

$$
[ B B ] = 2 \left[ \Phi_ {1, 1} \right] - 2 \left[ \Phi_ {1, 2} \right] + \left[ \Phi_ {1, 3} \right], \tag {6.45}
$$

$$
[ A B ] = [ \Phi_ {1, 1} ] - 2 [ \Phi_ {1, 2} ] + [ \Phi_ {1, 3} ].
$$

The other Schur indices of two line defects (e.g.  $\mathcal{I}_{A_1B_2}$ ) can be similarly shown to give the same definitions for [AA], [BB], [AB].

The fusion rule in the (2,7) Virasoro minimal model,

$$
\left[ \Phi_ {1, 2} \right] \times \left[ \Phi_ {1, 2} \right] = \left[ \Phi_ {1, 1} \right] + \left[ \Phi_ {1, 3} \right],
$$

$$
\left[ \Phi_ {1, 3} \right] \times \left[ \Phi_ {1, 3} \right] = \left[ \Phi_ {1, 1} \right] + \left[ \Phi_ {1, 2} \right] + \left[ \Phi_ {1, 3} \right], \tag {6.46}
$$

$$
\left[ \Phi_ {1, 2} \right] \times \left[ \Phi_ {1, 3} \right] = \left[ \Phi_ {1, 2} \right] + \left[ \Phi_ {1, 3} \right].
$$

We have omitted the trivial fusion rules between the identity  $\Phi_{1,1} = 1$  with others. It is straightforward to check that

$$
[ A A ] = [ A ] \times [ A ], \quad [ B B ] = [ B ] \times [ B ], \quad [ A B ] = [ A ] \times [ B ], \tag {6.47}
$$

where  $\times$  is the fusion product in the (2,7) Virasoro minimal model given in (6.46). To conclude, we find that the indices of products of  $4d$  line defects are reproduced by the Verlinde algebra of the (2,7) Virasoro minimal model.

# Acknowledgements

We thank Tomoyuki Arakawa, Chris Beem, Chih-Kai Chang, Heng-Yu Chen, Thomas Dumitrescu, Sarah Harrison, Leonardo Rastelli, Cumrun Vafa, Herman Verlinde, Masahito

Yamazaki for interesting discussions. CC is supported by a Schmidt fellowship at the Institute for Advanced Study and DOE grant de-sc0009988. The research of DG was supported by the Perimeter Institute for Theoretical Physics. Research at Perimeter Institute is supported by the Government of Canada through Industry Canada and by the Province of Ontario through the Ministry of Economic Development & Innovation. SHS would like to thank National Taiwan University, University of Amsterdam, Tata Institute of Fundamental Research, Perimeter Institute for Theoretical Physics for their hospitality during various stages of this work. We thank the 2015 Simons workshop in Mathematics and Physics and the Simons Center for Geometry and Physics for hospitality.

# A Supercharges of Line Defects and Chiral Algebras

In this appendix we will work out the supercharges shared by the line defects and the chiral algebra. We will see that both the full line defects and the half line defects share the same set of supercharges (A.17) with the chiral algebra plane.

We will follow the convention in [13] for the  $4d\mathcal{N} = 2$  superconformal algebra.  $A, B, \dots = 1, 2$  will denote the doublet index of  $SU(2)_R$ .  $\alpha, \beta, \dots = +, -$  and  $\dot{\alpha}, \dot{\beta}, \dots = \dot{+}, \dot{-}$  will denote the doublet indices of  $SU(2)_1 \times SU(2)_2 = SO(4)_{\mathrm{rotation}}$ . All the doublet indices will be raised and lowered by  $\epsilon^{12} = \epsilon_{21} = +1$ . The nonzero anticommutators between the sixteen fermionic generators  $\{Q_{\alpha}^{A}, \tilde{Q}_{A\dot{\alpha}}, S_{A}^{\alpha}, \tilde{S}^{A\dot{\alpha}}\}$  in the  $4d\mathcal{N} = 2$  superconformal algebra are

$$
\{Q _ {\alpha} ^ {A}, \tilde {Q} _ {B \dot {\beta}} \} = 2 \delta_ {B} ^ {A} \sigma_ {\alpha \dot {\beta}} ^ {\mu} P _ {\mu} = \delta_ {B} ^ {A} P _ {\alpha \dot {\beta}},
$$

$$
\{\tilde {S} ^ {A \dot {\alpha}}, S _ {B} ^ {\beta} \} = 2 \delta_ {B} ^ {A} \bar {\sigma} ^ {\mu \dot {\alpha} \beta} K _ {\mu} = \delta_ {B} ^ {A} K ^ {\dot {\alpha} \beta},
$$

$$
\left\{Q _ {\alpha} ^ {A}, S _ {B} ^ {\beta} \right\} = \frac {1}{2} \delta_ {B} ^ {A} \delta_ {\alpha} ^ {\beta} D + \delta_ {B} ^ {A} M _ {\alpha} ^ {\beta} - \delta_ {\alpha} ^ {\beta} R _ {B} ^ {A}, \tag {A.1}
$$

$$
\{\tilde {S} ^ {A \dot {\alpha}}, \tilde {Q} _ {B \dot {\beta}} \} = \frac {1}{2} \delta_ {B} ^ {A} \delta_ {\dot {\beta}} ^ {\dot {\alpha}} D + \delta_ {B} ^ {A} M _ {\dot {\beta}} ^ {\dot {\alpha}} + \delta_ {\dot {\beta}} ^ {\dot {\alpha}} R _ {B} ^ {A}.
$$

The  $SU(2)_R$  generators  $R^{\pm}, R$  and the  $U(1)_r$  generator  $r$  sit inside  $R_B^A$  as

$$
R _ {2} ^ {1} = R ^ {+}, \quad R _ {1} ^ {2} = R ^ {-}, \quad R _ {1} ^ {1} = \frac {1}{2} r + R, \quad R _ {2} ^ {2} = \frac {1}{2} r - R, \tag {A.2}
$$

where  $[R^{+}, R^{-}] = 2R$  and  $[R, R^{\pm}] = \pm R^{\pm}$ .

# A.1 Supercharges Preserved by Full Lines

The eight supercharges preserved by an infinitely extended line defect (a full line) pointing in the direction  $n^{\mu}$  in  $\mathbb{R}^4$  are [27]

$$
G _ {\alpha} ^ {A} \equiv \xi^ {- 1} Q _ {\alpha} ^ {A} + \xi n _ {\mu} \sigma_ {\alpha \dot {\alpha}} ^ {\mu} \tilde {Q} ^ {A \dot {\alpha}}, \tag {A.3}
$$

$$
H _ {A} ^ {\alpha} \equiv \xi S _ {A} ^ {\alpha} - \xi^ {- 1} n _ {\mu} \bar {\sigma} ^ {\mu \dot {\alpha} \alpha} \tilde {S} _ {A \dot {\alpha}},
$$

where  $\xi$  is a phase related to the  $u$  and  $\vartheta$  in the previous sections by  $u = \xi^{-2} = e^{i\vartheta}$ . Here  $\bar{\sigma}^{\mu \dot{\alpha}\alpha} = \epsilon^{\dot{\alpha}\dot{\beta}}\epsilon^{\alpha \beta}\sigma_{\dot{\beta}\beta}^{\mu}$ .

To make sure the above linear combinations are the correct supercharges preserved by the line defect, one has to check, for example, their anticommutators do not contain the  $U(1)_r$  generator  $r$ , the translations, nor the special conformal transformations along other directions than  $n_\mu$ . Let us check this for a few anticommutators.[20]

$$
\left\{G _ {\alpha} ^ {A}, G _ {\beta} ^ {B} \right\} = 4 \epsilon^ {A B} n _ {\mu} P _ {\nu} \epsilon^ {\dot {\alpha} \dot {\beta}} \sigma_ {[ \alpha | \dot {\alpha}} ^ {\mu} \sigma_ {\beta ] \dot {\beta}} ^ {\nu} = 4 \epsilon^ {A B} \epsilon_ {\alpha \beta} n ^ {\mu} P _ {\mu}, \tag {A.4}
$$

$$
\{H _ {A} ^ {\alpha}, H _ {B} ^ {\beta} \} = 4 \epsilon_ {A B} \epsilon_ {\dot {\alpha} \dot {\gamma}} \bar {\sigma} ^ {\mu \dot {\alpha} [ \beta} \bar {\sigma} ^ {\nu \dot {\gamma} | \alpha ]} n _ {\mu} K _ {\nu} = - 4 \epsilon_ {A B} \epsilon^ {\alpha \beta} n ^ {\mu} K _ {\mu},
$$

where we have used  $\epsilon^{\dot{\alpha}\dot{\beta}}\sigma_{[\alpha |\dot{\alpha}}^{\mu}\sigma_{\beta ]\dot{\beta}}^{\nu} = \delta^{\mu \nu}\epsilon_{\alpha \beta}$  and  $\epsilon_{\dot{\alpha};\bar{\gamma}}\bar{\sigma}^{\mu \dot{\alpha} [\beta}\bar{\sigma}^{\nu \dot{\gamma} | \alpha ]} = \delta^{\mu \nu}\epsilon_{\beta \alpha}$ . Also, let us check that there is no  $r$  in the anticommutator between  $G_{\alpha}^{A}$  and  $H_{A}^{\alpha}$ .

$$
\begin{array}{l} \left\{G _ {\alpha} ^ {A}, H _ {B} ^ {\beta} \right\} = \left\{Q _ {\alpha} ^ {A}, S _ {B} ^ {\beta} \right\} - \left(n _ {\mu} \sigma_ {\alpha \dot {\alpha}} ^ {\mu}\right) \left(n _ {\nu} \bar {\sigma} ^ {\nu \dot {\beta} \beta}\right) \epsilon^ {A C} \epsilon^ {\dot {\alpha} \dot {\gamma}} \epsilon_ {B D} \epsilon_ {\dot {\beta} \dot {\delta}} \left\{\tilde {Q} _ {C}; \tilde {S} ^ {D \dot {\delta}} \right\} \tag {A.5} \\ = - \delta_ {\alpha} ^ {\beta} \left(R _ {B} ^ {A} + \epsilon^ {A C} \epsilon_ {B D} R _ {C} ^ {D}\right) + \text {r o t a t i o n s a n d d i l a t i o n} \\ \end{array}
$$

where we have used  $\sigma_{\alpha \dot{\alpha}}^{(\mu}\bar{\sigma}^{\nu)\dot{\alpha}\beta} = -\delta^{\mu \nu}\delta_{\alpha}^{\beta}$ . Indeed, the righthand side does not contain the  $U(1)_r$  generator.

Finally, we would like to determine the supercharges preserved by the Schur operators in the presence of a full line defect. To do so, it is convenient to fix our conventions for the Pauli matrices to be

$$
\sigma_ {\alpha \dot {\beta}} ^ {1} = \left( \begin{array}{c c} 0 & 1 \\ 1 & 0 \end{array} \right), \quad \sigma_ {\alpha \dot {\beta}} ^ {2} = \left( \begin{array}{c c} 0 & - i \\ i & 0 \end{array} \right), \quad \sigma_ {\alpha \dot {\beta}} ^ {3} = \left( \begin{array}{c c} 1 & 0 \\ 0 & - 1 \end{array} \right), \quad \sigma_ {\alpha \dot {\beta}} ^ {4} = \left( \begin{array}{c c} i & 0 \\ 0 & i \end{array} \right). \tag {A.6}
$$

We have  $\bar{\sigma}^{\mu} = (-\sigma^{1}, - \sigma^{2}, - \sigma^{3},\sigma^{4})$  numerically. Further, we will choose our line defect to be along the 1-direction, i.e.,  $n_\mu = (1,0,0,0)$ , and  $\xi = 0$ .

The ordinary Schur index without insertions of line defects receives contributions only from operators that are annihilated by four supercharges, which can be chosen to be  $Q_{-}^{1}$ ,

$\tilde{Q}_{2\dot{-}}$ ,  $S_{1}^{-}$ ,  $\tilde{S}^{2 - }$  [2,3]. In the presence of a full line defect, two out of the four supercharges above are shared by the eight supercharges (A.3) that are preserved by a full line defect,[21]

$$
G _ {-} ^ {1} = Q _ {-} ^ {1} + \tilde {Q} _ {2 \cdot},
$$

$$
H _ {1} ^ {-} = S _ {1} ^ {-} + \tilde {S} ^ {2 \dot {\mathrm {一}}}. \tag {A.7}
$$

Their anticommutators are

$$
\left\{G _ {-} ^ {1}, G _ {-} ^ {1} \right\} = 0, \quad \left\{H _ {1} ^ {-}, H _ {1} ^ {-} \right\} = 0,
$$

$$
\left\{G _ {-} ^ {1}, H _ {1} ^ {-} \right\} = E - \left(j _ {1} + j _ {2}\right) - 2 R. \tag {A.8}
$$

Here  $E$ ,  $j_{1}$ ,  $j_{2}$  are the eigenvalue of the dilation charge  $D$ ,  $M_{+}^{+}$ ,  $M_{+}^{+}$ , respectively. Note that  $\mathcal{M} \equiv j_{1} + j_{2}$  is the rotation of the plane orthogonal to the line defect.

In summary, the Schur operators in the presence of a line defect are annihilated by the two supercharges (A.7), and they obey the following condition,

$$
\text {L i n e} \hat {L} _ {0} \equiv \frac {1}{2} (E - (j _ {1} + j _ {2})) - R = 0. \tag {A.9}
$$

Recall that the ordinary Schur operators without the line defect obey an additional condition  $\mathcal{Z} \equiv r + (j_1 - j_2) = 0$ .

More generally, we consider multiple full line defects  $L_{i}$  with phases

$$
\xi_ {i} = e ^ {- i \vartheta_ {i} / 2} \tag {A.10}
$$

and pointing in the directions  $(n_i)^\mu$ . To preserve supersymmetry, we have to put them on the 12-plane (the plane orthogonal to the chiral algebra plane) with orientations

$$
\left(n _ {i}\right) ^ {\mu} = \left(\cos \vartheta_ {i}, \sin \vartheta_ {i}, 0, 0\right). \tag {A.11}
$$

Note that  $(n_i)_{\mu}\sigma_{\alpha \dot{\alpha}}^{\mu} = \left( \begin{array}{cc}0 & e^{-i\vartheta_i}\\ e^{i\vartheta_i} & 0 \end{array} \right)$  and  $(n_i)_{\mu}\bar{\sigma}^{\mu \dot{\alpha}\alpha} = -\left( \begin{array}{cc}0 & e^{-i\vartheta_i}\\ e^{i\vartheta_i} & 0 \end{array} \right).$

The 4 supercharges shared by all the (full) line defects are

$$
G _ {-} ^ {A} = Q _ {-} ^ {A} + \tilde {Q} ^ {A \dot {+}},
$$

$$
H _ {A} ^ {-} = S _ {A} ^ {-} + \tilde {S} _ {A \dot {+}}. \tag {A.12}
$$

which in particular include two supercharges  $G_{-}^{1}$  and  $H_{1}^{-}$  (A.7) that are used to define the line defect Schur index.

# A.2 Supercharges Preserved by Half Lines

In this subsection we will check that half line defects preserve the same two supercharges (A.7) that are used to define the line defect Schur index. To simplify the notations, we define

$$
G _ {\alpha} \equiv G _ {\alpha} ^ {1}, \qquad H ^ {\alpha} \equiv H _ {1} ^ {\alpha}, \qquad \Delta \equiv E - 2 R. \tag {A.13}
$$

They satisfy the following hermicity conditions  $(G_{\alpha})^{\dagger} = H^{\alpha}$ ,  $\Delta^{\dagger} = \Delta$ .

The superalgebra preserved by a half line pointing along the 1-direction is

$$
\begin{array}{l} \{G _ {\alpha}, G _ {\beta} \} = \{H ^ {\alpha}, H ^ {\beta} \} = 0, \qquad \{G _ {\alpha}, H ^ {\beta} \} = M _ {\alpha} ^ {\beta} + \delta_ {\alpha} ^ {\beta} \Delta , \\ [ \Delta , G _ {\alpha} ] = - \frac {1}{2} G _ {\alpha}, \quad [ \Delta , H ^ {\alpha} ] = \frac {1}{2} H ^ {\alpha}, \tag {A.14} \\ [ M _ {\alpha} ^ {\beta}, G _ {\gamma} ] = \delta_ {\gamma} ^ {\beta} G _ {\alpha} - \frac {1}{2} \delta_ {\alpha} ^ {\beta} G _ {\gamma}, [ M _ {\alpha} ^ {\beta}, H ^ {\gamma} ] = - \delta_ {\alpha} ^ {\gamma} H ^ {\beta} + \frac {1}{2} \delta_ {\alpha} ^ {\beta} H ^ {\gamma}, \\ [ M _ {\alpha} ^ {\beta}, M _ {\gamma} ^ {\delta} ] = \delta_ {\gamma} ^ {\beta} M _ {\alpha} ^ {\delta} - \delta_ {\alpha} ^ {\delta} M _ {\gamma} ^ {\beta}, [ M _ {\alpha} ^ {\beta}, \Delta ] = 0. \\ \end{array}
$$

Here  $M_{\alpha}^{\beta}$  are the generators of the  $SO(3)$  rotating the  $\mathbb{R}^3$  transverse to the half line. Note that in contrast to the case of a full line defect, the translation  $P^1$  and special conformal transformation  $K^1$  are no longer symmetries of the configuration, hence do not show up in the above algebra. The preserved four supercharges are  $G_{\alpha}^{1}$  and  $H_{1}^{\alpha}$ , which include (A.7).

More generally, if we include multiple half lines with phases  $\xi_{i} = e^{-i\vartheta_{i} / 2}$  and orientations given as in (A.11), the only preserved rotation symmetry is the one rotating the 34-plane,  $M_{-}^{-}$  (whose eigenvalue is  $j_{1} + j_{2}$  in the notation before). The preserved superalgebra is

$$
\begin{array}{l} \{G _ {-}, G _ {-} \} = \{H ^ {-}, H ^ {-} \} = 0, \qquad \{G _ {-}, H ^ {-} \} = \Delta + M _ {-} ^ {-}, \\ [ \Delta , G _ {-} ] = - \frac {1}{2} G _ {-}, \quad [ \Delta , H ^ {-} ] = \frac {1}{2} H ^ {-}, \tag {A.15} \\ [ M _ {-} ^ {-}, G _ {-} ] = \frac {1}{2} G _ {-}, [ M _ {-} ^ {-}, H ^ {-} ] = - \frac {1}{2} H ^ {-}, [ M _ {-} ^ {-}, \Delta ] = 0, \\ \end{array}
$$

The preserved supercharges  $G_{-}^{1}$  and  $H_{1}^{-}$  (A.7) are precisely those used to define the line defect Schur index. Incidentally,  $2\hat{L}_0 = \Delta + M_{-}^{-}$  where  $\hat{L}_0 = \frac{1}{2} (E - 2R - (j_1 + j_2)) = 0$  is satisfied by all the line defect Schur operators (A.9).

# A.3 Supercharges Shared by the Chiral Algebra and Line Defects

In this subsection we will determine the supercharges shared by the chiral algebra plane and the (full or half) line defects. The chiral algebra operators live in the cohomology of

the following four supercharges [13]

$$
\mathbb {Q} _ {1} \equiv Q _ {-} ^ {1} + \tilde {S} ^ {2 -}, \quad \mathbb {Q} _ {2} \equiv S _ {1} ^ {-} - \tilde {Q} _ {2 -}, \tag {A.16}
$$

$$
\mathbb {Q} _ {1} ^ {\dagger} \equiv S _ {1} ^ {-} + \tilde {Q} _ {2 -}, \quad \mathbb {Q} _ {2} ^ {\dagger} \equiv Q _ {-} ^ {1} - \tilde {S} ^ {2 -}.
$$

We find two supercharges  $G_{-}^{1}$  and  $H_{1}^{-}$  (in the notation of (A.7)) that are preserved by both (full or half) line defects and the chiral algebra,

$$
G _ {-} ^ {1} = \frac {1}{2} \left(\mathbb {Q} _ {1} + \mathbb {Q} _ {1} ^ {\dagger} - \mathbb {Q} _ {2} + \mathbb {Q} _ {2} ^ {\dagger}\right) = Q _ {-} ^ {1} + \tilde {Q} _ {2 -}, \tag {A.17}
$$

$$
H _ {1} ^ {-} = \frac {1}{2} \left(\mathbb {Q} _ {1} + \mathbb {Q} _ {1} ^ {\dagger} + \mathbb {Q} _ {2} - \mathbb {Q} _ {2} ^ {\dagger}\right) = S _ {1} ^ {-} + \tilde {S} ^ {2 \dot {-}}, \tag {A.17}
$$

They satisfy the hermicity conditions  $(G_{-}^{1})^{\dagger} = H_{1}^{-}$ . These are precisely the two supercharges preserved by multiple half line defects lying on the 12-plane (A.15). Let us record the anticommutators of the supercharges  $G_{-}^{1}$  and  $H_{1}^{-}$  again,

$$
\left\{G _ {-} ^ {1}, G _ {-} ^ {1} \right\} = \left\{H _ {1} ^ {-}, H _ {1} ^ {-} \right\} = 0, \quad \left\{G _ {-} ^ {1}, H _ {1} ^ {-} \right\} = 2 \hat {L} _ {0}, \tag {A.18}
$$

where  $\hat{L}_0 = \frac{1}{2} (E - 2R - (j_1 + j_2)) = 0$  is the defining condition for the line defect Schur operators (A.9). Note that  $\hat{L}_0$  involves the rotation  $\mathcal{M} = j_{1} + j_{2}$  on the chiral algebra plane. Since  $\hat{L}_0$  shows up on the righthand side of the anticommutator of supercharges preserved by the line defect, the defect must lie on the plane transverse to the chiral algebra plane, with which it intersects at a point. In the convention in (A.16), we have chosen the chiral algebra plane to be the 34-plane where  $x^{1} = x^{2} = 0$  [13], and the line defects lie on the 12-plane where  $x^{3} = x^{4} = 0$ . See Figure 11 for the incidence geometry of line defects and the chiral algebra plane.

Notice that since  $G_{-}^{1}$  and  $H_{1}^{-}$  only involve of  $Q$ 's and  $S$ 's, respectively, the  $\hat{L}_{\pm}$  generators cannot be both exact under either of them. Thus it seems difficult to construct the chiral algebra for the defect operators by generalizing [13].

On the other hand, we can consider the following linear combinations of  $G_{-}^{1}$  and  $H_{1}^{-}$ ,

$$
\mathfrak {Q} _ {1} \equiv G ^ {1} _ {-} + H _ {1} ^ {-} = \mathbb {Q} _ {1} + \mathbb {Q} _ {1} ^ {\dagger}, \quad \mathfrak {Q} _ {2} \equiv - G ^ {1} _ {-} + H _ {1} ^ {-} = \mathbb {Q} _ {2} - \mathbb {Q} _ {2} ^ {\dagger}, \tag {A.19}
$$

such that  $\hat{L}_{\pm 1}$  are both  $\mathfrak{Q}_1$  -exact and  $\mathfrak{Q}_2$  -exact,

$$
\left\{\mathfrak {Q} _ {1}, \tilde {Q} _ {1 -} \right\} = \left\{\mathfrak {Q} _ {2}, - Q ^ {2} _ {-} \right\} = \hat {L} _ {- 1} = P _ {-} \cdot + R ^ {2} _ {1}, \tag {A.20}
$$

$$
\left\{\mathfrak {Q} _ {1}, S _ {2} ^ {-} \right\} = \left\{\mathfrak {Q} _ {2}, \tilde {S} ^ {1 -} \right\} = \hat {L} _ {+ 1} = K ^ {- -} - R _ {2} ^ {1}.
$$

However, the disadvantage of these combinations  $\mathfrak{Q}_i$  is that they are not nilpotent, but

square to  $\hat{L}_0$

$$
\left\{\mathfrak {Q} _ {1}, \mathfrak {Q} _ {1} \right\} = 4 \hat {L} _ {0}, \quad \left\{\mathfrak {Q} _ {2}, \mathfrak {Q} _ {2} \right\} = - 4 \hat {L} _ {0}. \tag {A.21}
$$

Note that  $\mathfrak{Q}_1$  and  $\mathfrak{Q}_2$  anticommute  $\{\mathfrak{Q}_1,\mathfrak{Q}_2\} = 0$

# B Framed Quivers for the Argyres-Douglas Theories

In this Appendix we derive the core charges and the framed BPS quivers for the generators of line defects in the  $A_{3}$  and  $A_{4}$  Argyres-Douglas theories, following closely the example of the  $A_{2}$  Argyres-Douglas theory in [65].

Given a fixed point on the moduli space and a choice of the half plane, the BPS quiver, if exists, is unique. The charges  $\{\gamma_i\}$  of the nodes of the quiver are called a seed. Associated to this seed is a cone  $\mathcal{C}$  on the charge lattice  $\Gamma$ ,

$$
\mathcal {C} = \left\{\sum_ {i} a _ {i} \gamma_ {i} \in \Gamma \mid a _ {i} \in \mathbb {R} _ {\geq 0} \right\}, \tag {B.1}
$$

generated by the seed with non-negative coefficients. We also define the dual cone  $\check{C}$  in the charge lattice  $\Gamma$  as

$$
\check {\mathcal {C}} = \left\{\check {\gamma} \in \Gamma \mid \langle \check {\gamma}, \gamma \rangle \geq 0 \right\}. \tag {B.2}
$$

One important feature of the dual cone is that the UV line defects in  $\check{\mathcal{C}}$  satisfy a universal OPE,

$$
L _ {\gamma_ {1}} L _ {\gamma_ {2}} = q ^ {\frac {1}{2} \langle \gamma_ {1}, \gamma_ {2} \rangle} L _ {\gamma_ {1} + \gamma_ {2}}, \quad \gamma_ {1}, \gamma_ {2} \in \check {\mathcal {C}}. \tag {B.3}
$$

where  $\gamma_{i}$ 's are the core charges of the defect  $L_{i}$ .

Starting from an initial seed, we can generate other seeds by mutation and obtain their associated dual cones. It follows that on the charge lattice for the defects, there are many distinct dual cones  $\check{C}$ , inside which the defect OPE takes the simple form above. For the Argyres-Douglas theories considered in this paper, the dual cones cover the full charge lattice. This simplifies the study of defect OPE significantly. In particular, the defect OPEs are completely encoded in the OPEs between those defects whose core charges lie at the boundaries of the dual cones. These distinguished defects will be called the generators of defects. In the rest of this Appendix, we will compute the core charges of these generators and their associated framed BPS quiver in the  $A_{3}$  and  $A_{4}$  Argyres-Douglas theories.

# B.1  $A_{3}$  Argyres-Douglas Theory

As shown in Figure 12, from the initial seed  $\{\gamma_1,\gamma_2,\gamma_3\}$ , we generate fourteen seeds by mutation. Each of the fourteen seeds is associated to a dual cone  $\tilde{C}$  in the space of line defects. Out of the fourteen seeds, eight of them give rise to degenerate dual cones that are of higher codimensions. The remaining six non-degenerate dual cones are

$$
\check {\mathcal {C}} _ {\{\gamma_ {1}, \gamma_ {2}, \gamma_ {3} \}} = \left\{a _ {1} \gamma_ {1} + a _ {2} \gamma_ {2} + a _ {3} \gamma_ {3} \Big | a _ {1} + a _ {3} \geq 0, a _ {2} \leq 0 \right\},
$$

$$
\check {\mathcal {C}} _ {\{\gamma_ {1}, - \gamma_ {2}, \gamma_ {3} \}} = \left\{a _ {1} \gamma_ {1} + a _ {2} \gamma_ {2} + a _ {3} \gamma_ {3} \mid 0 \geq a _ {1} + a _ {3} \geq a _ {2} \right\},
$$

$$
\check {\mathcal {C}} _ {\{- \gamma_ {1}, \gamma_ {1} + \gamma_ {2} + \gamma_ {3}, - \gamma_ {3} \}} = \left\{a _ {1} \gamma_ {1} + a _ {2} \gamma_ {2} + a _ {3} \gamma_ {3} \mid a _ {1} + a _ {3} \geq 0, a _ {2} \geq 0 \right\}, \tag {B.4}
$$

$$
\check {\mathcal {C}} _ {\{\gamma_ {2} + \gamma_ {3}, - \gamma_ {1} - \gamma_ {2} - \gamma_ {3}, \gamma_ {1} + \gamma_ {2} \}} = \left\{a _ {1} \gamma_ {1} + a _ {2} \gamma_ {2} + a _ {3} \gamma_ {3} \Big | a _ {1} + a _ {3} \leq 0, a _ {2} \geq 0 \right\},
$$

$$
\check {\mathcal {C}} _ {\{- \gamma_ {1}, - \gamma_ {2}, - \gamma_ {3} \}} = \left\{a _ {1} \gamma_ {1} + a _ {2} \gamma_ {2} + a _ {3} \gamma_ {3} \Big | a _ {2} \geq a _ {1} + a _ {3} \geq 2 a _ {2} \right\},
$$

$$
\check {\mathcal {C}} _ {\{- \gamma_ {1} - \gamma_ {2}, \gamma_ {2}, - \gamma_ {2} - \gamma_ {3} \}} = \left\{a _ {1} \gamma_ {1} + a _ {2} \gamma_ {2} + a _ {3} \gamma_ {3} \mid 0 \geq 2 a _ {2} \geq a _ {1} + a _ {3} \right\}.
$$

For example, to obtain the dual cone of the seed  $\{\gamma_1, -\gamma_2, \gamma_3\}$ , we first find those charges which have positive Dirac pairing with the seed, and then (right) mutate them back to the original seed to see how this dual cone is embedded in the charge lattice. Explicitly, we have

$$
\begin{array}{l} \check {\mathcal {C}} _ {\{\gamma_ {1}, - \gamma_ {2}, \gamma_ {3} \}} = \mu_ {R _ {2}} \left(\left\{a _ {1} \gamma_ {1} + a _ {2} \gamma_ {2} + a _ {3} \gamma_ {3} \mid a _ {1} + a _ {3} \leq 0, a _ {2} \leq 0 \right\}\right) \tag {B.5} \\ = \left\{a _ {1} \gamma_ {1} + a _ {2} \gamma_ {2} + a _ {3} \gamma_ {3} \mid 0 \geq a _ {1} + a _ {3} \geq a _ {2} \right\}. \\ \end{array}
$$

Note that since  $\gamma_{1} - \gamma_{3}$  is a flavor node, only the combination  $a_1 + a_3$  shows up but not  $a_1$  and  $a_3$  individually. We show the two-dimensional projection  $(a_{2},a_{1} + a_{3})$  of the geometry of the dual cones in Figure 13. The six boundaries of the dual cones are two-dimensional half-planes, each generated by a flavor charge  $\gamma_{1} - \gamma_{3}$  and an electromagnetic core charge  $A_{i}$ ,  $B_{i}$  with  $i = 1,2,3$ . The core charges, defined as the images of the RG map [65], are

$$
\begin{array}{l l} \mathbf {R G} (A _ {1}) = \gamma^ {\prime}, & \mathbf {R G} (A _ {2}) = - \gamma^ {\prime} - \gamma_ {2}, \\ \mathbf {R G} (B) & \mathbf {R G} (B) = - \gamma^ {\prime}, \end{array} \quad \begin{array}{l l} \mathbf {R G} (A _ {3}) = - \gamma^ {\prime}, \\ \mathbf {R G} (B) = - \gamma^ {\prime}. \end{array} \tag {B.6}
$$

$$
\mathbf {R G} (B _ {1}) = - \gamma_ {2}, \quad \mathbf {R G} (B _ {2}) = - 2 \gamma^ {\prime} - \gamma_ {2}, \quad \mathbf {R G} (B _ {3}) = \gamma_ {2},
$$

where  $\gamma'$  is any charge vector of the form  $x\gamma_{1} + (1 - x)\gamma_{3}$ . The associated framed BPS quivers for these six defects are given in Figure 14. By applying the mutation method to these framed quivers, we obtain the generating functions for these defects (5.11).

![](images/de0d000c87d63da2b86133a18b8c1d6ca31f1a7bb2c64f604001eabd5387b6be.jpg)  
Figure 12: The fourteen seeds of the  $A_{3}$  Argyres-Douglas theory.  $\mu_{L_i}$  denotes the left mutation with respect to the  $i$ -th node counting from the top.

![](images/89cc956692828a6ca657f0b57e1caf348c5e59b7db9aed508d15695f8c5bae6f.jpg)

![](images/342eb12cfb503fc9760b8dbf1805ed0219c50afa72be541305ada2e49ed06ee4.jpg)  
(a)  $A_{1}$

![](images/2996bf94f2bb32ecf1208df0e7183966261dadb0f216ac3d77c1d9b46ee395a8.jpg)  
Figure 13: The six dual cones of the  $A_{3}$  Argyres-Douglas theory. Here we show the two-dimensional projection  $(a_{2}, a_{1} + a_{3})$  of the three-dimensional space  $\Gamma \otimes_{\mathbb{Z}} \mathbb{R}$ , which is identified as  $\mathbb{R} \oplus \mathbb{R} \oplus \mathbb{R}$  by expressing the seed as  $\gamma = a_{1}\gamma_{1} + a_{2}\gamma_{2} + a_{3}\gamma_{3}$ . The black arrows are the projection of the boundary half-planes of the dual cones. The red dots are the generators  $A_{i}$ ,  $B_{i}$  for the line defects.  
(b)  $A_{2}$

![](images/16c1d2a6bff4139677a878de8ccad983e3ace5830a181294c616e523da37002a.jpg)  
(c)  $A_{3}$

![](images/69713a04780e5e8d1fe235b4cbd9312c1778ffb43ae78ddf03e333a887f92202.jpg)  
(d)  $B_{1}$

![](images/b0e4f450bdb3843d5b248fa4bdebfd4bee185e917b20b29ea1a4c85dedfb0f10.jpg)  
(e)  $B_{2}$

![](images/13b9de649959891a91fab9f1db819bf5994b294398ce37efc7e0c6b891ca6e94.jpg)  
(f)  $B_{3}$  
Figure 14: The framed quivers for the six generators  $A_{i}$  and  $B_{i}$  ( $i = 1,2,3$ ) of the  $A_{3}$  Argyres-Douglas theory. The core charges are labeled above the framed nodes (the square nodes).  $\gamma'$  is any charge vector of the form  $\gamma' = x\gamma_{1} + (1 - x)\gamma_{3}$  with real  $x$ .

# B.2  $A_{4}$  Argyres-Douglas Theory

For the  $A_4$  Argyres-Douglas theory, there are 42 seeds and their associated dual cones are listed below.

$$
\check {\mathcal {C}} _ {\{\gamma_ {1}, \gamma_ {2}, \gamma_ {3}, \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \Big | a _ {2} \leq 0, a _ {1} + a _ {3} \geq 0, a _ {2} + a _ {4} \leq 0, a _ {3} \geq 0 \right\},
$$

$$
\check {\mathcal {C}} _ {\{\gamma_ {1}, \gamma_ {2}, \gamma_ {3}, - \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \mid a _ {2} \leq 0, a _ {1} + a _ {3} \geq 0, a _ {3} \geq a _ {2} + a _ {4}, a _ {3} \leq 0 \right\},
$$

$$
\check {\mathcal {C}} _ {\{\gamma_ {1}, \gamma_ {2} + \gamma_ {3}, - \gamma_ {3}, \gamma_ {3} + \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \mid a _ {2} \leq 0, a _ {1} + a _ {3} \geq 0, a _ {2} + a _ {4} \geq 0, a _ {3} \geq 0 \right\},
$$

$$
\check {\mathcal {C}} _ {\{\gamma_ {1}, - \gamma_ {2}, \gamma_ {3}, \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \left| a _ {1} + a _ {3} \geq a _ {2}, a _ {1} + a _ {3} \leq 0, a _ {1} + a _ {3} \geq a _ {2} + a _ {4}, a _ {3} \geq 0 \right. \right\},
$$

$$
\check {\mathcal {C}} _ {\{- \gamma_ {1}, \gamma_ {1} + \gamma_ {2}, \gamma_ {3}, \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \mid a _ {2} \geq 0, a _ {1} + a _ {3} \geq 0, a _ {2} + a _ {4} \leq 0, a _ {3} \geq 0 \right\},
$$

$$
\check {\mathcal {C}} _ {\{\gamma_ {1}, - \gamma_ {2}, - \gamma_ {3}, \gamma_ {3} + \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \left| a _ {1} + a _ {3} \geq a _ {2}, a _ {2} + a _ {4} \leq 0, a _ {2} + a _ {4} \geq a _ {1} + a _ {3}, a _ {3} \geq 0 \right. \right\},
$$

$$
\check {\mathcal {C}} _ {\{- \gamma_ {1}, - \gamma_ {2}, \gamma_ {3}, \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \Big | a _ {2} \leq 0, a _ {1} + a _ {3} \leq a _ {2}, a _ {1} + a _ {3} \geq a _ {2} + a _ {4}, a _ {3} \geq 0 \right\},
$$

$$
\check {\mathcal {C}} _ {\{\gamma_ {1}, - \gamma_ {2}, \gamma_ {3}, - \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \left| a _ {1} + a _ {3} \geq a _ {2}, a _ {1} + a _ {3} \leq 0, a _ {1} + 2 a _ {3} \geq a _ {2} + a _ {4}, a _ {3} \leq 0 \right. \right\},
$$

$$
\check {\mathcal {C}} _ {\{- \gamma_ {1}, \gamma_ {1} + \gamma_ {2} + \gamma_ {3}, - \gamma_ {3}, \gamma_ {3} + \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \mid a _ {2} \geq 0, a _ {1} + a _ {3} \geq 0, a _ {2} + a _ {4} \geq 0, a _ {3} \geq 0 \right\},
$$

$$
\check {\mathcal {C}} _ {\{\gamma_ {2}, - \gamma_ {1} - \gamma_ {2}, \gamma_ {3}, \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \mid a _ {2} \geq 0, a _ {1} + a _ {3} \leq 0, a _ {1} + a _ {3} \geq a _ {2} + a _ {4}, a _ {3} \geq 0 \right\},
$$

$$
\check {\mathcal {C}} _ {\{- \gamma_ {1}, \gamma_ {1} + \gamma_ {2}, \gamma_ {3}, - \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \mid a _ {2} \geq 0, a _ {1} + a _ {3} \geq 0, a _ {3} \geq a _ {2} + a _ {4}, a _ {3} \leq 0 \right\},
$$

$$
\check {\mathcal {C}} _ {\{\gamma_ {1}, \gamma_ {2} + \gamma_ {3}, - \gamma_ {3}, - \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \Big | a _ {2} \leq 0, a _ {1} + a _ {3} \geq 0, a _ {3} \leq a _ {2} + a _ {4}, a _ {2} + a _ {4} \leq 0 \right\},
$$

$$
\check {\mathcal {C}} _ {\{\gamma_ {1}, \gamma_ {2} + \gamma_ {3}, \gamma_ {4}, - \gamma_ {3} - \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \Big | a _ {2} \leq 0, a _ {1} + a _ {3} \geq 0, a _ {2} + a _ {4} \geq 0, a _ {3} \leq 0 \right\},
$$

$$
\begin{array}{l} \check {\mathcal {C}} _ {\{\gamma_ {1}, - \gamma_ {2} - \gamma_ {3}, \gamma_ {2}, \gamma_ {3} + \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \Big | a _ {1} + a _ {3} \geq a _ {2}, a _ {1} + a _ {3} \leq 0, a _ {2} + a _ {4} \geq 0, a _ {3} \geq 0 \right\}, \\ \check {\mathcal {C}} _ {\{\gamma_ {1}, - \gamma_ {2}, - \gamma_ {3}, - \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \left| a _ {1} + a _ {3} \geq a _ {2}, a _ {3} \geq a _ {2} + a _ {4} \geq a _ {1} + 2 a _ {3}, a _ {1} + a _ {3} \geq a _ {2} + a _ {4} \right. \right\}, \\ \check {\mathcal {C}} _ {\{- \gamma_ {1}, - \gamma_ {2}, \gamma_ {3}, - \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \mid a _ {2} \leq 0, a _ {1} + a _ {3} \leq a _ {2}, a _ {1} + 2 a _ {3} \geq a _ {2} + a _ {4}, a _ {3} \leq 0 \right\}, \\ \check {\mathcal {C}} _ {\{- \gamma_ {1}, \gamma_ {1} + \gamma_ {2} + \gamma_ {3}, \gamma_ {4}, - \gamma_ {3} - \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \mid a _ {2} \geq 0, a _ {1} + a _ {3} \geq 0, a _ {2} + a _ {4} \geq 0, a _ {3} \leq 0 \right\}, \\ \check {\mathcal {C}} _ {\{\gamma_ {2} + \gamma_ {3}, - \gamma_ {1} - \gamma_ {2} - \gamma_ {3}, \gamma_ {1} + \gamma_ {2}, \gamma_ {3} + \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \Big | a _ {2} \geq 0, a _ {1} + a _ {3} \leq 0, a _ {2} + a _ {4} \geq 0, a _ {3} \geq 0 \right\}, \\ \check {\mathcal {C}} _ {\{\gamma_ {2} + \gamma_ {3}, - \gamma_ {1} - \gamma_ {2}, - \gamma_ {3}, \gamma_ {3} + \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \mid a _ {2} \geq 0, a _ {1} + a _ {3} \leq a _ {2} + a _ {4}, a _ {2} + a _ {4} \leq 0, a _ {3} \geq 0 \right\}, \\ \check {\mathcal {C}} _ {\{- \gamma_ {1}, \gamma_ {1} + \gamma_ {2} + \gamma_ {3}, - \gamma_ {3}, - \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \mid a _ {2} \geq 0, a _ {1} + a _ {3} \geq 0, a _ {2} + a _ {4} \leq 0, a _ {3} \leq a _ {2} + a _ {4} \right\}, \\ \check {\mathcal {C}} _ {\{\gamma_ {2}, - \gamma_ {1} - \gamma_ {2}, \gamma_ {3}, - \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \Big | a _ {2} \geq 0, a _ {1} + a _ {3} \leq 0, a _ {1} + 2 a _ {3} \geq a _ {2} + a _ {4}, a _ {3} \leq 0 \right\}, \\ \check {\mathcal {C}} _ {\{\gamma_ {1}, - \gamma_ {2}, - \gamma_ {3} - \gamma_ {4}, \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \left| a _ {1} + a _ {3} \geq a _ {2}, a _ {3} \geq a _ {2} + a _ {4}, a _ {1} + a _ {3} \leq a _ {2} + a _ {4}, a _ {3} \leq 0 \right. \right\}, \\ \check {\mathcal {C}} _ {\{\gamma_ {1}, - \gamma_ {2} - \gamma_ {3}, \gamma_ {2}, - \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \left| 0 \geq a _ {1} + a _ {3} \geq a _ {2}, a _ {3} \leq a _ {2} + a _ {4}, a _ {1} + a _ {3} \geq a _ {2} + a _ {4} \right. \right\}, \\ \check {\mathcal {C}} _ {\{\gamma_ {1}, - \gamma_ {2} - \gamma_ {3}, \gamma_ {2} + \gamma_ {3} + \gamma_ {4}, - \gamma_ {3} - \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \left| 0 \geq a _ {1} + a _ {3} \geq a _ {2}, a _ {2} + a _ {4} \geq 0, a _ {3} \leq 0 \right. \right\}, \\ \check {\mathcal {C}} _ {\{- \gamma_ {2} - \gamma_ {3}, - \gamma_ {1}, \gamma_ {1} + \gamma_ {2}, \gamma_ {3} + \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \mid a _ {2} \leq 0, a _ {1} + a _ {3} \leq a _ {2}, a _ {2} + a _ {4} \geq 0, a _ {3} \geq 0 \right\}, \\ \check {\mathcal {C}} _ {\{- \gamma_ {1}, - \gamma_ {2}, - \gamma_ {3}, \gamma_ {3} + \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \left| a _ {2} \geq a _ {1} + a _ {3} \geq 2 a _ {2} + a _ {4}, a _ {1} + a _ {3} \leq a _ {2} + a _ {4}, a _ {3} \geq 0 \right. \right\}, \\ \check {\mathcal {C}} _ {\{- \gamma_ {1}, - \gamma_ {2} - \gamma_ {3}, \gamma_ {1} + \gamma_ {2}, - \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \left| 0 \geq a _ {2} \geq a _ {1} + a _ {3} \geq a _ {2} + a _ {4} \geq a _ {3} \right. \right\}, \\ \check {\mathcal {C}} _ {\{\gamma_ {1}, \gamma_ {4}, - \gamma_ {2} - \gamma_ {3} - \gamma_ {4}, \gamma_ {2} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \Big | 0 \geq a _ {2} + a _ {4} \geq a _ {1} + a _ {3} \geq a _ {2}, a _ {3} \leq a _ {2} + a _ {4} \right\}, \\ \end{array}
$$

$$
\check {\mathcal {C}} _ {\{- \gamma_ {2} - \gamma_ {3}, - \gamma_ {1}, \sum_ {i = 1} ^ {4} \gamma_ {i}, - \gamma_ {3} - \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \Big | 0 \geq a _ {2} \geq a _ {1} + a _ {3}, a _ {2} + a _ {4} \geq 0, a _ {3} \leq 0 \right\},
$$

$$
\check {\mathcal {C}} _ {\{- \gamma_ {1}, - \gamma_ {2}, - \gamma_ {3} - \gamma_ {4}, \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \Big | a _ {1} + a _ {3} \leq a _ {2}, a _ {1} + 2 a _ {3} \geq 2 a _ {2} + a _ {4}, a _ {1} + a _ {3} \leq a _ {2} + a _ {4}, a _ {3} \leq 0 \right\},
$$

$$
\check {\mathcal {C}} _ {\{- \gamma_ {1}, - \gamma_ {2}, - \gamma_ {3}, - \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \Big | a _ {2} \geq a _ {1} + a _ {3} \geq a _ {2} + a _ {4}, a _ {2} + a _ {4} \geq 0 a _ {1} + 2 a _ {3} \geq 2 a _ {2} + a _ {4} \right\},
$$

$$
\check {\mathcal {C}} _ {\{\gamma_ {2} + \gamma_ {3}, - \gamma_ {1} - \gamma_ {2} - \gamma_ {3}, \sum_ {i = 1} ^ {4} \gamma_ {i}, - \gamma_ {3} - \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \Big | a _ {2} \geq 0, a _ {1} + a _ {3} \leq 0, a _ {2} + a _ {4} \geq 0, a _ {3} \leq 0 \right\},
$$

$$
\check {\mathcal {C}} _ {\{\gamma_ {2} + \gamma_ {3}, \gamma_ {4}, - \gamma_ {3} - \gamma_ {4}, - \gamma_ {1} - \gamma_ {2} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \mid a _ {2} \geq 0, a _ {3} \geq a _ {2} + a _ {4} \geq a _ {1} + a _ {3}, a _ {3} \leq 0 \right\},
$$

$$
\check {\mathcal {C}} _ {\{- \gamma_ {2} - \gamma_ {3}, - \gamma_ {1} - \gamma_ {2}, \gamma_ {2}, \gamma_ {3} + \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \Big | a _ {2} \leq 0, a _ {1} + a _ {3} \leq 2 a _ {2} + a _ {4}, a _ {2} + a _ {4} \leq 0, a _ {3} \geq 0 \right\},
$$

$$
\check {\mathcal {C}} _ {\{\gamma_ {2} + \gamma_ {3}, - \gamma_ {1} - \gamma_ {2} - \gamma_ {3}, \gamma_ {1} + \gamma_ {2}, - \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \Big | a _ {2} \geq 0, 0 \geq a _ {1} + a _ {3} \geq a _ {2} + a _ {4} \geq a _ {3} \right\},
$$

$$
\check {\mathcal {C}} _ {\{\gamma_ {2} + \gamma_ {3}, - \gamma_ {1} - \gamma_ {2}, - \gamma_ {3}, - \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \left| a _ {2} \geq 0, a _ {1} + a _ {3} \geq a _ {2} + a _ {4} \geq a _ {1} + 2 a _ {3}, a _ {3} \geq a _ {2} + a _ {4} \right. \right\},
$$

$$
\check {\mathcal {C}} _ {\{- \gamma_ {1}, \gamma_ {4}, - \gamma_ {2} - \gamma_ {3} - \gamma_ {4}, \gamma_ {1} + \gamma_ {2} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \Big | a _ {2} \geq a _ {1} + a _ {3} \geq 2 a _ {2} + a _ {4}, a _ {1} + a _ {3} \leq a _ {2} + a _ {4}, a _ {3} \leq a _ {2} + a _ {4} \right\},
$$

$$
\check {\mathcal {C}} _ {\{- \gamma_ {2} - \gamma_ {3}, \gamma_ {2} + \gamma_ {3} + \gamma_ {4}, - \sum_ {i = 1} ^ {4} \gamma_ {i}, \gamma_ {1} + \gamma_ {2} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \Big | a _ {2} \leq 0, a _ {1} + a _ {3} \leq 2 a _ {2} + a _ {4}, 0 \geq a _ {2} + a _ {4} \geq a _ {3} \right\},
$$

$$
\check {\mathcal {C}} _ {\{\gamma_ {2} + \gamma_ {3}, \gamma_ {4}, - \sum_ {i = 1} ^ {4} \gamma_ {i}, \gamma_ {1} + \gamma_ {2} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \Big | a _ {2} \geq 0, 0 \geq a _ {2} + a _ {4} \geq a _ {1} + a _ {3}, a _ {3} \leq a _ {2} + a _ {4}, a _ {3} \geq 0 \right\},
$$

$$
\check {\mathcal {C}} _ {\{- \gamma_ {2} - \gamma_ {3}, - \gamma_ {1} - \gamma_ {2}, \gamma_ {2} + \gamma_ {3} + \gamma_ {4}, - \gamma_ {3} - \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \Big | a _ {2} \leq 0, a _ {1} + a _ {3} \leq 2 a _ {2} + a _ {4}, 0 \geq a _ {3} \geq a _ {2} + a _ {4} \right\},
$$

$$
\check {\mathcal {C}} _ {\{\gamma_ {4}, - \gamma_ {2} - \gamma_ {3} - \gamma_ {4}, \gamma_ {2}, - \gamma_ {1} - \gamma_ {2} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \Big | a _ {3} \geq a _ {2} + a _ {4} \geq a _ {1} + a _ {3} \geq 2 a _ {2} + a _ {4}, a _ {1} + 2 a _ {3} \leq 2 a _ {2} + a _ {4} \right\},
$$

$$
\check {\mathcal {C}} _ {\{- \gamma_ {2} - \gamma_ {3}, - \gamma_ {1} - \gamma_ {2}, \gamma_ {2}, - \gamma_ {4} \}} = \left\{\sum_ {i = 1} ^ {4} a _ {i} \gamma_ {i} \Bigg | a _ {2} \leq 0, a _ {3} \geq a _ {2} + a _ {4}, a _ {1} + a _ {3} \geq a _ {2} + a _ {4}, a _ {1} + 2 a _ {3} \leq 2 a _ {2} + a _ {4} \right\}.
$$

Each of the 42 dual cones is generated by four boundary rays, and there are in total fourteen distinct boundary rays which are the generators of line defects in the  $A_4$  Argyres-Douglas theory. We will denote them by  $A_i$  and  $B_i$  with  $i = 1, \dots, 7$ . Their core charges are

$$
\mathbf {R G} (A _ {1}) = - \gamma_ {2} + \gamma_ {4}, \quad \mathbf {R G} (A _ {2}) = \gamma_ {2} - \gamma_ {4}, \quad \mathbf {R G} (A _ {3}) = - \gamma_ {1} - \gamma_ {2} + \gamma_ {4},
$$

$$
\mathbf {R G} (A _ {4}) = \gamma_ {1} - \gamma_ {3} - \gamma_ {4}, \quad \mathbf {R G} (A _ {5}) = - \gamma_ {1} + \gamma_ {3}, \quad \mathbf {R G} (A _ {6}) = \gamma_ {1} - \gamma_ {3},
$$

$$
\mathbf {R G} \left(A _ {7}\right) = - \gamma_ {1} - \gamma_ {4}, \tag {B.7}
$$

$$
\mathbf {R G} (B _ {1}) = - \gamma_ {1}, \qquad \mathbf {R G} (B _ {2}) = - \gamma_ {3} - \gamma_ {4}, \quad \mathbf {R G} (B _ {3}) = - \gamma_ {2} - \gamma_ {3},
$$

$$
\mathbf {R G} (B _ {4}) = - \gamma_ {1} - \gamma_ {2}, \quad \mathbf {R G} (B _ {5}) = - \gamma_ {4}, \quad \mathbf {R G} (B _ {6}) = \gamma_ {1},
$$

$$
\mathbf {R G} (B _ {\tau}) = \gamma_ {4}.
$$

Their associated framed quivers are shown in Figure 15. By applying the mutation method to these framed quivers, we obtain the generating functions for these defects (5.19).

# C Affine Characters of Kac-Moody Algebra at Negative Level

In this appendix we review a generalization of the Weyl-Kac formula, known as the Kazhdan-Lusztig conjecture [94], for affine characters of Kac-Moody algebra at negative levels, following [95]. We will compute the affine characters for several modules in  $\widehat{su(2)}_{-\frac{4}{3}}$  and  $\widehat{so(8)}_{-2}$ , which are the chiral algebras of the  $A_3$  Argyres-Douglas theory and the  $SU(2)$  with  $N_f = 4$  flavors theory, respectively.

# C.1 Generalities on Affine Lie Algebra

We begin by reviewing some basic facts about affine Lie algebra (see, for example, [93]). Let  $\mathfrak{g}$  be an affine Lie algebra associated with a finite dimensional simple Lie algebra  $\underline{\mathfrak{g}}$  of rank  $r$ . An affine weight  $\lambda$  of  $\mathfrak{g}$  will be denoted by

$$
\lambda = (\underline {{\lambda}}; k; n), \tag {C.1}
$$

where  $\underline{\lambda}$  is a weight of the finite dimensional Lie algebra  $\underline{\mathfrak{g}}$ .  $k$  is the level of the weight and  $n$  is the eigenvalue with respect to  $-L_0$ . The inner product between weights is given by

$$
\left(\lambda_ {1}, \lambda_ {2}\right) = \left(\underline {{\lambda}} _ {1}, \underline {{\lambda}} _ {2}\right) + k _ {1} n _ {2} + k _ {2} n _ {1}. \tag {C.2}
$$

![](images/e69713c9bcc352e354dc8d5436e618962ed92565ff6c298c8b2ae93766612280.jpg)  
(a)  $A_{1}$

![](images/cea69cf37c9a2df4d0cd9dfbd2efc2553523d770fbfe47b69fa0f300a2f9f18a.jpg)  
(b)  $A_{2}$

![](images/a5b234002a5c75557268280eabd181d8d2843dd91c2aeab56e9d4d02850a946e.jpg)  
(c)  $A_{3}$

![](images/f533b1b3ec9755401182eea1a2a28328325b4ffbbf252c502aed02590a6e872b.jpg)  
(d)  $A_4$

![](images/b7b1825b00858ccdbef7f57875500d9f97a51c223565d1f8272557caba1f3a6c.jpg)  
(e)  $A_{5}$

![](images/bbcfe0a05bddccebc10802ec109b11287e2bf7d4c58007d6917140c35904eb67.jpg)  
(f)  $A_{6}$

![](images/a27da47b057216d7d983bc56821f30a724ab8864f38b7a10da0caa75f3385a0b.jpg)  
(g)  $A_7$

![](images/dc36310047d4501d2e0dd4dc136e57e8fdb2961b0d128bd3b444115abaef286c.jpg)  
(h)  $B_{1}$

![](images/9514b9a7e4da91d44f975a0226ec8c6bd2cd4bd7d9f0a0c0110441263c2038da.jpg)  
(i)  $B_{2}$

![](images/308bf11ac888daded7fdfe4c3a3d05c77ff86e84a60c465e66f6cd2347373696.jpg)  
(j)  $B_{3}$

![](images/630e123909d0ae9a5f8d0ba922c97ce21ea6943f8989edd8292c0640f9f9923e.jpg)  
(k)  $B_{4}$

![](images/4e518445c7314e1049d5905773b2fd51d2827e142ed7010f01dd2e6cbf4877ae.jpg)  
(1)  $B_{5}$

![](images/f2c671f063ae8ec56d6ee8f83e315dce9005987cae3f19ca391e19d4423bd2d3.jpg)  
(m)  $B_{6}$

![](images/ba7cecfd651a66564f3e4a116cbb560b33976fdea536319df8854025bc4a5d06.jpg)  
(n)  $B_{7}$  
Figure 15: The framed quivers associated to the fourteen generators of line defects  $A_{i}$ ,  $B_{i}$  in the  $A_{4}$  Argyres-Douglas theory. The core charges are labeled above the framed nodes (the square nodes).

The simple roots of the affine Lie algebra  $\mathfrak{g}$  consist of

$$
\alpha_ {0} = (- \theta ; 0; 1) = - \theta + \delta , \tag {C.3}
$$

$$
\alpha_ {i} = (\underline {{\alpha}} _ {i}; 0; 0), \quad i = 1, \dots , r,
$$

where  $\underline{\alpha}_i$ 's are the simple roots of  $\underline{\mathfrak{g}}$ .  $\theta$  is the highest root of  $\underline{\mathfrak{g}}$  normalized such that  $|\theta|^2 = 2$ .  $\delta$  is defined as  $\delta = (0;0;1)$ . Since  $(\delta, \delta) = 0$ ,  $n\delta$  is called an imaginary root for all  $n$ . The other roots are said to be real.

The set of positive roots  $\Delta_{+}$  of the affine Lie algebra is

$$
\Delta_ {+} = \left\{\alpha + n \delta \mid n > 0, \alpha \in \underline {{\Delta}} \right\} \cup \left\{\alpha \mid \alpha \in \underline {{\Delta}} _ {+} \right\}, \tag {C.4}
$$

where  $\underline{\Delta}$  and  $\underline{\Delta}_{+}$  are the sets of roots and positive roots of the finite dimensional Lie algebra  $\mathfrak{g}$ , respectively. The set of real positive roots  $\Delta_{+}^{re}$  is defined as

$$
\Delta_ {+} ^ {r e} = \Delta_ {+} / \{n \delta \} = \{\alpha + n \delta \mid n > 0, \alpha \in \underline {{\Delta}}, \alpha \neq 0 \} \cup \{\alpha | \alpha \in \underline {{\Delta}} _ {+} \}. \tag {C.5}
$$

The Cartan matrix  $A_{ij}$  of the affine Lie algebra  $\mathfrak{g}$  is defined as  $A_{ij} = (\alpha_i, \alpha_j^\vee)$  with  $0 \leq i, j \leq r$ . Note in particular,  $A_{00} = 2$  and  $A_{0i} = -(\theta, \underline{\alpha}_i^\vee)$ .

The marks  $a_{i}$  and comarks  $a_i^\vee (i = 1,\dots ,r)$  for the finite Lie algebra  $\mathfrak{g}$  are defined as

$$
\theta = \sum_ {i = 1} ^ {r} a _ {i} \alpha_ {i} = \sum_ {i = 1} ^ {r} a _ {i} ^ {\vee} \alpha_ {i} ^ {\vee}. \tag {C.6}
$$

For the affine Lie algebra  $\mathfrak{g}$ , the mark and comark of the extra simple root  $\alpha_0$  is defined to be  $1$ ,  $a_0 = a_0^\vee = 1$ . The dual Coxeter number is defined as  $h^\vee = 1 + \sum_{i=1}^{r} a_i^\vee = \sum_{i=0}^{r} a_i^\vee$ .

We will sometimes label an affine weight  $\lambda$  by its Dynkin labels,

$$
\lambda = \left[ \lambda_ {0}, \lambda_ {1}, \dots , \lambda_ {r} \right], \tag {C.7}
$$

where  $\lambda_{i} = (\lambda, \alpha_{i}^{\vee})$ . Note that the Dynkin labels have  $r + 1$  components, one less than the  $(\underline{\lambda}, k, n)$  notation. In other words, the Dynkin labels do not completely specify an affine weight, but up to an imaginary root  $n\delta$ . Note that the level  $k$  of a weight is related to the Dynkin labels by

$$
k = \sum_ {i = 0} ^ {r} a _ {i} ^ {\vee} \lambda_ {i}. \tag {C.8}
$$

Finally, the affine Weyl vector is defined as  $\rho = [1,1,\dots ,1]$

# Weyl Reflections

Let  $W$  be the Weyl group of the affine Lie algebra  $\mathfrak{g}$ . The Weyl reflection  $s_{\alpha} \in W$  generated by real affine root  $\alpha$  is given by

$$
s _ {\alpha} = \lambda - (\lambda , \alpha^ {\vee}) \alpha . \tag {C.9}
$$

The shifted Weyl reflection  $\circ$  generated by a real affine root  $\alpha$  is defined as

$$
s _ {\alpha} \circ \lambda = s _ {\alpha} (\lambda + \rho) - \rho = \lambda - (\lambda + \rho , \alpha^ {\vee}) \alpha . \tag {C.10}
$$

The Weyl group  $W$  is generated by  $s_{\alpha_i}$ , where  $\alpha_i$ 's are the simple roots, with the following relations

$$
s _ {i} ^ {2} = 1,
$$

$$
s _ {i} s _ {j} = s _ {j} s _ {i}, \quad \text {i f} A _ {i j} = 0, \tag {C.11}
$$

$$
(s _ {i} s _ {j}) ^ {m _ {i j}} = 1, \quad \text {i f} i \neq j,
$$

where  $m_{ij} = 2,3,4,6,\infty$  if the number of lines joining the  $i$ -th and  $j$ -th node is  $0,1,2,3,4$ .

# The case of  $\widehat{su(2)}$

Let us collect some properties of the affine Lie algebra  $\widehat{su(2)}$ . The affine Dynkin diagram has two nodes connected by 4 lines. The Cartan matrix is

$$
A _ {i j} = \left( \begin{array}{c c} 2 & - 2 \\ - 2 & 2 \end{array} \right). \tag {C.12}
$$

The highest root of  $su(2)$  is  $\theta = [2]$ . Hence the marks and comarks of  $\widehat{su(2)}$  are  $a_{i} = a_{i}^{\vee} = (1,1)$ . The level of an affine weight  $\lambda$  is then given by

$$
k = \lambda_ {0} + \lambda_ {1}. \tag {C.13}
$$

The dual Coxeter number is  $h^{\vee} = 2$ .

The affine Weyl group of  $\widehat{su(2)}$  is generated by  $s_0, s_1$  satisfying  $(s_i)^2 = 1$ . The affine Weyl group elements are

$$
W = \left\{s _ {0}, s _ {1}, s _ {1} s _ {0}, s _ {0} s _ {1}, \dots \right\}. \tag {C.14}
$$

# The case of  $\widehat{so(8)}$

Let us collect some properties of the affine Lie algebra  $\widehat{so(8)}$ . The Cartan matrix is

$$
A _ {i j} = \left( \begin{array}{r r r r} 2 & - 1 & 0 & 0 \\ - 1 & 2 & - 1 & - 1 \\ 0 & - 1 & 2 & 0 \\ 0 & - 1 & 0 & 2 \end{array} \right), \tag {C.15}
$$

where the central node in the affine Dynkin diagram is  $\alpha_{2}$ . The highest root of  $so(8)$  is  $\theta = [0,1,0,0]$ . Hence the marks and comarks of  $\widehat{so(8)}$  are  $a_{i} = a_{i}^{\vee} = (1,1,2,1,1)$ . The level of an affine weight  $\lambda$  is then given by

$$
k = \lambda_ {0} + \lambda_ {1} + 2 \lambda_ {2} + \lambda_ {3} + \lambda_ {4}. \tag {C.16}
$$

The dual Coxeter number is  $h^{\vee} = 6$ .

The affine Weyl group of  $\widehat{so(8)}$  is generated by  $s_0, s_1, s_2, s_3, s_4$  satisfying

$$
\left(s _ {i}\right) ^ {2} = 1, \qquad i = 0, 1, \dots , 4,
$$

$$
s _ {i} s _ {j} = s _ {j} s _ {i}, \quad i, j = 0, 1, 3, 4, \tag {C.17}
$$

$$
\left(s _ {2} s _ {i}\right) ^ {3} = 1, \qquad i = 0, 1, 3, 4.
$$

# C.2 Affine Characters and the Kazhdan-Lusztig Polynomials

In this subsection we present the formula for affine characters following [95] (see also [13]). We will assume

$$
k + h ^ {\vee} > 0, \tag {C.18}
$$

which is indeed the case for  $\widehat{su(2)}_{-\frac{4}{3}}$  and  $\widehat{so(8)}_{-2}$ .

To every weight  $\lambda$ , we define a subset  $\Delta_{+, \lambda}^{re}$  of the real positive roots of  $\mathfrak{g}$  to be

$$
\Delta_ {+, \lambda} ^ {r e} = \left\{\alpha \in \Delta_ {+} ^ {r e} \mid (\lambda , \alpha^ {\vee}) \in \mathbb {Z} \right\}, \tag {C.19}
$$

and let  $W_{\lambda}$  be the subgroup of the affine Weyl group  $W$  generated by  $s_{\alpha}$  with  $\alpha \in \Delta_{+, \lambda}^{re}$ . In the case when  $\lambda$  is integral (i.e. all the Dynkin labels are integers),  $W_{\lambda} = W$ .

Let  $\lambda$  be the highest affine weight of a module, then in the orbit

$$
W _ {\lambda} \circ \lambda \tag {C.20}
$$

there is exactly one element  $\Lambda$  such that the Dynkin labels of  $\Lambda +\rho$  are all non-negative,

$$
(\Lambda + \rho , \alpha_ {i} ^ {\vee}) \geq 0, \tag {C.21}
$$

where  $\alpha_{i}$ 's are the simple roots of the affine Lie algebra  $\mathfrak{g}$ .

Let  $\operatorname{ch} M(\mu)$  be the character of the Verma module with highest weight  $\mu$ ,

$$
\operatorname {c h} M (\mu) = \frac {e ^ {\mu}}{\prod_ {\alpha \in \Delta_ {+}} (1 - e ^ {- \alpha}) ^ {\operatorname {m u l t} (\alpha)}}, \tag {C.22}
$$

where  $\Delta_{+}$  is the set of positive roots for the affine Lie algebra  $\mathfrak{g}$  and  $\mathrm{mult}(\alpha)$  is the multiplicity of the root  $\alpha$ . Suppose  $\mu = (\underline{\mu}; k; n)$ , then  $e^{\mu}$  is understood as

$$
e ^ {\mu} = q ^ {- n} \eta_ {-} ^ {\mu}, \tag {C.23}
$$

where we have defined a compact notation  $\eta_{-}^{\mu} = \prod_{i=1}^{r} \eta_{i}^{\frac{\mu}{\mu_{i}}}$ . Here  $\eta_{i}$  are the fugacities and  $\underline{\mu}_{i}$ 's are the Dynkin labels of the weight  $\underline{\mu}$  of the finite Lie algebra  $\underline{\mathfrak{g}}$ . The character  $\operatorname{ch} L(\lambda)$  of the irreducible module with the highest weight  $\lambda = w \circ \Lambda$  is then given by

$$
\operatorname {ch} L (w \circ \Lambda) = \sum_ {\substack {w ^ {\prime} \in W _ {\Lambda} / W _ {\Lambda} ^ {0} \\ w ^ {\prime} \geq w}} m _ {w, w ^ {\prime}} \operatorname {ch} M \left(w ^ {\prime} \circ \Lambda\right). \tag{C.24}
$$

Here  $W_{\Lambda}^{0}$  is the subgroup of  $W_{\Lambda}$  that leaves  $\Lambda$  invariant.

To define the order  $>$  on the coset  $W_{\Lambda} / W_{\Lambda}^{0}$ , we first define the Bruhat order on the Weyl group  $W_{\Lambda}$ . An arbitrary element  $w$  in  $W_{\Lambda}$  can be written as  $w = s_{i_1}\dots s_{i_k}$ . An expression of minimal length is called reduced.[22] Let  $w, w' \in W_{\Lambda}$ , then we write

$$
w <   w ^ {\prime} \tag {C.25}
$$

if the reduced expression for  $w$  can be obtained by dropping simple reflections from a reduced expression for  $w'$ . The resulting relation  $w \leq w'$  is called the Bruhat order.

The order on the coset space  $W_{\Lambda} / W_{\Lambda}^{0}$  is then defined as

$$
w \leq w ^ {\prime}, \quad \text {w i t h} w, w ^ {\prime} \in W _ {\Lambda} / W _ {\Lambda} ^ {0} \quad \text {i f f} \quad \underline {{w}} \leq \underline {{w}} ^ {\prime}, \quad \text {w i t h} \underline {{w}}, \underline {{w}} ^ {\prime} \in W _ {\Lambda}, \quad \tag {C.26}
$$

where  $\underline{w}$  is the minimal representative of  $w$  in the coset, defined by  $\ell(\underline{w}s) > \ell(\underline{w})$  for all  $s \in W_{\Lambda}^{0}$ . Here  $\ell$  is the length of a reduced expression of a Weyl group element  $s$ . The determination of the multiplicities  $m_{w,w'}$  is the content of the Kazhdan-Lusztig conjecture.

# C.2.1 The Kazhdan-Lusztig Conjecture

The Kazhdan-Lusztig conjecture states that the multiplicities  $m_{w, w'}$  are given by the inverse Kazhdan-Lusztig polynomials  $\tilde{Q}_{w, w'}^I(q)$  for the coset  $W_{\Lambda} / W_{\Lambda}^{0}$  evaluated at  $q = 1$ ,[23]

$$
m _ {w, w ^ {\prime}} = \tilde {Q} _ {w, w ^ {\prime}} ^ {I} (1). \tag {C.27}
$$

The inverse Kazhdan-Lusztig polynomials  $\tilde{Q}_{w,w^{\prime}}^{I}(q)$  for the coset  $W_{\Lambda} / W_{\Lambda}^{0}$  are in turn related to those  $Q_{w,w^{\prime}}(q)$  of  $W_{\Lambda}$  by

$$
\tilde {Q} _ {x, y} ^ {I} (q) = \sum_ {z \in [ y ]} Q _ {\bar {x}, z} (q) (- 1) ^ {\ell (\bar {x})} (- 1) ^ {\ell (z)}, \tag {C.28}
$$

where  $\bar{z}$  and  $\underline{z}$  are the maximal and minimal representative of the coset  $[z]$  of  $z$ . For rational  $k$  with  $k + h^{\vee} > 0$ , which is indeed the case for  $\widehat{su(2)}_{-\frac{4}{3}}$  and  $\widehat{so(8)}_{-2}$ , the Kazhdan-Lusztig conjecture has been proven in [96] (see also [97-100] for earlier works).

The inverse Kazhdan-Lusztig polynomials for the Weyl group  $W_{\Lambda}$  are determined using a recurrence relation (see, for example, [95]). In the case of  $\widehat{so(8)}_{-2}$ , we will use the C code Coxeter [101] to compute the inverse Kazhdan-Lusztig polynomials  $Q_{w,w^{\prime}}(q)$ .

# Weyl-Kac Character Formula

Consider the special case  $\lambda$  itself is a dominant weight, i.e.  $(\lambda, \alpha_i^\vee) \geq 0$  for all the simple roots  $\alpha_i$ . In this case  $\Lambda = \lambda$  and  $w = e$  is the identity element in the affine Weyl group  $W$ . Since  $\lambda$  is integral,  $W_\Lambda = W_\lambda = W$ . Further, the subgroup  $W_\Lambda^0$  is trivial and the coset  $W / W_\lambda^0$  is the full Weyl group  $W$ . In this case we have  $\tilde{Q}_{e,y} = (-1)^{\ell(y)}$ . Then the character (C.24) for the module with highest weight  $\lambda$  reduces to the familiar Weyl-Kac formula.

# C.3 Affine Characters of  $\widehat{su(2)}_{-\frac{4}{3}}$

In this subsection we apply the above formalism to compute the characters of admissible representations (6.29) of  $\widehat{su(2)}_{-\frac{4}{3}}$ . The vacuum character of  $\widehat{su(2)}_{-\frac{4}{3}}$  has been previously computed in [7] (see also [93]) so we will not repeat it here. We will explicitly compute the character  $\chi_1$  for the admissible representation with highest weight  $\Phi_1 = [-\frac{2}{3}, -\frac{2}{3}]$ . The character  $\chi_2$  for the other admissible representation  $\Phi_2 = [0, -\frac{4}{3}]$  is completely analogous.

The real positive roots associated to the highest weight  $\lambda = \left[-\frac{2}{3}, -\frac{2}{3}\right]$  is

$$
\Delta_ {+, \lambda} ^ {r e} = \left\{\left(\alpha_ {1}; 0; 3 m + 1\right) \mid m \geq 0 \right\} \cup \left\{\left(- \alpha_ {1}; 0; 3 m + 2\right) \mid m \geq 0 \right\}. \tag {C.29}
$$

One can easily check that  $\langle \lambda + \rho, \alpha^{\vee} \rangle > 0$  for all  $\alpha \in \Delta_{+, \lambda}^{re}$ , hence there is no need to perform a further Weyl reflection  $w$ . In the notations of the previous section, we have  $\Lambda = \lambda$  and  $w = e$ . Furthermore, the inverse Kazhdan-Lusztig polynomials are just signs in this case,  $\tilde{Q}_{e,y} = (-1)^{\ell(y)}$ . The character formula (C.24) reduces to

$$
\operatorname {c h} L (\lambda) = \sum_ {w ^ {\prime} \in W _ {\lambda}} (- 1) ^ {\ell \left(w ^ {\prime}\right)} \operatorname {c h} M \left(w ^ {\prime} \circ \lambda\right). \tag {C.30}
$$

This special case of the character formula is known as the Kac-Wakimoto formula [92].

Let us take a closer look into  $W_{\lambda}$ , which is the subgroup of the affine Weyl group that is generated by roots in  $\Delta_{+, \lambda}^{re}$ . The simple roots  $\tilde{\alpha}_{0}, \tilde{\alpha}_{1}$  of  $\Delta_{+, \lambda}^{re}$  are

$$
\tilde {\alpha} _ {0} = 2 \alpha_ {0} + \alpha_ {1}, \quad \tilde {\alpha} _ {1} = \alpha_ {0} + 2 \alpha_ {1}. \tag {C.31}
$$

$W_{\lambda}$  is then generated by the corresponding Weyl reflections  $\tilde{s}_0, \tilde{s}_1$

$$
\tilde {s} _ {0} = s _ {0} s _ {1} s _ {0}, \quad \tilde {s} _ {1} = s _ {1} s _ {0} s _ {1}, \tag {C.32}
$$

where  $s_0$  and  $s_1$  are the Weyl reflections generated by the simple roots  $\alpha_0$  and  $\alpha_1$  of  $\widehat{su(2)}$ , respectively. The elements in  $W_{\lambda}$  include  $e$ ,  $(\tilde{s}_0\tilde{s}_1)^{n - 1}\tilde{s}_0$ ,  $(\tilde{s}_1\tilde{s}_0)^{n - 1}\tilde{s}_1$ ,  $(\tilde{s}_1\tilde{s}_0)^n$ , and  $(\tilde{s}_0\tilde{s}_1)^n$  with  $n \geq 1$ . After working out the shifted Weyl reflection on the highest weight  $w \circ \lambda$ , we obtain the character  $\chi_1(q,z)$  for  $\lambda = [-\frac{2}{3}, -\frac{2}{3}]$ ,

$$
\chi_ {1} (q, z) = \frac {1 + \sum_ {n = 1} ^ {\infty} (- 1) ^ {n} \left(z ^ {- 2 n} q ^ {\frac {n}{2} (3 n - 1)} + z ^ {2 n} q ^ {\frac {n}{2} (3 n + 1)}\right)}{(1 - z ^ {- 2}) \prod_ {n = 1} ^ {\infty} (1 - q ^ {n}) (1 - z ^ {2} q ^ {n}) (1 - z ^ {- 2} q ^ {n})}. \tag {C.33}
$$

Similarly the character  $\chi_2(q,z)$  for  $[0, -\frac{4}{3}]$  is

$$
\chi_ {2} (q, z) = \frac {1 + \sum_ {n = 1} ^ {\infty} (- 1) ^ {n} \left(z ^ {2 n} q ^ {\frac {n}{2} (3 n - 1)} + z ^ {- 2 n} q ^ {\frac {n}{2} (3 n + 1)}\right)}{\left(1 - z ^ {- 2}\right) \prod_ {n = 1} ^ {\infty} \left(1 - q ^ {n}\right) \left(1 - z ^ {2} q ^ {n}\right) \left(1 - z ^ {- 2} q ^ {n}\right)}. \tag {C.34}
$$

# C.4 Affine Characters of  $\widehat{so(8)}_{-2}$

In this subsection we will record the answers of the affine characters for several highest weight modules in  $\widehat{so(8)}_{-2}$ . The line defect indices of the  $SU(2)$  with  $N_f = 4$  flavors theory turn out to be linear combinations of these affine characters. The computation of these

characters are done with the help of Mathematica and the C code Coxeter [101] which computes the Kazhdan-Lusztig polynomials efficiently. The vacuum character of  $\widehat{so(8)}_{-2}$  has been previously computed in [13].

Recall that to apply the Kazhdan-Lusztig formula for a given module with highest affine weight  $\lambda$ , we need to find an element  $w$  of the affine Weyl group  $W$  such that  $\Lambda = w^{-1} \circ \lambda$  has all affine Dynkin labels no smaller than -1. We list the affine Dynkin labels, dimensions<sup>26</sup>,  $\Lambda$ ,  $w$  of the highest weight modules that we will compute their characters below:

<table><tr><td>Affine Dynkin Label</td><td>Dimension</td><td>Λ</td><td>w</td></tr><tr><td>[-2,0,0,0,0]</td><td>1</td><td>[0,0,-1,0,0]</td><td>s0</td></tr><tr><td>[-3,1,0,0,0]</td><td>8v</td><td>[0,0,0,-1,-1]</td><td>s0s2</td></tr><tr><td>[-4,2,0,0,0]</td><td>35v</td><td>[0,0,-1,0,0]</td><td>s0s2s3s4</td></tr><tr><td>[-5,3,0,0,0]</td><td>112v</td><td>[-1,-1,0,0,0]</td><td>s0s2s3s4s2</td></tr><tr><td>[-6,4,0,0,0]</td><td>294v</td><td>[0,0,-1,0,0]</td><td>s0s2s3s4s2s1s0</td></tr><tr><td>[-7,5,0,0,0]</td><td>672v</td><td>[0,0,0,-1,-1]</td><td>s0s2s3s4s2s1s0s2</td></tr><tr><td>[-8,6,0,0,0]</td><td>1386v</td><td>[0,0,-1,0,0]</td><td>s0s2s3s4s2s1s0s2s3s4</td></tr><tr><td>[-9,7,0,0,0]</td><td>2640v</td><td>[-1,-1,0,0,0]</td><td>s0s2s3s4s2s1s0s2s3s4s2</td></tr><tr><td>[-10,8,0,0,0]</td><td>4719v</td><td>[0,0,-1,0,0]</td><td>s0s2s3s4s2s1s0s2s3s4s2s1s0</td></tr><tr><td>[-11,9,0,0,0]</td><td>8008v</td><td>[0,0,0,-1,-1]</td><td>s0s2s3s4s2s1s0s2s3s4s2s1s0s2</td></tr><tr><td>[-12,10,0,0,0]</td><td>13013v</td><td>[0,0,-1,0,0]</td><td>s0s2s3s4s2s1s0s2s3s4s2s1s0s2s3s4</td></tr></table>

We will denote the affine character of a highest weight module with affine Dynkin labels  $[a_0, a_1, a_2, a_3, a_4]$  by  $\chi_{[a_0, a_1, a_2, a_3, a_4]}$ . We record several affine characters of  $\widehat{so(8)}_{-2}$  with flavor fugacities set to be 1 below:

$$
\begin{array}{l} \chi_ {[ - 2, 0, 0, 0, 0 ]} = 1 + 2 8 q + 3 2 9 q ^ {2} + 2 6 3 2 q ^ {3} + 1 6 3 8 0 q ^ {4} + 8 5 7 6 4 q ^ {5} + 3 9 3 5 8 9 q ^ {6} \\ + 1 6 2 8 5 4 8 q ^ {7} + 6 1 9 0 5 2 7 q ^ {8} + 2 1 9 2 1 9 0 0 q ^ {9} + 7 3 0 7 0 2 9 1 q ^ {1 0} + 2 3 1 1 1 8 3 8 4 q ^ {1 1} \\ + 6 9 8 1 2 8 3 8 9 q ^ {1 2} + 2 0 2 4 4 3 3 4 6 0 q ^ {1 3} + 5 6 5 9 7 3 0 0 7 5 q ^ {1 4} + \mathcal {O} (q ^ {1 5}), \\ \end{array}
$$

$$
\begin{array}{l} \chi_ {[ - 3, 1, 0, 0, 0 ]} = 8 + 1 6 8 q + 1 9 0 4 q ^ {2} + 1 5 5 1 2 q ^ {3} + 1 0 1 6 9 6 q ^ {4} + 5 6 9 0 7 2 q ^ {5} + 2 8 1 7 6 4 0 q ^ {6} \\ + 1 2 6 4 2 0 1 6 q ^ {7} + 5 2 2 7 5 2 1 6 q ^ {8} + 2 0 1 7 1 6 0 3 2 q ^ {9} + 7 3 3 3 2 6 4 4 0 q ^ {1 0} \\ + 2 5 3 0 6 0 9 5 3 6 q ^ {1 1} + \mathcal {O} (q ^ {1 2}), \\ \end{array}
$$

$$
\begin{array}{l} \chi_ {[ - 4, 2, 0, 0, 0 ]} = 3 5 + 6 3 0 q + 6 5 2 4 q ^ {2} + 4 9 4 9 0 q ^ {3} + 3 0 5 7 9 5 q ^ {4} + 1 6 2 5 0 6 0 q ^ {5} + 7 6 8 3 5 5 0 q ^ {6} \\ + 3 3 0 5 8 9 5 6 q ^ {7} + 1 3 1 5 2 9 9 4 4 q ^ {8} + 4 8 9 7 0 0 5 1 2 q ^ {9} + 1 7 2 1 7 5 4 3 9 1 q ^ {1 0} + 5 7 5 7 9 3 7 5 2 8 q ^ {1 1} \\ + 1 8 4 2 1 7 0 6 9 2 4 q ^ {1 2} + 5 6 6 5 2 3 2 2 6 3 6 q ^ {1 3} + 1 6 8 1 2 8 8 6 3 1 9 6 q ^ {1 4} + \mathcal {O} (q ^ {1 5}), \\ \end{array}
$$

$$
\begin{array}{l} \chi_ {[ - 5, 3, 0, 0, 0 ]} = 1 1 2 + 1 8 4 0 q + 1 7 9 2 0 q ^ {2} + 1 3 0 4 8 0 q ^ {3} + 7 8 3 4 4 0 q ^ {4} + 4 0 8 0 2 7 2 q ^ {5} \\ + 1 9 0 2 1 2 9 6 q ^ {6} + 8 1 0 4 7 5 6 8 q ^ {7} + 3 2 0 3 9 0 9 4 4 q ^ {8} + 1 1 8 8 1 7 7 3 1 2 q ^ {9} \\ + 4 1 6 9 2 4 9 7 2 8 q ^ {1 0} + 1 3 9 3 6 1 9 8 3 0 4 q ^ {1 1} + \mathcal {O} (q ^ {1 2}), \\ \end{array}
$$

$$
\begin{array}{l} \chi_ {[ - 6, 4, 0, 0, 0 ]} = 2 9 4 + 4 5 5 7 q + 4 2 5 1 6 q ^ {2} + 2 9 9 1 0 3 q ^ {3} + 1 7 4 4 1 0 6 q ^ {4} + 8 8 5 2 9 6 3 q ^ {5} \\ + 4 0 3 2 6 1 3 2 q ^ {6} + 1 6 8 2 2 4 5 2 5 q ^ {7} + 6 5 2 0 8 9 8 7 2 q ^ {8} + 2 3 7 4 3 1 6 2 2 8 q ^ {9} + 8 1 8 8 5 3 2 2 9 6 q ^ {1 0} \\ + 2 6 9 2 6 3 0 1 2 0 6 q ^ {1 1} + 8 4 8 7 2 8 6 0 4 0 8 q ^ {1 2} + \mathcal {O} (q ^ {1 3}), \\ \end{array}
$$

$$
\begin{array}{l} \chi_ {[ - 7, 5, 0, 0, 0 ]} = 6 7 2 + 1 0 0 1 6 q + 9 0 6 0 8 q ^ {2} + 6 2 1 2 6 4 q ^ {3} + 3 5 4 7 0 4 0 q ^ {4} + 1 7 6 9 0 9 6 0 q ^ {5} \\ + 7 9 4 1 0 4 6 4 q ^ {6} + 3 2 7 2 1 2 7 0 4 q ^ {7} + 1 2 5 5 2 9 9 5 6 8 q ^ {8} + 4 5 3 0 9 1 0 7 2 0 q ^ {9} + \mathcal {O} (q ^ {1 0}), \\ \end{array}
$$

$$
\begin{array}{l} \chi_ {[ - 8, 6, 0, 0, 0 ]} = 1 3 8 6 + 2 0 0 9 7 q + 1 7 7 7 1 6 q ^ {2} + 1 1 9 4 9 6 3 q ^ {3} + 6 7 0 7 2 0 4 q ^ {4} + 3 2 9 4 6 0 5 3 q ^ {5} \\ + 1 4 5 8 5 3 4 9 8 q ^ {6} + 5 9 3 3 8 3 0 2 8 q ^ {7} + 2 2 4 9 6 0 9 6 5 6 q ^ {8} + 8 0 3 0 0 8 4 5 9 4 q ^ {9} \\ + 2 7 2 0 4 2 0 9 1 1 6 q ^ {1 0} + \mathcal {O} (q ^ {1 1}), \\ \end{array}
$$

$$
\begin{array}{l} \chi_ {[ - 9, 7, 0, 0, 0 ]} = 2 6 4 0 + 3 7 5 2 0 q + 3 2 6 1 4 4 q ^ {2} + 2 1 6 0 1 4 4 q ^ {3} + 1 1 9 6 2 8 3 2 q ^ {4} + 5 8 0 6 3 3 7 6 q ^ {5} \\ + 2 5 4 3 1 8 2 8 8 q ^ {6} + 1 0 2 4 8 2 1 1 3 6 q ^ {7} + \mathcal {O} (q ^ {8}), \\ \end{array}
$$

$$
\begin{array}{l} \chi_ {[ - 1 0, 8, 0, 0, 0 ]} = 4 7 1 9 + 6 6 0 6 6 q + 5 6 6 7 4 8 q ^ {2} + 3 7 0 9 5 2 4 q ^ {3} + 2 0 3 2 4 1 9 2 q ^ {4} + 9 7 6 8 5 6 7 2 q ^ {5} \\ + 4 2 4 0 2 1 3 3 2 q ^ {6} + 1 6 9 4 4 0 5 9 4 8 q ^ {7} + \mathcal {O} (q ^ {8}), \\ \end{array}
$$

$$
\chi_ {[ - 1 1, 9, 0, 0, 0 ]} = 8 0 0 8 + 1 1 0 8 2 4 q + 9 4 0 9 1 2 q ^ {2} + \mathcal {O} (q ^ {3}),
$$

$$
\chi_ {[ - 1 2, 1 0, 0, 0, 0 ]} = 1 3 0 1 3 + 1 7 8 4 6 4 q + \mathcal {O} (q ^ {2}).
$$

# References

[1] J. Kinney, J. M. Maldacena, S. Minwalla, and S. Raju, “An Index for 4 dimensional super conformal theories,” Commun. Math. Phys. 275 (2007) 209–254, hep-th/0510251.  
[2] A. Gadde, L. Rastelli, S. S. Razamat, and W. Yan, “The 4d Superconformal Index from q-deformed 2d Yang-Mills,” Phys.Rev.Lett. 106 (2011) 241602, 1104.3850.  
[3] A. Gadde, L. Rastelli, S. S. Razamat, and W. Yan, “Gauge Theories and Macdonald Polynomials,” Commun. Math. Phys. 319 (2013) 147–193, 1110.3740.  
[4] A. Gadde, L. Rastelli, S. S. Razamat, and W. Yan, “The Superconformal Index of the  $E_6$  SCFT,” JHEP 08 (2010) 107, 1003.4244.  
[5] D. Gaiotto, L. Rastelli, and S. S. Razamat, “Bootstrapping the superconformal index with surface defects,” JHEP 01 (2013) 022, 1207.3577.

[6] L. Rastelli and S. S. Razamat, “The Superconformal Index of Theories of Class  $S$ ," in New Dualities of Supersymmetric Gauge Theories, J. Teschner, ed., pp. 261-305. 2016. 1412.7131.  
[7] M. Buican and T. Nishinaka, “On the superconformal index of Argyres-Douglas theories,” J. Phys. A49 (2016), no. 1, 015401, 1505.05884.  
[8] C. Córdova and S.-H. Shao, “Schur Indices, BPS Particles, and Argyres-Douglas Theories,” JHEP 01 (2016) 040, 1506.00265.  
[9] A. Gadde, E. Pomoni, L. Rastelli, and S. S. Razamat, “S-duality and 2d Topological QFT,” JHEP 1003 (2010) 032, 0910.2225.  
[10] T. Kawano and N. Matsumiya, “5D SYM on 3D Sphere and 2D YM,” Phys.Lett. B716 (2012) 450–453, 1206.5966.  
[11] Y. Fukuda, T. Kawano, and N. Matsumiya, “5D SYM and 2D q-Deformed YM,” Nucl.Phys. B869 (2013) 493–522, 1210.2855.  
[12] J. Song, “Superconformal indices of generalized Argyres-Douglas theories from 2d TQFT,” JHEP 02 (2016) 045, 1509.06730.  
[13] C. Beem, M. Lemos, P. Liendo, W. Peelaers, L. Rastelli, et al., “Infinite Chiral Symmetry in Four Dimensions,” Commun. Math. Phys. 336 (2015), no. 3, 1359–1433, 1312.5344.  
[14] C. Beem, W. Peelaers, L. Rastelli, and B. C. van Rees, “Chiral algebras of class S,” JHEP 1505 (2015) 020, 1408.6522.  
[15] M. Lemos and W. Peelaers, “Chiral Algebras for Trinion Theories,” JHEP 1502 (2015) 113, 1411.3252.  
[16] M. Buican and T. Nishinaka, “Conformal Manifolds in Four Dimensions and Chiral Algebras,” 1603.00887.  
[17] D. Xie, W. Yan, and S.-T. Yau, “Chiral algebra of Argyres-Douglas theory from M5 brane,” 1604.02155.  
[18] S. Cecotti, J. Song, C. Vafa, and W. Yan, “Superconformal Index, BPS Monodromy and Chiral Algebras,” 1511.01516.  
[19] T. Arakawa and A. Moreau, “Joseph ideals and lisse minimal W-algebras,” 1506.00710.

[20] T. Arakawa, V. Futorny, and L. E. Ramirez, “Weight representations of admissible affine vertex algebras,” 1605.07580.  
[21] C. Beem and L. Rastelli, “Vertex operator algebras, Higgs branches, and modular differential equations,” to appear (2016).  
[22] A. Kapustin, “Wilson-’t Hooft operators in four-dimensional gauge theories and S-duality,” Phys. Rev. D74 (2006) 025005, hep-th/0501015.  
[23] S. Gukov and E. Witten, “Gauge Theory, Ramification, And The Geometric Langlands Program,” hep-th/0612073.  
[24] A. Kapustin, “Holomorphic reduction of  $\mathrm{N} = 2$  gauge theories, Wilson-'t Hooft operators, and S-duality,” hep-th/0612119.  
[25] A. Kapustin and N. Saulina, “The Algebra of Wilson-’t Hooft operators,” Nucl. Phys. B814 (2009) 327–365, 0710.2097.  
[26] N. Drukker, D. R. Morrison, and T. Okuda, “Loop operators and S-duality from curves on Riemann surfaces,” JHEP 09 (2009) 031, 0907.2593.  
[27] D. Gaiotto, G. W. Moore, and A. Neitzke, “Framed BPS States,” Adv. Theor. Math. Phys. 17 (2013), no. 2, 241–397, 1006.0146.  
[28] D. Xie, “Higher laminations, webs and  $\mathrm{N} = 2$  line operators,” 1304.2390.  
[29] O. Aharony, N. Seiberg, and Y. Tachikawa, “Reading between the lines of four-dimensional gauge theories,” JHEP 08 (2013) 115, 1305.0318.  
[30] D. Xie, “Aspects of line operators of class S theories,” 1312.3371.  
[31] I. Coman, M. Gabella, and J. Teschner, “Line operators in theories of class  $\mathcal{S}$ , quantized moduli space of flat connections, and Toda field theory,” JHEP 10 (2015) 143, 1505.05898.  
[32] O. DeWolfe, D. Z. Freedman, and H. Ooguri, “Holography and defect conformal field theories,” Phys. Rev. D66 (2002) 025009, hep-th/0111135.  
[33] D. Gaiotto and E. Witten, “Supersymmetric Boundary Conditions in N=4 Super Yang-Mills Theory,” J. Statist. Phys. 135 (2009) 789–855, 0804.2902.  
[34] D. Gaiotto and E. Witten, “S-Duality of Boundary Conditions In N=4 Super Yang-Mills Theory,” Adv. Theor. Math. Phys. 13 (2009), no. 3, 721-896, 0807.3720.  
[35] S. Cecotti, C. Cordova, and C. Vafa, “Braids, Walls, and Mirrors,” 1110.2115.

[36] T. Dimofte, D. Gaiotto, and S. Gukov, “3-Manifolds and 3d Indices,” Adv. Theor. Math. Phys. **17** (2013), no. 5, 975–1076, 1112.5179.  
[37] T. Dimofte and D. Gaiotto, “An E7 Surprise,” JHEP 10 (2012) 129, 1209.1404.  
[38] T. Dimofte, D. Gaiotto, and R. van der Veen, “RG Domain Walls and Hybrid Triangulations,” Adv. Theor. Math. Phys. 19 (2015) 137–276, 1304.6721.  
[39] Y. Ito, T. Okuda, and M. Taki, "Line operators on  $S^1 \times R^3$  and quantization of the Hitchin moduli space," JHEP 04 (2012) 010, 1111.4221. [Erratum: JHEP03,085(2016)].  
[40] D. Gang, E. Koh, and K. Lee, “Line Operator Index on  $S^1 \times S^3$ ,” JHEP 05 (2012) 007, 1201.5539.  
[41] C.-K. Chang, H.-Y. Chen, D. Jain, and N. Lee, “Connecting Localization and Wall-Crossing via D-Branes,” 1512.02645.  
[42] V. Pestun, “Localization of gauge theory on a four-sphere and supersymmetric Wilson loops,” Commun. Math. Phys. 313 (2012) 71–129, 0712.2824.  
[43] K. Hosomichi, S. Lee, and J. Park, “AGT on the S-duality Wall,” JHEP 12 (2010) 079, 1009.0340.  
[44] N. Drukker, J. Gomis, T. Okuda, and J. Teschner, “Gauge Theory Loop Operators and Liouville Theory,” JHEP 02 (2010) 057, 0909.1105.  
[45] N. Drukker, D. Gaiotto, and J. Gomis, “The Virtue of Defects in 4D Gauge Theories and 2D CFTs,” JHEP 06 (2011) 025, 1003.1112.  
[46] J. Gomis, T. Okuda, and V. Pestun, “Exact Results for ’t Hooft Loops in Gauge Theories on  $S^4$ ,” JHEP 05 (2012) 141, 1105.2568.  
[47] N. Hama and K. Hosomichi, “Seiberg-Witten Theories on Ellipsoids,” JHEP 09 (2012) 033, 1206.6359. [Addendum: JHEP10,051(2012)].  
[48] N. Seiberg and E. Witten, "Electric - magnetic duality, monopole condensation, and confinement in  $\mathrm{N} = 2$  supersymmetric Yang-Mills theory," Nucl. Phys. B426 (1994) 19-52, hep-th/9407087. [Erratum: Nucl. Phys.B430,485(1994)].  
[49] N. Seiberg and E. Witten, “Monopoles, duality and chiral symmetry breaking in  $\mathbf{N} = 2$  supersymmetric QCD,” Nucl. Phys. B431 (1994) 484–550, hep-th/9408099.  
[50] T. Dumitrescu, G. Festuccia, and M. Del Zotto, work in progress.

[51] S. Cecotti, A. Neitzke, and C. Vafa, “R-Twisting and 4d/2d Correspondences,” 1006.3435.  
[52] A. Iqbal and C. Vafa, “BPS Degeneracies and Superconformal Index in Diverse Dimensions,” Phys. Rev. D90 (2014), no. 10, 105031, 1210.3605.  
[53] S. Cecotti and C. Vafa, “On classification of  $\mathbf{N} = 2$  supersymmetric theories,” Commun. Math. Phys. 158 (1993) 569-644, hep-th/9211097.  
[54] D. Gaiotto, G. W. Moore, and E. Witten, “An Introduction To The Web-Based Formalism,” 1506.04086.  
[55] D. Gaiotto, G. W. Moore, and E. Witten, “Algebra of the Infrared: String Field Theoretic Structures in Massive  $\mathcal{N} = (2,2)$  Field Theory In Two Dimensions,” 1506.04087.  
[56] C. Cordova, D. Gaiotto, and S.-H. Shao, “Surface Defect Indices and 2d-4d BPS States,” 1703.02525.  
[57] C. Cordova, D. Gaiotto, and S.-H. Shao, “Surface Defects and Chiral Algebras,” 1704.01955.  
[58] C. Beem, W. Peelaers, and L. Rastelli, work in progress.  
[59] M. Kontsevich and Y. Soibelman, “Stability structures, motivic Donaldson-Thomas invariants and cluster transformations,” 0811.2435.  
[60] T. Dimofte, S. Gukov, and Y. Soibelman, “Quantum Wall Crossing in N=2 Gauge Theories,” Lett. Math. Phys. 95 (2011) 1–25, 0912.1346.  
[61] C. Papageorgakis, A. Pini, and D. Rodriguez-Gomez, “The NS limit of the 5D Superconformal Index,” 1602.02647.  
[62] S. Lee and P. Yi, “Framed BPS States, Moduli Dynamics, and Wall-Crossing,” JHEP 04 (2011) 098, 1102.1729.  
[63] W.-y. Chuang, D.-E. Diaconescu, J. Manschot, G. W. Moore, and Y. Soibelman, “Geometric engineering of (framed) BPS states,” Adv. Theor. Math. Phys. 18 (2014), no. 5, 1063–1231, 1301.3065.  
[64] M. Cirafici, “Line defects and (framed) BPS quivers,” JHEP 11 (2013) 141, 1307.7134.  
[65] C. Córdova and A. Neitzke, “Line Defects, Tropicalization, and Multi-Centered Quiver Quantum Mechanics,” JHEP 09 (2014) 099, 1308.6829.

[66] G. W. Moore, A. B. Royston, and D. V. d. Bleeken, “ $L^2$ -Kernels Of Dirac-Type Operators On Monopole Moduli Spaces,” 1512.08923.  
[67] G. W. Moore, A. B. Royston, and D. V. d. Bleeken, “Semiclassical framed BPS states,” 1512.08924.  
[68] M. Gabella, “Quantum Holonomies from Spectral Networks and Framed BPS States,” 1603.05258.  
[69] P. C. Argyres and M. R. Douglas, “New phenomena in SU(3) supersymmetric gauge theory,” Nucl. Phys. B448 (1995) 93–126, hep-th/9505062.  
[70] P. C. Argyres, M. R. Plesser, N. Seiberg, and E. Witten, “New  $\mathbf{N} = 2$  superconformal field theories in four-dimensions,” Nucl. Phys. B461 (1996) 71-84, hep-th/9511154.  
[71] T. Eguchi, K. Hori, K. Ito, and S.-K. Yang, “Study of N=2 superconformal field theories in four-dimensions,” Nucl. Phys. B471 (1996) 430-444, hep-th/9603002.  
[72] G. Bonelli, K. Maruyoshi, and A. Tanzini, “Wild Quiver Gauge Theories,” JHEP 1202 (2012) 031, 1112.1691.  
[73] D. Xie, “General Argyres-Douglas Theory,” JHEP 1301 (2013) 100, 1204.2270.  
[74] D. Xie, “Network, cluster coordinates and  $\mathbf{N} = 2$  theory II: Irregular singularity,” 1207.6112.  
[75] D. Xie and P. Zhao, “Central charges and RG flow of strongly-coupled N=2 theory,” JHEP 1303 (2013) 006, 1301.0210.  
[76] A. D. Shapere and C. Vafa, “BPS structure of Argyres-Douglas superconformal theories,” hep-th/9910182.  
[77] S. Cecotti and C. Vafa, “Classification of complete  $\mathrm{N} = 2$  supersymmetric theories in 4 dimensions,” Surveys in differential geometry 18 (2013) 1103.5832.  
[78] M. Alim, S. Cecotti, C. Córdova, S. Espahbodi, A. Rastogi, and C. Vafa, “BPS Quivers and Spectra of Complete N=2 Quantum Field Theories,” Commun. Math. Phys. 323 (2013) 1185–1227, 1109.4941.  
[79] M. Alim, S. Cecotti, C. Córdova, S. Espahbodi, A. Rastogi, et al., “ $\mathcal{N} = 2$  quantum field theories and their BPS quivers,” Adv. Theor. Math. Phys. 18 (2014) 27–127, 1112.3984.  
[80] M. Buican and T. Nishinaka, “Argyres-Douglas Theories, the Macdonald Index, and an RG Inequality,” JHEP 02 (2016) 159, 1509.05402.

[81] K. Maruyoshi and J. Song, “The Full Superconformal Index of the Argyres-Douglas Theory,” 1606.05632.  
[82] D. Gaiotto, G. W. Moore, and A. Neitzke, “Wall-crossing, Hitchin Systems, and the WKB Approximation,” 0907.3987.  
[83] A. Kapustin and E. Witten, "Electric-Magnetic Duality And The Geometric Langlands Program," Commun. Num. Theor. Phys. 1 (2007) 1-236, hep-th/0604151.  
[84] L. F. Alday, D. Gaiotto, S. Gukov, Y. Tachikawa, and H. Verlinde, “Loop and surface operators in  $\mathbf{N} = 2$  gauge theory and Liouville modular geometry,” JHEP 1001 (2010) 113, 0909.0945.  
[85] A. Braverman, M. Finkelberg, and H. Nakajima, “Coulomb branches of  $3d\mathcal{N} = 4$  quiver gauge theories and slices in the affine Grassmannian (with appendices by Alexander Braverman, Michael Finkelberg, Joel Kamnitzer, Ryosuke Kodera, Hiraku Nakajima, Ben Webster, and Alex Weekes),” 1604.03625.  
[86] M. Bullimore, T. Dimofte, and D. Gaiotto, “The Coulomb Branch of 3d  $\mathcal{N} = 4$  Theories,” 1503.04817.  
[87] D. Galakhov, P. Longhi, T. Mainiero, G. W. Moore, and A. Neitzke, “Wild Wall Crossing and BPS Giants,” JHEP 11 (2013) 046, 1305.5454.  
[88] D. Gaiotto, “Domain Walls for Two-Dimensional Renormalization Group Flows,” JHEP 12 (2012) 103, 1201.0767.  
[89] S. S. Razamat, “On a modular property of  $\mathrm{N} = 2$  superconformal theories in four dimensions,” JHEP 10 (2012) 191, 1208.5056.  
[90] E. P. Verlinde, “Fusion Rules and Modular Transformations in 2D Conformal Field Theory,” Nucl. Phys. B300 (1988) 360–376.  
[91] C. Beem, M. Lemos, P. Liendo, L. Rastelli, and B. C. van Rees, “The  $\mathcal{N} = 2$  superconformal bootstrap,” JHEP 03 (2016) 183, 1412.7541.  
[92] V. G. Kac and M. Wakimoto, “Modular invariant representations of infinite-dimensional Lie algebras and superalgebras,” Proceedings of the National Academy of Sciences 85 (1988), no. 14, 4956–4960.  
[93] P. D. Francesco, P. Mathieu, and D. Senechal, Conformal Field Theory. Graduate Texts in Contemporary Physics. Springer, 1997.

[94] D. Kazhdan and G. Lusztig, “Representations of Coxeter Groups and Hecke Algebras,” Inventiones Mathematicae 53 (1979), no. 2, 165–184.  
[95] K. De Vos and P. Van Driel, “The Kazhdan-Lusztig conjecture for W algebras,” J. Math. Phys. 37 (1996) 3587, hep-th/9508020.  
[96] M. Kashiwara and T. Tanisaki, “Kazhdan-Lusztig Conjecture for Affine Lie Algebras with Negative Level II: Nonintegral Case,” Duke Math. J. 84 (09, 1996) 771-813.  
[97] M. Kashiwara, “Kazhdan-Lusztig Conjecture for A Symmetrizable Kac-Moody Lie Algebra,” The Grothendieck Festschrift: A Collection of Articles Written in Honor of the 60th Birthday of Alexander Grothendieck (1990) 407–433.  
[98] L. Casian, “Kazhdan-Lusztig Multiplicity Formulas for Kac-Moody Algebras,” Comptes Rendus de Lacademie des Sciences Serie I-Mathematique 310 (1990), no. 6, 333-337.  
[99] M. Kashiwara and T. Tanisaki, “Kazhdan-Lusztig conjecture for symmetrizable Kac-Moody Lie algebras II,” Operator algebras, unitary representations, enveloping algebras, and invariant theory 2 (1979) 159–195.  
[100] M. Kashiwara and T. Tanisaki, “Kazhdan-Lusztig Conjecture for Affine Lie Algebras with Negative Level,” Duke Math. J. 77 (01, 1995) 21–62.  
[101] F. Du Cloux, “Computing Kazhdan-Lusztig Polynomials for Arbitrary Coxeter Groups,” Experimental Mathematics 11 (2002), no. 3, 371-381.