# TASK

Consider affine Lie algebra $(\widehat{D}_6)_{-4}$. Consider $\widehat \lambda = -4 \widehat \omega_0$. Compute the character using **Kazhan-Lusztig method**.


# FACT

- the character `ch` start with $1$
- The $q$ term coefficient must be the character of the adjoint representation of $D_6$. There should be $66$ states at this level.

# CONSTRAINTS

- use `pyw` to compute. Key class to use
  - `KazhdanLusztigCharacter`
  - `KazhdanLusztigPolynomials`
- **FORBIDDEN**: you are not allowed to use `Algebra.py` (deprecated implementation)
- **FORBIDDEN**: you are not allowed to use `IntegrableModuleCharacter`