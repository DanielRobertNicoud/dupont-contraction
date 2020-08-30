# Dupont contraction

This package provides tools to work and do explicit computations on Sullivan and Dupont forms, as well as calculating the action of the various maps involved in the Dupont contraction and the transferred structure from the Sullivan algebra (a commutative algebra) to the Dupont algebra (which receives the structure of a commutative algebra up to homotopy).

# Classes

## DupontForm

Implements the Dupont forms, given by the cellular cochains on the simplices. The Dupont forms in simplicial degree n form a chain complex with basis $\omega_{i_0|...|i_k}$ (of degree k), where $0\le i_0 < i_1 < ... < i_k \le n$, representing the sub-simplex of the n-simplex spanned by $i_0, ..., i_k$.

### Constructors

* `DupontForm(n, form)`: basic constructor; arguments are
    * `n`: simplicial dimension
    * `form`: the actual content of the form, given as a dictionary with keys of the form `i_0|...|i_k` (representing the basic form $\omega_{i_0|...|i_k}$) and corresponding elements given by the rational coefficient of the key form (written as string or integer). For example, the dictionary `{'0|2|3': '1/2', '1': '-1/3', '': 2}` will result in the form $\tfrac{1}{2}\omega_{0|2|3} - \tfrac{1}{3}\omega_{1}} + 2$.
* `DupontForm.zero(n)`: the zero form; arguments:
    * `n`: simplicial dimension.

### Operations

The following basic operations are supported:
* sum of Dupont forms (`+`, `sum`),
* comparison of Dupont forms (`==`).

### Class methods

* `i()` returns the image of the Dupont form in the Sullivan complex.

## SullivanForm

Implements the Sullivan forms, given by the polynomial differential forms on the simplices. In simplicial degree n they are given by the dg commutative algebra $\mathbb{Q}[t_0, ..., t_n, dt_1, ..., dt_n]/\sim$, where the equivalence relation is given by $t_0 + ... + t_n \sim 1$ and $dt_0 + ... + dt_n \sim 0$.

### Constructors

* `SullivanForm(n, form)`: basic constructor; arguments are
    * `n`: simplicial dimension
    * `form`: the actual content of the form, given as a dictionary with keys of the form `i_0|...|i_k` indicating a term of the form $p(t_0,\ldots,t_n)dt_{i_0}\ldots dt_{i_k}$ and associated element representing the polynomial. This is again a dictionary with keys of the form `k_0|...|k_n` (notice that $n+1$ terms need always be present) indicating the monomial $t_0^{k_0}\ldotst_n^{k_n}$ and the rational coefficient as associated element.<br>
    For example, `SullivanForm(2, {'': {'0|0|0': 1, '2|1|0': -1}, '1|2': {'1|1|0': '1/2', '1|0|1': 1}})` gives the form $1 - t_0^2 t_1 + (\tfrac{1}{2}t_1 t_2 + t_1 t_3)dt_1 dt_2$.
* `SullivanForm.zero(n)`: the zero form; arguments:
    * `n`: simplicial dimension

### Operations

The following basic operations are supported:
* sum of Sullivan forms (`+`, `sum`),
* multiplication of Sullivan forms and scalar times Sullivan form (`*`, `np.product`),
* comparison of Sullivan forms (`==`).

### Class methods

* `reduce(eliminate=0)`: using the algebraic relations, simplifies the Sullivan form by eliminating completely $t_{\text{eliminate}}$ from the expression.
* `p()`: projection from Sullivan forms to Dupont forms.
* `h()`: contraction of Sullivan forms.

# References

 1. L. Lunardon, <i>Some remarks on Dupont contraction</i>, [arXiv:1807.02517](https://arxiv.org/pdf/1807.02517.pdf).

# Credits

Developer and maintainer: Daniel Robert-Nicoud
