# Dupont contraction

[![DOI](https://zenodo.org/badge/289950438.svg)](https://zenodo.org/badge/latestdoi/289950438)

**Author:** Daniel Robert-Nicoud<br>
**Contact:** daniel.robertnicoud@gmail.com

**Table of contents**
* [Mathematical objects](#mathematicalobjects)
    * [Sullivan forms](#sullivanforms)
    * [Dupont forms and Dupont contraction](#dupontforms)
    * [Transferred algebraic structures](#transferredstructures)
* [Classes](#classes)
    * [DupontForms](#classdupont)
    * [SullivanForms](#classsullivan)
* [References](#references)

This package provides tools to work and do explicit computations on Sullivan and Dupont forms, as well as calculating the action of the various maps involved in the Dupont contraction and the transferred structure from the Sullivan algebra (a commutative algebra) to the Dupont algebra (which receives the structure of a commutative algebra up to homotopy).

Install this package from [PyPi](https://pypi.org/project/dupont-contraction/1.0.0/): `pip install dupont-contraction`<br>
Use this package in your code: `import dupontcontraction`

# Mathematical objects <a name="mathematicalobjects"></a>

We work in the homological convention over the field $\mathbb{Q}$ of rational numbers.

We give a brief recollection of the mathematical objects in play. For more
details, we direct the reader to the article of Lunardon.

## Sullivan forms <a name="sullivanforms"></a>

The algebra of Sullivan forms $\Omega_\bullet$ is the simplicial differential
graded (dg) commutative algebra of polynomial differential forms over the
standard simplices.

Explicitly, for $k\ge1$
$$\Omega_n = \mathbb{Q}[t_0, \ldots, t_n, dt_0,\ldots, dt_n]/\sim$$
with $t_0 + ... + t_n \sim 1$ and $dt_0 + ... + dt_n \sim 0$. Here the $t_i$
have degree $0$, the $dt_i$ degree $-1$, and $d(t_i) = dt_i$.

This object has its origin in the works of Sullivan and Bousfield-Guggenheim in
rational homotopy theory.

## Dupont forms and Dupont contraction <a name="dupontforms"></a>

Dupont defined a simplicial sub-complex $C_\bullet$ of the Sullivan forms as the chain
complex spanned by the forms
$$\omega_{i_0|\ldots|i_k} = \sum_{j=1}^n(-1)^jt_{i_j}dt_{i_0}\cdots\widehat{dt_{i_j}}\cdots dt_{i_k}$$
for $0\le i_0 < i_1 < \ldots < i_k \le n$.

This is nothing else than the (co)cellular complex of the $n$-simplex, with the
form $\omega_{i_0|\ldots|i_k}$ representing the sub-simplex spanned by the
vertices $i_0, \ldots, i_k$.

Dupont additionally showed that there exists a simplicial contraction from the
Sullivan forms to the Dupont forms, i.e. three simplicial maps
$i:C_\bullet\to\Omega_\bullet$, $p:\Omega_\bullet\to C_\bullet$ and
$h:\Omega_\bullet\to\Omega_\bullet$ satisfying certain relations, such as
$pi = \mathrm{id}$ and $ip - \mathrm{id} = dh + hd$.

## Transferred algebraic structures <a name="transferredstructures"></a>

Given a contraction from a dg commutative algebra to a sub-chain complex one can
transfer the comutative structure of the algebra to a structure on the sub-complex
that is commutative *up to homotopy*. This is called the *homotopy transfer theorem*
and we apply it in two different (but related) ways to obtain
* a structure of $\Omega\mathrm{BCom}$-algebra on the Dupont forms, and
* a structure of $\mathrm{Com}_\infty = \Omega\mathrm{Lie}^\vee$-algebra on the
Dupont forms.

By Cheng-Getzler, this second structure can be faithfully encoded by a structure
of associative algebra up to homotopy, with the full structure recovered by taking
permutations.

# Classes <a name="classes"></a>

## DupontForm <a name="classdupont"></a>

Implements the Dupont forms $C_\bullet$, given by the cellular cochains on the simplices.

### Constructors

* `DupontForm(n, form)`: basic constructor; arguments are
    * `n`: simplicial dimension
    * `form`: the actual content of the form, given as a dictionary with keys of the form `i_0|...|i_k` (representing the basic form $\omega_{i_0|...|i_k}$) and corresponding elements given by the rational coefficient of the key form (written as string or integer). For example, the dictionary `{'0|2|3': '1/2', '1': '-1/3', '': 2}` will result in the form $\tfrac{1}{2}\omega_{0|2|3} - \tfrac{1}{3}\omega_{1} + 2$.
* `DupontForm.zero(n)`: the zero form; arguments:
    * `n`: simplicial dimension.

### Operations

The following basic operations are supported:
* sum of Dupont forms (`+`, `sum`),
* comparison of Dupont forms (`==`).
* `tree_product(tree)`: transferred $\Omega\mathrm{BCom}$ structure on Dupont forms from the Sullivan forms via the Dupont contraction. The generating operations f this structure are indexed by rooted trees. This function does not encode the abstract operation but calculates it on arguments directly. Its argument represents a planar tree with Dupont forms at the leaves by writing it as a nested list of Dupont forms. For example, `[[w0, w01], w012, [w2, w12]]` for `w0` given by $\omega_0$ and so on, represents the tree with a 3-corolla at the root, with at its leaves the 2-corolla with $\omega_0$ and $\omega_{0|1}$ at the leaves, the Dupont form $\omega_{0|1|2}$, and the 2-corolla with $\omega_2$ and $\omega_{1|2}$ at the leaves.
* `a_infinity_product(*args)`: transferred $\mathrm{Com}_\infty$ structure on Dupont forms from the Sullivan forms via the Dupont contraction. The arguments `args` are Dupont forms.<br>Warning: This can be very slow for high arities as the function `SullivanForm.reduce()` needs to be called on complex Sullivan forms.

### Class methods

* `d()`: differential of a Dupont form. For `w` a Dupont form we have
`w.d() == w.i().d().p()`.
* `i()`: image of the Dupont form in the Sullivan complex.

## SullivanForm <a name="classsullivan"></a>

Implements the Sullivan forms, given by the polynomial differential forms on the simplices.

### Constructors

* `SullivanForm(n, form)`: basic constructor; arguments are
    * `n`: simplicial dimension
    * `form`: the actual content of the form, given as a dictionary with keys of the form `i_0|...|i_k` indicating a term of the form $p(t_0,\ldots,t_n)dt_{i_0}\ldots dt_{i_k}$ and associated element representing the polynomial. This is again a dictionary with keys of the form `k_0|...|k_n` (notice that $n+1$ terms need always be present) indicating the monomial $t_0^{k_0}\ldots t_n^{k_n}$ and the rational coefficient as associated element.<br>
    For example, `SullivanForm(2, {'': {'0|0|0': 1, '2|1|0': -1}, '1|2': {'1|1|0': '1/2', '1|0|1': 1}})` gives the form $1 - t_0^2 t_1 + (\tfrac{1}{2}t_1 t_2 + t_1 t_3)dt_1 dt_2$.
* `SullivanForm.zero(n)`: the zero form; arguments:
    * `n`: simplicial dimension

### Operations

The following basic operations are supported:
* sum of Sullivan forms (`+`, `sum`),
* multiplication of Sullivan forms and scalar times Sullivan form (`*`, `np.product`),
* comparison of Sullivan forms (`==`).

### Class methods

* `d()`: differential of a Sullivan form.
* `reduce(eliminate=0)`: using the algebraic relations, simplifies the Sullivan form by eliminating completely $t_{\text{eliminate}}$ from the expression.
* `p()`: projection from Sullivan forms to Dupont forms.
* `h()`: contraction of Sullivan forms.

# References <a name="references"></a>

1. X. Z. Cheng and E. Getzler. <i>Transferring homotopy commutative algebraic structures</i>. Journal of Pure and Applied Algebra, 212:2535–2542, 2008. [arXiv:math/0610912](https://arxiv.org/pdf/math/0610912.pdf).
2. L. Lunardon. <i>Some remarks on Dupont contraction</i>. [arXiv:1807.02517](https://arxiv.org/pdf/1807.02517.pdf).
