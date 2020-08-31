##################
Dupont contraction
##################

**Author:** Daniel Robert-Nicoud

**Contact:** daniel.robertnicoud@gmail.com

.. contents:: :local:
    :depth: 3

This package provides tools to work and do explicit computations on Sullivan and
Dupont forms, as well as calculating the action of the various maps involved in
the Dupont contraction and the transferred structure from the Sullivan algebra
(a commutative algebra) to the Dupont algebra (which receives the structure of a
commutative algebra up to homotopy).

For the discussion below, we adopt the homological convention and we work over
the field :math:`\mathbb{Q}` of rational numbers.

********************
Mathematical objects
********************

TBD

*******
Classes
*******

DupontForm
==========

.. code-block:: python

  dupontcontraction.dupontforms.dupontforms

Implements the Dupont forms :math:`C_\bullet`, given by the cellular cochains on
the simplices. The Dupont forms in simplicial degree :math:`n` form a chain complex with
basis :math:`\omega_{i_0|\ldots|i_k}\in C_\bullet` for (non-empty) sequences
:math:`0\le i_0 < i_1 < \ldots < i_k \le n`$`, representing the sub-simplex of
the :math:`n`-simplex spanned by :math:`i_0, \ldots, i_k`$`.  The form
:math:`\omega_{i_0|\ldots|i_k}` has degree :math:`k`.

Constructors
------------

- :code:`DupontForm(n, form)`: basic constructor; arguments are
    - :code:`n`: simplicial dimension
    - :code:`form`: the actual content of the form, given as a dictionary with keys of the form :code:`'i_0|...|i_k'` (representing the basic form :math:`\omega_{i_0|\ldots|i_k}``) and corresponding elements given by the rational coefficient of the key form (written as string or integer). For example, the dictionary :code:`{'0|2|3': '1/2', '1': '-1/3', '0': 2}` will result in the form :math:`\tfrac{1}{2}\omega_{0|2|3} - \tfrac{1}{3}\omega_{1} + 2\omega_0`. The empty string key :code:`''` corresponds to the sum :math:`\omega_0 + \cdots + \omega_n \sim 1`.
- :code:`DupontForm.zero(n)`: the zero form; arguments:
    - `n`: simplicial dimension.

Operations
----------

The following basic operations are supported:

- sum of Dupont forms (:code:`+`, :code:`sum`),
- comparison of Dupont forms (:code:`==`),
- multiplication by scalars (:code:`*`),
- :code:`tree_product(tree)`: transferred :math:`\Omega\mathrm{BCom}`` structure on Dupont forms from the Sullivan forms via the Dupont contraction. The generating operations of this structure are indexed by rooted trees. This function does not encode the abstract operation but calculates it on arguments directly. Its argument represents a planar tree with Dupont forms at the leaves by writing it as a nested list of Dupont forms. For example, :code:`[[w0, w01], w012, [w2, w12]]` for :code:`w0` given by :math:`\omega_0` and so on, represents the tree with a :math:`3`-corolla at the root, with at its leaves the :math:`2`-corolla with :math:`\omega_0` and :math:`\omega_{0|1}` at the leaves, the Dupont form :math:`\omega_{0|1|2}`, and the :math:`2`-corolla with :math:`\omega_2` and :math:`\omega_{1|2}` at the leaves.
- :code:`a_infinity_product(\*args)`: transferred :math:`\mathscr{C}_\infty` structure on Dupont forms from the Sullivan forms via the Dupont contraction. By Cheng-Getzler, equivalently we only give the transferred :math:`\mathscr{A}_\infty` structure (all other operations are obtained by applying permutations). The arguments :code:`args` are Dupont forms. Warning: This can be very slow for high (i.e. greater than :math:`6`) arities as the function :code:`SullivanForm.reduce()` needs to be called on complex Sullivan forms.

Class methods
-------------

- :code:`i()` returns the image of the Dupont form in the Sullivan complex.

SullivanForm
============

.. code-block:: python

  dupontcontraction.sullivanforms.sullivanforms

Implements the Sullivan forms, given by the polynomial differential forms on the
simplices. In simplicial degree :math:`n` they are given by the dg commutative
algebra :math:`\mathbb{Q}[t_0, ..., t_n, dt_1, ..., dt_n]/\sim`, where the
equivalence relation is given by :math:`t_0 + ... + t_n \sim 1` and
:math:`dt_0 + ... + dt_n \sim 0`.

Constructors
------------

- :code:`SullivanForm(n, form)`: basic constructor; arguments are
    - :code:`n`: simplicial dimension
    - :code:`form`: the actual content of the form, given as a dictionary with keys of the form :code:`i_0|...|i_k` indicating a term of the form :math:`p(t_0,\ldots,t_n)dt_{i_0}\ldots dt_{i_k}` and associated element representing the polynomial. This is again a dictionary with keys of the form :code:`k_0|...|k_n` (notice that :math:`n+1` terms need always be present) indicating the monomial :math:`t_0^{k_0}\cdots t_n^{k_n}` and the rational coefficient as associated element. For example, :code:`SullivanForm(2, {'': {'0|0|0': 1, '2|1|0': -1}, '1|2': {'1|1|0': '1/2', '1|0|1': 1}})` gives the form :math:`1 - t_0^2 t_1 + (\tfrac{1}{2}t_1 t_2 + t_1 t_3)dt_1 dt_2`.
- :code:`SullivanForm.zero(n)`: the zero form; arguments:
    - :code:`n`: simplicial dimension

Operations
----------

The following basic operations are supported:

- sum of Sullivan forms (:code:`+`, :code:`sum`),
- multiplication of Sullivan forms and scalar multiplication (:code:`*`, :code:`np.product`),
- comparison of Sullivan forms (:code:`==`).

Class methods
-------------

- :code:`reduce(eliminate=0)`: using the algebraic relations, simplifies the Sullivan form by eliminating completely :math:`t_{\text{eliminate}}` and :math:`dt_{\text{eliminate}}` from the expression.
- :code:`p()`: projection from Sullivan forms to Dupont forms.
- :code:`h()`: contraction of Sullivan forms.

**********
References
**********

- X. Z. Cheng and E. Getzler. *Transferring homotopy commutative algebraic structures*. Journal of Pure and Applied Algebra, 212:2535–2542, 2008. `arXiv:math/0610912 <https://arxiv.org/pdf/math/0610912.pdf>`_.
- L. Lunardon. *Some remarks on Dupont contraction*. `arXiv:1807.02517 <https://arxiv.org/pdf/1807.02517.pdf>`_.
