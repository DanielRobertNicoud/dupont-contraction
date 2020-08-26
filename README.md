# Dupont contraction

This package provides tools to work and do explicit computations on Sullivan and Dupont forms, as well as calculating the action of the various maps involved in the Dupont contraction and the transferred structure from the Sullivan algebra (a commutative algebra) to the Dupont algebra (which receives the structure of a commutative algebra up to homotopy).

# Classes

## DupontForm

### Constructors

* `DupontForm(n, form)`: basic constructor, arguments are
    * `n`: simplicial dimension
    * `form`: the actual content of the form, given as a dictionary with keys of the form `i_0|...|i_k` (representing the basic form $\omega_{i_0|...|i_k}$) and corresponding elements given by the rational coefficient of the key form (written as string or integer). For example, the dictionary `{'0|2|3': '1/2', '1': '-1/3', '': 2}` will result in the form 1/2 w_(0|2|3) - 1/3 w_(1) + 2.

### Operations

The following basic operations are supported:
* sum of Dupont forms.

### Class methods

* `i()` returns the image of the Dupont form in the Sullivan complex.

## SullivanForm

TBD

# References

 1. L. Lunardon, <i>Some remarks on Dupont contraction</i>, [arXiv:1807.02517](https://arxiv.org/pdf/1807.02517.pdf).

# Credits

Developer and maintainer: Daniel Robert-Nicoud
