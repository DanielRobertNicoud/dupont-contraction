"""
Generate all binary trees and their sign for the transferred Coo-structure on
Dupont forms.

The generatr only outputs planar trees, to obtain non-planar trees one can
apply a permutation to the leaves (but we will not need it).

Trees are given by nested lists, for example
[[1,2], [3,4]] 2-corolla with two 2-corollas at the leaves
[[[1,2],3],4] left comb of arity 4

The generator yields couples (sign, tree, arities), where the sign is the sign
of the tree in Delta(corolla of the corresponding arity) in S* (the dual of the
operad of endomorphism of the 1-dimensional vector space in degree 1).

The arities are [arity of the whole tree, arities of left subtree, arities of
right subtree]. For example for the tree [[1, [2, 3]], [4, 5]]] the arities
are [5, [3, [1], [2, [1, [1]]]], [2, [1], [1]]]

For the signs, one can find that if (-1)^e(t) is the sign of the tree t and
t = c2 o (t1, t2), then e(t) = n1 + 1 + e(t1) + e(t2) modulo 2.
"""

import itertools as it

def binary_tree_generator(arity, shift=0):
    
    if not isinstance(arity, int):
        raise TypeError('Invalid arity')
    if arity < 1:
        raise ValueError('Invalid arity')
        
    if arity == 1:
        yield (1, 0 + shift, [1])
        return
    
    if arity == 2:
        yield (1, [0 + shift, 1 + shift], [2, [1], [1]])
        return
    
    if arity > 2:
        for n1 in range(1, arity):
            n2 = arity - n1
            for (s1, t1, a1), (s2, t2, a2) in it.product(
                    binary_tree_generator(n1, shift=shift),
                    binary_tree_generator(n2, shift=n1 + shift)
                    ):
                yield (-s1*s2*(-1)**n1, [t1, t2], [arity, a1, a2])


def map_args(tree, args):
    """
    Given a tree and a list of arguments, produce the tree with the arguments
    instead of integers at the leaves of the tree.
    
    E.g. for tree = [[1, 2], [3, 4]] and args = [a, b, c, d] we get
    [[a, b], [c, d]].
    """
    (s, t, a) = tree
    if a[0] == 1:
        return args[0]
    return [map_args((1, t[0], a[1]), args[:a[1][0]]),
            map_args((1, t[1], a[2]), args[a[1][0]:])]
    
    