"""
Class for Dupont forms.
"""

import itertools as it
import fractions
from copy import deepcopy
import math
import numpy as np

import sys
import os
sys.path.append(
    os.path.dirname(
        os.path.dirname(
            os.path.dirname(os.path.abspath(__file__))
        )
    )
)

import dupontcontraction.sullivanforms.sullivanforms as sf
import dupontcontraction.dupontforms.binary_tree_generator as btg

class DupontForm:
    
    def __init__(self, n, form):
        """
        Dupont forms.
        
        Attributes
        ----------
        n : int
            Simplicial dimension.
        is_zero : bool
            Indicates if the form is zero or not.
        form : dict
            Content of the form.

        Parameters
        ----------
        n : int
            Simplicial dimension.
        form : dict or string
            For a dict argument:
                Keys are of the form i_0|...|i_k indicating the basic form
                \omega_{i_0...i_k}, and the associated element is the
                coefficient of that form (a rational number).
                Note that 1 = \omega_0 + ... + \omega_n

        Raises
        ------
        TypeError
            For mis-specified parameters.
        Exception
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        if not isinstance(n, int):
            raise TypeError('n is not an integer')
        if not (isinstance(form, dict) or isinstance(form, str)):
            raise TypeError('invalid form')
        
        self.n = n
        self.is_zero = False
        
        if isinstance(form, str):
            raise NotImplementedError('string argument for form is not '
                                      'implemented yet')
        
        if isinstance(form, dict):
            out_form = {}
            
            # transform 1 in the sum of the degree 0 basic forms
            if '' in form:
                c = fractions.Fraction(form[''])
                for v in range(n + 1):
                    if str(v) in form:
                        form[str(v)] = fractions.Fraction(form[str(v)]) + c
                    else:
                        form[str(v)] = c
            
            for w in form:
                if w == '':
                    continue
                else:
                    w_split = [int(i) for i in w.split('|')]
                    
                # check validity
                for i in w_split:
                    if i < 0 or i > n:
                        raise TypeError('invalid form')
                    
                # check for invalid (zero) forms
                # notice that degree = len(w_split) - 1
                if len(w_split) > n + 1:
                    continue
                
                # sign of permutation
                sign = 1
                for (i, j) in it.combinations(w_split, 2):
                    if i == j:
                        sign = 0
                        break
                    if i > j:
                        sign *= -1
                
                if sign == 0:
                    continue
                
                # add valid forms to out
                w_split.sort()
                w_out = '|'.join([str(i) for i in w_split])
                out_form[w_out] = sign*fractions.Fraction(form[w])
                if out_form[w_out] == 0:
                    del out_form[w_out]
            
            # check for zero forms
            if not out_form:
                self.is_zero = True
            
            self.form = out_form
        return
    
    
    def __eq__(self, other):
        """
        Check equality of Dupont forms by comparing coefficients on the basis.
        """        
        for w, c in self.form.items():
            if w not in other.form or c != other.form[w]:
                return False
        for w, c in other.form.items():
            if w not in self.form or c != self.form[w]:
                return False
        return True
    
    
    def zero(n):
        """
        Zero Dupont form of simplicial degree n.
        """
        return DupontForm(n, {'': 0})
        
    
    def __repr__(self):
        """
        Represent Dupont form in LaTeX code.

        Returns
        -------
        str
            LaTeX string.

        """
        if self.is_zero:
            return '0'
        
        r = ''
        for k, w in enumerate(self.form):
            c = self.form[w]
            
            if k > 0:
                r += ' '
            
            if c < 0:
                r += '-'
                c = -c
            elif k > 0:
                r += '+'
            
            if k > 0:
                r += ' '
            
            if c == 1:
                if w == '':
                    if c.denominator == 1:
                        r += f"{c.numerator}"
                    else:
                        r += f"\\frac{{{c.numerator}}}{{{c.denominator}}}"
                else:
                    pass
            elif c.denominator == 1:
                r += f"{c.numerator}"
            else:
                r += f"\\frac{{{c.numerator}}}{{{c.denominator}}}"
            
            if w == '':
                continue
            else:
                r += f"\omega_{{{w}}}"
            
        return r
                
    
    def __add__(self, other):
        """
        Sum of two Dupont forms.

        Parameters
        ----------
        sf : DupontForm

        Raises
        ------
        TypeError
            If simplicial dimensions are incompatible.

        Returns
        -------
        DupontForm
            Sum.

        """
        # check same simplicial dimension
        if self.n != other.n:
            raise TypeError('Sullivan forms need to have the same simplicial'
                             ' dimension to be added together.')
        
        # case where one of the forms is zero
        if self.is_zero:
            return deepcopy(other)
        if other.is_zero:
            return deepcopy(self)
        
        # actual sum
        out_form = deepcopy(self.form)
        for w, c in other.form.items():
            if w in out_form:
                out_form[w] += c
                if out_form[w] == 0:
                    del out_form[w]
            else:
                out_form[w] = c
        
        return DupontForm(self.n, out_form)
    
    
    def __rmul__(self, other):
        """
        Scalar multiplicaiton of Dupont forms.
        """
        try:
            other = fractions.Fraction(other)
        except:
            raise TypeError('Invalid scalar multiplication.')
        
        return DupontForm(
            self.n,
            {w: other*c for w, c in self.form.items()}
        )
    
    
    def __radd__(self, other):
        """
        Sum of scalar with form.
        """
        return DupontForm(self.n, {'': other}) + self
        
        
    def i(self):
        """
        Map i of the contraction, returns a SullivanForm.
        """
        out_n = self.n
        out_form = sf.SullivanForm(out_n, {})
        
        for w, c in self.form.items():
            if w == '':
                aux_sf = sf.SullivanForm(
                    out_n,
                    {'': {'|'.join(['0' for i in range(out_n + 1)]): c}}
                )
            else:
                w_split = w.split('|')
                # Sullivan form associated to w
                aux_sf = {}
                for k in range(len(w_split)):
                    sign = (-1)**k
                    dt = '|'.join(w_split[:k] + w_split[k+1:])
                    
                    p = ['0' for i in range(out_n + 1)]
                    p[int(w_split[k])] = '1'
                    p = '|'.join(p)
                    
                    aux_sf[dt] = {p: sign*c}
                
                aux_sf = math.factorial(len(w_split) - 1) * \
                    sf.SullivanForm(out_n, aux_sf)
            
            out_form += aux_sf
        
        return out_form
    
    
    def tree_product(tree):
        """
        Compute the basic operation of the Cobar(Bar(Com))-algebra structure on
        Dupont forms obtained by HTT from Sullivan forms via the Dupont contra-
        ction. These basic operations are indexed by rooted trees.
        

        Parameters
        ----------
        tree : nested lists of Dupont forms.
            For example [[a, b], [c, d, e]] with a,...,e Dupont forms. Every
            vertex needs to be at least binary.

        Returns
        -------
        Dupont form
            The product of the transferred structure.

        """
        return DupontForm._tree_product(tree)
    
    def _tree_product(tree, root=True):
        
        if len(tree) < 2:
            raise TypeError('Invalid tree.')
        
        aux_tree = []
        for sub_tree in tree:
            # if it is a Dupont form
            if isinstance(sub_tree, DupontForm):
                # any appearance of zero makes the whole thing zero
                if sub_tree.is_zero:
                    return DupontForm.zero(sub_tree.n)
                aux_tree.append(sub_tree.i())
            # otherwise iterate
            else:
                aux_tree.append(DupontForm._tree_product(sub_tree, root=False))
                # any appearance of zero makes the whole thing zero
                if aux_tree[-1].is_zero:
                    return DupontForm.zero(aux_tree[-1].n)
        
        # now we have just a list of Sullivan form, take the product and apply
        # h or p depending if we are at the root or not
        if root:
            return np.product(aux_tree).p()
        else:
            return np.product(aux_tree).h()
    
    
    def a_infinity_product(*args):
        
        for duf in args:
            if not isinstance(duf, DupontForm):
                raise TypeError('All elements need to be Dupont forms.')
        
        out_n = args[0].n
        for duf in args:
            if duf.n != out_n:
                raise ValueError('All Dupont forms need to have the same '
                                 'simplicial dimension.')
        
        arity = len(args)
        
        out_form = DupontForm.zero(out_n)
        for tree in btg.binary_tree_generator(arity):
            sign = tree[0]
            tree = btg.map_args(tree, args)
            
            form = sign*DupontForm.tree_product(tree)
            
            out_form += form
        
        return out_form


if __name__ == '__main__':
    # playground
    w0 = DupontForm(3, {'0': 1})
    w1 = DupontForm(3, {'1': 1})
    w2 = DupontForm(3, {'2': 1})
    w3 = DupontForm(3, {'3': 1})
    w01 = DupontForm(3, {'0|1': 1})
    w02 = DupontForm(3, {'0|2': 1})
    w03 = DupontForm(3, {'0|3': 1})
    w12 = DupontForm(3, {'1|2': 1})
    w13 = DupontForm(3, {'1|3': 1})
    w23 = DupontForm(3, {'2|3': 1})
    w012 = DupontForm(3, {'0|1|2': 1})
    w013 = DupontForm(3, {'0|1|3': 1})
    w023 = DupontForm(3, {'0|2|3': 1})
    w123 = DupontForm(3, {'1|2|3': 1})
    w0123 = DupontForm(3, {'0|1|2|3': 1})