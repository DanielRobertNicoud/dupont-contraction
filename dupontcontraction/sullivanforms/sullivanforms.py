"""
Class for Sullivan forms.

This class implements polynomial differential forms on the simplices.

The base field is Q (the rational numbers).
"""

import itertools as it
import fractions
from copy import deepcopy
import numpy as np
import math

import sys
import os
sys.path.append(
    os.path.dirname(
        os.path.dirname(
            os.path.dirname(os.path.abspath(__file__))
        )
    )
)

from dupontcontraction.sullivanforms import auxiliary_functions as af
from dupontcontraction.dupontforms import dupontforms as duf

class SullivanForm:
    
    def __init__(self, n, form):
        """
        Sullivan forms.
        
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
                Keys are of the form i_0|...|i_k indicating a component of
                type polynomial*dt_{i_0}...dt_{i_k}
                For each key, the element is a dictionary expressing the
                polynomial. Each of the keys of this secondary dictionary is
                of the form k_0|...|k_n where n is the dimension indicating a
                monomial coeff*t_0^{k_0}...t_n^{k_n}, and the associated
                element is the coefficient

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
            raise Exception('string argument for form is not implemented yet')
        
        if isinstance(form, dict):
            out_form = {}
            
            # Keys denote the dt_i and need to be of the form
            # i_0|i_1|...|i_k with i_j between 0 and n
            # Repetitions make the form be zero
            for key in form:
                if key:  # non-empty string
                    split_key = [int(i) for i in key.split('|')]
                    
                    # check validity
                    for i in split_key:
                        if i < 0 or i > n:
                            raise TypeError('invalid form')
                    
                    # if we have too many dt_i we get zero
                    if len(split_key) > n:
                        continue
                    
                    # sign of permutation to order the dt_i
                    sign = 1
                    for (x, y) in it.combinations(split_key, 2):
                        if x == y:
                            sign = 0
                            break
                        if x > y:
                            sign *= -1
                    
                    if sign == 0:
                        continue
                    
                    # sort
                    split_key.sort()
                else:  # empty string
                    split_key = []
                    sign = 1
                
                # Element associated to the key are again dictionaries with
                # keys k_0|k_1|...|k_n (corresponding to the monomial
                # t_0^k_0*...*t_n^k_n) linked to the coefficient (string or
                # fraction.Fraction)
                monomials = form[key]
                out_monomials = {}
                
                for m in monomials:
                    split_m = [int(k) for k in m.split('|')]
                    
                    # checks
                    if len(split_m) != n + 1:
                        raise TypeError('invalid form')
                    for k in split_m:
                        if k < 0:
                            raise TypeError('invalid form')
                    
                    # coefficient
                    coeff = sign*fractions.Fraction(monomials[m])
                    if coeff == 0:
                        continue
                    out_monomials[m] = coeff
                    
                # if non-empty list of monomials, save to output
                if out_monomials:
                    out_key = '|'.join([str(i) for i in split_key])
                    out_form[out_key] = out_monomials
            
            if out_form:
                self.form = out_form
            else:
                self.is_zero = True
                self.form = {}
            
            return
    
    def zero(n):
        """
        The zero form.
        """
        return SullivanForm(n, {})
        
        
    def __repr__(self):
        """
        Represent Sullivan form in LaTeX code.

        Returns
        -------
        str
            LaTeX string.

        """
        
        if self.is_zero:
            return '0'
        
        r = ''
        for k, ds in enumerate(self.form):
            monomials = self.form[ds]
            n_monomials = len(monomials)
            
            # polynomial
            p = ''
            for i, m in enumerate(monomials):
                coeff = monomials[m]
                
                if i > 0:
                    p += ' '
                
                if coeff < 0:
                    c = -coeff
                    p += '-'
                elif i == 0:
                    c = coeff
                else:
                    c = coeff
                    p += '+'
                
                if i > 0:
                    p += ' '
                    
                if c == 1:
                    if m == '|'.join(['0' for j in range(self.n + 1)]) \
                        and (ds == '' or len(monomials) > 1):
                        
                        p += '1'
                    else:
                        pass
                        
                else:
                    if c.denominator != 1:
                        p += f"\\frac{{{c.numerator}}}{{{c.denominator}}}"
                    else:
                        p += f"{c.numerator}"
                    
                for j, e in enumerate(m.split('|')):
                    if e == '1':
                        p += f"t_{{{j}}}"
                    elif e != '0':
                        p += f"t_{{{j}}}^{{{e}}}"
                        
            if len(p) > 0:
                if n_monomials > 1:
                    if k > 0:
                        r += ' + '
                    r += f"\\left({p}\\right)"
                else:
                    if k > 0 and p[0] == '-':
                        r += f" - {p[1:]}"
                    else:
                        if k > 0:
                            r += ' + '
                        r += f"{p}"
            if ds:
                for i in ds.split('|'):
                    r += f"dt_{{{i}}}"
                
        return r
    
    def __add__(self, sf):
        """
        Sum of two Sullivan forms.

        Parameters
        ----------
        sf : SullivanForm

        Raises
        ------
        TypeError
            If simplicial dimensions are incompatible.

        Returns
        -------
        SullivanForm
            Sum.

        """
        # check same simplicial dimension
        if self.n != sf.n:
            raise TypeError('Sullivan forms need to have the same simplicial'
                             ' dimension to be added together.')
        
        # case where one of the forms is zero
        if self.is_zero:
            return deepcopy(sf)
        if sf.is_zero:
            return deepcopy(self)
        
        # actual addition
        n_out = self.n
        form_out = deepcopy(self.form)
        for ds, p in sf.form.items():
            if ds in form_out:
                form_out[ds] = af._add_polynomials(form_out[ds], p)
                if not form_out[ds]:  # zero polynomial
                    del form_out[ds]
            else:
                form_out[ds] = p
        
        return SullivanForm(n_out, form_out)
    
    
    def __neg__(self):
        """
        Negation of a form (unary minus).

        Returns
        -------
        SullivanForm
            Negation.

        """
        return SullivanForm(
            self.n,
            {dt: {m: -c for m, c in p.items()} for dt, p in self.form.items()}
        )
    
    
    def __sub__(self, other):
        return self + (-other)
    
    
    def __mul__(self, sf):
        """
        Product of two Sullivan forms.

        Parameters
        ----------
        sf : SullivanForm

        Raises
        ------
        TypeError
            If simplicial dimensions are incompatible.

        Returns
        -------
        SullivanForm
            Product.

        """
        # check same simplicial dimension
        if self.n != sf.n:
            raise TypeError('Sullivan forms need to have the same simplicial'
                             ' dimension to be multiplied together.')
        
        # case where one of the forms is zero
        if self.is_zero or sf.is_zero:
            return SullivanForm(self.n, {})
        
        # actual multiplication
        n_out = self.n
        form_out = {}
        # pairs of dt_i combinations
        for dt1, dt2 in it.product(self.form, sf.form):
            dt1_split, dt2_split = dt1.split('|'), dt2.split('|')
            if dt1 == '':
                dt1_split = []
            if dt2 == '':
                dt2_split = []
            
            # if dt1 and dt2 have common element, we get zero
            if not set(dt1_split).isdisjoint(dt2_split):
                continue
            # otherwise the resulting dt is the union of the two and the
            # polynomial is the product of polynomials
            dt = '|'.join(dt1_split + dt2_split)
            p = af._multiply_polynomials(self.form[dt1], sf.form[dt2])
            
            if dt in form_out:
                form_out[dt] = af._add_polynomials(form_out[dt], p)
            else:
                form_out[dt] = p
            
            if not form_out[dt]:
                del form_out[dt]
        
        return SullivanForm(n_out, form_out)
    
    
    def __rmul__(self, other):
        """
        Scalar multiplication of Sullivan form.

        Parameters
        ----------
        other : int or Fraction
            Scalar.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return SullivanForm(
            self.n,
            {dt: {m: other*c for m, c in p.items()} \
             for dt, p in self.form.items()}
        )
    
    
    def pullback(self, f):
        """
        Pullback of Sullivan form of simplicial degree n by injective,
        increasing map f\in I(m, n), encoded as a subset of [0,...,n].

        Parameters
        ----------
        f : list
            Subset of [0,...,n].

        Returns
        -------
        SullivanForm
            Pullback along f.
        """
        # check that f is valid
        if len([i for i in f if i < 0 or i > self.n]) > 0:
            raise TypeError('Invalid f.')
        
        # f as map from image to range
        f = {j: i for i, j in enumerate(f)}
        
        out_n = len(f) - 1 # this is m (the simplicial dimension of the output)
        out_form = {}
        for dt, p in self.form.items():
            # pullback of form
            if dt == '':
                split_dt = []
            else:
                split_dt = [int(i) for i in dt.split('|')]
                
                # zero form
                if len(split_dt) > out_n:
                    continue
                # also if we have forms not in the sub-simplex
                if len([i for i in split_dt if i not in f]) > 0:
                    continue
            
            # pullback of polynomials
            out_p = {}
            for m, c in p.items():
                split_m = [int(i) for i in m.split('|')]
                
                # zero polynomial on restriction
                if len([i for i, e in enumerate(split_m) \
                        if e > 0 and i not in f]) > 0:
                    continue
                
                out_m = [e for i, e in enumerate(split_m) if i in f]
                out_m = '|'.join([str(i) for i in out_m])
                
                out_p[out_m] = c
            
            # map dt to range of f
            out_form['|'.join([str(f[i]) for i in split_dt])] = out_p
            
        return SullivanForm(out_n, out_form)
            
        
        raise Exception('to be implemented')
        
        return SullivanForm(m, out_form)
    
    
    def p(self):
        
        # initialize output form as 0
        out_n = self.n
        out_form = duf.DupontForm(out_n, {'': 0})
        
        # for each dt, we integrate and add the result to the output form
        for dt, p in self.form.items():            
            # in case there is no dt (ie: just a polynomial): evaluation at
            # the vertices
            if dt == '':
                aux_duf = [duf.DupontForm(out_n, {'': 0})]
                for m, c in p.items():
                    split_m = [int(e) for e in m.split('|')]
                    
                    # loop through vertices of the simplex
                    n_nonzero, v_nonzero = 0, -1
                    for v, e in enumerate(split_m):
                        if e > 0:
                            n_nonzero += 1
                            v_nonzero = v
                    # if more than one are non-zero, then we get zero
                    if n_nonzero != 1:
                        continue
                    else:
                        aux_duf.append(
                            duf.DupontForm(
                                out_n,
                                {str(v_nonzero): c}
                            )
                        )
                        
                aux_duf = sum(aux_duf)
                out_form += aux_duf
            # otherwise, for proper forms
            else:
                aux_duf = [duf.DupontForm(out_n, {'': 0})]
                
                split_dt = [int(i) for i in dt.split('|')]

                # we consider all the sub-simplices of dimension where dt could
                # integrate to something
                for f in it.combinations(range(out_n+1), len(split_dt) + 1):
                    # dt needs to cover all coordinates except one
                    if len([i for i in f if i not in split_dt]) != 1:
                        continue
                    
                    # f as map from image to range
                    f = {j: i for i, j in enumerate(f)}
                    
                    # pullback form
                    pbf = SullivanForm(out_n, {dt: p}).pullback(f)
                    
                    # skip if pullback is zero
                    if pbf.is_zero:
                        continue
                    
                    # integrate the polynomial
                    int_poly = 0
                    # pullback of dt
                    pb_dt = '|'.join([str(f[i]) for i in split_dt])
                    sign = [(-1)**i for i in range(pbf.n + 1) \
                            if i not in [f[j] for j in split_dt]][0]
                    for m, c in pbf.form[pb_dt].items():
                        split_m = [int(e) for e in m.split('|')]
                        int_poly += sign * c * \
                            np.product([math.factorial(e) for e in split_m])/ \
                            math.factorial(sum(split_m) + pbf.n)
                    
                    # add resulting form
                    aux_duf.append(
                        duf.DupontForm(
                            out_n,
                            {'|'.join([str(i) for i in f]): int_poly}
                        )
                    )
                    
                aux_duf = sum(aux_duf)
                out_form += aux_duf            
        
        return out_form
    
    def d(self):
        """
        Differential.
        """
        out_n = self.n
        out_form = SullivanForm.zero(out_n)
        
        for dt, p in self.form.items():
            for m, c in p.items():
                # d(monomial)
                split_m = [int(e) for e in m.split('|')]
                for i, e in enumerate(split_m):
                    if e == 0:
                        continue
                    
                    aux_split_m = deepcopy(split_m)
                    aux_split_m[i] -= 1
                    
                    aux_m = '|'.join([str(x) for x in aux_split_m])
                    aux_c = c * e
                    if dt == '':
                        aux_dt = str(i)
                    else:
                        aux_dt = f"{i}|{dt}"
                    
                    # add term to final form
                    out_form += SullivanForm(
                        out_n,
                        {aux_dt: {aux_m: aux_c}}
                    )
        
        return out_form
    
    
    def apply_permutation(self, permutation):
        """
        Action of permutation on a form (mapping t_i to t_permutation(i)).
        
        Permutation given as a dictionary, an be only partial in which case it
        is automatically completed with identity in all missing values.
        
        Eg: {1: 2, 2: 1} is the switch (1 2).

        Parameters
        ----------
        permutation : dict
            Permutation to apply.

        Returns
        -------
        SullivanForm
            Permuted form.

        """
        
        # complete permutation
        for i in range(self.n + 1):
            if i not in permutation:
                permutation[i] = i
        
        # inverse permutation
        permutation_inv = {j: i for i, j in permutation.items()}
        
        
        out_n = self.n
        out_form = {}
        
        for dt, p in self.form.items():
            
            # permute dt
            if dt == '':
                dt_split = []
            else:
                dt_split = [int(i) for i in dt.split('|')]
            aux_dt = [str(permutation[i]) for i in dt_split]
            aux_dt = '|'.join(aux_dt)
            
            aux_p = {}
            for m, c in p.items():
                # permute t
                m_split = [int(e) for e in m.split('|')]
                aux_m = [str(m_split[permutation_inv[i]]) \
                         for i in range(out_n + 1)]
                aux_m = '|'.join(aux_m)
                
                aux_p[aux_m] = c
            
            out_form[aux_dt] = aux_p
        
        return SullivanForm(out_n, out_form)
    
    
    def reduce(self, eliminate=0):
        """
        Simplify the form by eliminating all occurrences of t_[eliminate].

        Parameters
        ----------
        eliminate : int, optional
            The t to eliminate from the expression. The default is 0.

        Returns
        -------
        SullivanForm
            Simplified form.

        """
        
        out_n = self.n
        
        # first reduce the dt
        temp_form = SullivanForm.zero(out_n)
        for dt, p in self.form.items():
            if dt == '':
                dt_split = []
            else:
                dt_split = [int(i) for i in dt.split('|')]
                        
            if eliminate not in dt_split:
                temp_form += SullivanForm(out_n, {dt: p})
            else:
                i_elim = dt_split.index(eliminate)
                
                aux_form = {}
                for i in range(out_n + 1):
                    if i in dt_split:  # either eliminate or already occurring
                        continue
                    
                    aux_dt = deepcopy(dt_split)
                    aux_dt[i_elim] = i
                    aux_dt = '|'.join([str(j) for j in aux_dt])
                    aux_form[aux_dt] = p
                    
                # notice the minus sign
                aux_form = -SullivanForm(out_n, aux_form)
                
                temp_form += aux_form
        
        # then reduce the polynomials
        
        # probably not optimal, but the easiest and clearest way of
        # implementing it
        
        # t_eliminate = 1 - t_0 - ... - t_n () only t_eliminate not appearing
        # on the right-hand side)
        replacement_poly = {
            '|'.join([str(int(j == i)) for j in range(out_n + 1)]): -1 \
                for i in range(out_n + 1) if i != eliminate
        }
        replacement_poly['|'.join(['0']*(out_n + 1))] = 1
        replacement_poly = SullivanForm(out_n, {'': replacement_poly})
        
        # replace occurrences of t_eliminate with the replacement polynomial
        out_form = SullivanForm.zero(out_n)
        for dt, p in temp_form.form.items():
            for m, c in p.items():
                m_split = [int(e) for e in m.split('|')]
                
                # exponent of t_eliminate
                exponent_elim = m_split[eliminate]
                
                if exponent_elim == 0:
                    out_form += SullivanForm(out_n, {dt: {m: c}})
                else:
                    aux_m = m_split
                    aux_m[eliminate] = 0
                    aux_m = '|'.join([str(i) for i in aux_m])
                    aux_form = SullivanForm(out_n, {dt: {aux_m: c}})
                    
                    for _ in range(exponent_elim):
                        aux_form *= replacement_poly
                        
                    out_form += aux_form
        
        return out_form
    
    
    def hj(self, j):
        """
        Implement the auxiliary function h_j (Lunardon, p.7) which will then
        be used for the contraction.

        Parameters
        ----------
        j : int
            DESCRIPTION.

        Returns
        -------
        SullivanForm.
        """
        
        if j != 0:
            # reduce to the case j = 0
            out_form = self.apply_permutation({0: j, j: 0})
            out_form = out_form.hj(0)
            return out_form.apply_permutation({0: j, j: 0})
        
        # case j = 0
        out_n = self.n
        # first eliminate t_0 from the expression
        aux_form = self.reduce()
        # there is an easy formula for a form written as
        # x1^k1...xn^kn dt_c1 dt_c2 ... dt_cm with 1 < c1 < c2 < ... < cm,
        # see Lunardon, p. 12
        out_form = SullivanForm.zero(out_n)
        for dt, p in aux_form.form.items():
            if dt == '':
                # in this case we get 0 for the monomial
                continue
            else:
                dt_split = [int(i) for i in dt.split('|')]
                
            for m, c in p.items():
                m_split = [int(e) for e in m.split('|')]
                
                aux_c = c / (sum(m_split) + len(dt_split))
                
                for i, k in enumerate(dt_split):
                    aux_dt = deepcopy(dt_split)
                    aux_dt.remove(k)
                    aux_dt = '|'.join([str(j) for j in aux_dt])
                    
                    aux_m = deepcopy(m_split)
                    aux_m[k] += 1
                    aux_m = '|'.join([str(j) for j in aux_m])
                    
                    # notice we start with c_1 and not c_0 in the formula in
                    # Lunardon, so we take (-1)^(i-1) here
                    out_form += SullivanForm(
                        out_n,
                        {aux_dt: {aux_m: -((-1)**i) * aux_c}}
                    )
        
        return out_form
    
    
    def hf(self, f):
        """
        Implement the auxiliary function h_f (Lunardon, p.7) which will then
        be used for the contraction.

        Parameters
        ----------
        f : list or tuple
            Increasing, elements between 0 and n.

        Returns
        -------
        SullivanForm.
        """
        # check that f is valid
        last_i = -1
        for i in f:
            if i <= last_i or i < 0 or i > self.n or not isinstance(i, int):
                raise TypeError('invalid f')
        
        out_form = deepcopy(self)
        for j in f:
            out_form = out_form.hj(j)
        
        return out_form
    
    
    def h(self):
        """
        The contraction h.

        Returns
        -------
        SullivanForm.
        """
        out_form = SullivanForm.zero(self.n)
        for k in range(1, self.n + 2):
            for f in it.combinations(range(0, self.n + 1), k):
                wf = duf.DupontForm(self.n, {'|'.join([str(i) for i in f]): 1})
                wf = wf.i()
                
                out_form += wf * self.hf(f)
        
        return out_form
        
                
            
if __name__ == '__main__':
    n = 3
    form = {'1|3|2': {'0|1|1|0': '-3/7', '2|1|0|0': '1'},
            '': {'1|1|0|0': '14'}}
    sf = SullivanForm(n, form)
    print(sf, '\n')
    
    form2 = {'': {'0|0|0|0': '1', '1|1|0|0': '-1'},
             '1|2|3': {'0|0|0|1': '1'}}
    sf2 = SullivanForm(n, form2)
    print(sf2, '\n')
    
    one = SullivanForm(n, {'': {'0|0|0|0': '1'}})
    print(one, '\n')
    
    print(sf + sf2, '\n')
    print(sf * sf2, '\n')
    print(sf * one, '\n')
    print(-one * one, '\n')
    print(-sf, '\n')
    print(3*sf, '\n')
    print(fractions.Fraction('1/2')*sf, '\n')
    
    form3 = {'': {'1|0|0|0': 2, '0|1|0|0': '3/4'}}
    sf3 = SullivanForm(n, form3)
    print(sf3.p(), '\n')
    
    form_p = {'1|2': {'1|2|0|0': '1/2', '0|1|0|0': '1'}}
    sfp = SullivanForm(n, form_p)
    print(sfp, '\n')
    print(sfp.p(), '\n')
