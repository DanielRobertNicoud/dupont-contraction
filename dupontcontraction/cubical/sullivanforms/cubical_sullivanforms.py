"""
Class for Sullivan forms in the cubical context.

The base field is Q (the rational numbers).
"""

import fractions
import itertools as it
import sys
import os


sys.path.append(
    os.path.dirname(
        os.path.dirname(
            os.path.dirname(os.path.abspath(__file__))
        )
    )
)

from dupontcontraction.cubical.sullivanforms import cubical_auxiliary_functions as caf

class SullivanForm:
    
    def __init__(self, n, form):
        """
        Cubical Sullivan forms.
        
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
                Keys are of the form i_1|...|i_k indicating a component of
                type polynomial*dx_{i_1}...dx_{i_k}
                For each key, the element is a dictionary expressing the
                polynomial. Each of the keys of this secondary dictionary is
                of the form k_1|...|k_n where n is the dimension indicating a
                monomial coeff*x_1^{k_1}...x_n^{k_n}, and the associated
                element is the coefficient

        Raises
        ------
        TypeError
            For mis-specified parameters.

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
            
            # Keys denote the dx_i and need to be of the form
            # i_1|i_1|...|i_k with i_j between 1 and n
            # Repetitions make the form be zero
            for key in form:
                if key:  # non-empty string
                    split_key = [int(i) for i in key.split('|')]
                    
                    # check validity
                    for i in split_key:
                        if i < 1 or i > n:
                            raise TypeError('invalid form')
                    
                    # if we have too many dx_i we get zero
                    if len(split_key) > n:
                        continue
                    
                    # sign of permutation to order the dx_i
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
                # keys k_1|k_1|...|k_n (corresponding to the monomial
                # x_1^k_1*...*t_n^k_n) linked to the coefficient (string or
                # fraction.Fraction)
                monomials = form[key]
                out_monomials = {}
                
                for m in monomials:
                    split_m = [int(k) for k in m.split('|')]
                    
                    # checks
                    if len(split_m) != n:
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
    
    
    def copy(self):
        """
        Copy form.

        Returns
        -------
        SullivanForm.
            A copy of self

        """
        return SullivanForm(
            self.n,
            {dt: {m: c for m, c in p.items()} for dt, p in self.form.items()}
        )
        
        
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
                        p += f"x_{{{j + 1}}}"
                    elif e != '0':
                        p += f"x_{{{j + 1}}}^{{{e}}}"
                        
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
                if k > 0 and len(p) == 0:
                    r += ' + '
                for i in ds.split('|'):
                    r += f"dx_{{{i}}}"
                
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
            return sf.copy()
        if sf.is_zero:
            return self.copy()
        
        # actual addition
        n_out = self.n
        form_out = self.copy().form
        for ds, p in sf.form.items():
            if ds in form_out:
                form_out[ds] = caf._add_polynomials(form_out[ds], p)
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
            p = caf._multiply_polynomials(self.form[dt1], sf.form[dt2])
            
            if dt in form_out:
                form_out[dt] = caf._add_polynomials(form_out[dt], p)
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
    
    
    def __eq__(self, other):
        """
        Check for equality of Sullivan forms. Since we work in free algebras in
        the cubical case, we only need to compare the coefficients of the
        monomials.
        """
        
        for dt, p in self.form.items():
            if dt not in other.form:
                return False
            for m, c in p.items():
                if m not in other.form[dt] or c != other.form[dt][m]:
                    return False
        
        for dt, p in other.form.items():
            if dt not in self.form:
                return False
            for m, c in p.items():
                if m not in self.form[dt] or c != self.form[dt][m]:
                    return False
        
        return True
    
    def d(self):
        """
        Differential.
        """
        out_n = self.n
        out_form = SullivanForm.zero(out_n)
        
        for dx, p in self.form.items():
            for m, c in p.items():
                # d(monomial)
                split_m = [int(e) for e in m.split('|')]
                for i, e in enumerate(split_m):
                    if e == 0:
                        continue
                    
                    aux_split_m = [*split_m]
                    aux_split_m[i] -= 1
                    
                    aux_m = '|'.join([str(x) for x in aux_split_m])
                    aux_c = c * e
                    if dx == '':
                        aux_dx = str(i + 1)
                    else:
                        aux_dx = f"{i + 1}|{dx}"
                    
                    # add term to final form
                    out_form += SullivanForm(
                        out_n,
                        {aux_dx: {aux_m: aux_c}}
                    )
        
        return out_form


if __name__ == '__main__':
    sf1 = SullivanForm(2,
                       {'2|1': {'0|3': 1}, '': {'1|0': -2}, '1': {'0|0': -1}})
    sf2 = SullivanForm(2, {'1|2': {'0|2': -3, '0|3':2}})
    print(sf1)
    print(sf2)
    print(sf1 + sf2)
    print(-sf1)
    print(sf1 * sf1)
    print(sf1.d())
    