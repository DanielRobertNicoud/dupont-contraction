"""
Class for Dupont forms.
"""

import itertools as it
import fractions
from copy import deepcopy

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
            for w in form:
                if w == '':
                    w_split = []
                else:
                    w_split = [int(i) for i in w.split('|')]
                    
                # check validity
                for i in w_split:
                    if i < 0 or i > n:
                        raise TypeError('invalid form')
                    
                # check for invalid (zero) forms
                if len(w_split) >= n + 1:
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
            
            # check for zero forms
            if not out_form:
                self.is_zero = True
            
            self.form = out_form
        return
    
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
            print(k, w)
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
                
    
    def __add__(self, df):
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
        if self.n != df.n:
            raise TypeError('Sullivan forms need to have the same simplicial'
                             ' dimension to be added together.')
        
        # case where one of the forms is zero
        if self.is_zero:
            return deepcopy(df)
        if df.is_zero:
            return deepcopy(self)
        
        # actual sum
        out_form = deepcopy(self.form)
        for w, c in df.form.items():
            if w in out_form:
                out_form[w] += c
                if out_form[w] == 0:
                    del out_form[w]
            else:
                out_form[w] = c
        
        return DupontForm(self.n, out_form)
        
        
        

if __name__ == '__main__':
    n = 3
    form = {'0|2|1': '3/4',
            '0|1': '-1'}
    df = DupontForm(n, form)
    print(df, '\n')
    
    form2 = {'': '3'}
    df2 = DupontForm(n, form2)
    print(df2, '\n')
    
    print(df + df2, '\n')
