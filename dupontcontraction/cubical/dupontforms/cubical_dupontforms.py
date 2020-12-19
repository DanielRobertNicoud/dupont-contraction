"""
Class for cubical Dupont forms.
"""

import sys
import os
import fractions
from copy import deepcopy

sys.path.append(
    os.path.dirname(
        os.path.dirname(
            os.path.dirname(os.path.abspath(__file__))
        )
    )
)

from cubical.signed_ordered_set import SignedOrderedSet
import cubical.sullivanforms.cubical_sullivanforms as csf


class DupontForm:
    
    def __init__(self, n, form):
        """
        

        Parameters
        ----------
        n : int
            Cubical dimension.
        form : dict or string
            For a dict argument:
                Keys are of the form 'i_1|...|i_k,j_1...j_{n-k} indicating the
                basic form \omega_{I, J} where I is the set {i_1,...,i_k}
                (ordered) and j_1,...,j_{n-k} are the images of the complement
                of I (ordered). For example, for n=4 the key '2|4,01' denotes
                I={2,4}, J(1) = 0, J(3) = 1. The associated element is the
                coefficient of that form (a rational number).
        """
        
        if not isinstance(n, int):
            raise TypeError('n is not an integer')
        if not (isinstance(form, dict) or isinstance(form, str)):
            raise TypeError('invalid form')
        
        self.n = n
        
        if form == '0':
            self.is_zero = True
            self.form = dict()
            return
        
        self.is_zero = False
        
        if isinstance(form, str):
            raise NotImplementedError('string argument for form is not '
                                      'implemented yet')
        
        if isinstance(form, dict):
            out_form = dict()
            for k, c in form.items():
                I, J = k.split(',')
                
                if I != '':
                    I = [int(i) for i in I.split('|')]
                else:
                    I = []
                    
                for i in I:
                    if i <= 0 or i > n or not isinstance(i, int):
                        raise TypeError('Invalid form')
                I = SignedOrderedSet(I)
                
                sign = I.sign()
                degree = len(I)
                
                if len(J) != n - degree:
                    raise TypeError('Invalid form')
                
                for j in J:
                    if j not in '01':
                        raise TypeError('Invalid form')
                        
                key = f"{str(I)},{J}"
                coeff = sign*fractions.Fraction(c)
                
                if key in out_form:
                    out_form[key] = coeff + out_form[key]
                else:
                    out_form[key] = coeff
                
                if out_form[key] == 0:
                    del out_form[key]
            
            self.form = out_form
            
            # check for zero forms
            if not out_form:
                self.is_zero = True
        
        
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
        return DupontForm(n, '0')
        
    
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
        for i, (k, c) in enumerate(self.form.items()):
            
            if c > 0 and i > 0:
                r += ' + '
            elif c < 0 and i > 0:
                r += ' - '
            elif c < 0 and i == 0:
                r += '-'
            
            if c.denominator == 1:
                if abs(c.numerator) != 1:
                    r += str(abs(c.numerator))
            else:
                r += f"\\frac{{{abs(c.numerator)}}}{{{c.denominator}}}"
            
            if k[0] != ',' and k[-1] != ',':
                r += f"\\omega_{{{k}}}"
            elif k[-1] != ',':
                r += f"\\omega_{{\\emptyset{k}}}"
            else:
                r += f"\\omega_{{{k}\\emptyset}}"
        
        return r
    
    
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
    
    def d(self):
        """
        Differential of Dupont form.
        """
        out_n = self.n
        out_form = DupontForm.zero(out_n)
        
        for k, c in self.form.items():
            I, J = k.split(',')
            I = [int(i) for i in I.split('|')]
            J = [int(j) for j in list(J)]
            for cnt, (i, j) in  enumerate(
                    zip([x for x in range(1, out_n+1) if x not in I], J)):
                    
                sign = (-1)**(j + (i - 1 - cnt))
                out_I = '|'.join([str(y) for y in 
                                  I[:i - 1 - cnt] + [i] + I[i - 1 - cnt:]])
                out_J = ''.join([str(z) for z in J[:cnt] + J[cnt + 1:]])
                
                out_form += DupontForm(out_n, {f"{out_I},{out_J}": sign * c})
        
        return out_form
        
        
    def i(self):
        """
        Map i of the contraction, returns a SullivanForm.
        """
        out_n = self.n
        out_form = csf.SullivanForm(out_n, dict())
        
        for k, c in self.form.items():
            I, J = k.split(',')
            dx = csf.SullivanForm(out_n, {I: {'|'.join(['0'] * out_n): c}})
            
            I = I.split('|')
            
            J = list(J)
            J.reverse()
            for i in range(1, out_n + 1):
                if str(i) not in I:
                    j = J.pop()
                    
                    if j == '0':
                        dx *= csf.SullivanForm(out_n, f"1 - x_{i}")
                    else:
                        dx *= csf.SullivanForm(out_n, f"x_{i}")
            
            out_form += dx
        
        return out_form


if __name__ == '__main__':
    df1 = DupontForm(3, {'1|2,0': 1})
    df2 = DupontForm(3, {'1|3,0': 1})
    df3 = DupontForm(3, {'3,01': 1})
    
    print((df1.i() * df2.i()).p())
    print(df1.i() * df3.i())
    print((df1.i() * df3.i()).p())
    