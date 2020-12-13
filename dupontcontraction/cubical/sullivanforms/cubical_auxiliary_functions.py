import itertools as it

def _add_polynomials(p1, p2):
    """
    Add two polynomials (auxiliary function).

    Parameters
    ----------
    p1 : dict
        Polynomial as in the forms.
    p2 : dict
        Polynomial as in the forms.

    Returns
    -------
    p_out : dict
        Sum of the polynomials.

    """        
    p_out = p1.copy()
    
    for m, coeff in p2.items():
        if m in p_out:
            p_out[m] += coeff
            if p_out[m] == 0:
                del p_out[m]
        else:
            p_out[m] = coeff
            
    return p_out

def _multiply_polynomials(p1, p2):
    """
    Multiply two polynomials (auxiliary function).

    Parameters
    ----------
    p1 : dict
        Polynomial as in the forms.
    p2 : dict
        Polynomial as in the forms.

    Returns
    -------
    p_out : dict
        Product of the polynomials.

    """
    p_out = {}
    
    for m1, m2 in it.product(p1, p2):
        c1, c2 = p1[m1], p2[m2]
        
        m_out = '|'.join(
            [str(int(k1) + int(k2)) for k1, k2 \
             in zip(m1.split('|'), m2.split('|'))]
        )
        
        c_out = c1*c2
        
        if m_out in p_out:
            p_out[m_out] += c_out
            if p_out[m_out] == 0:
                del p_out[m_out]
        else:
            p_out[m_out] = c_out
    
    return p_out