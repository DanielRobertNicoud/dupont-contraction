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