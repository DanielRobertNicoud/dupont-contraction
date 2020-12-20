import itertools as it
import pandas as pd

from dupontcontraction.cubical import DupontForm, SullivanForm

# Sullivan forms
sf1 = SullivanForm(3, '3*x_1*x_2^2*x_3*dx_1*dx_2')
sf2 = SullivanForm(3, '-x_1*dx_3 + 1/2*x_2*dx_2')
print(sf1)
print(sf2)
print(sf1*sf2)

# Dupont forms
# generate basis for n=3
n = 3

basis = dict()
for degree in range(n + 1):
    set_I = list(it.combinations(list(range(1, n + 1)), degree))
    set_J = list(it.product([0, 1], repeat=n-degree))
    for I, J in it.product(set_I, set_J):
        I = '|'.join([str(i) for i in I])
        J = ''.join([str(j) for j in J])
        
        k = f"{I},{J}"
        
        basis[k] = DupontForm(n, {k: 1})

# binary products
f1, f2, bp = [], [], []
for form1, form2 in it.product(basis.values(), repeat=2):
    f1.append(form1)
    f2.append(form2)
    bp.append(DupontForm.tree_product([form1, form2]))
    
binary_products = pd.DataFrame(data={
    'form_1': f1, 'form_2': f2, '\\ell_2': bp
    })