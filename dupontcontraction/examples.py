import pandas as pd

from sullivanforms import sullivanforms as sf
from dupontforms import dupontforms as duf

# =============================================================================
# Example 1: calculate p(i(w_{0|1|2}) * t_{1})
# =============================================================================
print('Example 1: calculate p(i(w_{0|1|2}) * t_{1})\n----------')

w012 = duf.DupontForm(2, {'0|1|2': 1})
t1 = sf.SullivanForm(2, {'': {'0|1|0': 1}})

df_ex1 = pd.DataFrame(
    data=[
        ['w012', repr(w012)],
        ['t1', repr(t1)],
        ['i(w012)', repr(w012.i())],
        ['i(w012)*t1', repr(w012.i() * t1)],
        ['p(i(w012)*t1)', repr((w012.i() * t1).p())]
    ],
    columns=['expression', 'result']
)
    
print(df_ex1)