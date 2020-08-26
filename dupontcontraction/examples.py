import pandas as pd
import itertools as it

from sullivanforms import sullivanforms as sf
from dupontforms import dupontforms as duf

if __name__ == '__main__':
    
    n_ex = 0
    
    # =============================================================================
    # Example: basic Dupont forms and their expression as Sullivan forms
    # =============================================================================
    n_ex += 1
    print(f"Example {n_ex}: basic Dupont forms and their expression as Sullivan "
          "forms\n---------" + "-"*len(str(n_ex)))
    
    n = 2
    data = [['1', duf.DupontForm(n, {'': 1})]]
    for p in range(1, n+2):
        for ind in it.combinations(range(0, n+1), p):
            ind = '|'.join([str(i) for i in ind])
            name = f"w_{{{ind}}}"
            
            data.append([
                name,
                duf.DupontForm(n, {ind: 1})
            ])
    
    df_ex = pd.DataFrame(data=data,
                         columns=['expression', 'dupont'])
    
    df_ex['sullivan'] = df_ex.apply(lambda row: row.dupont.i(), axis=1)
    
    print(df_ex)
    
    # =============================================================================
    # Example: check that pi = id on low order forms
    # =============================================================================
    n_ex += 1
    print(f"\n\nExample {n_ex}: check that pi = id on low order forms\n---------" +
          "-"*len(str(n_ex)))
    
    df_ex['expression'] = 'pi(' + df_ex['expression'] + ')'
    df_ex['pi(dupont)'] = df_ex.apply(lambda row: row.dupont.i().p(), axis=1)
    
    print(df_ex[['expression', 'dupont', 'pi(dupont)']])
    
    # =============================================================================
    # Example: differentials of basic forms
    # =============================================================================
    n_ex += 1
    print(f"\n\nExample {n_ex}: differentials of basic forms\n---------" +
          "-"*len(str(n_ex)))
    
    df_ex['d(form)'] = df_ex.apply(lambda row: row.sullivan.d(), axis=1)
    
    print(df_ex[['dupont', 'sullivan', 'd(form)']])
    
    # =============================================================================
    # Example: hj, hf, and h
    # =============================================================================
    n_ex += 1
    print(f"\n\nExample {n_ex}: hj\n---------" +
          "-"*len(str(n_ex)))
    
    form = sf.SullivanForm(
        3,
        {'1|2': {'0|2|2|3': 1}}
    )
    
    df_ex = pd.DataFrame(
        data=[
            ['original form', repr(form)],
            ['hj(form), j=0', repr(form.hj(0))],
            ['hj(form), j=1', repr(form.hj(1))],
            ['hf(form), f=[0,1]', repr(form.hf([0, 1]))]
        ],
        columns=['expression', 'result']
    )
        
    print(df_ex)
    
    # =============================================================================
    # Example: calculate p(i(w_{0|1|2}) * t_{1}) and h(i(w_{0|1|2}) * t_{1})
    # =============================================================================
    n_ex += 1
    print(f"\n\nExample {n_ex}: calculate p(i(w_{0|1|2}) * t_{1})\n---------" +
          "-"*len(str(n_ex)))
    
    w012 = duf.DupontForm(2, {'0|1|2': 1})
    t1 = sf.SullivanForm(2, {'': {'0|1|0': 1}})
    
    df_ex = pd.DataFrame(
        data=[
            ['w012', repr(w012)],
            ['t1', repr(t1)],
            ['i(w012)', repr(w012.i())],
            ['i(w012)*t1', repr(w012.i() * t1)],
            ['p(i(w012)*t1)', repr((w012.i() * t1).p())],
            ['h(i(w012)*t1)', repr((w012.i() * t1).h())]
        ],
        columns=['expression', 'result']
    )
        
    print(df_ex)
    
    # =============================================================================
    #     Example: check that $ip - 1 = hd + dh$
    # =============================================================================
    n_ex += 1
    print(f"\n\nExample {n_ex}: check that $ip - 1 = hd + dh$\n---------" +
          "-"*len(str(n_ex)))
    
    f1 = sf.SullivanForm(3, {'': {'0|1|1|0': '1'}})
    f2 = sf.SullivanForm(3, {'1': {'0|1|1|0': '1'}})
    f3 = sf.SullivanForm(3, {'1|2': {'0|1|1|0': '1'}})
    f4 = sf.SullivanForm(3, {'1|2': {'2|0|0|0': '1', '1|1|0|0': '1/2'}})
    f5 = sf.SullivanForm(3, {'1|2': {'0|1|1|0': '1'}, '1': {'0|0|0|1': '-1'}})
    
    df_ex = pd.DataFrame(
        data=[
            [f1],
            [f2],
            [f3],
            [f4],
            [f5]
        ],
        columns=['sullivan']
    )
    
    df_ex['ip - 1'] = df_ex.apply(lambda row: (row.sullivan.p().i() - row.sullivan).reduce(), axis=1)
    df_ex['dh + hd'] = df_ex.apply(lambda row: (row.sullivan.d().h() + row.sullivan.h().d()).reduce(), axis=1)
    
    df_ex = df_ex.applymap(repr)
        
    print(df_ex)
