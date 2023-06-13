import numpy as np
from test_struc import *


def interpole_three(values: list | np.ndarray, 
                    n: float, 
                    xvalues: list | np.ndarray = [], 
                    centre: float | None = None, 
                    val: str | None = None
                    ) -> float | tuple[float, dict]:
    dim = len(values)
    sample = np.copy(values)
    if dim > 3:
        if centre == None: raise Exception('!Error in central value!\n You have to pass the central value for interpolation')
        if len(xvalues) == 0: raise Exception('!Error in `xvalues`!\n You have to pass x values for interpolation')
        diff = np.diff(sample)
        for i in range(2):
            diff = np.diff(diff)
        maxdiff = diff.max()
        if np.abs(maxdiff) > 2: raise Exception('!Error in interpolation!\nThe third difference is not negligible')
        else:
            idx = (np.abs(xvalues - centre)).argmin()
            sample = sample[idx-1:idx+2] 
    elif dim < 3: raise Exception('!Error in length!\nInterpolation needs three points at least')
    
    centre_val = sample[1]
    a, b = np.diff(sample)
    c = b - a
    result = centre_val + n * (a + b + n*c) / 2
    
    if val != None:
        dval = {}

        if val == 'ym' or val == 'all':
            dval['ym'] = centre_val - (a+b)**2/(8*c)
            dval['nm'] = - (a+b) / (2*c)
        if val == 'n0' or val == 'all':
            n0 = 0
            cnt = 0
            while True:
                dn0 = -(2*centre_val + n0*(a+b+c*n0))/(a+b+2*c*n0)
                n0 += dn0
                cnt += 1
                if np.abs(dn0) < 1e-6:
                    break
                if cnt > 20: 
                    print('Algorithm did not converge for n0!')
                    break
            dval['n0'] = n0 if cnt <= 20 else None    

        result = (result, dval)

    return result    

if __name__ == '__main__':
    
    starting_test('TEST FOR INTERPOLATION')
    
    try:
        x,y,results = get_test_data('test_interpole-01_data.csv')
        print('> Find distance Mars-Earth at 1992 November 8 at 4h21m TD')
        y0 = results[0]

        n = 4.35/24
        centre = 8 + n

        yint = interpole_three(y,n,x,centre)

        if np.abs(yint - y0) > 1e-5:        
            print('DATA01: NOT MATCHED!')
            print(f'  Result from interpolation:\t{yint}')
        else: 
            print('DATA01: OK!')
        
        x,y,results = get_test_data('test_interpole-02_data.csv')
        print('> Find the extremum value of distance Sun-Mars for data in May 1992')
        ym0, nm0 = results[:2]

        n = 0.25 

        yint, dvals = interpole_three(y,n,val='ym')
        ym = dvals['ym']
        nm = dvals['nm']
        if np.abs(ym - ym0) > 1e-5 and np.abs(nm - nm0) > 1e-4:        
            print('DATA02: NOT MATCHED!')
            print(f'  Result from interpolation:\t{ym} and {nm}')
        else: 
            print('DATA02: OK!')
        
        x,y,results = get_test_data('test_interpole-03_data.csv')
        print('> Find n at which function becomes zero')
        n00 = results[0]

        n = 1

        yint, dvals = interpole_three(y,n,val='n0')

        n0 = dvals['n0']
        
        if np.abs(n0 - n00) > 1e-5:        
            print('DATA03: NOT MATCHED!')
            print(f'  Result from interpolation:\t{n0}')
        else: 
            print('DATA03: OK!')

        ending_test()
    except:
        test_error()