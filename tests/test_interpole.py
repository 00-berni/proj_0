import numpy as np
from test_struc import *


def interpole_three(values: list | np.ndarray, centre: float | None = None, xvalues: list | np.ndarray | None = None, val: str | None = None,**kargs) -> float | tuple[float, dict]:
    if 'warning' in kargs.keys():
        warning = kargs['warning']
    else:
        warning = True
    # storing dimension and values
    dim = len(values)
    sample = np.copy(values)
    # condition for more than three values
    if dim > 3:
        # central value is necessary
        if centre == None: raise Exception('Error in central value!\n You have to pass the central value for interpolation')
        # also x values are necessary
        if len(xvalues) == 0: raise Exception('Error in `xval`!\n You have to pass x values for interpolation')
        # checking the magnitude of the third difference
        diff = np.diff(np.diff(np.diff(sample)))
        maxdiff = diff.max()
        if np.abs(maxdiff) > 3: 
            raise Exception(f'Error in interpolation!\nThe third difference is not negligible: {maxdiff}')
        # computing the index of the central value
        idx = (np.abs(xvalues - centre)).argmin()
        # correction for extrema
        if idx == dim-1: idx -= 1
        elif idx == 0: idx += 1
        # cutting the sample to have only three
        # centered in central value
        sample = sample[idx-1:idx+2]
        xvalues = xvalues[idx-1:idx+2]
    # checking at least three values are passed
    elif dim < 3: raise Exception('!Error in length!\nInterpolation needs three points at least')
    # getting the central values
    centre_val = sample[1]
    # computing differences
    a, b = np.diff(sample)
    c = b - a
    if centre != None:
        x2 = xvalues[1]
        # computing the interpolation factor
        n = (centre - x2)/np.abs(x2-xvalues[0])
        # computing the interpolated value
        result = centre_val + n/2 * (a + b + n*c) 
    else:
        if warning: print(': No interpolation was done: parameter `centre` set to `None`')
        result = None

    # condition to compute other quantities
    if val != None:
        # defining the dictionary for additional values
        dval = {}
        # condition to compute extrema through iterpolation
        if val == 'ym' or val == 'all':
            dval['ym'] = centre_val - (a+b)**2/(8*c)
            dval['nm'] = - (a+b) / (2*c)
        # condition to compute where function is null through iterpolation
        if val == 'n0' or val == 'all':
            # recursive method  
            n0 = 0
            cnt = 0     # control variable
            while True:
                dn0 = -(2*centre_val + n0*(a+b+c*n0))/(a+b+2*c*n0)
                n0 += dn0
                cnt += 1
                if np.abs(dn0) < 1e-6:
                    break
                # control condition to prevent a Stack Overflow error
                if cnt > 20 and warning:
                    n0 = None 
                    print('Algorithm did not converge for n0!')
                    break
            dval['n0'] = n0 
        # appending the dictionary to the result
        result = (result, dval)
    return result    

if __name__ == '__main__':
    
    starting_test('TEST FOR INTERPOLATION')
    
    try:
        x,y,results = get_test_data('test_interpole-01_data.csv')
        print('> Find distance Mars-Earth at 1992 November 8 at 4h21m TD')
        y0 = results[0]

        centre = 8 + (4+21/60)/24

        yint = interpole_three(y,centre=centre,xvalues=x)

        if np.abs(yint - y0) > 1e-5:        
            print('DATA01: NOT MATCHED!')
            print(f'  Result from interpolation:\t{yint}')
        else: 
            print('DATA01: OK!')
        
        x,y,results = get_test_data('test_interpole-02_data.csv')
        print('> Find the extremum value of distance Sun-Mars for data in May 1992')
        ym0, nm0 = results[:2]


        yint, dvals = interpole_three(y,val='ym',warning=False)
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

        yint, dvals = interpole_three(y,val='n0',warning=False)

        n0 = dvals['n0']
        
        if np.abs(n0 - n00) > 1e-5:        
            print('DATA03: NOT MATCHED!')
            print(f'  Result from interpolation:\t{n0}')
        else: 
            print('DATA03: OK!')

        ending_test()
    except:
        test_error()