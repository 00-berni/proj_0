## required packages
import os.path as ph
import numpy as np

PROJECT_FOLDER = ph.split(ph.dirname(ph.realpath(__file__)))[0]
DATA_FOLDER = ph.join(PROJECT_FOLDER,'data')


def get_data(filename: str, unpack: bool = True, delimiter: str = ',', message: bool = False) -> np.ndarray:
    DATA_FILE = ph.join(DATA_FOLDER,filename)
    if message:
        print('> Open file: ' + DATA_FILE)
    ext = filename[-3:]
    if ext == 'csv':
        from pandas import read_csv
        data = read_csv(DATA_FILE, delimiter=delimiter).to_numpy().transpose()
    elif ext == 'txt':
        data = np.loadtxt(DATA_FILE,unpack=unpack)
    else: raise Exception(f'.{ext} is not allowed!\nOnly .txt or .csv files')
    return data

def interpole_three(values: list | np.ndarray, centre: float | None = None, xvalues: list | np.ndarray | None = None, val: str | None = None,**kargs) -> float | tuple[float, dict]:
    """_summary_

    :param values: _description_
    :type values: list | np.ndarray
    :param centre: _description_, defaults to None
    :type centre: float | None, optional
    :param xvalues: _description_, defaults to None
    :type xvalues: list | np.ndarray | None, optional
    :param val: _description_, defaults to None
    :type val: str | None, optional
    :raises Exception: _description_
    :raises Exception: _description_
    :raises Exception: _description_
    :raises Exception: _description_
    :return: _description_
    :rtype: float | tuple[float, dict]
    """
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