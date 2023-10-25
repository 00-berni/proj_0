## required packages
import os.path as ph
import numpy as np

PROJECT_FOLDER = ph.split(ph.dirname(ph.realpath(__file__)))[0]
DATA_FOLDER = ph.join(PROJECT_FOLDER,'data')
RESULTS_FOLDER = ph.join(PROJECT_FOLDER,'results')

def get_data(filename: str, rows: int | np.ndarray | slice = slice(None), cols: int | np.ndarray | slice = slice(None), unpack: bool = True, delimiter: str = ',', message: bool = False) -> np.ndarray:
    """Extracting data from `.csv` or `.txt` files.

    Function takes the name of a data file in `/data` directory and extracts data.

    :param filename: name of the data filter
    :type filename: str
    :param unpack: parameter of `np.loadtxt()`, defaults to `True`
    :type unpack: bool, optional
    :param delimiter: only for `.csv` files, the delimiter between data, defaults to `','`
    :type delimiter: str, optional
    :param message: control parameter, if it is `True` update message will be printed, defaults to `False`
    :type message: bool, optional
    :raises Exception: only `.txt` and `.csv` files are allowed
    
    :return: extracted data
    :rtype: np.ndarray
    """
    # file path
    DATA_FILE = ph.join(DATA_FOLDER,filename)
    if message:
        print('> Open file: ' + DATA_FILE)
    # getting extension of the file 
    ext = filename[-3:]
    if ext == 'csv':
        from pandas import read_csv
        data = read_csv(DATA_FILE, delimiter=delimiter).to_numpy().transpose()
    elif ext == 'txt':
        data = np.loadtxt(DATA_FILE,unpack=unpack)
    else: raise Exception(f'.{ext} is not allowed!\nOnly .txt or .csv files')
    data = data[rows,cols]
    return data

def import_data(filename: str, sel: int | np.ndarray | slice = slice(None), delimiter: str =',', sep: str = ':') -> tuple[np.ndarray]:
    """Extracting target data.

    Functions extracts data for target object(s) as name (`names`), right ascension (`ras`), 
    declination (`decs`) and proper motion (`muas`, `muds`).

    For multiple objects it is possible to select just some, for which data are extracted, 
    through `sel` parameter.

    The standard format for angles (times) is a string:
    `'degrees:primes:seconds.seconds'` (`'hh:mm:ss.ss'`)

      The parameter `sep` is the character separates different units (in this case `':'`) 

    :param filename: name of the data file
    :type filename: str
    :param sel: parameter to select object for which data are extracted, defaults to `slice(None)`
    :type sel: int | np.ndarray | slice, optional
    :param delimiter: only for `.csv` files, the delimiter between data, defaults to `','`
    :type delimiter: str, optional
    :param sep: character to separate different unit, defaults to `':'`
    :type sep: str, optional
    
    :return: extracted data
    :rtype: tuple[np.ndarray]
    """
    def extract_values(data_val: str, valtype: str = 'ang') -> tuple[tuple[list[str | list[float]]]]:
        if valtype == 'ang':
            return [data_val[0], [float(val) for val in data_val[1:].split(sep)]]
        elif valtype == 'time':
            return [float(val) for val in data_val.split(sep)]
        else: raise Exception(f'Error in `valtype`! `{valtype}` is not allowed, only `ang` or `time`')

    if type(sel) == int: sel = slice(sel,sel+1) 
    # extracting data
    names, ras, decs, muas, muds, epochs, obs_names, lats, lons, hs, obs_dates, obs_times = get_data(filename, cols=sel, delimiter=delimiter)
    
    # defining an empty list for proper motion data
    prmts = []
    for i in range(len(names)):
        # object info
        alpha = ras[i]
        delta = decs[i]
        mua = muas[i]
        mud = muds[i]

        # from `'+hh:mm:ss.ss'` to `['+', [hh, mm, ss]]` 
        ras[i] = extract_values(alpha)
        # from `'+dd:pp:ss.ss'` to `['+', [dd, pp, ss]]` 
        decs[i] = extract_values(delta)
        
        mua = float(mua) if mua != 'None' else None
        mud = float(mud) if mud != 'None' else None
        prmts += [[mua,mud]]
        

        # observatory info
        lat = lats[i]
        lon = lons[i]
        hs[i] = float(hs[i])

        lat = '+' + lat[1:] if lat[0] == 'N' else '-' + lat[1:]
        lon = '+' + lon[1:] if lon[0] == 'W' else '-' + lon[1:]
        lats[i] = extract_values(lat)
        lons[i] = extract_values(lon)

        # observation time info
        obs_date = obs_dates[i]
        obs_time = obs_times[i]  

        obs_dates[i] = extract_values(obs_date, valtype='time')
        obs_times[i] = extract_values(obs_time, valtype='time')

    obj  = (names, ras, decs, prmts, epochs)
    obs  = (obs_names, lats, lons, hs)
    date = (obs_dates, obs_times) 

    return obj, obs, date

def interpole_three(values: list | np.ndarray, centre: float | None = None, xvalues: list | np.ndarray | None = None, val: str | None = None, **kargs) -> float | tuple[float, dict]:
    """Computing interpolated value from three data points.

    If more than three points are passed, function computes third difference(s) in order to 
    check whether interpolation is possible or not.

    It possible to interpolate also extrema or zero points through the parameter `val`:

      * `val = None`: no additional result is returned
      * `val = 'ym'`: the extrema values and the corresponding interpolation parameter are computed
      * `val = 'n0'`: the interpolation parameter for zero point is computed
      * `val = 'all'`: both are computed

    For `val is not None` a dictionary with selected additional values are returned, in addition to the interpolated one. 

    To disable the warning messages one can use the parameter `warning = False` of kargs.

    :param values: y values from which interpolating
    :type values: list | np.ndarray
    :param centre: x value for which interpolating, defaults to `None`
    :type centre: float | None, optional
    :param xvalues: x values corrisponding to y values, defaults to `None`
    :type xvalues: list | np.ndarray | None, optional
    :param val: for additional results, defaults to `None`
    :type val: str | None, optional
    :param kargs: additional parameters
    :type kargs: Any 
    
    :return: the interpolated value and if `val is not None` selected additional results
    :rtype: float | tuple[float, dict]
    """
    # condition for warning messages
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
        # finding the index of the x value that is nearest to the `centre`
        idx = (np.abs(xvalues - centre)).argmin()
        # check for edges
        if idx == dim-1: idx -= 1
        elif idx == 0: idx += 1
        # cutting the sample to have only three elements
        # centered in central value
        sample = sample[idx-1:idx+2]
        xvalues = xvalues[idx-1:idx+2]
    # checking at least three values are passed
    elif dim < 3: 
        raise Exception('!Error in length!\nInterpolation needs three points at least')
    # getting the central value
    centre_val = sample[1]
    # computing differences
    a, b = np.diff(sample)
    c = b - a
    # central value is necessary
    if centre != None:
        # getting the x central value
        x2 = xvalues[1]
        # computing the interpolation factor
        n = (centre - x2)/np.abs(x2-xvalues[0])
        # computing the interpolated value
        result = centre_val + n/2 * (a + b + n*c) 
    else:
        if warning: 
            print(': No interpolation was done: parameter `centre` set to `None`')
        result = None

    # condition to compute additional quantities
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
                # condition to stop the loop
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
