import numpy as np
from visibility.stuff import get_data, interpole_three



def time_correction(year: float,**kargs) -> float:
    """Computing the difference between TD and UT

    For years between 1620 and 1998 the 
    difference is computed through 
    interpolation from a data table (see
    `/data/dT_data.csv`)
    
    For years out that range an 
    approximated formula is used, instead

    :param year: year for which computing TD-UT
    :type year: float
    
    :return: TD - UT
    :rtype: float
    """
    # interpolation from data table
    if 1620 <= year <= 1998:
        # name of data file
        filename = 'dT_data.csv'
        # extracting data
        y0, dt0 = get_data(filename)
        # value is already in data
        if year in y0:
            dt = dt0[year == y0][0]
        # computing interpolation
        else:
            dt = interpole_three(dt0,centre=year,xvalues=y0,kargs=kargs)
    # approximated formula
    else:        
        t = (year - 2000)/100
        if year < 948:
            dt = 2177 + 497*t + 44.1*t**2
        else:
            dt =  102 + 102*t + 25.3*t**2
        if 2000 <= year <= 2100: dt += 0.37*(year-2100)
    return dt

class Time():

    @staticmethod
    def seconds(hh,mm,ss) -> float:
        return ss + mm*60 + hh*3600 

    def __init__(self, value: list | float | np.ndarray | None = None, timetype: str = 'UT') -> None:
        if value is None: 
            value = 0.
        elif type(value) == list:
            value = Time.seconds(*value)
        elif type(value) == np.ndarray:
            value = np.copy(value)
        self.val = value
        self.tytime = timetype
    
    def change_time_type(self, year: int, tytime: str, out: str = 'sec'):
        if type(self.val) == np.ndarray:
            time = np.copy(self.val)
        else:
            time = self.val
        dt = time_correction(year)
        if tytime == 'TD':   time += dt
        elif tytime == 'UT': time -= dt
        if out == 'sec':
            return time
        elif out == 'Time':
            return Time(time,timetype=tytime)

    def minute(self):
        return self.val / 60

    def hour(self):
        return self.val / 3600

    def hms_form(self) -> list:
        time = self.val
        ss = round(time % 60,4)
        mm = int((time // 60) % 60)
        hh = int(time // 3600)
        return [hh,mm,ss]

    def str_time(self) -> str:
        hh,mm,ss = self.hms_form()
        strtime = f'{hh}h'
        if ss != 0.:
            form = lambda x: f'{x:.0f}' if x % 1 == 0 else f'{x:.4f}'
            strtime += f'{mm}m'+form(ss)+'s'
        elif mm != 0.:
            strtime += f'{mm}m'
        tytime = self.tytime
        return strtime + ' ' + tytime

    def copy(self):
        return Time(self.val,self.tytime)

    def __add__(self, time):
        sumtime = self.val + time.val
        return Time(sumtime,self.tytime)
    
    def __sub__(self, time):
        subtime = self.val - time.val
        return Time(subtime,self.tytime)
              


class Date():

    # epochs
    J2000 = 2451545.
    B1950 = 2433282.4235
    B1900 = 2415020.3135



    MONTHS = {
        1 : 'January',
        2 : 'February',
        3 : 'March',
        4 : 'April',
        5 : 'May',
        6 : 'June',
        7 : 'July',
        8 : 'August',
        9 : 'September',
        10 : 'October',
        11 : 'November',
        12 : 'December'
    }


    def __init__(self, date: list , time: list | float | np.ndarray | Time | None = None, timetype: str = 'UT', epoch: str = 'J2000.0') -> None:
        if type(time) == Time:
            timetype = time.tytime
            time = time.val
        self.date = [*date]
        self.time = Time(time,timetype)
        self.epoch = epoch

    def daydec(self):
        time = self.time.val
        return self.date[-1] + time/86400

    def std_format(self):
        year, month, day = self.date 
        time = self.time.val
        if self.date[-1] % 1 != 0:        
            time += (day % 1)*86400
            day -= (day % 1)
        if time > 86400:
            day += time // 86400
            time -= time // 86400 * 86400
        return Date([year,month,day],time,self.time.tytime,epoch=self.epoch)

    def str_date(self, sel: str = 'all', epoch: bool = False) -> str:
        date = self.std_format()
        year, month, day = date.date
        time = date.time
        strdate = f'{year:.0f} {Date.MONTHS[month]} {day:.0f}'
        strtime = time.str_time()

        if sel == 'date': strresult = strdate
        elif sel == 'time': strresult = strtime
        elif sel == 'all': strresult = strdate + ' ' + strtime

        if epoch: strresult += ' ' + self.epoch
        return strresult
    
    def copy(self):
        return Date([*self.date],self.time.val,self.time.tytime,self.epoch)
