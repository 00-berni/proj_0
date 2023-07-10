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

    DAYSEC = 86400.

    @staticmethod
    def seconds(hh: int, mm: int, ss: float) -> float:
        return ss + mm*60 + hh*3600 
    
    @staticmethod
    def hms_form(time: float) -> list:
        ss = np.round(time % 60, 4)
        mm = int((time // 60) % 60)
        hh = int(time // 3600)
        if ss == 60: 
            ss -= 60
            mm += 1
        if mm == 60:
            mm -= 1
            hh += 1
        return [hh,mm,ss]
    
    @staticmethod
    def str_time(time: float, timetype: str) -> str:
        hh,mm,ss = Time.hms_form(time)
        strtime = f'{hh}h'
        if ss != 0.:
            form = lambda x: f'{x:.0f}' if x % 1 == 0 else f'{x:.4f}'
            strtime += f'{mm}m'+form(ss)+'s'
        elif mm != 0.:
            strtime += f'{mm}m'
        return strtime + ' ' + timetype

    def __init__(self, value: list | float | np.ndarray | None = None, timetype: str = 'UT') -> None:
        if value is None: 
            value = 0.0
        elif type(value) == list:
            value = Time.seconds(*value)
        elif type(value) == np.ndarray:
            value = np.copy(value)
        self.val = value
        self.tytime = timetype

    def copy(self):
        return Time(self.val,self.tytime)
    
    def change_time_type(self, year: int, tytime: str, out: str = 'sec'):
        time = self.copy().val
        if tytime != self.tytime:
            if type(year) == np.ndarray:  
                dt = np.array([time_correction(yy) for yy in year])
            else:
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

    def print_time(self, enum: bool = True) -> list[str] | str:
        time = self.copy()
        if type(time.val) != np.ndarray:
            time.val = np.array([time.val])
        sec = time.val
        str_list = []
        for i in range(len(sec)):
            str_res = ''
            if len(sec) > 1 and enum:
                str_res += f'time {i}:\n'
            str_res += Time.str_time(sec[i], time.tytime)
            str_list += [str_res]
        if len(sec) == 1:
            str_list = str_res
        return str_list
            
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

    @staticmethod
    def julian_day(year: int | np.ndarray, month: int | np.ndarray, day: float | int | np.ndarray, time: float | np.ndarray, calendar: str = 'Gregorian', MJD: bool = False) -> float:
        """Computing the Julian Day

        One can compute the JD from a date of Gregorian or Julian calendar and also the 
        Modified JD

        From J.,Meeus, *Astronomical Algorithms*, pg. 61. 
        Follow the same notation.

        Same from M., Gallaway, *An Introduction to Observational Astrophysics*, 
        pg. 41.

        Different in *Explatatory supplement to the Astronomical Almanac*, pg. 604

        :param year: year
        :type year: int
        :param month: month
        :type month: int
        :param day: day 
        :type day: float
        :param calendar: kind of chosen calendar, defaults to 'Gregorian'
        :type calendar: str, optional
        :param MJD: set `True` to have the Modified JD, defaults to False
        :type MJD: bool, optional

        :return: Julian Day
        :rtype: float

        :raise: 
        """
        if type(year) != np.ndarray:
            year = np.array([year])
            month = np.array([month])
            day = np.array([day])
            time = np.array([time])
        day = day.astype(float)
        day += time/Time.DAYSEC
        year  = np.where(month <= 2, year-1, year)
        month = np.where(month <= 2, month+12, month)

        JD = np.floor(365.25 * (year + 4716)) + np.floor(30.6001 * (month + 1)) + day - 1524.5

        if calendar == 'Gregorian':
            A = np.floor(year / 100)
            B = (2 - A + np.floor(A / 4)) 
            JD += B
        elif calendar != 'Julian':
            raise Exception("!ERROR in chosen calendar!\nFunction accepts only the string:\n\t- 'Gregorian'\n\t- 'Julian'")
        
        if MJD:
            JD -= 2400000.5

        if len(JD) == 1:
            JD = JD[0]
        return JD

    @staticmethod
    def calendar_date(JD: float | np.ndarray, timetype: str = 'UT', epoch: str ='J2000.0'):
        """Computing the calendar date from JD

        Different in *Explatatory supplement to the Astronomical Almanac*, pg. 605  

        :param JD: Julian Day
        :type JD: float

        :return: calendar date in [year, month, day.decimal]
        :rtype: list
        """
        if type(JD) != np.ndarray:
            JD = np.array([JD])
        
        JD = np.copy(JD)        
        JD += 0.5
        F, Z = np.modf(JD)
        alpha = np.floor((Z-1867216.25)/36524.25)
        A = np.where(Z < 2299161, Z, Z + 1 + alpha - np.floor(alpha / 4))
        B = A + 1524
        C = np.floor((B-122.1)/365.25)
        D = np.floor(365.25 * C)
        E = np.floor((B-D) / 30.6001)
        
        day = B - D - np.floor(30.6001 * E) + F
        month = np.where(E < 14, E-1, E-13).astype(int)
        year = np.where(month > 2, C-4716, C-4715).astype(int)        
        time = (day % 1) 
        day -= time
        if len(year) == 1:
            year  = year[0]
            month = month[0]
            day   = day[0]
            time  = time[0]
        return [year,month,day], Time(time*Time.DAYSEC,timetype=timetype)

    @staticmethod
    def std_format(day: float, time: float) -> str:
        if day % 1 != 0:        
            time += (day % 1)*86400
            day -= (day % 1)
        if time > Time.DAYSEC:
            day += time // Time.DAYSEC
            time -= time // Time.DAYSEC * Time.DAYSEC
        return day, time

    @staticmethod
    def str_date(year: int, month: int, day: int | float) -> str:
        strdate = f'{year:.0f} {Date.MONTHS[month]} {day:.0f}'
        return strdate

    def __init__(self, date: list | None = None, time: list | float | np.ndarray | Time | None = None, jd: float | None = None, timetype: str = 'UT', calendar: str = 'Gregorian', epoch: str = 'J2000.0') -> None:
        if type(time) != Time:
            time = Time(time,timetype)
        
        if (jd is None) and (date is None):
            raise Exception('error')
        elif jd is None:
            if type(time.val) == np.ndarray and type(date[0]) != np.ndarray:
                date = [np.array([dd]*len(time.val)) for dd in date]
            jd = Date.julian_day(*date,time.val,calendar=calendar)
        date, time = Date.calendar_date(jd,timetype=timetype)
        
        self.jd = jd
        self.date = [*date]
        self.time = time.copy()
        self.calendar = calendar
        self.epoch = epoch

    def copy(self):
        return Date(jd=self.jd,calendar=self.calendar,timetype=self.time.tytime,epoch=self.epoch)

    def daydec(self):
        time = self.time.val
        return self.date[-1] + time/Time.DAYSEC
  
    def change_time_type(self, tytime: str):
        date = self.copy()
        year = date.date[0]
        time = date.time.change_time_type(year,tytime,out='Time')
        return Date(date.date,time)

    def print_date(self, sel: str = 'all', eph: bool = False) -> list[str] | str:
        if sel == 'jd':
            jd = self.jd 
            if type(jd) != np.ndarray:
                return f'{jd}'
            else: 
                return [f'{jdi}' for jdi in jd]
        else:
            date = self.copy()
            year, month, day = date.date
            time = date.time
            if type(year) != np.ndarray:
                year = np.array([year])
                month = np.array([month])
                day = np.array([day])
                time.val = np.array([time.val])
            enum = False
            str_list = []
            for i in range(len(year)):
                day[i], time.val[i] = Date.std_format(day[i], time.val[i])
                if sel == 'time':
                    enum = True
                else:
                    str_res = ''
                    if sel == 'all' and len(year) > 1:
                        str_res += f'Date {i}:\n'
                    str_res += Date.str_date(year[i],month[i],day[i])
                    str_list += [str_res]
            if sel == 'date':
                if len(year) == 1:
                    str_list = str_res
                return str_list
            else:
                str_time = time.print_time(enum=enum)
                if sel == 'time':
                    return str_time
                elif sel == 'all':
                    if type(str_time) == list:
                        for i in range(len(str_time)):
                            str_list[i] = str_list[i] + ' ' + str_time[i]
                            if eph:
                                str_list[i] = str_list[i] + ' ' + date.epoch
                    else:
                        if len(year) == 1:
                            str_list = str_res
                        str_list += ' ' + str_time
                        if eph:
                            str_list += ' ' + date.epoch
                    return str_list

    def __add__(self,day):
        if isinstance(day,(float,int)):
            jd = self.jd + day
            return Date(jd=jd,timetype=self.time.tytime,calendar=self.calendar,epoch=self.epoch)

    def __sub__(self,day):
        if isinstance(day,(float,int)):
            jd = self.jd - day
            return Date(jd=jd,timetype=self.time.tytime,calendar=self.calendar,epoch=self.epoch)

if __name__ == '__main__':
    print()
    date = Date([1957,10,4.81])
    print(f"Julian Day of {date.print_date()}:\t{date.print_date(sel='jd')}")
    
    print()
    start = Date.J2000
    jd = np.linspace(start, start+1000, 50)
    date = Date(jd=jd)
    for ele in date.print_date()[:5]:
        print(ele)