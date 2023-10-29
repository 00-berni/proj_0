import numpy as np
from visibility.stuff import get_data, interpole_three


def time_correction(year: float | int, **kargs) -> float:
    """Computing the difference between TD and UT

    For years between 1620 and 1998 the 
    difference is computed through 
    interpolation from a data table (see
    `/data/dT_data.csv`)
    
    For years out that range an 
    approximated formula is used, instead

    :param year: year for which computing TD-UT
    :type year: float | int
    
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
    """Class to store and manipolate time data

    The class takes value of time in seconds (or in [hh,mm,ss.ss])
    and the name of the kind of time to store them

    The attributes of the class are:

    :ivar val: the value of the angle in radiants
    :vartype val: float | np.ndarray | None
    :ivar tytime: the value of the angle in deg
    :vartype tytime: str    
    """

    #: the number of seconds in a day
    DAYSEC = 86400.

    @staticmethod
    def seconds(hh: int, mm: int, ss: float) -> float:
        """Function to convert time from [hh,mm,ss.ss] to seconds

        :param hh: hours
        :type hh: int
        :param mm: minutes
        :type mm: int
        :param ss: seconds
        :type ss: float
        
        :return: value in seconds only
        :rtype: float
        """
        return ss + mm*60 + hh*3600 
    
    @staticmethod
    def hms_form(time: float) -> list:
        """Function to convert time from seconds to [hh, mm (<60), ss.ss (<60)]

        :param time: value in seconds only
        :type time: float
        
        :return: value in [hh, mm (<60), ss.ss (<60)]
        :rtype: list
        """
        # getting seconds < 60
        # seconds are appoximated 
        ss = np.round(time % 60, 4)
        # getting minutes < 60
        mm = int((time // 60) % 60)
        # getting hours
        hh = int(time // 3600)
        # check
        if ss == 60: 
            ss -= 60
            mm += 1
        if mm == 60:
            mm -= 1
            hh += 1
        return [hh,mm,ss]
    
    @staticmethod
    def str_time(time: float, timetype: str, timezone: int = 0, dl_save: bool = False, sep: list[str] = [':',':','']) -> str:
        """Function to get a string to print time value

        The format of the string is `'hh:mm:ss.ssss [timetype]'`.

        If the time is local civil one the time zone is printed, 
        instead of time type, like `'hh:mm:ss.ssss ([time zone])'`.
        If it is present the information about daylight saving
        is also added: `'hh:mm:ss.ssss ([time zone], dl.)'`

        :param time: value in seconds only
        :type time: float
        :param timetype: kind of time 
        :type timetype: str
        :param timezone: the time zone, defaults to `0`
        :type timezone: int, optional
        :param dl_save: if it is `True` daylight saving is taken into account, defaults to `False`
        :type dl_save: bool, optional
        
        :return: the string with time value and informations
        :rtype: str
        """
        # getting time in [hh,mm,ss.ss] format
        hh,mm,ss = Time.hms_form(time)
        # adding a 0 before the number if it is < 10
        strzero = '0' if hh < 10 else ''
        # starting to built the string
        strtime = strzero + f'{hh}' + sep[0]
        # adding a 0 before the number if it is < 10
        strzero = '0' if mm < 10 else ''
        # updating the string
        strtime += strzero + f'{mm}' + sep[1]
        # adding a 0 before the number if it is < 10
        strzero = '0' if ss < 10 else ''
        # defing a method to print decimals of seconds only if they are
        form = lambda x: f'{x:.0f}' if x % 1 == 0 else f'{x:.4f}'
        # updating the string
        strtime += strzero + form(ss) + sep[2]
        # condition for local civil time
        if timezone != 0 or timetype == 'local':
            timetype = '('
            if timezone > 0: timetype += '+'
            timetype += f'{timezone}'
            # condition for daylight saving
            if dl_save: timetype += ', dl.s.'
            timetype += ')'
        return strtime + ' ' + timetype

    def __init__(self, value: list | float | np.ndarray | None = None, timetype: str = 'UT') -> None:
        """Constructor of the class

        Function takes a time value in seconds or in [hh,mm,ss.ss] and the kind of time and
        stores them.

        :param value: time value in seconds, defaults to `None`
        :type value: list | float | np.ndarray | None, optional
        :param timetype: kind of time, defaults to `'UT'`
        :type timetype: str, optional
        """
        # for None value returning a zero value
        if value is None: 
            value = 0.0
        # condition for time in [hh,mm,ss.ss]
        elif type(value) == list:
            # converting in seconds only
            value = Time.seconds(*value)
        # condition for array
        elif type(value) == np.ndarray:
            value = np.copy(value)
        self.val = value
        self.tytime = timetype

    def copy(self):
        """Function to make an exact copy of a `Time` object

        :return: the exact copy of the `Time` object
        :rtype: Time
        """
        time = self.val
        timetype = self.tytime
        return Time(time,timetype)
    
    def change_time_type(self, year: int | float | np.ndarray, tytime: str, out: str = 'sec'):
        """Converting UT in TD and vice versa.

        Function calls `time_correction()` function to converting time 
        from `self.tytime` to `tytime`. If they are the same nothing 
        happends.

        It is possible to choose the output of the function:

          * `out = 'sec'`: the computed value in seconds is returned
          * `out = 'Time'`: the computed value is returned as `Time` object

        :param year: year for which time correction is computed 
        :type year: int
        :param tytime: the wanted type of time
        :type tytime: str
        :param out: parameter to choose the output, defaults to `'sec'`
        :type out: str, optional
        
        :return: computed value either in seconds or in `Time` class object
        :rtype: float | np.ndarray | Time

        :raises Exception: only UT and TD are allowed as tipies
        """
        # getting the value
        time = self.copy().val
        # condition to make the correction
        if tytime != self.tytime:
            # computing the time correction
            # condition for array type
            if type(year) == np.ndarray:  
                dt = np.array([time_correction(yy) for yy in year])
            # condition for float
            else:
                dt = time_correction(year)
            # convertion
            if tytime == 'TD':   time += dt
            elif tytime == 'UT': time -= dt
            else: raise Exception('Error in time type!\nConversion is possible only for TD and UT, not for '+tytime)
        # condition for the output
        if out == 'sec':
            return time
        elif out == 'Time':
            return Time(time,timetype=tytime)

    def minute(self) -> float | np.ndarray:
        """Function to get time value in minutes.

        :return: value in minutes
        :rtype: float | np.ndarray
        """
        return self.val / 60

    def hour(self) -> float | np.ndarray:
        """Function to get time value in hours

        :return: value in hours
        :rtype: float | np.ndarray
        """
        return self.val / 3600

    def local_to_ut(self, timezone: int, dl_save: bool):
        """Function to convert local civil time in UT.

        Function takes information about time zone and
        daylight saving as input to return a `Time`
        object of time value in UT.

        :param timezone: time zone value
        :type timezone: int
        :param dl_save: if it is `True` daylight saving is taken into account
        :type dl_save: bool

        :return: time in UT
        :rtype: Time
        """
        hour = self.copy().hour()
        # condition for daylight saving
        if dl_save:
            hour -= 1
        hour -= timezone
        return Time(hour*3600,'UT')

    def print_time(self, enum: bool = True) -> list[str] | str:
        """Function to print time value.

        The parameter `enum` allows to print 
        the index for array values.

        :param enum: parameter for array, defaults to True
        :type enum: bool, optional
        
        :return: string (or list of string) with time value(s)
        :rtype: list[str] | str
        """
        time = self.copy()
        # condition to generalize the method for not-array type
        if type(time.val) != np.ndarray:
            time.val = np.array([time.val])
        sec = time.val
        # defing the list of strings
        str_list = []
        for i in range(len(sec)):
            str_res = ''
            if len(sec) > 1 and enum:
                # for array of values printing the index
                str_res += f'time {i}:\n'
            str_res += Time.str_time(sec[i], time.tytime)
            str_list += [str_res]
        # only a string for not-array type
        if len(sec) == 1:
            str_list = str_res
        return str_list

    def __add__(self, time):
        """Function to sum times

        If the adding value is not a `Time` object,
        it is considered a time in seconds.

        :param time: the value to sum
        :type time: Time | int | float | np.ndarray
        
        :return: the sum
        :rtype: Time
        """
        # condition for not-`Time` object
        if isinstance(time,(float,int,np.ndarray)):
            sumtime = self.val + time
        else:
            sumtime = self.val + time.val
        return Time(sumtime,self.tytime)
    
    def __sub__(self, time):
        """Function to subctraction times

        If the subctracting value is not a `Time` 
        object, it is considered a time in seconds.

        :param time: the value to subtract
        :type time: Time | int | float | np.ndarray
        
        :return: the subtraction
        :rtype: Time
        """
        # condition for not-`Time` object
        if isinstance(time,(float,int,np.ndarray)):
            subtime = self.val - time
        else:
            subtime = self.val - time.val
        return Time(subtime,self.tytime)
              
    def __mul__(self,val: float | int | np.ndarray):
        """Function to implement the time-number product

        :param val: a number
        :type val: float | int | np.ndarray

        :return: time-number product
        :rtype: Time
        """
        time = self.val * val
        return Time(time,timetype=self.tytime)

    def __neg__(self):
        return self * -1



class Date():
    """Class to store and handle dates and times

    Date, time, corresponding Julian day and calendar 
    type are stored in this class object.

    One can pass the Julian day or the date and time to the 
    constructor of the class; for the latter the stored 
    date is not the input value, but the computed one from 
    the corresponding Julian day. This procedure was 
    implemented to avoid possible typing errors in date or 
    time (such as wrong day number for a month, leap years,
    more than 24 hours for time).

    If date and time are local civil ones they are 
    converted in UT.

    The attributes of the class are:

    :ivar jd: Julian day
    :vartype jd: float | np.ndarray
    :ivar date: date value [year, month, day]
    :vartype date: list[int] | list[np.ndarray]
    :ivar time: time value
    :vartype time: Time
    :ivar calendar: calendar type
    :vartype calendar: str    
    """
    #: corresponding julian day for each epoch
    J2000 = 2451545.
    B1950 = 2433282.4235
    B1900 = 2415020.3135

    #: a dictionary to pass from the number of a month to its name
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
    def julian_day(year: int | np.ndarray, month: int | np.ndarray, day: float | int | np.ndarray, time: float | np.ndarray, calendar: str = 'Gregorian', MJD: bool = False) -> float | np.ndarray:
        """Computing the Julian day.

        One can compute the JD from a date of either the Gregorian or Julian calendar 
        and also the Modified JD (MJD) through the parameter `MJD`.


        :param year: year value
        :type year: int | np.ndarray
        :param month: month number
        :type month: int | np.ndarray
        :param day: day value
        :type day: float | int | np.ndarray
        :param time: time value in seconds
        :type time: float | np.ndarray
        :param calendar: calendar type, defaults to `'Gregorian'`
        :type calendar: str, optional
        :param MJD: if `True` the MJD is computed, defaults to `False`
        :type MJD: bool, optional
        
        :return: corresponding Julian day
        :rtype: float | np.ndarray
        
        :raises Exception: only Julian and Gregorian calendar are allowed
        """
        # condition to generalize the method for not-array type
        if type(year) != np.ndarray:
            year = np.array([year])
            month = np.array([month])
            day = np.array([day])
            time = np.array([time])
        
        # checking if there some dates that are
        # previous than 1582 October 15
        tmp = np.where(year < 1582)[0]
        if len(tmp) != 0:
            tmp = np.where(month[tmp] < 10)[0]
            if len(tmp) != 0:
                tmp = np.where(day[tmp] < 15)[0]
        # if there are, indecies are stored
        if len(tmp) != 0: 
            idx = np.delete(np.arange(len(year)),tmp)
        else:
            idx = slice(None)
        del tmp

        # converting time in decimal of day
        day = day.astype(float)
        day += time/Time.DAYSEC
        # condition for month number = 1 or 2
        year  = np.where(month <= 2, year-1, year)
        month = np.where(month <= 2, month+12, month)
        # computing JD
        JD = np.floor(365.25 * (year + 4716)) + np.floor(30.6001 * (month + 1)) + day - 1524.5
        # condition for Gregorian calendar
        if calendar == 'Gregorian':            
            A = np.floor(year[idx] / 100)
            B = (2 - A + np.floor(A / 4)) 
            JD[idx] += B
        # only Julian and Gregorian calendar are allowed
        elif calendar != 'Julian':
            raise Exception("!ERROR in chosen calendar!\nFunction accepts only the string:\n\t- 'Gregorian'\n\t- 'Julian'")
        # condition to compute MJD
        if MJD:
            JD -= 2400000.5
        # condition for not-array type
        if len(JD) == 1:
            JD = JD[0]
        return JD

    @staticmethod
    def calendar_date(JD: float | np.ndarray, timetype: str = 'UT') -> tuple[list[int], Time]:
        """Computing the calendar date from JD.

        :param JD: Julian day value
        :type JD: float | np.ndarray
        :param timetype: kind of time, defaults to `'UT'`
        :type timetype: str, optional
        
        :return: list date and time as `Time` object
        :rtype: tuple[list[int], Time]
        """
        # condition to generalize the method for not-array type
        if type(JD) != np.ndarray:
            JD = np.array([JD])
        JD = np.copy(JD)        
        JD += 0.5
        # taking the int and decimal parts
        F, Z = np.modf(JD)
        alpha = np.floor((Z-1867216.25)/36524.25)
        A = np.where(Z < 2299161, Z, Z + 1 + alpha - np.floor(alpha / 4))
        B = A + 1524
        C = np.floor((B-122.1)/365.25)
        D = np.floor(365.25 * C)
        E = np.floor((B-D) / 30.6001)
        # computing day
        day = B - D - np.floor(30.6001 * E) + F
        # computing month
        month = np.where(E < 14, E-1, E-13).astype(int)
        # computing year
        year = np.where(month > 2, C-4716, C-4715).astype(int)        
        # computing time from decimals of day
        time = (day % 1) 
        day -= time
        # condition for not-array type
        if len(year) == 1:
            year  = year[0]
            month = month[0]
            day   = day[0]
            time  = time[0]
        return [year,month,day], Time(time*Time.DAYSEC,timetype=timetype)

    @staticmethod
    def std_format(day: float, time: float) -> tuple[int,float]:
        """Function to check time and day values.

        Function checks that day has no decimals 
        and time value is not greater than the 
        number of seconds in a day and it 
        corrects them. 

        :param day: day
        :type day: float
        :param time: time in seconds
        :type time: float
        
        :return: day without decimals and updated time
        :rtype: tuple[int,float]
        """
        # condition for decimals of day
        if day % 1 != 0:        
            time += (day % 1)*86400
            day -= (day % 1)
        # condition for time value greater than seconds in a day
        if time >= Time.DAYSEC:
            day += time // Time.DAYSEC
            time %= Time.DAYSEC
            # time -= time // Time.DAYSEC * Time.DAYSEC
        return day, time

    @staticmethod
    def str_date(year: int, month: int, day: int | float) -> str:
        """Function to get a string to print a date

        The format of the output is `'[year] [month name] [day]'`
        
        :param year: year value
        :type year: int
        :param month: month number
        :type month: int
        :param day: day number
        :type day: int | float

        :return: the string of the date
        :rtype: str
        """
        strdate = f'{day:.0f} {Date.MONTHS[month]} {year:.0f}'
        return strdate



    def __init__(self, date: list[int | float | np.ndarray] | None = None, time: list | float | np.ndarray | Time | None = None, jd: float | np.ndarray | None = None, timetype: str = 'UT', timezone: int = 0, dl_save: bool = False, calendar: str = 'Gregorian') -> None:
        """Constructor of the class

        There are two possible entries:

          1. *date* and *time*: 
            
              function computes the corresponding JD, from which new date and time values
              are computed. This procedure avoid possible typing errors in the input values.

              The allowed format for date is a list [year, month, day], while for time is 
             `Time` object, the value in seconds or a list [hh, mm, ss.ss].

              If date and time values are civil local ones, they are converted in UT.
              Informations about time zone and daylight saving are necessary.

          2. *Julian day*:

              function computes the corrisponding date and time for the input JD.    

        :param date: list date, defaults to `None`
        :type date: list | None, optional
        :param time: time value, defaults to `None`
        :type time: list | float | np.ndarray | Time | None, optional
        :param jd: julian day, defaults to `None`
        :type jd: float | None, optional
        :param timetype: kind of time, defaults to `'UT'`
        :type timetype: str, optional
        :param timezone: only for local civil date, time zone value, defaults to `0`
        :type timezone: int, optional
        :param dl_save: only for local civil date, if `True` daylight saving is taken into account, defaults to `False`
        :type dl_save: bool, optional
        :param calendar: calendar type, defaults to `'Gregorian'`
        :type calendar: str, optional

        :raises Exception: one must pass either date and time or julian day
        """
        # condition for not-Time object
        if type(time) != Time:
            time = Time(time,timetype)
        # condition for Time object
        else:
            timetype = time.tytime
        # converting local civil date in UT
        if timetype == 'local' and (timezone != 0 or dl_save):
            time = time.local_to_ut(timezone,dl_save)
            timetype = 'UT'    
        # checking that one entry is given at least
        if (jd is None) and (date is None):
            raise Exception('Error in entries!\nYou have to pass either date and time or julian day at least')
        # if jd is not passed
        elif jd is None:
            # if time variable is an array referred to the same date
            if type(time.val) == np.ndarray and type(date[0]) != np.ndarray:
                date = [np.array([dd]*len(time.val)) for dd in date]
            # computing the jd
            jd = Date.julian_day(*date,time.val,calendar=calendar)
        # computing date and time from the jd
        date, time = Date.calendar_date(jd,timetype=timetype)
        
        self.jd = jd
        self.date = [*date]
        self.time = time.copy()
        self.calendar = calendar
        self.timezone = timezone
        self.dls = dl_save

    def copy(self) -> object:
        """Function to make an exact copy of the date

        :return: exact copy of the date 
        :rtype: Date
        """
        return Date(jd=self.jd,calendar=self.calendar,timetype=self.time.tytime,timezone=self.timezone,dl_save=self.dls)

    def leapyear(self) -> bool | np.ndarray:
        year = self.date[0]
        if self.calendar == 'Julian': 
            condition = (year % 4)==0
        elif self.calendar == 'Gregorian': 
            condition = ((year % 4)==0) * ((year % 400)==0)
        # if type(year) == np.ndarray:
            # print('cond',condition[:5])
            # print('len(year)',len(year))
        return condition

    def daydec(self) -> float | np.ndarray:
        """Function to get time in decimals of day

        :return: the day with decimals
        :rtype: float | np.ndarray
        """
        time = self.time.val
        return self.date[-1] + time/Time.DAYSEC
  
    def daynumber(self) -> float | np.ndarray:
        year, month, days = self.copy().date
        leap = self.leapyear()
        if type(month) != np.ndarray:
            year  = np.array([year])
            month = np.array([month])
            days  = np.array([days])
            leap  = np.array([leap])
        month = np.where(month > 2, month+1, month-1)
        leap = np.where(leap, [62]*len(leap), [63]*len(leap))
        # print('len(leap)',len(leap))
        # print('len(month)',len(month))
        daynum = np.where(month > 2, np.trunc(month*30.6) - leap, np.trunc((month*leap)/2))
        days = self.daydec()
        daynum += days
        if len(daynum) == 1: daynum = daynum[0]
        return daynum

    def yeardec(self) -> float | np.ndarray:
        """Function to get date in decimals of year

        Function converts month and day in
        decimals of year and adds them to
        the current year

        :return: the year with decimals
        :rtype: float | np.ndarray
        """
        date = self.copy()
        daynum = date.daynumber()
        year = date.date[0]
        # checking the calendar type
        if date.calendar == 'Julian':
            if type(year) == np.ndarray:
                dpy = np.where(self.leapyear(),366,365)
            else:
                dpy = 366 if self.leapyear() else 365
        elif date.calendar == 'Gregorian':
            dpy = 365.2425
        return year + daynum/dpy

    def change_time_type(self, tytime: str):
        """Function to convert UT to TD and vice versa

        The function is the same as 
        `Time.change_time_type()` function, but it 
        returns a `Date` object.

        :param tytime: the wanted type of time
        :type tytime: str
        
        :return: converted date
        :rtype: Date
        """
        date = self.copy()
        # getting the date (in decimals of year) for which correction will be estimated
        year = date.yeardec()
        time = date.time.change_time_type(year,tytime,out='Time')
        return Date(date.date,time,calendar=date.calendar)

    def ut_to_local(self, out: str = 'date') -> tuple[list[int],Time] | str:
        """Function to convert UT in local civil date and time
        
        One can choose to return the converted value in two possible forms
        through the parameter `out`: 
        
          * `out = 'date'`: a date list [year, month, day] and `Time` object 
          * `out = 'str`: a string to print
          * `out = 'all'`: as `'date'` but the value is printed

        :param timezone: time zone value
        :type timezone: int
        :param dl_save: if `True` daylight saving is taken into account
        :type dl_save: bool
        :param out: parameter to choose the output, defaults to 'date'
        :type out: str, optional

        :return: values as date list and time or a string with them
        :rtype: tuple[list[int],Time] | str
        """
        date = self.copy()
        timezone = date.timezone
        dl_save = date.dls
        time = date.time.hour()
        # conversion
        time += timezone
        if dl_save:
            time += 1
        # computing the local julian day
        ljd = Date.julian_day(*date.date,0.,calendar=date.calendar) + time/24
        # storing the results
        results = Date.calendar_date(ljd,timetype='local')
        if out == 'date':
            return results
        # making the string with values
        date, time = results
        if type(ljd) == np.ndarray:
            strdate = []
            for i in range(len(ljd)):
                strdate += [Date.str_date(date[0][i],date[1][i],date[2][i])+' '+ Time.str_time(time.val[i],'local',timezone,dl_save)]
        else:
            strdate = Date.str_date(*date)+' '+ Time.str_time(time.val,'local',timezone,dl_save)
        if out == 'str':
            return strdate
        elif out == 'all':
            if type(strdate) == list:
                for value in strdate:
                    print(value)
            else:
               print(strdate)
            return results

    def print_date(self, sel: str = 'all', local: bool = True) -> list[str] | str:
        """Function to print date values


        One can select to print both date and time or just one 
        through the `sel` parameter:

            * `sel = 'date'`: print only the date
            * `sel = 'time'`: print only the time
            * `sel = 'all'`: print both
            * `sel = 'jd'`: print only the corresponding JD        

        :param sel: selection parameter, defaults to `'all'`
        :type sel: str, optional

        :return: string (or list of strings) with date value(s)
        :rtype: list[str] | str
        """
        # only julian day
        if sel == 'jd':
            jd = self.jd 
            if type(jd) != np.ndarray:
                return f'{jd}'
            else: 
                return [f'{jdi}' for jdi in jd]
        # date and time
        else:
            date = self.copy()
            year, month, day = date.date
            time = date.time
            # condition to generalize the method for not-array type
            if type(year) != np.ndarray:
                year = np.array([year])
                month = np.array([month])
                day = np.array([day])
                time.val = np.array([time.val])
            # this is the parameter for time enumeration
            # it is initialiazed to `False`
            enum = False
            # defining list for string
            str_list = []
            for i in range(len(year)):
                day[i], time.val[i] = Date.std_format(day[i], time.val[i])
                # condition for only time value
                if sel == 'time':
                    enum = True
                else:
                    str_res = ''
                    if sel == 'all' and len(year) > 1:
                        # for array indecies are printed
                        str_res += f'Date {i}:\n'
                    str_res += Date.str_date(year[i],month[i],day[i])
                    str_list += [str_res]
            # condition for only date
            if sel == 'date':
                # only a string for not-array type
                if len(year) == 1:
                    str_list = str_res
                return str_list
            else:
                str_time = time.print_time(enum=enum)
                # condition for only time value
                if sel == 'time':
                    return str_time
                elif sel == 'all':
                    local = local if (date.timezone != 0 or date.dls) else False 
                    if local:
                        localtime = self.ut_to_local('str')
                    # condition for list of string
                    if type(str_time) == list:
                        for i in range(len(str_time)):
                            str_list[i] = str_list[i] + ' ' + str_time[i]
                            if local: str_list[i] += ' ==> ' + localtime[i]
                                
                    # condition for string
                    else:
                        if len(year) == 1:
                            str_list = str_res
                        str_list += ' ' + str_time
                        if local: str_list += ' ==> ' + localtime
                    return str_list

    def __add__(self,day):
        """Funtion to sum value to date

        If the value is a number, it is added to the 
        Julian day of the Date.
        If it is a `Date` or `Time` object, the 
        time part of each one are summed.

        :param day: value to sum
        :type day: float | int | np.ndarray | Date | Time
        
        :return: the sum as a `Date` object
        :rtype: Date
        """
        # condition for number
        if isinstance(day,(float,int,np.ndarray)):
            # summing to julian day
            jd = self.jd + day
            return Date(jd=jd,timetype=self.time.tytime,calendar=self.calendar)
        # condition for `Date` or `Time` object
        else:
            # condition for `Date` object
            if type(day) == Date:
                # taking time part
                day = day.time.copy()
            sumdate = self.time + day
            return Date(self.date,sumdate,calendar=self.calendar) 

    def __sub__(self,day):
        """Funtion to subtract value to date

        If the value is a number, it is 
        subctracted to the Julian day of the Date.
        If it is a `Date` or `Time` object, the 
        time part of each one are subtracted.

        :param day: value to subtract
        :type day: float | int | np.ndarray | Date | Time
        
        :return: the sum as a `Date` object
        :rtype: Date
        """
        # condition for number
        if isinstance(day,(float,int,np.ndarray)):
            # summing to julian day
            jd = self.jd - day
            return Date(jd=jd,timetype=self.time.tytime,calendar=self.calendar)
        # condition for `Date` or `Time` object
        else:
            # condition for `Date` object
            if type(day) == Date:
                # taking time part
                day = day.time.copy()
            subdate = self.time - day
            return Date(self.date,subdate,calendar=self.calendar) 

    def __mul__(self, val: float | int | np.ndarray):
        """Function to implement the date-number product

        :param val: a number
        :type val: float | int | np.ndarray
        :return: date-number product 
        :rtype: Date
        """
        jd = self.jd*val
        return Date(jd=jd,timetype=self.time.tytime,calendar=self.calendar)
