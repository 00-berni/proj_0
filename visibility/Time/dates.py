import numpy as np
from visibility.Time.Tclasses import *
from visibility.Angles import Angles, HAngles, ArrAngle

    

def julian_day(date: Date, calendar: str = 'Gregorian', MJD: bool = False) -> float:
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
    year, month, day = date.date
    day = date.daydec()
    if month <= 2:
        year -= 1
        month += 12
    JD = np.floor(365.25 * (year + 4716)) + np.floor(30.6001 * (month + 1)) + day - 1524.5
    if calendar == 'Gregorian':
        A = np.floor(year / 100)
        B = (2 - A + np.floor(A / 4)) 
        JD += B
    elif calendar != 'Julian':
        raise Exception("!ERROR in chosen calendar!\nFunction accepts only the string:\n\t- 'Gregorian'\n\t- 'Julian'")
    if MJD:
        JD -= 2400000.5
    return JD


def calendar_date(JD: float, timetype: str = 'UT', epoch: str ='J2000.0') -> Date:
    """Computing the calendar date from JD

    Different in *Explatatory supplement to the Astronomical Almanac*, pg. 605  

    :param JD: Julian Day
    :type JD: float

    :return: calendar date in [year, month, day.decimal]
    :rtype: list
    """
    JD += 0.5
    F, Z = np.modf(JD)
    if Z < 2299161:
        A = Z
    else:
        alpha = np.floor((Z-1867216.25)/36524.25)
        A = Z + 1 + alpha - np.floor(alpha / 4)
    B = A + 1524
    C = np.floor((B-122.1)/365.25)
    D = np.floor(365.25 * C)
    E = np.floor((B-D) / 30.6001)
    
    day = B - D - np.floor(30.6001 * E) + F
    month = E-1 if E < 14 else E-13
    year = C-4716 if month > 2 else C-4715
    return Date([year,month,day],timetype=timetype,epoch=epoch)


def mean_Green_HA(date: Date) -> HAngles:
    JD = julian_day(date)
    T = (JD - 2451545) / 36525
    theta0 = 280.46061837 + 360.98564736629 * (JD - 2451545) + 3.87933e-4 * T**2 - T**3 / 3871e4
    theta0 = theta0 / 15
    if type(theta0) != np.ndarray:
        #! Da capire e rivedere questa parte
        if np.abs(theta0) > 24:
            numbdays = np.trunc(theta0 / 24).astype(int)
            theta0 -= numbdays*24
        if theta0 < 0:
            theta0 += 24
        #!
        return HAngles(theta0,'hms')
    else:
        theta0 = np.where(np.abs(theta0) > 24, theta0 - np.trunc(theta0/24).astype(int)*24, theta0)
        theta0 = np.where(theta0 < 0, theta0 + 24, theta0)
        return ArrAngle(theta0,'hms')


def local_ST(date: Date, long: Angles) -> HAngles | ArrAngle:
    GST = mean_Green_HA(date)
    return GST - long


