import numpy as np
from test_struc import *


def julian_day(year: int, month: int, day: float, calendar: str = 'Gregorian', MJD: bool = False) -> float:
    """Computing the Julian Day

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
    """
    if month <= 2:
        year -= 1
        month += 12
    A = np.floor(year / 100)
    if calendar == 'Gregorian':
        B = (2 - A + np.floor(A / 4)) 
    elif calendar == 'Julian':
        B = 0
    else:
        raise Exception("!ERROR in chosen calendar!\nFunction accepts only the string:\n\t- 'Gregorian'\n\t- 'Julian'")
    JD = np.floor(365.25 * (year + 4716)) + np.floor(30.6001 * (month + 1)) + day + B - 1524.5
    if MJD:
        JD -= 2400000.5
    return JD

def calendar_date(JD: float) -> list:
    """_summary_

    Different in *Explatatory supplement to the Astronomical Almanac*, pg. 605  

    :param JD: _description_
    :type JD: float
    :return: _description_
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
    return [int(year), int(month), day]


if __name__ == '__main__':

    starting_test('TEST FOR JULIAN DAY')

    try:
        years, months, days, cal, JD0 = get_test_data('test_JD_data.csv')

        JD = np.array([JulianDay(years[i],months[i],days[i],cal[i]) for i in range(len(years))])
        print('> Compute the JD')
        print('JD0')
        for jd, jd0 in zip(JD,JD0):
            if jd == jd0:
                print(f'{jd0}: OK!')
            else:
                print(f'{jd0}: NO MATCH --> DIFF = {jd-jd0}\n\tJD = {jd}')

        print('> All results matched')

        ending_test()
    except:
        test_error()