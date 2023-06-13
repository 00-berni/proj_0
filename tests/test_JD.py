import numpy as np
from test_struc import *


def JulianDay(year: int, month: int, day: float, calendar: str = 'Gregorian') -> float:
    """Computing the Julian Day

    From J.,Meeus, *Astronomical Algorithms*, pg. 61. 
    Follow the same notation.

    Same from M., Gallaway, *An Introduction to Observational Astrophysics*, 
    pg. 41.

    Different in *Explanatory supplement to the Astronomical Almanac*, pg. 604

    :param year: year
    :type year: int
    :param month: month
    :type month: int
    :param day: day 
    :type day: float
    :param calendar: kind of chosen calendar, defaults to 'Gregorian'
    :type calendar: str, optional

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
    return JD


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