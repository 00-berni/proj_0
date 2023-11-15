import numpy as np
# import os.path as ph
# PROJECT_FOLDER = ph.split(ph.split(ph.dirname(ph.realpath(__file__)))[0])[0]
# import sys
# sys.path.append(PROJECT_FOLDER)

from .Tclasses import *
from visibility.Angles import Angles, HAngles

def nutation_corr(date: Date) -> tuple[Angles, Angles]:
    """Computing nutation correction

    The accuracy of this formula is 
    0.5'' for celestial longitude and
    0.1'' for obliquity.

    :param date: selected date against which to compute correction
    :type date: Date
    
    :return: correction for celestial longitude and for obliquity
    :rtype: tuple[HAngles, HAngles]
    """
    # computing the Julian Ephemeris Day 
    date = date.change_time_type('TD')
    T = (date.jd - Date.J2000)/36525
    # mean longitude of the Sun
    L  = Angles.deg_to_rad(280.4665 + 36000.7698*T)
    # mean longitude of the Moon
    l  = Angles.deg_to_rad(218.3165 + 481267.8813*T)
    # mean longitude of the ascending node of the Moon's mean orbit
    Om = Angles.deg_to_rad(125.04452 - 1934.136261*T)
    # correction in longitude
    Dlon = -17.20*np.sin(Om) - 1.32*np.sin(2*L) - 0.23*np.sin(2*l) + 0.21*np.sin(2*Om) 
    # correction in obliquity
    Deps =   9.20*np.cos(Om) + 0.57*np.cos(2*L) + 0.10*np.cos(2*l) - 0.09*np.cos(2*Om) 
    return Angles(Dlon/3600,'deg'), Angles(Deps/3600,'deg')

def mean_obliquity(date: Date) -> Angles:
    """Computing mean obliquity

    The error remains less than
    1'' in a range of 2000 years
    from the epoch J2000.0 

    :param date: selected date against which to compute obliquity
    :type date: Date
    :param epoch: standard epoch, defaults to `'J2000.0'`
    :type epoch: str, optional

    :return: mean obliquity
    :rtype: HAngles
    """
    # computing the Julian Ephemeris Day 
    date = date.change_time_type('TD')
    T = (date.jd - Date.J2000)/36525
    # computing the mean obliquity
    eps = 23.4392911*3600 - 46.815*T - 59e-5*T**2 + 1.813e-3*T**3
    return Angles(eps/3600,'deg')


def Green_ST(date: Date, nut_corr: bool = False, epoch: str = 'J2000.0') -> HAngles:
    """Computing the ST at Greenwich

    It possible to compute the AST through the
    parameter `nut_corr`.

    :param date: date for which to compute ST
    :type date: Date
    :param nut_corr: if `True` AST is computed, defaults to `False`
    :type nut_corr: bool, optional
    :param epoch: standard epoch, defaults to `'J2000.0'`
    :type epoch: str, optional
    
    :return: ST at Greenwich
    :rtype: HAngles
    """
    JD = date.jd
    T = (JD - Date.J2000) / 36525
    theta0 = 280.46061837 + 360.98564736629 * (JD - 2451545) + 3.87933e-4 * T**2 - T**3 / 3871e4
    # converting deg in hours
    theta0 = theta0 / 15
    # condition to generalize the method for not-array type
    if type(theta0) != np.ndarray:
        theta0 = np.array([theta0])
    # check
    theta0 = np.where(np.abs(theta0) > 24, theta0 - np.trunc(theta0/24).astype(int)*24, theta0)
    theta0 = np.where(theta0 < 0, theta0 + 24, theta0)
    # condition for not-array type
    if len(theta0) == 1:
        theta0 = theta0[0]
    # condition to compute AST
    if nut_corr:
        Dlon, De = nutation_corr(date)
        e = mean_obliquity(date) + De
        theta0 += Dlon.deg*np.cos(e.rad) / 15
    return HAngles(theta0,'hms')


def local_ST(date: Date, long: Angles, nut_corr: bool = False, epoch: str = 'J2000.0') -> HAngles:
    """Computing the LST

    It possible to compute the AST through the
    parameter `nut_corr`.    

    Terrestrial longitude is taken positive 
    west and negative east from Greenwich
    
    :param date: date for which to compute ST
    :type date: Date
    :param long: Earth longitude for which to compute ST
    :type long: Angles
    :param nut_corr: if `True` AST is computed, defaults to `False`
    :type nut_corr: bool, optional
    :param epoch: standard epoch, defaults to `'J2000.0'`
    :type epoch: str, optional

    :return: LST
    :rtype: HAngles
    """
    # ST at Greenwich
    GST = Green_ST(date, nut_corr=nut_corr,epoch=epoch)
    return GST - long
