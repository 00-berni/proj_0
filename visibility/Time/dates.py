import numpy as np
# import os.path as ph
# PROJECT_FOLDER = ph.split(ph.split(ph.dirname(ph.realpath(__file__)))[0])[0]
# import sys
# sys.path.append(PROJECT_FOLDER)

from visibility.Time.Tclasses import *
from visibility.Angles import Angles, HAngles

def mean_Green_HA(date: Date) -> HAngles:
    JD = date.jd
    T = (JD - 2451545) / 36525
    theta0 = 280.46061837 + 360.98564736629 * (JD - 2451545) + 3.87933e-4 * T**2 - T**3 / 3871e4
    theta0 = theta0 / 15
    if type(theta0) != np.ndarray:
        theta0 = np.array([theta0])
    #! Da capire e rivedere questa parte
    theta0 = np.where(np.abs(theta0) > 24, theta0 - np.trunc(theta0/24).astype(int)*24, theta0)
    theta0 = np.where(theta0 < 0, theta0 + 24, theta0)
    #!
    if len(theta0) == 1:
        theta0 = theta0[0]
    return HAngles(theta0,'hms')


def local_ST(date: Date, long: Angles) -> HAngles:
    GST = mean_Green_HA(date)
    return GST - long


if __name__ == '__main__':
    date = Date([2000,11,28])
    long = Angles(43.231,'deg')
    print(date.print_date())
    lst = local_ST(date,long)
    print(lst.print_angle('hms'))   
    print() 
    jd = np.linspace(Date.J2000+500.5,Date.J2000+505,200)
    date = Date(jd=jd)
    long = Angles(43.231,'deg')
    lst = local_ST(date,long)
    for i in range(5):
        print(date.print_date()[i])
        print(lst.print_angle('hms',unit=True)[i])
        print('- -')
    print()
    long = Angles(np.linspace(43.22,120.32,len(jd)),'deg')
    lst = local_ST(date,long)
    for i in range(5):
        print(date.print_date()[i])
        print(long.print_angle('deg')[i])
        print(lst.print_angle('hms',unit=True)[i])
        print('- -')

