import numpy as np
from test_angle_class import *

class Equatorial():
    def __init__(self, alpha: HAngles, delta: HAngles) -> None:
        self.alpha = alpha
        self.delta = delta

class Ecliptical():
    def __init__(self, lon: HAngles, lat: HAngles) -> None:
        self.lon = lon
        self.lat = lat

class AltAz():
    def __init__(self, alt: HAngles, az: HAngles) -> None:
        self.alt = alt
        self.az = az

def equat_to_eclipt(eq_coor: Equatorial, e: Angles) -> Ecliptical:
    alpha = eq_coor.alpha
    delta = eq_coor.delta

    lon = np.arctan( (np.sin(alpha.rad)*np.cos(e.rad) + np.tan(delta.rad)*np.sin(e.rad)) / np.cos(alpha.rad) )
    lat = np.arcsin( np.sin(delta.rad)*np.cos(e.rad) - np.cos(delta.rad)*np.sin(e.rad)*np.sin(alpha.rad) )

    lon = HAngles(lon,'rad')
    if lon.rad < 0: lon = lon + FLAT
    lat = HAngles(lat,'rad',lim=90)
    return Ecliptical(lon,lat)

def eclipt_to_equat(ec_coor: Ecliptical, e: Angles) -> Equatorial:
    lon = ec_coor.lon
    lat = ec_coor.lat

    alpha = np.arctan( (np.sin(lon.rad)*np.cos(e.rad) - np.tan(lat.rad)*np.sin(e.rad)) / np.cos(lon.rad) ) 
    delta = np.arcsin( np.sin(lat.rad)*np.cos(e.rad) + np.cos(lat.rad)*np.sin(e.rad)*np.sin(lon.rad) )

    alpha = HAngles(alpha,'rad')
    if alpha.rad < 0: alpha = alpha + FLAT
    delta = HAngles(delta,'rad',lim=90)
    return Equatorial(alpha,delta)

def equat_to_altaz(eq_coor: Equatorial, LST: HAngles, lat: Angles) -> AltAz:
    alpha = eq_coor.alpha
    delta = eq_coor.delta

    LHA = LST - alpha
    
    az = np.arctan(np.sin(LHA.rad) / (np.cos(LHA.rad)*np.sin(lat.rad) - np.tan(delta.rad)*np.cos(lat.rad)))
    alt = np.arcsin(np.sin(lat.rad)*np.sin(delta.rad) + np.cos(lat.rad)*np.cos(delta.rad)*np.cos(LHA.rad))

    alt = HAngles(alt,'rad',lim=90)
    az = HAngles(az,'rad')
    return AltAz(alt,az)


def altaz_to_equat(loc_coor: AltAz, LST: HAngles, lat: Angles) -> Equatorial:
    alt = loc_coor.alt
    az = loc_coor.az

    
    LHA = np.arctan(np.sin(az.rad) / (np.cos(az.rad)*np.sin(lat.rad) + np.tan(alt.rad)*np.cos(lat.rad)))
    delta = np.arcsin(np.sin(lat.rad)*np.sin(alt.rad) - np.cos(lat.rad)*np.cos(alt.rad)*np.cos(az.rad))

    alpha = LST - LHA
    
    alpha = HAngles(alpha,'rad')
    delta = HAngles(delta,'rad')
    return Equatorial(alpha,delta)


if __name__ == '__main__':
    SEP = lambda obj : '------' + obj + '------\n'
    print(SEP('TEST COORDINATE SYSTEM'))

    try:    
        alpha = [7,45,18.946]
        delta = [28,1,34.26]
        epsilon = 23.4392911

        alpha = HAngles(HAngles.decimal(alpha),'hms')
        delta = HAngles(HAngles.decimal(delta),'deg')
        epsilon = Angles(epsilon,'deg')

        print('TEST: FROM EQUATORIAL TO ECLIPTIC')
        print('Epoch 2000')
        print('alpha:\t'+alpha.print_angle('rad'))
        print('delta:\t'+delta.print_angle('rad'))


        eq_coor = Equatorial(alpha,delta)

        ec_coor = equat_to_eclipt(eq_coor,epsilon)

        lon = ec_coor.lon
        lat = ec_coor.lat

        print('Ecliptic coordinates')
        print('lambda:\t'+lon.print_angle('deg')+f'--> {HAngles.decimal(lon.deg[1])}')
        print('beta:\t'+lat.print_angle('deg')+f'--> {HAngles.decimal(lat.deg[1])}')

        print()

        print('Try to recover equatorial from results')
        print('\nTEST: FROM ECLIPTIC TO EQUATORIAL')

        test_data = [alpha.rad,delta.rad]

        eq_coor = eclipt_to_equat(ec_coor,epsilon)

        alpha = eq_coor.alpha
        delta = eq_coor.delta
        

        if alpha.rad == test_data[0] and delta.rad/pi == test_data[1]/pi:
            print('Equatorial coordinates')
            print('alpha:\t'+alpha.print_angle('hms'))
            print('delta:\t'+delta.print_angle('deg'))

            print('\n!TEST COMPLETE!')
            print(SEP('----------'))
        else:
            print('Equatorial coordinates')
            print('alpha:\t'+alpha.print_angle('rad'))
            print('delta:\t'+delta.print_angle('rad'))

            # print('\n!TEST FAILD!')
            raise Exception('Results are not the same!')           
    except:
        print('\n!TEST FAILD!')
        raise