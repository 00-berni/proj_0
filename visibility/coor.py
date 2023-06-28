import numpy as np
from visibility.angles import Angles, HAngles, ArrAngle


class Equatorial():

    def __init__(self, alpha: HAngles | float | list, delta: HAngles | float | list) -> None:
        if type(alpha) != HAngles:
            alpha = HAngles(alpha,'hms')
        if type(delta) != HAngles:
            delta = HAngles(delta,'deg',lim=90)
        self.alpha = alpha.copy()
        self.delta = delta.copy()
    
    def copy(self):
        return Equatorial(self.alpha,self.delta)


    def print_values(self, sel: str = 'all'):
        alpha_str = 'alpha:\t'+self.alpha.str_angle('hms',unit=True)
        delta_str = 'delta:\t'+self.delta.str_angle('deg',unit=True)
        if sel == 'alpha': return alpha_str
        elif sel == 'delta': return delta_str
        elif sel == 'all': return 'Equatorial Coordinates:\n' + alpha_str + '\n' + delta_str
        else: raise Exception(f'!Error: `{sel}` is not allowed!\nSee the docstring of `print_values()` function')


class Ecliptical():
    def __init__(self, lon: Angles | float | list, lat: Angles | float | list) -> None:
        if type(lon) != Angles:
            lon = Angles(lon,'deg')
        if type(lat) != Angles:
            lat = Angles(lat,'deg',lim=90)
        self.lon = lon.copy()
        self.lat = lat.copy()

    def copy(self):
        return Ecliptical(self.lon, self.lat)

    def print_values(self, sel: str = 'all'):
        lon_str = 'lambda:\t'+self.lon.str_angle('deg',unit=True)
        lat_str = 'beta:\t'+self.lat.str_angle('deg',unit=True)
        if sel == 'lon': return lon_str
        elif sel == 'lat': return lat_str
        elif sel == 'all': return 'Ecliptical Coordinates:\n' + lon_str + '\n' + lat_str
        else: raise Exception(f'!Error: `{sel}` is not allowed!\nSee the docstring of `print_values()` function')


class AltAz():
    def __init__(self, alt: Angles | float | list, az: Angles | float | list) -> None:
        if type(alt) != Angles:
            alt = Angles(alt,'deg',lim=90)
        if type(az) != Angles:
            az = Angles(az,'deg')
        self.alt = alt.copy()
        self.az = az.copy()
    
    def copy(self):
        return AltAz(self.alt,self.az)

    def print_values(self, sel: str = 'all'):
        alt_str = 'alt:\t'+self.alt.str_angle('deg',unit=True)
        az_str = 'az:\t'+self.az.str_angle('deg',unit=True)
        if sel == 'alt': return alt_str
        elif sel == 'az': return az_str
        elif sel == 'all': return 'Local Coordinates:\n' + alt_str + '\n' + az_str
        else: raise Exception(f'!Error: `{sel}` is not allowed!\nSee the docstring of `print_values()` function')


class GeoPos():
    def __init__(self,lon: Angles | float | list, lat: Angles | float | list, h: float) -> None:
        if type(lon) != Angles:
            lon = Angles(lon,'deg')
        if type(lat) != Angles:
            lat = Angles(lat,'deg',lim=90)
        self.lon = lon.copy()
        self.lat = lat.copy()
        self.h = h

    def copy(self):
        return GeoPos(self.lon,self.lat,self.h)


def equat_to_eclipt(eq_coor: Equatorial, e: Angles) -> Ecliptical:
    alpha = eq_coor.alpha
    delta = eq_coor.delta

    lon = np.arctan2( (np.sin(alpha.rad)*np.cos(e.rad) + np.tan(delta.rad)*np.sin(e.rad)), np.cos(alpha.rad) )
    lat = np.arcsin( np.sin(delta.rad)*np.cos(e.rad) - np.cos(delta.rad)*np.sin(e.rad)*np.sin(alpha.rad) )

    lon = HAngles(lon,'rad')
    lat = HAngles(lat,'rad',lim=90)
    return Ecliptical(lon,lat)

def eclipt_to_equat(ec_coor: Ecliptical, e: Angles) -> Equatorial:
    lon = ec_coor.lon
    lat = ec_coor.lat

    alpha = np.arctan2( (np.sin(lon.rad)*np.cos(e.rad) - np.tan(lat.rad)*np.sin(e.rad)), np.cos(lon.rad) ) 
    delta = np.arcsin( np.sin(lat.rad)*np.cos(e.rad) + np.cos(lat.rad)*np.sin(e.rad)*np.sin(lon.rad) )

    alpha = HAngles(alpha,'rad')
    delta = HAngles(delta,'rad',lim=90)
    return Equatorial(alpha,delta)

def equat_to_altaz(eq_coor: Equatorial, LST: HAngles, lat: Angles) -> AltAz:
    delta = eq_coor.delta
    LHA = eq_coor.ha(LST)
    
    az = np.arctan2(np.sin(LHA.rad), (np.cos(LHA.rad)*np.sin(lat.rad) - np.tan(delta.rad)*np.cos(lat.rad)))
    alt = np.arcsin(np.sin(lat.rad)*np.sin(delta.rad) + np.cos(lat.rad)*np.cos(delta.rad)*np.cos(LHA.rad))

    alt = HAngles(alt,'rad',lim=90)
    az = HAngles(az,'rad')
    return AltAz(alt,az)


def altaz_to_equat(loc_coor: AltAz, LST: HAngles, lat: Angles) -> Equatorial:
    alt = loc_coor.alt
    az = loc_coor.az
    
    LHA = np.arctan2(np.sin(az.rad), (np.cos(az.rad)*np.sin(lat.rad) + np.tan(alt.rad)*np.cos(lat.rad)))
    delta = np.arcsin(np.sin(lat.rad)*np.sin(alt.rad) - np.cos(lat.rad)*np.cos(alt.rad)*np.cos(az.rad))

    alpha = LST - LHA
    
    alpha = HAngles(alpha,'rad')
    delta = HAngles(delta,'rad')
    return Equatorial(alpha,delta)

