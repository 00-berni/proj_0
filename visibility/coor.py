import numpy as np
from visibility.Angles import Angles, HAngles


class Equatorial():

    def __init__(self, alpha: HAngles | float | list | np.ndarray | None = None, delta: HAngles | float | list | np.ndarray | None = None) -> None:
        if type(alpha) != HAngles:
            alpha = HAngles(alpha,'hms')
        if type(delta) != HAngles:
            delta = HAngles(delta,'deg',lim=90)
        self.alpha = alpha.copy()
        self.delta = delta.copy()
    
    def copy(self):
        return Equatorial(self.alpha,self.delta)
    
    def print_values(self, sel: str = 'all') -> list[str] | str:
        title_str = 'Equatorial Coordinates'
        alpha_str = self.alpha.print_angle('hms',unit=True)
        delta_str = self.delta.print_angle('deg',unit=True)
        if type(alpha_str) != list:
            alpha_str = [alpha_str]
            delta_str = [delta_str]
        title_list = []    
        res_list = []    
        for i in range(len(alpha_str)):
            alpha_str[i] = 'alpha:\t' + alpha_str[i]
            delta_str[i] = 'delta:\t' + delta_str[i]
            if sel == 'all':
                title_list += [title_str]
                if len(alpha_str) > 1:
                    title_list[i] +=  f' {i}'
                title_list[i] += ':\n'
                res_list += [title_list[i] + alpha_str[i] + '\n' + delta_str[i]]
            elif sel == 'alpha':
                res_list += [alpha_str[i]]
            elif sel == 'delta':
                res_list += [delta_str[i]]
            else: raise Exception(f'!Error: `{sel}` is not allowed!\nSee the docstring of `print_values()` function')
        if len(res_list) == 1:
            res_list = res_list[0] 
        return res_list

class Ecliptical():
    def __init__(self, lon: Angles | float | list | np.ndarray | None = None, lat: Angles | float | list | np.ndarray | None = None) -> None:
        if type(lon) != Angles:
            lon = Angles(lon,'deg')
        if type(lat) != Angles:
            lat = Angles(lat,'deg',lim=90)
        self.lon = lon.copy()
        self.lat = lat.copy()

    def copy(self):
        return Ecliptical(self.lon, self.lat)

    def print_values(self, sel: str = 'all') -> list[str] | str:
        title_str = 'Ecliptical Coordinates'
        lon_str = self.alt.print_angle('deg',unit=True)
        lat_str = self.az.print_angle('deg',unit=True)
        if type(lon_str) != list:
            lon_str = [lon_str]
            lat_str = [lat_str]
        title_list = []    
        res_list = []    
        for i in range(len(lon_str)):
            lon_str[i] = 'lon:\t' + lon_str[i]
            lat_str[i] = 'lat:\t' + lat_str[i]
            if sel == 'all':
                title_list += [title_str]
                if len(lon_str) > 1:
                    title_list[i] +=  f' {i}'
                title_list[i] += ':\n'
                res_list += [title_list[i] + lon_str[i] + '\n' + lat_str[i]]
            elif sel == 'lon':
                res_list += [lon_str[i]]
            elif sel == 'lat':
                res_list += [lat_str[i]]
            else: raise Exception(f'!Error: `{sel}` is not allowed!\nSee the docstring of `print_values()` function')
        if len(res_list) == 1:
            res_list = res_list[0] 
        return res_list



class AltAz():
    def __init__(self, alt: Angles | float | list | np.ndarray | None = None, az: Angles | float | list | np.ndarray | None = None) -> None:
        if type(alt) != Angles:
            alt = Angles(alt,'deg',lim=90)
        if type(az) != Angles:
            az = Angles(az,'deg')
        self.alt = alt.copy()
        self.az = az.copy()
    
    def copy(self):
        return AltAz(self.alt,self.az)

    def print_values(self, sel: str = 'all') -> list[str] | str:
        title_str = 'Local Coordinates'
        alt_str = self.alt.print_angle('deg',unit=True)
        az_str = self.az.print_angle('deg',unit=True)
        if type(alt_str) != list:
            alt_str = [alt_str]
            az_str = [az_str]
        title_list = []    
        res_list = []    
        for i in range(len(alt_str)):
            alt_str[i] = 'alt:\t' + alt_str[i]
            az_str[i] = 'az:\t' + az_str[i]
            if sel == 'all':
                title_list += [title_str]
                if len(alt_str) > 1:
                    title_list[i] +=  f' {i}'
                title_list[i] += ':\n'
                res_list += [title_list[i] + alt_str[i] + '\n' + az_str[i]]
            elif sel == 'alt':
                res_list += [alt_str[i]]
            elif sel == 'az':
                res_list += [az_str[i]]
            else: raise Exception(f'!Error: `{sel}` is not allowed!\nSee the docstring of `print_values()` function')
        if len(res_list) == 1:
            res_list = res_list[0] 
        return res_list



class GeoPos():
    def __init__(self,lon: Angles | float | list | None = None, lat: Angles | float | list | None = None, h: float | int = 0., name: str = '') -> None:
        if type(lon) != Angles:
            lon = Angles(lon,'deg')
        if type(lat) != Angles:
            lat = Angles(lat,'deg',lim=90)
        self.name = name 
        self.lon = lon.copy()
        self.lat = lat.copy()
        self.h = h

    def copy(self):
        return GeoPos(self.lon,self.lat,self.h)
    
    def place_info(self):
        return f"{self.name}\nlon:\t{self.lon.print_angle('deg',unit=True)}\nlat:\t{self.lat.print_angle('deg',unit=True)}\nheight:\t{self.h:.0f} m"


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

def equat_to_altaz(eq_coor: Equatorial, LHA: HAngles, lat: Angles) -> AltAz:
    delta = eq_coor.delta
    
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

