import numpy as np
from visibility.Angles import Angles, HAngles


class Equatorial():
    """Class for equatorial celestial coordinates.

    The attributes of the class are:

    :ivar alpha: right ascension
    :vartype alpha: HAngles
    :ivar delta: declination
    :vartype delta: HAngles

    .. note:: It is possible to generate an empty object 
    using the command `Equatorial()`.
    """

    def __init__(self, alpha: HAngles | float | list | np.ndarray | None = None, delta: HAngles | float | list | np.ndarray | None = None) -> None:
        """Constructor of the class

        Function takes the angle values for right ascension and 
        declination and stores them.

        :param alpha: right ascension, defaults to `None`
        :type alpha: HAngles | float | list | np.ndarray | None, optional
        :param delta: declination, defaults to `None`
        :type delta: HAngles | float | list | np.ndarray | None, optional
        """
        # conditions for not-HAngles values
        if type(alpha) != HAngles:
            alpha = HAngles(alpha,'hms')
        if type(delta) != Angles:
            delta = Angles(delta,'deg',lim=90)

        self.alpha = alpha.copy()
        self.delta = delta.copy()
    
    def copy(self):
        """Function to make an exact copy of the coordinates

        :return: the exact copy of the coordinates object
        :rtype: Equatorial
        """
        return Equatorial(self.alpha,self.delta)

    def coor(self):
        return self.alpha, self.delta
    
    def print_values(self, sel: str = 'all', eph: str | None = None) -> list[str] | str:
        """Function to print the values of equatorial coordinates.

        One can select the output through `sel` parameter:
        
          * `sel = 'alpha'`: only right ascension is printed
          * `sel = 'delta'`: only declination is printed
          * `sel = 'all'`: both

        The parameter `eph` is used to print the standard epoch at 
        which coordinates are measured or calculated.
          
        :param sel: selection parameter, defaults to `'all'`
        :type sel: str, optional
        :param eph: standard epoch of coordinates, defaults to `''`
        :type eph: str, optional

        :return: string (or list of string) with coordinates values 
        :rtype: list[str] | str
        
        :raises Exception: only values in docstring are allowed for `sel`
        """
        # title
        title_str = 'Equatorial Coordinates'
        if eph is not None: title_str += ' for epoch '+ eph 
        # getting string for angles
        alpha_str = self.alpha.print_angle('hms',unit=True)
        delta_str = self.delta.print_angle('deg',unit=True)
        # generalizing the method for single value
        if type(alpha_str) != list:
            alpha_str = [alpha_str]
            delta_str = [delta_str]
        # defing list for titles
        title_list = []
        # defing list for values    
        res_list = []    
        for i in range(len(alpha_str)):
            alpha_str[i] = 'alpha:\t' + alpha_str[i]
            delta_str[i] = 'delta:\t' + delta_str[i]
            # condition for all values
            if sel == 'all':
                title_list += [title_str]
                if len(alpha_str) > 1:
                    # printing the index for array 
                    title_list[i] +=  f' {i}'
                title_list[i] += ':\n'
                res_list += [title_list[i] + alpha_str[i] + '\n' + delta_str[i]]
            elif sel == 'alpha':
                res_list += [alpha_str[i]]
            elif sel == 'delta':
                res_list += [delta_str[i]]
            else: raise Exception(f'!Error: `{sel}` is not allowed!\nSee the docstring of `print_values()` function')
        # only a string for not-array value
        if len(res_list) == 1:
            res_list = res_list[0] 
        return res_list

class Ecliptical():
    """Class for ecliptical celestial coordinates.

    The attributes of the class are:

    :ivar lon: ecliptical longitude
    :vartype lon: Angles
    :ivar lat: ecliptical latitude
    :vartype lat: Angles

    .. note:: It is possible to generate an empty object 
    using the command `Ecliptical()`.
    """
    def __init__(self, lon: Angles | float | list | np.ndarray | None = None, lat: Angles | float | list | np.ndarray | None = None) -> None:
        """Constructor of the class

        Function takes the angle values for ecliptical longitude 
        and latitude and stores them.

        :param lon: ecliptical longitude, defaults to `None`
        :type lon: Angles | float | list | np.ndarray | None, optional
        :param lat: ecliptical latitude, defaults to `None`
        :type lat: Angles | float | list | np.ndarray | None, optional
        """
        # conditions for not-Angles values
        if type(lon) != Angles:
            lon = Angles(lon,'deg')
        if type(lat) != Angles:
            lat = Angles(lat,'deg',lim=90)
        
        self.lon = lon.copy()
        self.lat = lat.copy()

    def coor(self):
        return self.lon, self.lat

    def copy(self):
        """Function to make an exact copy of the coordinates

        :return: the exact copy of the coordinates object
        :rtype: Ecliptical
        """
        return Ecliptical(self.lon, self.lat)

    def print_values(self, sel: str = 'all', eph: str | None = None) -> list[str] | str:
        """Function to print the values of ecliptical coordinates.

        One can select the output through `sel` parameter:
        
          * `sel = 'lon'`: only ecliptical longitude is printed
          * `sel = 'lat'`: only ecliptical latitude is printed
          * `sel = 'all'`: both

        The parameter `eph` is used to print the standard epoch at 
        which coordinates are measured or calculated.

        :param sel: selection parameter, defaults to `'all'`
        :type sel: str, optional
        :param eph: standard epoch of coordinates, defaults to `''`
        :type eph: str, optional

        :return: string (or list of string) with coordinates values 
        :rtype: list[str] | str
        
        :raises Exception: only values in docstring are allowed for `sel`
        """
        # title
        title_str = 'Ecliptical Coordinates'
        if eph is not None: title_str += ' for epoch '+ eph 
        # getting string for angles
        lon_str = self.alt.print_angle('deg',unit=True)
        lat_str = self.az.print_angle('deg',unit=True)
        # generalizing the method for single value
        if type(lon_str) != list:
            lon_str = [lon_str]
            lat_str = [lat_str]
        # defing list for titles
        title_list = []    
        # defing list for values    
        res_list = []    
        for i in range(len(lon_str)):
            lon_str[i] = 'lon:\t' + lon_str[i]
            lat_str[i] = 'lat:\t' + lat_str[i]
            # condition for all values
            if sel == 'all':
                title_list += [title_str]
                if len(lon_str) > 1:
                    # printing the index for array 
                    title_list[i] +=  f' {i}'
                title_list[i] += ':\n'
                res_list += [title_list[i] + lon_str[i] + '\n' + lat_str[i]]
            elif sel == 'lon':
                res_list += [lon_str[i]]
            elif sel == 'lat':
                res_list += [lat_str[i]]
            else: raise Exception(f'!Error: `{sel}` is not allowed!\nSee the docstring of `print_values()` function')
        # only a string for not-array value
        if len(res_list) == 1:
            res_list = res_list[0] 
        return res_list



class AltAz():
    """Class for local celestial coordinates.

    The attributes of the class are:

    :ivar alt: altitude
    :vartype alt: Angles
    :ivar az: azimuth
    :vartype az: Angles

    .. note:: It is possible to generate an empty object 
    using the command `AltAz()`.
    """
    def __init__(self, alt: Angles | float | list | np.ndarray | None = None, az: Angles | float | list | np.ndarray | None = None) -> None:
        """Constructor of the class

        Function takes the angle values for altitude and 
        azimuth and stores them.

        :param alt: altitude, defaults to `None`
        :type alt: Angles | float | list | np.ndarray | None, optional
        :param az: azimuth, defaults to `None`
        :type az: Angles | float | list | np.ndarray | None, optional
        """
        # conditions for not-Angles values
        if type(alt) != Angles:
            alt = Angles(alt,'deg',lim=90)
        if type(az) != Angles:
            az = Angles(az,'deg')
       
        self.alt = alt.copy()
        self.az = az.copy()
    
    def copy(self):
        """Function to make an exact copy of the coordinates

        :return: the exact copy of the coordinates object
        :rtype: AltAz
        """
        return AltAz(self.alt,self.az)

    def coor(self):
        return self.alt, self.az

    def print_values(self, sel: str = 'all', eph: str | None = None) -> list[str] | str:
        """Function to print the values of local coordinates.

        One can select the output through `sel` parameter:
        
          * `sel = 'alt'`: only altitude is printed
          * `sel = 'az'`: only azimuth is printed
          * `sel = 'all'`: both

        The parameter `eph` is used to print the standard epoch at 
        which coordinates are measured or calculated.
          
        :param sel: selection parameter, defaults to `'all'`
        :type sel: str, optional
        :param eph: standard epoch of coordinates, defaults to `''`
        :type eph: str, optional

        :return: string (or list of string) with coordinates values 
        :rtype: list[str] | str
        
        :raises Exception: only values in docstring are allowed for `sel`
        """
        # title
        title_str = 'Local Coordinates'
        if eph is not None: title_str += ' for epoch '+ eph 
        # getting string for angles
        alt_str = self.alt.print_angle('deg',unit=True)
        az_str = self.az.print_angle('deg',unit=True)
        # generalizing the method for single value
        if type(alt_str) != list:
            alt_str = [alt_str]
            az_str = [az_str]
        # defing list for titles
        title_list = []    
        # defing list for values    
        res_list = []    
        for i in range(len(alt_str)):
            alt_str[i] = 'alt:\t' + alt_str[i]
            az_str[i] = 'az:\t' + az_str[i]
            # condition for all values
            if sel == 'all':
                title_list += [title_str]
                if len(alt_str) > 1:
                    # printing the index for array 
                    title_list[i] +=  f' {i}'
                title_list[i] += ':\n'
                res_list += [title_list[i] + alt_str[i] + '\n' + az_str[i]]
            elif sel == 'alt':
                res_list += [alt_str[i]]
            elif sel == 'az':
                res_list += [az_str[i]]
            else: raise Exception(f'!Error: `{sel}` is not allowed!\nSee the docstring of `print_values()` function')
        # only a string for not-array value
        if len(res_list) == 1:
            res_list = res_list[0] 
        return res_list

class GeoPos():
    """Class for terrestrial coordinates for observatory.

    The attributes of the class are:

    :ivar name: name of the location
    :vartype name: str
    :ivar lon: terrestrial longitude
    :vartype lon: Angles
    :ivar lat: terrestrial latitude
    :vartype lat: Angles
    :ivar h: height a.s.l.
    :vartype lat: float | int

    .. note:: It is possible to generate an empty object 
    using the command `GeoPos()`.
    """
    def __init__(self,lon: Angles | float | list | None = None, lat: Angles | float | list | None = None, h: float | int = 0., name: str = '') -> None:
        """Constructor of the class

        Function takes the angle values for terrestrial longitude and 
        latitude and height a.s.l., then it stores them.

        :param lon: terrestrial longitude, defaults to `None`
        :type lon: Angles | float | list | None, optional
        :param lat: terrestrial latitude, defaults to `None`
        :type lat: Angles | float | list | None, optional
        :param h: height a.s.l., defaults to `0.`
        :type h: float | int, optional
        :param name: name of the location, defaults to `''`
        :type name: str, optional
        """
        # conditions for not-Angles values
        if type(lon) != Angles:
            lon = Angles(lon,'deg')
        if type(lat) != Angles:
            lat = Angles(lat,'deg',lim=90)
        
        self.name = name 
        self.lon = lon.copy()
        self.lat = lat.copy()
        self.h = h

    def copy(self):
        """Function to make an exact copy of the coordinates

        :return: the exact copy of the coordinates object
        :rtype: GeoPos
        """
        return GeoPos(self.lon,self.lat,self.h)
    
    def coor(self):
        return self.lon, self.lat, self.h

    def place_info(self, plot: bool = False) -> str:
        """Function to print informations about the location.

        :return: string with location information
        :rtype: str
        """
        strlon = self.lon.print_angle('deg',unit=True)
        angsign = strlon[0]
        strlon = 'W ' + strlon[1:] if angsign == '+' else 'E ' + strlon[1:]
        strlat = self.lat.print_angle('deg',unit=True)
        angsign = strlat[0]
        strlat = 'N ' + strlat[1:] if angsign == '+' else 'S ' + strlat[1:]
        name = self.name
        if plot:
            if name != '': name = f'({name})'
            return f"{name} {strlat}, {strlon}, {self.h:.0f} m a.s.l."
        else:
            if name != '': name = name + '\n'
            return f"{name}lon:\t{strlon}\nlat:\t{strlat}\nheight:\t{self.h:.0f} m a.s.l."


def equat_to_eclipt(eq_coor: Equatorial, e: Angles) -> Ecliptical:
    """Converting equatorial in ecliptical coordinates

    :param eq_coor: equatorial coordinates
    :type eq_coor: Equatorial
    :param e: obliquity of ecliptic
    :type e: Angles

    :return: ecliptical coordinates
    :rtype: Ecliptical
    """
    alpha = eq_coor.alpha
    delta = eq_coor.delta

    lon = np.arctan2( (np.sin(alpha.rad)*np.cos(e.rad) + np.tan(delta.rad)*np.sin(e.rad)), np.cos(alpha.rad) )
    lat = np.arcsin( np.sin(delta.rad)*np.cos(e.rad) - np.cos(delta.rad)*np.sin(e.rad)*np.sin(alpha.rad) )

    lon = Angles(lon,'rad')
    lat = Angles(lat,'rad',lim=90)
    return Ecliptical(lon,lat)

def eclipt_to_equat(ec_coor: Ecliptical, e: Angles) -> Equatorial:
    """Converting ecliptical in equatorial coordinates

    :param ec_coor: ecliptical coordinates
    :type ec_coor: Ecliptical
    :param e: obliquity of ecliptic
    :type e: Angles

    :return: equatorial coordinates
    :rtype: Equatorial
    """
    lon = ec_coor.lon
    lat = ec_coor.lat

    alpha = np.arctan2( (np.sin(lon.rad)*np.cos(e.rad) - np.tan(lat.rad)*np.sin(e.rad)), np.cos(lon.rad) ) 
    delta = np.arcsin( np.sin(lat.rad)*np.cos(e.rad) + np.cos(lat.rad)*np.sin(e.rad)*np.sin(lon.rad) )

    alpha = HAngles(alpha,'rad')
    delta = Angles(delta,'rad',lim=90)
    return Equatorial(alpha,delta)

def equat_to_altaz(eq_coor: Equatorial, HA: HAngles, lat: Angles) -> AltAz:
    """Converting equatorial in local coordinates

    :param eq_coor: equatorial coordinates
    :type eq_coor: Equatorial
    :param HA: local HA
    :type HA: HAngles
    :param lat: terrestrial latitude
    :type lat: Angles

    :return: local coordinates
    :rtype: AltAz
    """
    delta = eq_coor.delta
    
    az  = np.arctan2(np.sin(HA.rad), (np.cos(HA.rad)*np.sin(lat.rad) - np.tan(delta.rad)*np.cos(lat.rad)))
    alt = np.arcsin(np.sin(lat.rad)*np.sin(delta.rad) + np.cos(lat.rad)*np.cos(delta.rad)*np.cos(HA.rad))

    alt = Angles(alt,'rad',lim=90)
    az  = Angles(az,'rad')
    return AltAz(alt,az)


def altaz_to_equat(loc_coor: AltAz, LST: HAngles, lat: Angles) -> Equatorial:
    """Converting local in equatorial coordinates

    :param loc_coor: local coordinates
    :type loc_coor: AltAz
    :param LST: local sidereal time
    :type LST: HAngles
    :param lat: terrestrial latitude
    :type lat: Angles

    :return: equatorial coordinates
    :rtype: Equatorial
    """
    alt = loc_coor.alt
    az = loc_coor.az
    
    HA = np.arctan2(np.sin(az.rad), (np.cos(az.rad)*np.sin(lat.rad) + np.tan(alt.rad)*np.cos(lat.rad)))
    delta = np.arcsin(np.sin(lat.rad)*np.sin(alt.rad) - np.cos(lat.rad)*np.cos(alt.rad)*np.cos(az.rad))

    alpha = LST - HA
    
    alpha = HAngles(alpha,'rad')
    delta = Angles(delta,'rad')
    return Equatorial(alpha,delta)
