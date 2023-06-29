import numpy as np
from numpy import pi

## angle classes
class Angles():
    """Class to handle with angles in degrees or radiants.
    
    The class takes the value of an angle in deg or rad,
    converts one in other and stored the values.

    The attributes of the class are:

    :ivar rad: the value of the angle in radiants
    :vartype rad: float | None
    :ivar deg: the value of the angle in deg
    :vartype deg: float    
    :ivar lim: the angle is in a range set by [0,`lim`] (see `std_format()` docstring)
    :vartype lim: int

    .. note::
        It is possible to generate a null angle, that is an empty `Angles` variable
        (`self.deg = None` and `self.rad = None`), through the command `Angles(None)`. 

    :Example:

    >>> ang = Angles()
    """

    #: a dictionary to pass from string symbol to integer and vice versa
    strsign = { '+' :  1,
                '-' : -1,
                 1  : '+',
                -1  : '-' }

    @staticmethod
    def decimal(ang: np.ndarray | list) -> float:
        """Function to convert a list angle
        [val, val, val.val] in a float 
        val.val 

        :param ang: angle array value in deg or hms
        :type ang: np.ndarray | list
        
        :return: angle float value
        :rtype: float
        """
        if sum(ang) == 0: return 0
        else: return ang[0] + ang[1]/60 + ang[2]/3600
    
    @staticmethod
    def deg_to_rad(deg: float) -> float:
        """Function to convert deg in rad

        :param deg: angle in deg
        :type deg: float

        :return: value in rad
        :rtype: float
        """
        # computing radiants
        rad = deg * pi / 180
        return rad

    @staticmethod
    def rad_to_deg(rad: float) -> float:
        """Function to convert rad in deg

        :param rad: value in rad
        :type rad: float
        
        :return: angle in deg
        :rtype: float
        """
        # computing degrees
        deg = rad * 180 / pi
        return deg

    def __init__(self, ang: float | list | None, unit: str, lim: int = 360) -> None:
        """Constructor of the class

        The function takes a value (`ang`) and the corrisponding 
        unit (`unit`) as input, computes and stores the angle
        values in deg and rad. 

        :param ang: angle value
        :type ang: float | list | None
        :param unit: unit of angle value, like 'deg' or 'rad'
        :type unit: str
        :param lim: angle edge (see `std_format()` docstring), defaults to 360
        :type lim: int, optional
        """
        # setting lim
        self.lim = lim
        # condition for a null angle
        if ang is None:
            self.deg = None
            self.rad = None
        # value in radiants
        elif unit == 'rad':
            self.rad = ang
            self.deg = Angles.rad_to_deg(ang)
        # value in degrees
        elif unit == 'deg':
            # list format
            if type(ang) == list:
                # getting the sign
                sign = Angles.strsign[ang[0]]
                self.deg = Angles.decimal(ang[1])*sign
            else:
                self.deg = ang
            self.rad = Angles.deg_to_rad(self.deg)

    def std_format(self) -> list:
        """Function to adjust format of a list angle in deg or hms.

        The wanted format is: [`sign`, [deg (<`lim`), pp (<60), ss.ss (<60)]]
          
        :return: correct format for the angle value
        :rtype: list
        """
        # getting value in deg
        ang = self.deg
        if ang == 0.:
            return ['+',np.zeros(3,dtype=int)]
        else:
            # taking the lim
            lim = self.lim
            # taking the sign as a str
            sign = Angles.strsign[np.sign(ang)]
            # taking the absolute value
            ang = abs(ang)
            # edge condition
            if ang > lim:
                ang -= lim * (ang // lim)
            # degrees
            val = [np.trunc(ang).astype(int)]
            # primes
            val += [np.trunc(ang % 1 * 60).astype(int)]
            # seconds
            val += [round(ang % 1 * 60 % 1 * 60,4)]
            return [sign, np.array(val)]

    def str_angle(self,sel: str = 'all', unit: bool = False) -> str:
        """Function to print the value of an angle
        
        One can select to print the angle in all 
        units or just one through the `sel`
        parameter:

            * `sel = 'all'`: print value in all units
            * `sel = 'deg'`: print value in deg
            * `sel = 'rad'`: print value in rad

        :param sel: to select in which unit printing the angle, defaults to 'all'
        :type sel: str, optional

        :return: the string with the values
        :rtype: str

        :raise: 
        """
        # adjusting the format in deg
        deg = self.std_format()
        # collecting string for each unit
        dd, mm, ss = deg[1]
        deg_str = deg[0]+f'[{dd:.0f}, {mm:.0f}, {ss:.4f}]'
        rad_str = f'{self.rad/pi} * pi'
        
        if unit:
            deg_str += ' deg'
            rad_str += ' rad'

        if sel == 'deg': return deg_str
        elif sel == 'rad': return rad_str
        elif sel == 'all': return 'deg =\t' + deg_str + '\nrad =\t' + rad_str
        else: raise Exception(f"-> Error in `sel` parameter!\nsel = {sel} is not allowed\nRead the documentation for correct values of the parameter")

    def copy(self):
        return Angles(self.deg,'deg',self.lim)

    def __add__(self, angle):
        """Function to sum two angles

        :param angle: second angle
        :type angle: Angles

        :return: sum of the angles
        :rtype: Angles
        """
        if type(angle) == float or type(angle) == int:
            angle = Angles(angle,'deg',self.lim)
        # check for edges
        if self.lim != angle.lim:
            print(f'\n!warning: you are summing angles with different limits: {self.lim} and {angle.lim}!\nThe limit of the sum is taken equal to that of ang1\n')
        # computing the sum in rad
        sumrad = self.rad + angle.rad
        return Angles(ang=sumrad,unit='rad',lim=self.lim)

    def __sub__(self, angle):
        """Function to subtract two angles

        :param angle: second angle
        :type angle: Angles

        :return: subtraction of the angles
        :rtype: Angles
        """
        if type(angle) == float or type(angle) == int:
            angle = Angles(angle,'deg',self.lim)
        if self.lim != angle.lim:
            print(f'\n!warning: you are subtracting angles with different limits: {self.lim} and {angle.lim}!\nThe limit of the sum is taken equal to that of ang1\n')
        subrad = self.rad - angle.rad
        return Angles(ang=subrad,unit='rad',lim=self.lim)

    def __mul__(self,val: float | int):
        """Function to implement the angle-number product

        :param val: a number
        :type val: float | int

        :return: angle-number product
        :rtype: Angles
        """
        return Angles(self.rad*val,'rad',lim=self.lim)
    
    def __neg__(self):
        return self * -1



class HAngles(Angles):
    """This is essentially the same class as :class: `Angles`, keeping in account 
    the description of an angle in [hours, minutes, seconds] format.

    The attributes are:

    :ivar rad: from :class: `Angles`; the value of the angle in rad 
    :vartype rad: float
    :ivar deg: from :class: `Angles`; the value of the angle in deg 
    :vartype deg: float    
    :ivar deg: the value of the angle in hours
    :vartype deg: float  
    :ivar lim: from :class: `Angles`; the angle is in a range set by [0,`lim`] (see `std_format()` docstring)
    :vartype lim: int

    .. note::
        It is possible to generate a null angle, that is an empty `HAngles` variable
        (`self.deg = None`, `self.rad = None` and `self.hms = None`), through the 
        command `HAngles(None)`. 
    """

    @staticmethod
    def deg_to_hms(deg: float) -> float:
        """Function to convert deg in hms

        1 hour = 15 degrees

        :param deg: angle in deg
        :type deg: float

        :return: angle in hms
        :rtype: float
        """
        hms = deg / 15
        return hms
    
    @staticmethod
    def hms_to_deg(hms: float) -> float:
        """Function to convert hms in deg

        1 hour = 15 degrees

        :param hms: angle in hms
        :type hms: float

        :return: angle in deg
        :rtype: float
        """
        deg = hms * 15
        return deg
    
    @staticmethod
    def rad_to_hms(rad: float) -> float:
        """Function to convert rad in hms

        It pass through the `rad_to_deg()` and
        `deg_to_hms()` functions

        :param rad: angle value in rad
        :type rad: float

        :return: angle in hms
        :rtype: float
        """
        deg = Angles.rad_to_deg(rad)
        hms = HAngles.deg_to_hms(deg)
        return hms
    
    @staticmethod
    def hms_to_rad(hms: float) -> float:
        """Function to convert hms in rad

        It pass through the `hms_to_deg()` and
        `deg_to_rad()` functions

        :param rad: angle in hms
        :type rad: float

        :return: angle value in rad 
        :rtype: float
        """
        deg = HAngles.hms_to_deg(hms)
        rad = Angles.deg_to_rad(deg)
        return rad
    
    def __init__(self, ang: float | list | None, unit: str, lim: int = 360):
        """Constructor of the class (inherited from :class: `Angles`)

        The function takes a value (`ang`) and the corrisponding 
        unit (`unit`) as input, computes and stores the angle
        values in deg, rad and hms. 

        :param ang: angle value
        :type ang: float | list
        :param unit: unit of angle value, like 'deg' or 'rad'
        :type unit: str
        :param lim: angle edge (see `std_format()` docstring), defaults to 360
        :type lim: int, optional
        """
        # `Angles.__init__()` function
        super().__init__(ang, unit, lim)
        # condition for a null angle 
        if ang is None:
            self.hms = None
        # angle in hms
        elif unit == 'hms':
            # list format
            if type(ang) == list:
                # taking the sign
                sign = Angles.strsign[ang[0]]
                self.hms = Angles.decimal(ang[1])*sign
            else:
                self.hms = ang
            self.deg = HAngles.hms_to_deg(self.hms)
            self.rad = HAngles.hms_to_rad(self.hms)
        # angle in rad or deg
        else:
            self.hms = HAngles.deg_to_hms(self.deg)

    def std_format(self, unit: str = 'deg') -> list:
        """Function to adjust format of a list angle in deg or hms.

        The wanted format is: [`sign`, [deg (<`lim`), pp (<60), ss.ss (<60)]]

        :param unit: 'deg' or 'hms' to convert, defaults to 'deg'
        :type unit: str, optional

        :return: correct format for the angle value
        :rtype: list
        """
        # deg 
        if unit == 'deg' or self.rad == 0.:
            return super().std_format()
        # hms
        elif unit == 'hms':
            # getting value in hms
            ang = self.hms
            # converting the edge
            lim = self.lim / 15
            # taking the sign as a str
            sign = Angles.strsign[np.sign(ang)]
            # taking the absolute value
            ang = abs(ang)
            # edge condition
            if ang > lim:
                ang -= lim * (ang // lim)
            # hours
            val = [np.trunc(ang).astype(int)]
            # minutes
            val += [np.trunc(ang % 1 * 60).astype(int)]
            # seconds
            val += [ang % 1 * 60 % 1 * 60]
            return [sign, np.array(val)]

    
    def str_angle(self, sel: str = 'all', unit: bool = False) -> str:
        """Function to print the value of an angle
        (inherited from :class: `Angles`)

        One can select to print the angle in all 
        units or just one through the `sel`
        parameter:

            - `sel = 'all'`: print value in all units
            - `sel = 'deg'`: print value in deg
            - `sel = 'rad'`: print value in rad
            - `sel = 'hms'`: print value in hms

        :param sel: to select in which unit printing the angle, defaults to 'all'
        :type sel: str, optional

        :return: the string with the values
        :rtype: str
        """
        # adjusting the format in hms
        hms = self.std_format('hms')
        # collecting string for each unit
        hh, mm, ss = hms[1]
        hms_str = hms[0]+f'[{hh:.0f}, {mm:.0f}, {ss:.4f}]'

        if unit:
            hms_str += ' hms'
        
        if sel == 'hms': return hms_str
        elif sel == 'all': return  super().str_angle(sel=sel,unit=unit) + '\nhms = \t' + hms_str
        else: return super().str_angle(sel=sel,unit=unit)

    def copy(self):
        return HAngles(self.deg,'deg',self.lim)
    
    def __add__(self, angle):
        """Function to sum two angles

        :param angle: second angle
        :type angle: Angles

        :return: sum of the angles
        :rtype: Angles
        """
        if type(angle) == float or type(angle) == int:
            angle = HAngles(angle,'deg',self.lim)
        # check for edges
        if self.lim != angle.lim:
            print('-> Warning: you are summing angles with different limits: {self.lim} and {angle.lim}!\nThe limit of the sum is taken equal to that of ang1\n')
        # computing the sum
        sumrad = self.rad + angle.rad
        return HAngles(ang=sumrad,unit='rad',lim=self.lim)

    def __sub__(self, angle):
        """Function to subtract two angles

        :param angle: second angle
        :type angle: HAngles

        :return: subtraction of the angles
        :rtype: HAngles
        """        
        if type(angle) == float or type(angle) == int:
            angle = HAngles(angle,'deg',self.lim)
        # check for edges
        if self.lim != angle.lim:
            print(f'-> Warning: you are subtracting angles with different limits: {self.lim} and {angle.lim}!\nThe limit of the sum is taken equal to that of ang1\n')
        # computing the subtaction
        subrad = self.rad - angle.rad
        return HAngles(ang=subrad,unit='rad',lim=self.lim)

    def __mul__(self,val: float | int):
        """Function to implement the angle-number product

        :param val: a number
        :type val: float | int

        :return: angle-number product
        :rtype: HAngles
        """
        return HAngles(self.rad*val,'rad',lim=self.lim)
    
    def __neg__(self):
        return self * -1



# defining some fixed angles
RIGHT = Angles(90.,'deg')
FLAT = Angles(180.,'deg')


class ArrAngle():
    """This class is used to have an array of :class: `Hangles`.

    Each attribute (deg, rad, hms) is an array of angle values 
    in that unit

    The attributes are:

    :ivar deg: each element of the array in deg
    :vartype deg: np.ndarray[float]
    :ivar rad: each element of the array in rad
    :vartype rad: np.ndarray[float]
    :ivar hms: each element of the array in hms
    :vartype hms: np.ndarray[float]

    .. note::
        It is possible to generate a null arr, that is an empty `ArrAngle` variable
        (`self.deg = np.array([])`, `self.rad = np.array([])` and 
        `self.hms = np.array([])`), through the command `ArrAngle(None)`.  
    """
    def __init__(self, arrang: np.ndarray | list | None, unit: str, lim: int = 360) -> None:
        """Array to store values of an array of :class: `Hangles`

        The function takes an array of values (`arrang`) in a certain unit (`unit`) 
        and after proper convertions builts and stores an array for each unit.

        The allowed `unit` values are the same of :class: `HAngles`: 
        'deg', 'rad', 'hms'

        :param arrang: array or list of values
        :type arrang: np.ndarray | list | None
        :param unit: unit of the values
        :type unit: str
        :param lim: edge of the deg values, defaults to 360
        :type lim: int, optional
        """
        # condition for null array
        if arrang is None or len(unit) == 0:
            self.deg = np.array([]) 
            self.rad = np.array([])
            self.hms = np.array([]) 
            self.lim = lim
        else:
            # builting an array of HAngles
            allvalues = np.array([HAngles(ang,unit=unit,lim=lim) for ang in arrang])
            # storing values for each unit as array
            self.deg = np.array([ val.deg for val in allvalues]) 
            self.rad = np.array([ val.rad for val in allvalues])
            self.hms = np.array([ val.hms for val in allvalues]) 
            self.lim = lim

    def copy(self):
        return ArrAngle([*self.deg],'deg',self.lim)

    def __add__(self, angle):
        if isinstance(angle, (float,int)):
            angle = HAngles(angle, 'deg', self.lim)
        sumangle = self.rad + angle.rad
        return ArrAngle(sumangle,'rad',lim=self.lim)
    
    def __sub__(self, angle):
        if isinstance(angle, (float,int)):
            angle = HAngles(angle, 'deg', self.lim)
        sumangle = self.rad - angle.rad
        return ArrAngle(sumangle,'rad',lim=self.lim)
    
    def __mul__(self, val: float | int):
        mulangle = ArrAngle(None,'')
        mulangle.deg = np.copy(self.deg) * val        
        mulangle.rad = np.copy(self.rad) * val        
        mulangle.hms = np.copy(self.hms) * val        
        mulangle.lim = self.lim
        return mulangle

    def __neg__(self):
        return self * -1        

