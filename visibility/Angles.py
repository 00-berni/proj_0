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

    @staticmethod
    def std_format(deg, lim) -> list:
        """Function to adjust format of a list angle in deg or hms.

        The wanted format is: [`sign`, [deg (<`lim`), pp (<60), ss.ss (<60)]]
          
        :return: correct format for the angle value
        :rtype: list
        """
        # getting value in deg
        ang = deg
        if ang == 0.:
            return ['+',np.zeros(3,dtype=int)]
        else:
            # taking the lim
            lim = lim
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
    
    @staticmethod
    def str_angle(deg: float, rad: float, lim: int, sel: str = 'all', unit: bool = False) -> str:
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
        deg = Angles.std_format(deg,lim)
        # collecting string for each unit
        dd, mm, ss = deg[1]
        deg_str = deg[0]+f'[{dd:.0f}, {mm:.0f}, {ss:.4f}]'
        rad_str = f'{rad/pi} * pi'
        
        if unit:
            deg_str += ' deg'
            rad_str += ' rad'

        if sel == 'deg': return deg_str
        elif sel == 'rad': return rad_str
        elif sel == 'all': return 'deg =\t' + deg_str + '\nrad =\t' + rad_str
        else: raise Exception(f"-> Error in `sel` parameter!\nsel = {sel} is not allowed\nRead the documentation for correct values of the parameter")

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

    def copy(self):
        return Angles(self.deg,'deg',self.lim)

    def print_angle(self, sel: str = 'all', unit: bool = False) -> list[str] | str:
        ang = self.copy()
        if type(ang.deg) != np.ndarray:
            ang.deg = [ang.deg]
            ang.rad = [ang.rad]
        deg = ang.deg
        rad = ang.rad
        str_list = []
        for i in range(len(deg)):
            str_res = ''
            if sel == 'all' and len(deg) > 1:
                str_res += f'ang {i}:\n'
            str_res += Angles.str_angle(deg[i], rad[i], ang.lim, sel=sel, unit=unit)
            str_list += [str_res]
        
        if len(deg) == 1:
            str_list = str_res
        return str_list


    def __add__(self, angle):
        """Function to sum two angles

        :param angle: second angle
        :type angle: Angles

        :return: sum of the angles
        :rtype: Angles
        """
        if isinstance(angle, (int,float)):
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
        if isinstance(angle, (int,float)):
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

    @staticmethod
    def str_angle(deg: float, rad: float, hms: float, lim: int, sel: str = 'all', unit: bool = False) -> str:
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
        hms = Angles.std_format(hms,lim/15)
        # collecting string for each unit
        hh, mm, ss = hms[1]
        hms_str = hms[0]+f'[{hh:.0f}, {mm:.0f}, {ss:.4f}]'

        if unit:
            hms_str += ' hms'
        
        if sel == 'hms': return hms_str
        elif sel == 'all': return  Angles.str_angle(deg, rad, lim, sel=sel, unit=unit) + '\nhms = \t' + hms_str
        else: return Angles.str_angle(deg, rad, lim, sel=sel, unit=unit)
    
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

    def copy(self):
        return HAngles(self.deg,'deg',self.lim)
    
    def print_angle(self, sel: str = 'all', unit: bool = False) -> list[str] | str:
        ang = self.copy()
        if type(ang.deg) != np.ndarray:
            ang.deg = [ang.deg]
            ang.rad = [ang.rad]
            ang.hms = [ang.hms]
        deg = ang.deg
        rad = ang.rad
        hms = ang.hms
        str_list = []
        for i in range(len(deg)):
            str_res = ''
            if sel == 'all' and len(deg) > 1:
                str_res += f'ang {i}:\n'
            str_res += HAngles.str_angle(deg[i], rad[i], hms[i], ang.lim, sel=sel, unit=unit)
            str_list += [str_res]
        if len(deg) == 1:
            str_list = str_res
        return str_list
    
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
ROUND = Angles(360.,'deg')


if __name__ == '__main__':
    ang = np.linspace(10,120,120)
    ang = Angles(ang,'deg')
    strlist = ang.print_angle()
    for elem in strlist[:5]:
        print('\n'+elem)
    print('- - -')
    ang = 12.5
    ang = Angles(ang,'deg')
    print(ang.print_angle())
    print('- - -')
    ang = np.linspace(10,120,120)
    ang = HAngles(ang,'deg')
    strlist = ang.print_angle()
    for elem in strlist[:5]:
        print('\n'+elem)
    print('- - -')
    ang = 12.345
    ang = HAngles(ang,'deg')
    print(ang.print_angle())
