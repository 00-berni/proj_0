import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from test_struc import *


## angle classes
class Angles():
    """Class to handle with angles in degrees or radiants.
    
    The class takes the value of an angle in deg or rad,
    converts one in other and stored the values.

    The format for values in def is a `list` of length 2:
    
      - the first element is a `str` that keeps the sign of the angle

      - the second one is a `numpy.ndarray` of length 3: [deg, pp, ss.ss] 

    The attributes are:

    :ivar rad: the value of the angle in radiants
    :vartype rad: float
    :ivar deg: the value of the angle in deg. The format is [sign, [deg,pp,ss.ss]]
    :vartype deg: list[str, np.ndarray]    
    :ivar lim: the angle is in a range set by [0,`lim`] (see `std_format()` docstring)
    :vartype lim: int
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
        return ang[0] + ang[1]/60 + ang[2]/3600

    @staticmethod
    def std_format(sign: int | str, ang: np.ndarray | list, lim: int | None = 360) -> list:
        """Function to adjust format of a list angle in deg or hms.

        The wanted format is:

          - [`sign`, [deg (<`lim`), pp (<60), ss.ss (<60)]]

          - [`sign`, [hh, mm (<60), ss.ss (<60)]]

        It is possible to change the range of the angle: [0, `lim`]
          
        :param sign: sign of the angle
        :type sign: int | str
        :param ang: angle value
        :type ang: np.ndarray | list
        :param lim: the angle edge , defaults to 360
        :type lim: int | None, optional
        
        :return: correct format for the angle value
        :rtype: list
        """
        # taking the sign as a str
        if type(sign) != str:
            sign = Angles.strsign[sign]
        # copying the angle values
        ang = np.copy(ang)
        # checking for decimal in degrees and primes
        for i in range(2):
            diff = ang[i] % 1
            if diff != 0:
                ang[i] = int(ang[i])
                ang[i+1] += diff*60
        # rounding the seconds 
        ang[2] = round(ang[2], 4)
        # getting the wanted format (see the docstring)
        if ang[2] >= 60 :
            tmp = ang[2]
            ang[2] = tmp % 60
            ang[1] += tmp // 60
            del tmp
        if ang[1] >= 60 :
            tmp = ang[1]
            ang[1] = tmp % 60
            ang[0] += tmp // 60
            del tmp
        # condition only for deg    
        if lim != None:
            if ang[0] > lim :
                ang[0] -= lim * np.trunc(ang[0] / lim).astype(int)
        # collecting the values
        ang = [sign,ang]
        return ang
    
    @staticmethod
    def deg_to_rad(deg: list) -> float:
        """Function to convert deg in rad

        :param deg: list angle in deg
        :type deg: list

        :return: value in radiants
        :rtype: float
        """
        # getting the sign as int
        sign = Angles.strsign[deg[0]]
        # extracting the array values
        ang = deg[1]
        # converting in float
        dec_deg = Angles.decimal(ang)
        # computing radiants
        rad = sign * dec_deg * pi / 180
        return rad

    @staticmethod
    def rad_to_deg(rad: float, lim: int = 360) -> list:
        """Function to convert rad in deg

        :param rad: value in radiants
        :type rad: float
        :param lim: edge of the angle range, defaults to 360
        :type lim: int, optional
        
        :return: list angle in deg
        :rtype: list
        """
        # getting the sign
        sign = np.sign(rad)
        # computing degrees
        ang = np.abs(rad) * 180 / pi
        # adjusting the format
        deg = Angles.std_format(sign,[ang,0,0],lim=lim)
        return deg

    def __init__(self, ang: float | list, unit: str, lim: int = 360) -> None:
        """Constructor of the class

        The function takes a value (`ang`) and the corrisponding 
        unit (`unit`) as input, computes and stores the angle
        values in deg and rad. 

        :param ang: angle value
        :type ang: float | list
        :param unit: unit of angle value, like 'deg' or 'rad'
        :type unit: str
        :param lim: angle edge (see `std_format()` docstring), defaults to 360
        :type lim: int, optional
        """
        # setting lim
        self.lim = lim
        # value in radiants
        if unit == 'rad':
            self.rad = ang
            self.deg = Angles.rad_to_deg(ang,lim=lim)
        # value in degrees
        elif unit == 'deg':
            # `ang` is a float
            if type(ang) != list:
                # getting the sign
                sign = np.sign(ang)
                self.deg = Angles.std_format(sign,[np.abs(ang),0,0],lim=lim)
            else:
                self.deg = [ang[0], np.copy(ang[1])]
            self.rad = Angles.deg_to_rad(self.deg)
            
    def print_angle(self,sel: str = 'all') -> str:
        """Function to print the value of an angle
        in all units

        One can select to print the angle in all 
        units or just one through the `sel`
        parameter:

            - `sel = 'all'`: print value in all units
            - `sel = 'deg'`: print value in deg
            - `sel = 'rad'`: print value in rad

        :param sel: to select in which unit printing the angle, defaults to 'all'
        :type sel: str, optional

        :return: the string with the values
        :rtype: str
        """
        deg_str = 'deg =\t'+self.deg[0]+f'{self.deg[1]}\n'
        rad_str = 'rad =\t'+f'{self.rad/pi} * pi\n'
       
        if sel == 'deg': return deg_str
        elif sel == 'rad': return rad_str
        elif sel == 'all': return deg_str + rad_str
        else: raise Exception(f"!Error in `sel` parameter!\nsel = {sel} is not allowed\nRead the documentation for correct values of the parameter")

    def __add__(self, angle):
        """Function to sum two angles

        :param angle: second angle
        :type angle: Angles

        :return: sum of the angles
        :rtype: Angles
        """
        # check for edges
        if self.lim != angle.lim:
            print('\n!warning: you are summing angles with different limits!\nThe limit of the sum is taken equal to that of ang1\n')
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
        if self.lim != angle.lim:
            print('\n!warning: you are summing angles with different limits!\nThe limit of the sum is taken equal to that of ang1\n')
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



class HAngles(Angles):
    """This is essentially the same class as :class: `Angles`, keeping in account 
    the description of an angle in [hours, minutes, seconds] format.

    The attributes are:

    :ivar rad: from :class: `Angles`; the value of the angle in radiants; 
    :vartype rad: float
    :ivar deg: from :class: `Angles`; the value of the angle in deg. The format is [sign, [deg,pp,ss.ss]] 
    :vartype deg: list[str, np.ndarray]    
    :ivar deg: the value of the angle in hms. The format is [sign, [hh,mm,ss.ss]]
    :vartype deg: list[str, np.ndarray]    
    :ivar lim: from :class: `Angles`; the angle is in a range set by [0,`lim`] (see `std_format()` docstring)
    :vartype lim: int
    """

    @staticmethod
    def deg_to_hms(deg: list) -> list:
        """Function to convert deg in hms

        1 hour = 15 degrees

        :param deg: list angle in deg
        :type deg: list

        :return: list angle in hms
        :rtype: list
        """
        # getting the sign
        sign = deg[0]
        ang = np.copy(deg[1])
        # converting in float and computing the hours
        ang = Angles.decimal(ang) / 15
        hms = Angles.std_format(sign,[ang,0,0],lim=None)
        return hms
    
    @staticmethod
    def hms_to_deg(hms: list, lim: int = 360) -> list:
        """Function to convert hms in deg

        1 hour = 15 degrees

        :param hms: list angle in hms
        :type hms: list
        :param lim: edge of the angle range, defaults to 360
        :type lim: int, optional

        :return: list angle in deg
        :rtype: list
        """
        # getting the sign
        sign = hms[0]
        ang = np.copy(hms[1])
        # converting in float and computing the degrees
        ang = Angles.decimal(ang) * 15
        deg = Angles.std_format(sign,[ang,0,0],lim=lim)
        return deg
    
    @staticmethod
    def rad_to_hms(rad: float, lim: int = 360) -> list:
        """Function to convert rad in hms

        It pass through the `rad_to_deg()` and
        `deg_to_hms()` functions

        :param rad: angle value in rad
        :type rad: float
        :param lim: edge of the angle range, defaults to 360
        :type lim: int, optional

        :return: list angle in hms
        :rtype: list
        """
        deg = Angles.rad_to_deg(rad,lim=lim)
        hms = HAngles.deg_to_hms(deg)
        return hms
    
    @staticmethod
    def hms_to_rad(hms: list, lim: int = 360) -> float:
        """Function to convert hms in rad

        It pass through the `hms_to_deg()` and
        `deg_to_rad()` functions

        :param rad: list angle in hms
        :type rad: list
        :param lim: edge of the angle range, defaults to 360
        :type lim: int, optional

        :return: angle value in rad 
        :rtype: float
        """
        deg = HAngles.hms_to_deg(hms,lim=lim)
        rad = Angles.deg_to_rad(deg)
        return rad
    
    def __init__(self, ang: float | list, unit: str, lim: int = 360):
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
        # angle in hms
        if unit == 'hms':
            # `ang` is a float
            if type(ang) != list:
                sign = np.sign(ang)
                self.hms = Angles.std_format(sign,[np.abs(ang),0,0],lim=None)
            else:
                self.hms = [ang[0], np.copy(ang[1])]
            self.deg = HAngles.hms_to_deg(self.hms,lim=lim)
            self.rad = HAngles.hms_to_rad(self.hms)
        # angle in rad or deg
        else:
            self.hms = HAngles.deg_to_hms(self.deg)
    
    def print_angle(self, sel: str = 'all'):
        """Function to print the value of an angle
        in all units (inherited from :class: `Angles`)

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
        hms_str = 'hms =\t'+self.hms[0]+f'{self.hms[1]}\n'
        
        if sel == 'hms': return hms_str
        elif sel == 'all': return  super().print_angle(sel=sel) + hms_str
        else: return super().print_angle(sel=sel)

    def __add__(self, angle):
        """Function to sum two angles

        :param angle: second angle
        :type angle: Angles

        :return: sum of the angles
        :rtype: Angles
        """
        # check for edges
        if self.lim != angle.lim:
            print('\n!warning: you are summing angles with different limits!\nThe limit of the sum is taken equal to that of ang1\n')
        # computing 
        sumrad = self.rad + angle.rad
        return HAngles(ang=sumrad,unit='rad',lim=self.lim)

    def __sub__(self, angle):
        """Function to subtract two angles

        :param angle: second angle
        :type angle: HAngles

        :return: subtraction of the angles
        :rtype: HAngles
        """        
        if self.lim != angle.lim:
            print('\n!warning: you are summing angles with different limits!\nThe limit of the sum is taken equal to that of ang1\n')
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


# defining some angles
RIGHT = Angles(90.,'deg')
FLAT = Angles(180.,'deg')
        

if __name__ == "__main__":

    starting_test('TEST ANGLE CLASSES')

    try:

        ang1 = HAngles(['+',[12,5,7]],'hms')
        ang2 = HAngles(30,'rad')
        print('> Define two angles')
        print('ang1\n' + ang1.print_angle())
        print('ang2\n' + ang2.print_angle())

        sumang = ang1 + ang2
        subang = ang1 - ang2
        mulang1 = ang1 * 2
        mulang2 = ang2 * -1

        print('> Compute operations')
        print('ang1 + ang2\n' + sumang.print_angle())
        print('ang1 - ang2\n' + subang.print_angle())
        print('ang1 * 2\n' + mulang1.print_angle())
        print('ang2 * -1\n' + mulang2.print_angle())

        ending_test()
    except:
        test_error()