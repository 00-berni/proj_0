import numpy as np
from numpy import pi
import matplotlib.pyplot as plt


# class angunit():
#     def __init__(self, ang: float | list[str,np.ndarray]) -> None:
#         if type(ang) != list:
#             sign = '+' if ang >= 0 else '-'
#             ang = [sign, np.array([np.abs(ang),0,0])]
#         self.sign = ang[0]
#         self.ang = ang[1]

#     def get_sign(self) -> int:
#         return 1 if self.sign == '+' else -1

#     def convert_sign(self,sign: int):
#         self.sign = '+' if sign > 0 else '-'

#     def print_ang(self):
#         return self.sign + f'{self.ang}'

#     def clear_ang(self):
#         self.sign = '+'
#         self.ang = np.zeros(3)
#         return self


## angle classes
class Angles():
    """Class to handle with angles in degrees or radiants.
    
    The class keeps the value of an angle in deg and rad
    and gives the functions to convert the one in the other.

    :param deg: the amplitude of the angle in [°,',''], defaults to np.zeros(3)
    :type deg: numpy_array, optional
    :param rad: the amplitude of the angle in radiants, defaults to 0.0
    :type rad: float, optional
    """

    strsign = { '+' :  1,
                '-' : -1,
                 1  : '+',
                -1  : '-' }

    @staticmethod
    def decimal(ang: np.ndarray) -> float:
        return ang[0] + ang[1]/60 + ang[2]/3600

    @staticmethod
    def std_format(sign: int | str, ang: np.ndarray | list, lim: int | None = 360) -> list:
        if type(sign) != str:
            sign = Angles.strsign[sign]
        ang = np.copy(ang)
        # checking for decimal
        for i in range(2):
            diff = ang[i] % 1
            if diff != 0:
                ang[i] = int(ang[i])
                ang[i+1] += diff*60
        ang[2] = round(ang[2], 3)
        # getting the [deg,pp,ss.ss] format
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
        if lim != None:
            if ang[0] > lim :
                ang[0] -= lim * np.trunc(ang[0] / lim).astype(int)
        ang = [sign,ang]
        return ang
    
    @staticmethod
    def deg_to_rad(deg: list) -> float:
        sign = Angles.strsign[deg[0]]
        ang = deg[1]
        dec_deg = Angles.decimal(ang)
        rad = sign * dec_deg * pi / 180
        return rad

    @staticmethod
    def rad_to_deg(rad: float, lim: int = 360) -> list:
        sign = np.sign(rad)
        ang = np.abs(rad) * 180 / pi
        deg = Angles.std_format(sign,[ang,0,0],lim=lim)
        return deg

    def __init__(self, ang: float | list, unit: str, lim: int = 360):
        self.lim = lim
        if unit == 'rad':
            self.rad = ang
            self.deg = Angles.rad_to_deg(ang,lim=lim)
        elif unit == 'deg':
            if type(ang) != list:
                sign = np.sign(ang)
                self.deg = Angles.std_format(sign,[np.abs(ang),0,0],lim=lim)
            else:
                self.deg = [ang[0], np.copy(ang[1])]
            self.rad = Angles.deg_to_rad(self.deg)
            
    def print_angle(self):
        deg_str = 'deg =\t'+self.deg[0]+f'{self.deg[1]}\n'
        rad_str = 'rad =\t'+f'{self.rad/pi} * pi\n'
        return deg_str + rad_str

    def __add__(self, angle):
        if self.lim != angle.lim:
            print('\n!warning: you are summing angles with different limits!\nThe limit of the sum is taken equal to that of ang1\n')
        sumrad = self.rad + angle.rad
        return Angles(ang=sumrad,unit='rad',lim=self.lim)

    def __mul__(self,val: float | int):
        return Angles(self.rad*val,'rad',lim=self.lim)

    def __sub__(self, angle):
        if self.lim != angle.lim:
            print('\n!warning: you are summing angles with different limits!\nThe limit of the sum is taken equal to that of ang1\n')
        subrad = self.rad - angle.rad
        return Angles(ang=subrad,unit='rad',lim=self.lim)


class HAngles(Angles):
    """This is essentially the same class as :class: `Angles`, keeping in account the description of an angle in [hours,minutes,seconds] format.

    :param deg: parameter from :class: `Angles`
        the amplitude of the angle in [degrees, primes, seconds], defaults to np.zeros(3)
    :type deg: numpy_array, optional
    :param rad: parameter from :class: `Angles`
        the amplitude of the angle in radiants, defaults to 0
    :type rad: float, optional
    :param hms: the amplitude of the angle in [h, m, s], defaults to np.zeros(3)
    :type deg: numpy_array, optional
    """

    @staticmethod
    def deg_to_hms(deg: list) -> list:
        sign = deg[0]
        ang = np.copy(deg[1])
        ang = Angles.decimal(ang) / 15
        hms = Angles.std_format(sign,[ang,0,0],lim=None)
        return hms
    
    @staticmethod
    def hms_to_deg(hms: list, lim: int = 360) -> list:
        sign = hms[0]
        ang = np.copy(hms[1])
        ang = Angles.decimal(ang) * 15
        deg = Angles.std_format(sign,[ang,0,0],lim=lim)
        return deg
    
    @staticmethod
    def rad_to_hms(rad: float, lim: int = 360) -> list:
        deg = Angles.rad_to_deg(rad,lim=lim)
        hms = HAngles.deg_to_hms(deg)
        return hms
    
    @staticmethod
    def hms_to_rad(hms: list, lim: int = 360) -> float:
        deg = HAngles.hms_to_deg(hms,lim=lim)
        rad = Angles.deg_to_rad(deg)
        return rad
    
    def __init__(self, ang: float | list, unit: str, lim: int = 360):
        super().__init__(ang, unit, lim)
        if unit == 'hms':
            if type(ang) != list:
                sign = np.sign(ang)
                self.hms = Angles.std_format(sign,[np.abs(ang),0,0],lim=None)
            else:
                self.hms = [ang[0], np.copy(ang[1])]
            self.deg = HAngles.hms_to_deg(self.hms,lim=lim)
            self.rad = HAngles.hms_to_rad(self.hms)
        else:
            self.hms = HAngles.deg_to_hms(self.deg)
    
    def print_angle(self):
        ang_str = super().print_angle()
        hms_str = 'hms =\t'+self.hms[0]+f'{self.hms[1]}\n'
        return ang_str + hms_str

    def __add__(self, angle):
        if self.lim != angle.lim:
            print('\n!warning: you are summing angles with different limits!\nThe limit of the sum is taken equal to that of ang1\n')
        sumrad = self.rad + angle.rad
        return HAngles(ang=sumrad,unit='rad',lim=self.lim)

    def __mul__(self,val: float | int):
        return HAngles(self.rad*val,'rad',lim=self.lim)

    def __sub__(self, angle):
        if self.lim != angle.lim:
            print('\n!warning: you are summing angles with different limits!\nThe limit of the sum is taken equal to that of ang1\n')
        subrad = self.rad - angle.rad
        return HAngles(ang=subrad,unit='rad',lim=self.lim)
        

if __name__ == "__main__":
    SEP = '-------'
    print(SEP+'TEST ANGLE CLASSES'+SEP+'\n')
    try:

        ang1 = HAngles(['+',[12,5,7]],'hms')
        ang2 = HAngles(30,'rad')
        print('Start')
        print('ang1\n' + ang1.print_angle())
        print('ang2\n' + ang2.print_angle())

        sumang = ang1 + ang2
        subang = ang1 - ang2
        mulang1 = ang1 * 2
        mulang2 = ang2 * -1


        print('ang1 + ang2\n' + sumang.print_angle())
        print('ang1 - ang2\n' + subang.print_angle())
        print('ang1 * 2\n' + mulang1.print_angle())
        print('ang2 * -1\n' + mulang2.print_angle())

        print('\nTEST COMPLETE!\n'+SEP)
    except:
        print('TEST FAILD!\n'+SEP)
        raise