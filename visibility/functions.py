import numpy as np
from numpy import pi
import matplotlib.pyplot as plt



class Angles():
    """Class to handle with angles in degrees or radiants.
    The class keeps the value of an angle in deg and rad and gives the functions to convert the one in the other.

    :param deg: the amplitude of the angle in [degrees, primes, seconds], defaults to np.zeros(3)
    :type deg: numpy_array, optional
    :param rad: the amplitude of the angle in radiants, defaults to 0
    :type rad: float, optional
    """
    def __init__(self, deg = np.zeros(3), rad = 0.0) -> None:
        """The constructor of the class"""
        self.deg = deg
        self.rad = rad

    def deg_check(self) -> None:
        """The shape of the angle value in deg is [int (<360), int (<60), float (<60)]. 
        This function checks and adjusts the format of the angle value.

        :rtype: None
        """
        for i in range(2):
            diff = self.deg[i] % 1
            if diff != 0:
                self.deg[i] = int(self.deg[i])
                self.deg[i+1] += diff*60
        self.deg[2] = round(self.deg[2], 3)
        if self.deg[2] >= 60 :
            tmp = self.deg[2]
            self.deg[2] = tmp % 60
            self.deg[1] += tmp // 60
        if self.deg[1] >= 60 :
            tmp = self.deg[1]
            self.deg[1] = tmp % 60
            self.deg[0] += tmp // 60
        if self.deg[0] > 360 :
            self.deg[0] -= 360

    def deg_to_rad(self) -> None:
        """The function converts deg value in rad one of an angle

        :rtype: None
        """
        ang = self.deg[0] + self.deg[1]*60 + self.deg[2]*3600
        self.rad = ang / 180 * pi

    def rad_to_deg(self) -> None:
        """The function converts rad value of an angle in deg one

        :rtype: None
        """
        ang = self.rad / pi * 180
        # if(ang > 360):
        #     ang - 360
        self.deg = np.array([ang, 0, 0])
        self.deg_check()

    def angle_init(self) -> None:
        """The function checks deg value of an angle and assigns the respective rad one 

        :rtype: None
        """
        self.deg_check()
        self.deg_to_rad()


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
    def __init__(self, deg = np.zeros(3), rad = 0, hms = np.zeros(3)) -> None:
        """The constructor function"""
        super().__init__(deg,rad)
        self.hms = hms

    def hms_check(self) -> None:
        """The shape of the angle value in hms is [int, int (<60), float (<60)]. 
        This function checks and adjusts the format of the angle value.
        """
        for i in range(2):
            diff = self.hms[i] % 1
            if diff != 0:
                self.hms[i] = int(self.hms[i])
                self.hms[i+1] += diff*60
        self.hms[2] = round(self.hms[2], 3)
        if self.hms[2] >= 60 :
            tmp = self.hms[2]
            self.hms[2] = tmp % 60
            self.hms[1] += tmp // 60
        if self.hms[1] >= 60 :
            tmp = self.hms[1]
            self.hms[1] = tmp % 60
            self.hms[0] += tmp // 60

    def deg_to_hms(self) -> None:
        """The function converts deg value of an angle in hms one

        :rtype: None
        """
        self.hms = self.deg / 15
        self.hms_check()

    def hms_to_deg(self) -> None:
        """The function converts hms value of an angle in deg one

        :rtype: None
        """
        self.deg = self.hms * 15
        self.deg_check()

    def angle_init(self) -> None:
        """The function checks if angle is given in deg or hms format then assigns the respective values for the other formats 

        :rtype: None
        """
        if all(self.hms == np.zeros(3)) :
            self.deg_check()
            self.deg_to_hms()
        else:
            self.hms_check()
            self.hms_to_deg()
        self.deg_to_rad()


# # # # #

if __name__ == "__main__":

    delta = Angles(deg=np.array([16.2, 58, 0]))
    delta.deg_check()
    delta.deg_to_rad()
    print(delta.deg)
    print(delta.rad)
    delta.rad = 4*pi
    delta.rad_to_deg()
    print(type(delta))



