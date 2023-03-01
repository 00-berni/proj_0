import numpy as np
from numpy import pi
import matplotlib.pyplot as plt



class Angles():
    """Class to handle with angles in degrees or radiants.

    :param deg: the amplitude of the angle in [°,',''], defaults to np.zeros(3)
    :type deg: numpy_array
    :param rad: the amplitude of the angle in radiants, defaults to 0
    :type rad: float
    """
    def __init__(self, deg = np.zeros(3), rad = 0) -> None:
    
        self.deg = deg
        self.rad = rad

    def deg_check(self) -> None:
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
        ang = self.deg[0] + self.deg[1]*60 + self.deg[2]*3600
        self.rad = ang / 180 * pi

    def rad_to_deg(self) -> None:
        ang = self.rad / pi * 180
        # if(ang > 360):
        #     ang - 360
        self.deg = np.array([ang, 0, 0])
        self.deg_check()

    def angle_init(self) -> None:
        self.deg_check()
        self.deg_to_rad()


class HAngles(Angles):

    def __init__(self, deg = np.zeros(3), rad = 0, hms = np.zeros(3)) -> None:
        super().__init__(deg,rad)
        self.hms = hms

    def hms_check(self) -> None:
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
        self.hms = self.deg / 15
        self.hms_check()

    def hms_to_deg(self) -> None:
        self.deg = self.hms * 15
        self.deg_check()

    def angle_init(self) -> None:
        if all(self.hms == np.zeros(3)) :
            self.deg_to_hms()
        else:
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
    print(delta.deg)



