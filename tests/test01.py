import numpy as np
from numpy import pi
import matplotlib.pyplot as plt



class Angles():

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
        # if self.deg[0] > 360 :
        #     self.deg[0] -= 360

    def deg_to_rad(self) -> None:
        ang = self.deg[0] + self.deg[1]*60 + self.deg[2]*3600
        self.rad = ang / 180 * pi

    def rad_to_deg(self) -> None:
        ang = self.rad / pi * 180
        # if(ang > 360):
        #     ang - 360
        self.deg[0] = ang
        self.deg[1:] = np.zeros(2)

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

    def hms_to_deg(self) -> None:
        self.deg = self.hms * 15

    def angle_init(self) -> None:
        if all(self.hms == np.zeros(3)) :
            self.deg_to_hms()
            self.hms_check()
        else:
            self.hms_to_deg()
            self.deg_check()
        self.deg_to_rad()



if __name__ == "__main__":

    HA = np.array([20.12, 15, 63.12])
    HA = HAngles(hms=HA)
    HA.hms_to_deg()
    HA.hms_check()
    print(HA.deg)
    HA.deg_check()

    print(HA.hms, HA.deg)

    delta = Angles(deg=np.array([16.2, 58, 0]))
    delta.deg_check()

    print(delta.deg)