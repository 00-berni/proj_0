import numpy as np
from visibility.stuff import get_data, interpole_three
from visibility.angles import Angles, HAngles, ArrAngle
from visibility.coor import Equatorial, GeoPos
from visibility.Time import Date, julian_day, local_ST


class StarObj():
    def __init__(self, name: str, coor: Equatorial | list[HAngles], prmt: list[float] | None) -> None:
        self.name = name
        if type(coor) == list:
            ra, dec = coor
            coor = Equatorial(ra,dec)
        self.coor = coor
        self.mua = prmt[0] if prmt is not None else None
        self.mud = prmt[1] if prmt is not None else None
    
    def copy(self):
        prmt = [self.mua, self.mud] if self.mua is not None else None
        return StarObj(self.name,self.coor.copy(),prmt)

    def getcoor(self):
        starobj = self.copy()
        return starobj.coor.alpha, starobj.coor.delta

    def getprmt(self):
        starobj = self.copy()
        return [starobj.mua, starobj.mud]

    def obj_info(self):
        print('Target: ' + self.name)
        print(self.coor.print_values())
        print('Proper motion')
        if self.mua is not None:
            print(f'mu_alpha =\t{self.mua} as/yr\nmu_delta =\t{self.mud} as/yr')

    def coor_in_date(self, date):
        return precession_corr(date,self)
    
    def lha(self, date: Date, lon: Angles):
        ra, _ = self.coor_in_date(date)
        LST = local_ST(date,lon)
        return LST - ra        


def precession_corr(date: Date, obj: StarObj):
    ra0, dec0  = obj.getcoor()
    tmp_date = date.copy()

    J2000 = Date.J2000
    B1950 = Date.B1950
    B1900 = Date.B1900

    if tmp_date.time.tytime == 'UT':
        tmp_date.time = tmp_date.time.change_time_type(tmp_date.date[0],'TD','Time')

    jd = julian_day(tmp_date)

    del tmp_date

    epoch = date.epoch
    if epoch[0] == 'J': 
        if epoch[:5] == 'J2000': jd0 = J2000  
        else: J2000 - (2000 - float(epoch[1:]))*365.25  
    elif epoch[0] == 'B':
        jd0 = B1900 if epoch == 'B1900' else B1950

    t = (jd - jd0) / 36525
    if obj.mua is not None:
        ra0  = ra0  + HAngles(obj.mua/3600,'deg') * (t*100)
        dec0 = dec0 + HAngles(obj.mud/3600,'deg',lim=90) * (t*100)

    xi = 2306.2181 * t + 0.30188 * t**2 + 0.017998 * t**3    
    z  = 2306.2181 * t + 1.09468 * t**2 + 0.018203 * t**3    
    th = 2004.3109 * t - 0.42665 * t**2 - 0.041833 * t**3
    if epoch != 'J2000':
        T = (jd0 - J2000) / 36525        
        xi += ( (1.39656 * T - 1.39e-4 * T**2)*t - 3.44e-4 * T *t**2 )
        z  += ( (1.39656 * T - 1.39e-4 * T**2)*t + 0.66e-4 * T *t**2 )
        th += (-(0.85330 * T + 2.17e-4 * T**2)*t + 2.17e-4 * T *t**2 )
    
    xi = HAngles(xi/3600,'deg')
    z  = HAngles( z/3600,'deg')
    th = HAngles(th/3600,'deg')

    A = np.cos(dec0.rad) * np.sin((ra0 + xi).rad)
    B = np.cos(th.rad) * np.cos(dec0.rad) * np.cos((ra0 + xi).rad) - np.sin(th.rad) * np.sin(dec0.rad)

    ra = np.arctan2(A,B) + z.rad

    C = np.sin(th.rad) * np.cos(dec0.rad) * np.cos((ra0 + xi).rad) + np.cos(th.rad) * np.sin(dec0.rad)

    D = np.sqrt(A**2 + B**2)

    if type(jd) == np.ndarray:
        dec = np.where(dec0.deg < 80,  np.arcsin(C), np.arccos(D))
    else:
        dec = np.arcsin(C)  if dec0.deg < 80 else np.arccos(D)
    
    return HAngles(ra,'rad'), HAngles(dec,'rad',lim=90)


def std_atm(height: float) -> tuple[float,float]:
    # name of data file
    filename = 'atm_data.csv'
    # extracting data
    h,T,p = get_data(filename)
    if height in h:
        idx = np.where(height == h)[0]
        temp = T[idx]
        pres = p[idx]
    else:
        temp = interpole_three(T,height,h)
        pres = interpole_three(p,height,h)
    return temp, pres

def refraction_corr(alt: Angles | ArrAngle, height: float):
    rcorr = 1.02/np.tan( Angles.deg_to_rad(alt.deg + 10.3/(alt.deg + 5.11)) )        
    if height != 0.:
        T, P = std_atm(height)
        T0, P0 = std_atm(0.)
        rcorr *= P/P0 * T0/T 
    if type(alt) == Angles:
        if alt.deg == 90:
            rcorr += 0.0019279
        rcorr = Angles(rcorr/60,'deg')
    elif type(alt) == ArrAngle:
        idx = np.where(alt.deg == 90)[0]
        if len(idx) != 0:
            rcorr[idx] += 0.0019279
        rcorr = ArrAngle(rcorr/60,'deg')
    return rcorr


def trajectory(date: Date, obs_pos: GeoPos, obj: StarObj, numpoint: int = 1000) -> ArrAngle:
    date = date.copy()
    date.time.val += np.linspace(0,24,numpoint)*3600 
    
    lat = obs_pos.lat
    lon = obs_pos.lon
    height = obs_pos.h

    ra, dec = obj.coor_in_date(date)
    print(HAngles(ra.hms[0],'hms').str_angle('hms'))
    # ra = obj.coor.alpha
    # dec = obj.coor.delta
    # print(obj.coor.alpha.str_angle('hms'))
    
    HA = obj.lha(date,lon)

    phi = lat.rad
    delta = dec.rad

    alt = np.arcsin( np.sin(phi)*np.sin(delta) + np.cos(phi)*np.cos(delta) * np.cos(HA.rad)  )

    alt = ArrAngle(alt,'rad',lim=90)
    alt = alt + refraction_corr(alt, height)
    return alt, date.time.hour()
