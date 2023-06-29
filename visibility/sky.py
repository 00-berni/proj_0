import numpy as np
from visibility.stuff import get_data, interpole_three
from visibility.Angles import Angles, HAngles, ArrAngle
from visibility.coor import Equatorial, GeoPos
from visibility.Time.Tclasses import * 
from visibility.Time.dates import julian_day, local_ST, mean_Green_HA


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
        if len(idx) == 1: idx = idx[0]
        temp = T[idx]
        pres = p[idx]
    else:
        temp = interpole_three(T,height,h)
        pres = interpole_three(p,height,h)
    return temp, pres

def refraction_corr(alt: Angles | ArrAngle, height: float, alt0: bool = False) -> HAngles:
    if alt0:
        a,b,c =  1, 7.31, 4.4
        cor90 = 1.3515e-3
    else:
        a,b,c =  1.02, 10.3, 5.11
        cor90 = 1.9279e-3
    rcorr = a/np.tan( Angles.deg_to_rad(alt.deg + b/(alt.deg + c)) )        
    if height != 0.:
        T, P = std_atm(height)
        T0, P0 = std_atm(0.)
        rcorr *= P/P0 * T0/T 
    if type(alt) == Angles:
        if alt.deg == 90:
            rcorr += cor90
        rcorr = Angles(rcorr/60,'deg',lim=90)
    elif type(alt) == ArrAngle:
        idx = np.where(alt.deg == 90)[0]
        if len(idx) != 0:
            rcorr[idx] += cor90
        rcorr = ArrAngle(rcorr/60,'deg',lim=90)
    if alt0: rcorr = -rcorr
    return rcorr


def compute_alt(date: Date, obs_pos: GeoPos, obj: StarObj) -> HAngles:
    date = date.copy()
    
    lat = obs_pos.lat
    lon = obs_pos.lon
    height = obs_pos.h

    _, dec = obj.coor_in_date(date)
       
    HA = obj.lha(date,lon)

    phi = lat.rad
    delta = dec.rad

    alt = np.arcsin( np.sin(phi)*np.sin(delta) + np.cos(phi)*np.cos(delta) * np.cos(HA.rad)  )
    if type(alt) == np.ndarray:
        alt = ArrAngle(alt,'rad',lim=90)
    else:
        alt = Angles(alt,'rad',lim=90)
    return alt + refraction_corr(alt, height)



def trajectory(date: Date, obs_pos: GeoPos, obj: StarObj, numpoint: int = 1000) -> ArrAngle:
    date = date.copy()
    date.time.val += np.linspace(0,24,numpoint)*3600 

    alt = compute_alt(date,obs_pos,obj)

    # alt = ArrAngle(alt,'rad',lim=90)
    return alt, date.time.hour()



def tran_ris_set(date: Date, obs_pos: GeoPos, obj: StarObj, results: bool = False, iter: int = 3) -> Time:
        date = date.copy()
        lon = obs_pos.lon
        lat = obs_pos.lat
        h0 = refraction_corr(Angles(0.,'deg',lim=90),obs_pos.h, alt0=True)

        tmpdate = Date(date.date)
        Dt = time_correction(tmpdate.date[0])
        GHA = mean_Green_HA(tmpdate)
        tmpdate.time.tytime = 'TD'
        day = tmpdate.date[-1]
        a2, d2 = precession_corr(tmpdate,obj)
        tmpdate.date[-1] = day - 1
        a1, d1 = precession_corr(tmpdate,obj)
        tmpdate.date[-1] = day + 1
        a3, d3 = precession_corr(tmpdate,obj)
        del tmpdate

        cosH0 = (np.sin(h0.rad) - np.sin(lat.rad)*np.sin(d2.rad)) / (np.cos(lat.rad)*np.cos(d2.rad))
        
        # transit
        mt = (a2 + lon - GHA).deg / 360        
        for k in range(iter):
            if abs(mt) > 1:
                mt -= np.sign(mt)
            if mt < 0:
                mt += 1
            LST = GHA.deg + 360.985647*mt
            if LST > 360:
                LST -= 360*(LST//360)
            LST = HAngles(LST,'deg')
            n = mt + Dt/86400
            a = interpole_three([a1.deg,a2.deg,a3.deg],n+day,[day-1,day,day+1])
            a = HAngles(a,'deg')
            H = LST - lon - a
            Dmt = - H.deg/360
            mt += Dmt
            if abs(mt) > 1:
                mt -= np.sign(mt)
            if mt < 0:
                mt += 1
        time = date.time
        if mt*24 < time.hour():
            mt += 1
        m = Time(mt*86400)

        if results:
            print()
            transit = Date(date.date,m,epoch=date.epoch)

            if transit.time.val > 86400:
                days = transit.time.val // 86400
                transit.date[-1] += days
                transit.time.val -= days * 86400
            elif transit.time.val < time.val:
                transit.date[-1] += 1
            print('transit:\t' + transit.str_date())

        # rising and setting
        if abs(cosH0) <= 1:
            H0 = HAngles(np.arccos(cosH0),'rad',lim=180)
            mr = mt - H0.deg/360
            ms = mt + H0.deg/360
            m = np.array([mr,ms])
            for k in range(iter):
                LST = []
                for i in range(2):
                    if abs(m[i]) > 1:
                        m[i] -= np.sign(m[i])
                    if m[i] < 0:
                        m[i] += 1
                    LST += [GHA.deg + 360.985647*m[i]]
                    if LST[i] > 360:
                        LST[i] -= 360*(LST[i]//360)
                LST = HAngles(np.array(LST),'deg')
                n = m + Dt/86400
                a = np.array([interpole_three([a1.deg,a2.deg,a3.deg],ni+day,[day-1,day,day+1]) for ni in n])
                a = HAngles(a,'deg')
                d = np.array([interpole_three([d1.deg,d2.deg,d3.deg],ni+day,[day-1,day,day+1]) for ni in n])
                d = HAngles(d,'deg',lim=90) 
                H = LST - lon - a
                h = np.arcsin(np.sin(lat.rad)*np.sin(d.rad) + np.cos(lat.rad)*np.cos(d.rad)*np.cos(H.rad))
                h = Angles(h,'rad',lim=90)
                Dm = (h-h0).deg / (360 * (np.cos(d.rad)*np.cos(lat.rad)*np.sin(H.rad)))
                m += Dm
                for i in range(2):
                    if abs(m[i]) > 1:
                        m[i] -= np.sign(m[i])
                    if m[i] < 0:
                        m[i] += 1
            m = Time(m*86400)
            m.val = np.where(m.hour() < time.hour(), m.val + 86400, m.val)
            if results:
                mr, ms = m.val
                rising  = Date(date.date,Time(mr),epoch=date.epoch)
                setting = Date(date.date,Time(ms),epoch=date.epoch)

                event = [rising,setting]
                names = ['rising ','setting']
                for i in range(2):
                    if event[i].time.val > 86400:
                        days = event[i].time.val // 86400
                        event[i].date[-1] += days
                        event[i].time.val -= days * 86400
                    elif event[i].time.val < time.val:
                        event[i].date[-1] += 1
                    print(names[i] + ':\t' + event[i].str_date())
            m.val = np.append(mt*86400,m.val)
        return m
