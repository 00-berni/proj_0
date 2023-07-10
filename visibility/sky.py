import numpy as np
from visibility.stuff import get_data, interpole_three
from visibility.Angles import Angles, HAngles
from visibility.coor import Equatorial, GeoPos
from visibility.Time.Tclasses import * 
from visibility.Time.dates import local_ST, mean_Green_HA


def precession_corr(date: Date, ra0: HAngles, dec0: HAngles, prmt: list[float] | None = None, sel: str = 'all') -> HAngles | tuple[HAngles, HAngles]:
    tmp_date = date.copy()

    J2000 = Date.J2000
    B1950 = Date.B1950
    B1900 = Date.B1900

    if tmp_date.time.tytime == 'UT':
        tmp_date = tmp_date.change_time_type('TD')

    jd = tmp_date.jd

    del tmp_date

    epoch = date.epoch
    if epoch[0] == 'J': 
        if epoch[:5] == 'J2000': jd0 = J2000  
        else: J2000 - (2000 - float(epoch[1:]))*365.25  
    elif epoch[0] == 'B':
        jd0 = B1900 if epoch == 'B1900' else B1950

    t = (jd - jd0) / 36525
    if prmt is not None:
        mua, mud = prmt
        ra0  = ra0  + HAngles(mua/3600,'deg') * (t*100)
        dec0 = dec0 + HAngles(mud/3600,'deg',lim=90) * (t*100)

    # all quantities are in arcsec
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

    if sel == 'all' or sel == 'ra':
        ra = np.arctan2(A,B) + z.rad

    if sel == 'all' or sel == 'dec':
        C = np.sin(th.rad) * np.cos(dec0.rad) * np.cos((ra0 + xi).rad) + np.cos(th.rad) * np.sin(dec0.rad)

        D = np.sqrt(A**2 + B**2)

        if type(jd) == np.ndarray:
            dec = np.where(dec0.deg < 80,  np.arcsin(C), np.arccos(D))
        else:
            dec = np.arcsin(C)  if dec0.deg < 80 else np.arccos(D)
    
    if sel == 'ra':    return HAngles(ra,'rad')
    elif sel == 'dec': return HAngles(dec,'rad',lim=90)
    elif sel == 'all': return HAngles(ra,'rad'), HAngles(dec,'rad',lim=90)

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

def refraction_corr(alt: Angles, height: float, alt0: bool = False) -> Angles:
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
    if type(alt.deg) != np.ndarray and alt.deg == 90:
        rcorr += cor90
    elif type(alt.deg) == np.ndarray:
        rcorr = np.where(alt.deg == 90, rcorr + cor90, rcorr)

    if alt0: rcorr = -rcorr
    return Angles(rcorr/60,'deg',lim=90)


class StarObj():
    def __init__(self, name: str, coor: Equatorial | list[HAngles | float | list], prmt: list[float] | None = None) -> None:
        if type(coor) == list:
            ra, dec = coor
            coor = Equatorial(ra,dec)
        
        self.name = name
        self.coor = coor.copy()
        self.mua = prmt[0] if prmt is not None else None
        self.mud = prmt[1] if prmt is not None else None
    
    def copy(self):
        prmt = [self.mua, self.mud] if self.mua is not None else None
        return StarObj(self.name,self.coor,prmt)

    def obj_info(self) -> None:
        print('Target: ' + self.name)
        print(self.coor.print_values())
        print('Proper motion')
        if self.mua is not None:
            print(f'mu_alpha =\t{self.mua} as/yr\nmu_delta =\t{self.mud} as/yr')
        print()

    def getcoor(self, sel: str = 'all') -> HAngles | tuple[HAngles]:
        starobj = self.copy()
        if   sel == 'ra' : return starobj.coor.alpha
        elif sel == 'dec': return starobj.coor.delta
        elif sel == 'all': return starobj.coor.alpha, starobj.coor.delta
        else: raise Exception(f'Wrong selection!\n`{sel}` is not allowed.')

    def getprmt(self, sel: str = 'all') -> float | list[float] | None:
        starobj = self.copy()
        if self.mua is None and self.mud is None: return None 
        elif sel == 'mua': return starobj.mua
        elif sel == 'mud': return starobj.mud
        elif sel == 'all': return [starobj.mua, starobj.mud]
        else: raise Exception(f'Wrong selection!\n`{sel}` is not allowed.')

    def coor_in_date(self, date: Date, sel: str = 'all') -> HAngles | tuple[HAngles]:
        ra, dec = self.getcoor()
        prmt = self.getprmt()
        return precession_corr(date,ra,dec,prmt,sel=sel)
    
    def lha(self, date: Date, lon: Angles):
        ra = self.coor_in_date(date, sel='ra')
        LST = local_ST(date,lon)
        return LST - ra        

    # def objcoor(self, sel: str = 'all', date: Date | None = None, lon: Angles | None = None) -> HAngles | tuple[HAngles]:
        
    #     if   sel == 'ra' or sel == 'dec' or sel == 'all' : return self.getcoor(sel=sel)
    #     elif sel == 'ha' and date is not None and lon is not None: 
    #         starobj = self.copy() 
    #         return starobj.lha(date, lon)
    #     elif sel == 'ha' and (date is None or lon is None):
    #         raise Exception('Miss anrguments `date` and/or `lon`!')  
    #     else: raise Exception(f'Wrong selection!\n`{sel}` is not allowed.')




def compute_alt(date: Date, obs_pos: GeoPos, obj: StarObj, refcor: bool = False) -> Angles:
    date = date.copy()
    
    lat = obs_pos.lat
    lon = obs_pos.lon
    height = obs_pos.h

    dec = obj.coor_in_date(date,sel='dec')
       
    HA = obj.lha(date,lon).rad

    phi = lat.rad
    delta = dec.rad

    alt = np.arcsin( np.sin(phi)*np.sin(delta) + np.cos(phi)*np.cos(delta) * np.cos(HA)  )
    alt = Angles(alt,'rad',lim=90)
    if refcor:
        alt += refraction_corr(alt, height)
    return alt 



def trajectory(date: Date, obs_pos: GeoPos, obj: StarObj, numpoint: int = 1000) -> tuple[Angles, Date]:

    dayrange = date.jd + np.linspace(0,1,numpoint)
    dayrange = Date(jd=dayrange,timetype=date.time.tytime,calendar=date.calendar,epoch=date.epoch)

    alt = compute_alt(dayrange,obs_pos,obj,refcor=True)

    # alt = ArrAngle(alt,'rad',lim=90)
    return alt, dayrange.jd



def tran_ris_set(date: Date, obs_pos: GeoPos, obj: StarObj, results: bool = False, iter: int = 3) -> Time:
        date = date.copy()
        lon = obs_pos.lon
        lat = obs_pos.lat
        h0 = refraction_corr(Angles(0.,'deg',lim=90),obs_pos.h, alt0=True)

        tmpdate = Date(date.date,0.)
        Dt = time_correction(tmpdate.date[0])
        GHA = mean_Green_HA(tmpdate)
        
        tmpdate = Date(date.date,0.,timetype='TD')
        jd = tmpdate.jd
        a1, d1 = obj.coor_in_date(tmpdate-1)
        a2, d2 = obj.coor_in_date(tmpdate)
        a3, d3 = obj.coor_in_date(tmpdate+1)
        del tmpdate

        cosH0 = (np.sin(h0.rad) - np.sin(lat.rad)*np.sin(d2.rad)) / (np.cos(lat.rad)*np.cos(d2.rad))
        
        # transit
        mt = (a2 + lon - GHA).deg / 360        
        for k in range(iter):
            if abs(mt) > 1:
                mt -= np.sign(mt)
            elif mt < 0:
                mt += 1
            LST = GHA.deg + 360.985647*mt
            if LST > 360:
                LST -= 360*(LST//360)
            LST = HAngles(LST,'deg')
            n = mt + Dt/Time.DAYSEC
            a = interpole_three([a1.deg,a2.deg,a3.deg],n+jd,[jd-1,jd,jd+1])
            a = HAngles(a,'deg')
            H = LST - lon - a
            Dmt = - H.deg/360
            mt += Dmt
            if abs(mt) > 1:
                mt -= np.sign(mt)
            elif mt < 0:
                mt += 1

        time = date.time
        if mt*24 < time.hour():
            mt += 1
        m = Time(mt*Time.DAYSEC)

        if results:
            print()
            transit = Date(date.date,m,calendar=date.calendar,epoch=date.epoch)

            # if transit.time.val > Time.DAYSEC:
            #     days = transit.time.val // Time.DAYSEC
            #     transit = transit + days
            # elif transit.time.val < time.val:
            #     transit = transit + 1
            print('transit:\t' + transit.print_date())

        # rising and setting
        if abs(cosH0) <= 1:
            H0 = HAngles(np.arccos(cosH0),'rad',lim=180)
            mr = mt - H0.deg/360
            ms = mt + H0.deg/360
            m = np.array([mr,ms])
            for k in range(iter):
                m = np.where(abs(m) > 1, m-np.sign(m), m)
                m = np.where(m < 0, m+1, m)
                LST = GHA.deg + 360.985647*m
                LST = np.where(LST > 360, LST - 360*(LST//360), LST)               
                LST = HAngles(LST,'deg')
                n = m + Dt/Time.DAYSEC
                a = np.array([interpole_three([a1.deg,a2.deg,a3.deg],ni+jd,[jd-1,jd,jd+1]) for ni in n])
                a = HAngles(a,'deg')
                d = np.array([interpole_three([d1.deg,d2.deg,d3.deg],ni+jd,[jd-1,jd,jd+1]) for ni in n])
                d = HAngles(d,'deg',lim=90) 
                H = LST - lon - a
                h = np.arcsin(np.sin(lat.rad)*np.sin(d.rad) + np.cos(lat.rad)*np.cos(d.rad)*np.cos(H.rad))
                h = Angles(h,'rad',lim=90)
                Dm = (h-h0).deg / (360 * (np.cos(d.rad)*np.cos(lat.rad)*np.sin(H.rad)))
                m += Dm
                m = np.where(abs(m) > 1, m-np.sign(m), m)
                m = np.where(m < 0, m+1, m)
            
            m = np.where(m*24 < time.hour(), m+1, m)
            m = Time(m*Time.DAYSEC)
            
            if results:
                mr, ms = m.val
                rising  = Date(date.date,Time(mr),calendar=date.calendar,epoch=date.epoch)
                setting = Date(date.date,Time(ms),calendar=date.calendar,epoch=date.epoch)

                event = [rising,setting]
                names = ['rising ','setting']
                for i in range(2):
                    # if event[i].time.val > Time.DAYSEC:
                    #     days = event[i].time.val // Time.DAYSEC
                    #     event[i] = event[i] + days
                    # elif event[i].time.val < time.val:
                    #     event[i] = event[i] + 1
                    print(names[i] + ':\t' + event[i].print_date())
            m.val = np.append(mt*Time.DAYSEC,m.val)
        return m
