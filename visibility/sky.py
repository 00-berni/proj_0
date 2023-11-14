import numpy as np
import matplotlib.pyplot as plt
from .stuff import get_data, import_data, interpole_three
from .Angles import Angles, HAngles, FLAT, ROUND
from .coor import Equatorial, Ecliptical, GeoPos, eclipt_to_equat, equat_to_eclipt, equat_to_altaz
from .Time.Tclasses import * 
from .Time.dates import local_ST, Green_ST

def ang_sep(coor1: Equatorial, coor2: Equatorial) -> Angles:
    ra1, dec1 = coor1.coor()
    ra2, dec2 = coor2.coor()
    ra1, dec1 = ra1.rad, dec1.rad
    ra2, dec2 = ra2.rad, dec2.rad

    z = np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra2-ra1)
    d = np.arctan(np.sqrt(1-z**2)/z)
    d = Angles(d,'rad')
    # if type(z) == np.ndarray:
    #     d = d + np.where(z<0,90,0).astype(float)
    # elif z < 0: 
    #     d = d + 90
    return d 


## Corrections 
def nutation_corr(date: Date) -> tuple[HAngles, HAngles]:
    """Computing nutation correction

    The accuracy of this formula is 
    0.5'' for celestial longitude and
    0.1'' for obliquity.

    :param date: selected date against which to compute correction
    :type date: Date
    
    :return: correction for celestial longitude and for obliquity
    :rtype: tuple[HAngles, HAngles]
    """
    # computing the Julian Ephemeris Day 
    date = date.change_time_type('TD')
    T = (date.jd - Date.J2000)/36525
    # mean longitude of the Sun
    L  = Angles.deg_to_rad(280.4665 + 36000.7698*T)
    # mean longitude of the Moon
    l  = Angles.deg_to_rad(218.3165 + 481267.8813*T)
    # mean longitude of the ascending node of the Moon's mean orbit
    Om = Angles.deg_to_rad(125.04452 - 1934.136261*T)
    # correction in longitude
    Dlon = -17.20*np.sin(Om) - 1.32*np.sin(2*L) - 0.23*np.sin(2*l) + 0.21*np.sin(2*Om) 
    # correction in obliquity
    Deps =   9.20*np.cos(Om) + 0.57*np.cos(2*L) + 0.10*np.cos(2*l) - 0.09*np.cos(2*Om) 
    return Angles(Dlon/3600,'deg'), Angles(Deps/3600,'deg')

def mean_obliquity(date: Date) -> Angles:
    """Computing mean obliquity

    The error remains less than
    1'' in a range of 2000 years
    from the epoch J2000.0 

    :param date: selected date against which to compute obliquity
    :type date: Date
    :param epoch: standard epoch, defaults to `'J2000.0'`
    :type epoch: str, optional

    :return: mean obliquity
    :rtype: HAngles
    """
    # computing the Julian Ephemeris Day 
    date = date.change_time_type('TD')
    T = (date.jd - Date.J2000)/36525
    # computing the mean obliquity
    eps = 23.4392911*3600 - 46.815*T - 59e-5*T**2 + 1.813e-3*T**3
    return Angles(eps/3600,'deg')


def precession_corr(date: Date, ra0: HAngles, dec0: HAngles, epoch: str = 'J2000.0', prmt: list[float] | None = None, sel: str = 'all') -> HAngles | Angles | tuple[HAngles, Angles]:
    """Computing correction to right ascension and declination due to precession

    One can choose to compute correction only for one value or both throught 
    `sel` parameter:

      * `sel = 'ra'`: only corrected right ascension is computed and returned
      * `sel = 'dec'`: only corrected declination is computed and returned
      * `sel = 'all'`: both
    
    :param date: date for which to compute the correction
    :type date: Date
    :param ra0: right ascension at the standard epoch
    :type ra0: HAngles
    :param dec0: declination at the standard epoch
    :type dec0: HAngles
    :param prmt: information about proper motion, defaults to None
    :type prmt: list[float] | None, optional
    :param sel: selection parameter, defaults to 'all'
    :type sel: str, optional
    
    :return: corrected right ascension and/or declination
    :rtype: HAngles | tuple[HAngles, HAngles]
    """
    tmp_date = date.copy()

    J2000 = Date.J2000

    if tmp_date.time.tytime == 'UT':
        tmp_date = tmp_date.change_time_type('TD')

    jd = tmp_date.jd

    del tmp_date

    if epoch[0] == 'J': 
        if epoch[:5] == 'J2000': jd0 = J2000  
        else: jd0 = J2000 - (2000 - float(epoch[1:]))*365.25  

    t = (jd - jd0) / 36525
    if prmt is not None:
        mua, mud = prmt
        ra0  = ra0  + HAngles(mua/3600,'deg') * (t*100)
        dec0 = dec0 + HAngles(mud/3600,'deg',lim=90) * (t*100)

    # all quantities are in arcsec
    xi = 2306.2181 * t + 0.30188 * t**2 + 0.017998 * t**3    
    z  = 2306.2181 * t + 1.09468 * t**2 + 0.018203 * t**3    
    th = 2004.3109 * t - 0.42665 * t**2 - 0.041833 * t**3
    if epoch[:5] != 'J2000':
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
    elif sel == 'dec': return Angles(dec,'rad',lim=90)
    elif sel == 'all': return HAngles(ra,'rad'), Angles(dec,'rad',lim=90)


def std_atm(height: float) -> tuple[float,float]:
    """Extrapolated temperature and pression

    Function takes the value of the 
    observatory height and iterpolates to
    get the corresponding values of 
    atmospheric temperature and pressure    

    :param height: heght of the observatory
    :type height: float
    :return: atmospheric temperature and pressure 
    :rtype: tuple[float,float]
    """
    # name of data table file
    filename = 'atm_data.csv'
    # extracting data
    h,T,p = get_data(filename)
    # condition for tabulated height
    if height in h:
        idx = np.where(height == h)[0]
        if len(idx) == 1: idx = idx[0]
        temp = T[idx]
        pres = p[idx]
    # condition for interpolation
    else:
        temp = interpole_three(T,height,h)
        pres = interpole_three(p,height,h)
    return temp, pres

def refraction_corr(alt: Angles, height: float, alt0: bool = False) -> Angles:
    """Computing correction due to atmospheric refraction

    The used formula is an empirical one.

    Function allows to compute correction to get apparent altitude
    from true one and vice versa, through the parameter `alt0`:
    if it is `True` the parameter `alt` is considered the apparent
    altitude and correction is negative.

    In both cases the wanted corrected altitude is `alt + rcorr`.

    :param alt: altitude for which to compute correction
    :type alt: Angles
    :param height: location height
    :type height: float
    :param alt0: if `True`, `alt` is apparent, defaults to `False`
    :type alt0: bool, optional
 
    :return: refraction correction
    :rtype: Angles
    """
    # setting the parameters of the formula
    # deepending on `alt0` value
    if alt0:
        a,b,c =  1, 7.31, 4.4
        cor90 = 1.3515e-3
    else:
        a,b,c =  1.02, 10.3, 5.11
        cor90 = 1.9279e-3
    # computing correction
    rcorr = a/np.tan( Angles.deg_to_rad(alt.deg + b/(alt.deg + c)) )        
    # taking into account height a.s.l. 
    if height != 0.:
        T, P = std_atm(height)
        T0, P0 = std_atm(0.)
        rcorr *= P/P0 * T0/T 
    # correcting the formula for alt = 90 deg
    if type(alt.deg) != np.ndarray and alt.deg == 90:
        rcorr += cor90
    elif type(alt.deg) == np.ndarray:
        rcorr = np.where(alt.deg == 90, rcorr + cor90, rcorr)
    # condition for true altitude computation
    if alt0: rcorr = -rcorr
    return Angles(rcorr/60,'deg',lim=90)

## Sun
class Sun():
    """_summary_

    :raises Exception: _description_
    :raises Exception: _description_
    :raises Exception: _description_
    :raises Exception: _description_
    :raises Exception: _description_
    :return: _description_
    :rtype: _type_
    """

    AU = 1.49597870707e8 # km

    def __init__(self):
        """_summary_
        """
        # mean longitude referred to the mean equinox 
        self.L = Angles(280.46646,'deg')
        # mean anomaly
        self.M = Angles(357.52911,'deg')
        # eccentricity
        self.ecc = 1.6708634e-2
        # epoch
        self.epoch = 'J2000.0'
    
    def mean_coor(self, date: Date, sel: str = 'LM') -> tuple[Angles, Angles] | Angles:
        """_summary_

        :param date: _description_
        :type date: Date
        :param sel: _description_, defaults to 'LM'
        :type sel: str, optional
        :raises Exception: _description_
        :return: _description_
        :rtype: tuple[Angles, Angles] | Angles
        """
        # list to store results
        results = []
        # mean longitude for epoch J2000
        L = self.L.deg
        # mean anomaly for epoch J2000
        M = self.M.deg
        # 
        date = date.copy()
        date = date.change_time_type('TD')
        T = (date.jd - Date.J2000)/36525
        # computing mean longitude
        if 'L' in sel:
            L += 36000.76983*T + 3.032e-4*(T**2)
            # for array type
            if type(L) == np.ndarray:
                L = np.where(abs(L) > 360, L % 360, L)
                # L = np.where(abs(L) > 360, L-(L//360)*360, L)
                L = np.where(L < 0, L+360, L)
            else:
                if abs(L) > 360: L -= (L//360)*360
                if L < 0: L += 360
            # print('L',L)
            # storing the result
            results += [Angles(L,'deg')]
        # computing mean anomaly
        if 'M' in sel:
            M += 35999.05029*T - 1.537e-4*(T**2) 
            # for array type
            if type(M) == np.ndarray:
                M = np.where(abs(M) > 360, M % 360, M)
                # M = np.where(abs(M) > 360, M-(M//360)*360, M)
                M = np.where(M < 0, M+360, M)            
            else:
                if abs(M) > 360: M -= (M//360)*360
                if M < 0: M += 360
            # print('M',M)
            # storing the result
            results += [Angles(M,'deg')]
        if len(results) == 1: results = results[0]
        # only 
        elif len(results) == 0: raise Exception(f"Error in selection! '{sel}' is not allowed")
        return results
    
    def ecc_in_date(self, date: Date) -> float | np.ndarray:
        """_summary_

        :param date: _description_
        :type date: Date
        :return: _description_
        :rtype: float | np.ndarray
        """
        date = date.copy()
        date = date.change_time_type('TD')
        T = (date.jd - Date.J2000)/36525
        return self.ecc - (4.2037e-5*T + 1.267e-7*(T**2))

    def orbit_in_date(self, date: Date, sel: str = 'all', true_val: bool = True) -> list[Angles | float] | Angles | float:
        date = date.copy()
        results = []
        L, M = self.mean_coor(date)
        if true_val:
            date = date.change_time_type('TD')
            T = (date.jd - Date.J2000)/36525
            C =  (1.914602 - 4.817e-3*T - 1.4e-5*(T**2))*np.sin(M.rad)  
            C += (1.9993e-2 - 1.01e-4*T) * np.sin(2*M.rad)
            C += 2.89e-4 * np.sin(3*M.rad)
            C = Angles(C,'deg')
            # print('C',C.deg)
            # true long referred to the mean equinox
            L = L + C
            # true anomaly
            M = M + C
        if ('L' in sel) or (sel == 'all'):
            if type(L.deg) == np.ndarray:
                L = L - np.where(abs(L) > 360, (L.deg//360)*360, 0).astype(float)
                L = L + np.where(L < 0, 360, 0).astype(float)
            else:
                if abs(L) > 360: L = L % 360 #L -= (L.deg//360)*360
                if L < 0: L += 360
            # print('l',L.deg,L.print_angle('deg',True))
            results += [L]
        if ('M' in sel) or (sel == 'all'):
            if type(M.deg) == np.ndarray:
                M = M - np.where(abs(M) > 360, (M//360)*360, 0).astype(float)
                M = M + np.where(M < 0, 360, 0).astype(float)
            else:
                if abs(M) > 360: M = M % 360 #M -= (M//360)*360
                if M < 0: M += 360
            # print('nu',M.deg)
            results += [M]
        if ('e' in sel) or (sel == 'all'):
            ecc = self.ecc_in_date(date)
            # print('ecc',ecc)
            results += [ecc]
        if len(results) == 1: results = results[0]
        elif len(results) == 0: raise Exception(f"Error in selection! '{sel}' is not allowed")        
        return results

    def app_lon(self,date: Date) -> Angles:
        date = date.copy()
        date = date.change_time_type('TD')
        # apparent longitude referred to the true equinox
        L = self.orbit_in_date(date,'L') 
        T = (date.jd - Date.J2000)/36525
        Om = 125.04 - 1934.136*T
        # print('Om',Om)
        l = L.deg - 5.69e-3 - 4.78e-3*np.sin(Angles(Om,'deg').rad)
        return Angles(l,'deg')

    def equat_coor(self, date: Date, sel:str = 'all', app_val: bool = True) -> tuple[HAngles, Angles] | HAngles | Angles:
        e = mean_obliquity(date)
        # print('e0',e.deg)
        if app_val: 
            L = self.app_lon(date)
            _, De = nutation_corr(date)
            e += De
            # print('lon',L.deg,L.print_angle('deg',True))
            # print('e',e.deg)
        else: 
            L = self.orbit_in_date(date,'L')
        
        coor = Ecliptical(L, 0.)
        coor = eclipt_to_equat(coor,e)
        if type(L.deg) == np.ndarray: 
           coor.alpha = coor.alpha + np.where(coor.alpha < 0, 360, 0).astype(float)
        elif coor.alpha < 0: 
            coor.alpha = coor.alpha + 360
        if   sel == 'ra' : return coor.alpha
        elif sel == 'dec': return coor.delta
        elif sel == 'all': return coor.alpha, coor.delta
        else: raise Exception(f"Error in selection! '{sel}' is not allowed")        

    def distance(self, date: Date) -> float | np.ndarray:
        nu, ecc = self.orbit_in_date(date,'Me')
        r = 1.000001018 * ((1-ecc**2)/(1+ecc*np.cos(nu.rad)))
        return r
    
    def ang_diameter(self, date: Date) -> Angles:
        nu, ecc = self.orbit_in_date(date,'Me')
        D = 0.533130 * ((1+ecc*np.cos(nu.rad))/(1-ecc**2))
        return Angles(D,'deg')

    def aberration_corr(self, date: Date, obj_coor: Equatorial, sel: str = 'all') -> HAngles | Angles | tuple[HAngles, Angles]:
        date = date.copy().change_time_type('TD')
        k = Angles(['+',[0,0,20.49552]],'deg')
        e = mean_obliquity(date)
        # _, De = nutation_corr(date)
        # e = e + De
        ecl_coor = equat_to_eclipt(obj_coor,e)
        lon, lat = ecl_coor.lon, ecl_coor.lat
        T = (date.jd - Date.J2000)/36525
        # print('T',T)
        peri = 102.93735 + 1.71946*T + 4.6e-4*(T**2)
        # print('pi',peri)
        peri = Angles(peri,'deg')
        L, ecc = self.orbit_in_date(date,'Le')
        # print('L',L.deg)
        # print('ecc',ecc)    
        Dlon = -k*(np.cos(L.rad-lon.rad)-ecc*np.cos(peri.rad-lon.rad))/np.cos(lat.rad)
        Dlat = -k*np.sin(lat.rad)*(np.sin(L.rad-lon.rad)-ecc*np.sin(peri.rad-lon.rad))
        Dlat.lim = 90
        lon = lon + Dlon
        lat = lat + Dlat
        corr_coor = eclipt_to_equat(Ecliptical(lon,lat),e)
        if sel != 'dec':
            Dra = corr_coor.alpha - obj_coor.alpha
            if sel == 'ra': return Dra
        if sel != 'ra':
            Ddec = corr_coor.delta - obj_coor.delta
            Ddec.lim = 90
            if sel == 'dec': return Ddec
        if sel != 'all': raise Exception(f'Error in selection!\n`{sel}` is not allowed')
        else: return Dra, Ddec

    def corr_coor(self, date: Date, sel: str = 'all', app_val: bool = True) -> HAngles | Angles | tuple[HAngles, Angles]:
        alpha, delta = self.equat_coor(date,app_val=app_val)
        coor = Equatorial(alpha, delta)
        aberr = self.aberration_corr(date,coor,sel)
        e = mean_obliquity(date)
        Dlon, De = nutation_corr(date)

        # condition for right ascension only
        if sel == 'ra':
            Dra = Dlon*(np.cos(e.rad)+np.sin(e.rad)*np.tan(delta.rad)) - De*(np.cos(alpha.rad)*np.tan(delta.rad))
            ra = alpha + Dra + aberr
            if type(ra.deg) == np.ndarray:
                ra = ra - np.where(ra < 0, (ra.deg//360)*360, 0).astype(float)
            elif ra < 0:
                ra = ra % 360
                # ra = ra - ((ra.deg//360)*360)
            return ra
        # condition for declination only
        elif sel == 'dec':
            Ddec = Dlon*(np.sin(e.rad)*np.cos(alpha.rad))+ De*np.sin(alpha.rad)
            Ddec.lim = 90
            dec = delta + Ddec + aberr
            return dec
        # condition for both results
        elif sel == 'all':
            Dra = Dlon*(np.cos(e.rad)+np.sin(e.rad)*np.tan(delta.rad)) - De*(np.cos(alpha.rad)*np.tan(delta.rad))
            Ddec = Dlon*(np.sin(e.rad)*np.cos(alpha.rad))+ De*np.sin(alpha.rad)
            Ddec.lim = 90
            Dra2, Ddec2 = aberr 
            ra = alpha + Dra + Dra2
            dec = delta + Ddec + Ddec2
            if type(ra.deg) == np.ndarray:
                ra = ra - np.where(ra < 0, (ra.deg//360)*360, 0).astype(float)
            elif ra < 0:
                ra = ra % 360
                # ra = ra - ((ra.deg//360)*360)
            return ra, dec
        else: raise Exception(f'Error in selection!\n`{sel}` is not allowed')

    def rise_set(self, date: Date, obs_pos: GeoPos, results: bool = True, iter: int = 3) -> Time | bool:
        date = date.copy()
        lon = obs_pos.lon
        lat = obs_pos.lat
        height = obs_pos.h
        # computing the GST for the meridian in date at 0h UT
        tmpdate = Date(date.date,0.)
        Dt = time_correction(tmpdate.date[0])
        GST0 = Green_ST(tmpdate,True)
        # computing r.a. and dec. on nearest days
        # (previous and consecutive one) at 0h TD     
        tmpdate = Date(date.date,0.,timetype='TD')
        jd = tmpdate.jd
        a1, d1 = self.corr_coor(tmpdate-1)
        a2, d2 = self.corr_coor(tmpdate)
        a3, d3 = self.corr_coor(tmpdate+1)
        # true altitude for a zero apparent altitude  
        h0 = refraction_corr(Angles(0.,'deg',lim=90),height, alt0=True) - self.ang_diameter(tmpdate).deg/2
        del tmpdate
        # transit computation in decimals of a day
        mt = (a2 + lon - GST0).deg / 360 
        # checking the value       
        if abs(mt) > 1:
            mt -= np.sign(mt)
        elif mt < 0:
            mt += 1
        # rising and setting
        # computing approximated value for HA
        cosH0 = (np.sin(h0.rad) - np.sin(lat.rad)*np.sin(d2.rad)) / (np.cos(lat.rad)*np.cos(d2.rad))
        # checking the presence of rising and setting
        if abs(cosH0) <= 1:
            H0 = HAngles(np.arccos(cosH0),'rad',lim=180)
            # time of rising and setting in decimals of a day
            mr = mt - H0.deg/360
            ms = mt + H0.deg/360
            # generalizing the method 
            m = np.array([mr,ms])
            # checking the value       
            m = np.where(abs(m) > 1, m-np.sign(m), m)
            m = np.where(m < 0, m+1, m)
            # routine
            for k in range(iter):
                # computing the corresponding GST
                GST = GST0.deg + 360.985647*m
                # checking the value       
                GST = np.where(GST > 360, GST - 360*(GST/360).astype(int), GST)               
                GST = HAngles(GST,'deg')
                # converting from TD to UT
                n = m + Dt/Time.DAYSEC
                # computing the corresponging r.a.
                a = np.array([interpole_three([a1.deg,a2.deg,a3.deg],ni+jd,[jd-1,jd,jd+1]) for ni in n])
                a = HAngles(a,'deg')
                # computing the corresponging dec.
                d = np.array([interpole_three([d1.deg,d2.deg,d3.deg],ni+jd,[jd-1,jd,jd+1]) for ni in n])
                d = HAngles(d,'deg',lim=90) 
                # computing the HA
                H = GST - lon - a
                # computing corresponding altitude
                h = np.arcsin(np.sin(lat.rad)*np.sin(d.rad) + np.cos(lat.rad)*np.cos(d.rad)*np.cos(H.rad))
                h = Angles(h,'rad',lim=90)
                # correcting the previous computation
                Dm = (h-h0).deg / (360 * (np.cos(d.rad)*np.cos(lat.rad)*np.sin(H.rad)))
                m += Dm
                # checking the value       
                m = np.where(abs(m) > 1, m-np.sign(m), m)
                m = np.where(m < 0, m+1, m)
            
            dm = Date(date.date,m*Time.DAYSEC).jd
            if dm[0] < date.jd:
                m[0] += 1
            elif dm[0] > date.jd+1:
                m[0] -= 1
            if dm[1] < date.jd:
                m[1] += 1
            elif dm[1] > date.jd+1:
                m[1] -= 1
            
            # if mt < m[0]:
            #     m[0] -= 1
            # if mt > m[1]:
            #     m[1] += 1
            # checking the value       
            # m = np.where(m*24 < time.hour(), m+1, m)
            m = Time(m*Time.DAYSEC)    
            # condition to print the results
            if results:
                mr, ms = m.val
                rising  = Date(date.date,Time(mr),calendar=date.calendar,timezone=date.timezone,dl_save=date.dls)
                setting = Date(date.date,Time(ms),calendar=date.calendar,timezone=date.timezone,dl_save=date.dls)

                event = [rising,setting]
                names = ['rising ','setting']
                for i in range(2):
                    print(names[i] + ':\t' + event[i].print_date())
        else:
            m = False if cosH0 > 0 else True
            if results:
                res_str = 'sets' if m else 'rises'
                print('Sun never ' + res_str)
        return m
        
    def twilight(self, date: Date, obs_pos: GeoPos, tw_type: str = 'astro',results: bool = True, add_res: bool = False) -> list[Time | bool | None | float]:
        lat = obs_pos.lat
        lon = obs_pos.lon
        height = obs_pos.h
        secinday = Time.DAYSEC
        tw = {
            'astro': -18,
            'civil':  -6,
            'nauti': -12
        }
        m = self.rise_set(date,obs_pos,results=results)
        if type(m) == bool:
            if m:
                res_val = [m]
                if add_res:
                    res_val += [None,m]
            else:
                tmpdate = Date(date.date,0.,timetype='TD')
                dec = self.corr_coor(tmpdate,'dec')
                h = Angles(tw[tw_type],'deg',90)
                h0 = h + refraction_corr(h,height,True)-self.ang_diameter(tmpdate).deg/2
                del tmpdate
                cosH1 = (np.sin(h0.rad) - np.sin(lat.rad)*np.sin(dec.rad)) / (np.cos(lat.rad)*np.cos(dec.rad))
                if cosH1 > 1:
                    res_val = [m]
                    if add_res:
                        res_val += [None,m]
                else:                
                    tmpdate = Date(date.date,0.,timetype='UT')
                    GST0 = Green_ST(tmpdate,True)
                    tmpdate = Date(date.date,0.,timetype='TD')
                    a2 = self.corr_coor(tmpdate,'ra')
                    del tmpdate
                    # transit computation in decimals of a day
                    mt = (a2 + lon - GST0).deg / 360 
                    H1 = np.arccos(cosH1)
                    # print(HAngles.rad_to_hms(H1))
                    # t = HAngles.rad_to_hms(np.array([H1,-H1])/360)
                    # t *= 0.9973
                    t = mt - np.array([H1,-H1])/ROUND.rad
                    t *= 24
                    t = np.where(abs(t)>24, t % 24, t)
                    t = np.where(t<0, t + 24, t)
                    t = np.where(Date.julian_day(*date.date,t*3600) < date.jd, t+24, t)
                    t = np.where(Date.julian_day(*date.date,t*3600) > date.jd+1, t-24, t)
                    mtw = t*3600
                    # mtw = np.where(mtw < date.time.val, mtw + secinday, mtw)
                    if results:
                        print(f'Twilight -> {(t[1]-t[0]):.3f} h') 
                        mor = Date(date.date,Time(mtw[0]),calendar=date.calendar,timezone=date.timezone,dl_save=date.dls) 
                        eve = Date(date.date,Time(mtw[1]),calendar=date.calendar,timezone=date.timezone,dl_save=date.dls) 
                        print('start:\t'+mor.print_date())
                        print('stop:\t'+eve.print_date())
                    res_val = [Time(mtw)]
                    if add_res:
                        res_val += [t[1]-t[0], m]

        else:
            rise, set = m.val
            tmpdate = Date(date.date,0.,timetype='TD')
            dec = self.corr_coor(tmpdate,'dec')
            h = Angles(tw[tw_type],'deg',90)
            h0 = h + refraction_corr(h,height,True)-self.ang_diameter(tmpdate).deg/2
            del tmpdate
            H  = np.arccos(-np.tan(lat.rad)*np.tan(dec.rad))
            H1 = np.arccos((np.sin(h0.rad) - np.sin(lat.rad)*np.sin(dec.rad)) / (np.cos(lat.rad)*np.cos(dec.rad)))
            t  = HAngles.rad_to_hms(H1-H)
            t *= 0.9973
            if abs(t) > 24: t %= 24 #t -= (t//24)*24
            if t < 0 : t += 24
            tt = t*3600
            mtw = rise - tt
            # if Date.julian_day(*date.date,mtw) < date.jd: mtw += secinday
            # if Date.julian_day(*date.date,mtw) > date.jd+1: mtw -= secinday
            etw = set + tt
            # if etw < date.time.val: etw += secinday
            # if Date.julian_day(*date.date,etw) < date.jd: etw += secinday
            # if Date.julian_day(*date.date,etw) > date.jd+1: etw -= secinday
            if results:
                print(f'Twilight -> {t:.3f} h') 
                mor = Date(date.date,Time(mtw),calendar=date.calendar,timezone=date.timezone,dl_save=date.dls) 
                eve = Date(date.date,Time(etw),calendar=date.calendar,timezone=date.timezone,dl_save=date.dls)
                print('morning:\t'+mor.print_date())
                print('evening:\t'+eve.print_date())
            twm = Time(np.array([mtw,etw]))
            res_val = [twm]
            if add_res:
                res_val += [t,m]
        return res_val

## MOON
class Moon():
    @staticmethod
    def sun_mean_anom(date: Date) -> Angles:
        date = date.copy().change_time_type('TD')
        T = (date.jd - Date.J2000) / 36525
        M = Angles(357.5291092,'deg') + 35999.052909*T - 1.536e-4*(T**2) + (T**3)/2449e4
        if type(M.deg) == np.ndarray:
            M = M - np.where(M > 360, 360*(M.deg//360), 0).astype(float)
            M = M + np.where(M < 0, 360, 0).astype(float)
        else:
            if abs(M) > 360: M = M % 360
            if M < 0: M = M + 360
        return M

    def __init__(self) -> None:
        self.L = Angles(218.3164477,'deg')
        self.D = Angles(297.8501921,'deg')
        self.M = Angles(134.9633964,'deg')
        self.F = Angles( 93.272095,'deg')
        self.epoch = 'J2000'
        self.dist = 385000.56

    def mean_in_date(self, date: Date, sel: str = 'all') -> list[Angles] | Angles:
        results = []
        date = date.copy().change_time_type('TD')
        T = (date.jd - Date.J2000) / 36525
        if ('L' in sel) or (sel == 'all'):
            L = self.L + 481267.88123421*T - 1.5786e-3*(T**2) + (T**3)/538841 - (T**4)/65194e3
            if type(L.deg) == np.ndarray:
                L = L - np.where(L > 360, 360*(L.deg//360), 0).astype(float)
                L = L + np.where(L < 0, 360, 0).astype(float)
            else:
                if abs(L) > 360: L = L % 360
                if L < 0: L = L + 360
            results += [L]
        if ('D' in sel) or (sel == 'all'):
            D = self.D + 445267.1114034*T  - 1.8819e-3*(T**2) + (T**3)/545868 - (T**4)/113065e3
            if type(D.deg) == np.ndarray:
                D = D - np.where(D > 360, 360*(D.deg//360), 0).astype(float)
                D = D + np.where(D < 0, 360, 0).astype(float)
            else:
                if abs(D) > 360: D = D % 360
                if D < 0: D = D + 360
            results += [D]
        if ('M' in sel) or (sel == 'all'):
            M = self.M + 477198.8675055*T  - 8.7414e-3*(T**2) + (T**3)/69699  - (T**4)/14712e3
            if type(M.deg) == np.ndarray:
                M = M - np.where(M > 360, 360*(M.deg//360), 0).astype(float)
                M = M + np.where(M < 0, 360, 0).astype(float)
            else:
                if abs(M) > 360: M = M % 360
                if M < 0: M = M + 360
            results += [M]
        if ('F' in sel) or (sel == 'all'):
            F = self.F + 483202.0175233*T  - 3.6539e-3*(T**2) - (T**3)/3526e3 + (T**4)/86331e4
            if type(F.deg) == np.ndarray:
                F = F - np.where(F > 360, 360*(F.deg//360), 0).astype(float)
                F = F + np.where(F < 0, 360, 0).astype(float)
            else:
                if abs(F) > 360: F = F % 360
                if F < 0: F = F + 360
            results += [F]
        if len(results) == 1: results = results[0]
        elif len(results) == 0: raise Exception(f'Error in selction!\n`{sel}` is not allowed')
        return results
    
    def get_coor(self, date: Date, sel: str = 'all') -> list[Angles | float | np.ndarray] | Angles | float | np.ndarray:
        results = []
        date = date.copy().change_time_type('TD')
        # print('JD',date.jd)
        T = (date.jd - Date.J2000) / 36525
        # print('T',T)
        L,D,M,F = self.mean_in_date(date)
        Msun = Moon().sun_mean_anom(date)
        E = 1 - 2.516e-3*T - 7.4e-6*(T**2)
        # print("L'",L.deg)
        # print("D",D.deg)
        # print("M",Msun.deg)
        # print("M'",M.deg)
        # print("F",F.deg)
        # print('E',E)
        A1 = 119.75 + 131.849*T
        if type(A1) == np.ndarray:
            A1 = np.where(abs(A1) > 360, A1 % 360, A1)
            A1 = np.where(A1 < 0, A1 + 360, A1)
        else:
            if abs(A1) > 360: A1 %= 360
            if A1 < 0: A1 += 360 
        # print('A1',A1)
        A1 = Angles.deg_to_rad(A1)

        L = L.rad
        D = D.rad
        M = M.rad
        F = F.rad
        Msun = Msun.rad

        if type(T) != np.ndarray:
            T = np.array([T])
            L = np.array([L])
            D = np.array([D])
            M = np.array([M])
            F = np.array([F])
            Msun = np.array([Msun])
            E = np.array([E])
            A1 = np.array([A1])

        alon = np.array([])
        alat = np.array([])
        adist = np.array([])
        
        if ('lon' in sel) or ('dist' in sel) or (sel == 'all'):
            cd1,cmsun1,cm1,cf1,lon1,dist1 = get_data('moon1_data.csv')
        if ('lat' in sel) or (sel == 'all'):
            cd2,cmsun2,cm2,cf2,lat2 = get_data('moon2_data.csv')
        
        for (t,l,d,m,f,msun,e,a1) in zip(T,L,D,M,F,Msun,E,A1):
            if ('lon' in sel) or (sel == 'all'):        
                A2 = 53.09 + 479264.29*t
                if abs(A2) > 360: A2 %= 360
                if A2 < 0: A2 += 360 
                # print('A2',A2)
                A2 = Angles.deg_to_rad(A2)
                
                lon = lon1*np.sin(cd1*d + cmsun1*msun + cm1*m + cf1*f)
                lon *= np.where(abs(cmsun1)==1,e,1)*np.where(abs(cmsun1)==2,e**2,1)
                lon = lon.sum()
                lon += 3958*np.sin(a1) + 1962*np.sin(l-f) + 318*np.sin(A2)
                # print('Sl',lon)

                lon = Angles.rad_to_deg(l) + lon*1e-6
                alon = np.append(alon,lon)

            if ('dist' in sel) or (sel == 'all'):
                dist = dist1*np.cos(cd1*d + cmsun1*msun + cm1*m + cf1*f)
                dist *= np.where(abs(cmsun1)==1,e,1)*np.where(abs(cmsun1)==2,e**2,1)
                dist = dist.sum()
                # print('Sr',dist)

                dist = self.dist + dist*1e-3
                adist = np.append(adist,dist)


            if ('lat' in sel) or (sel == 'all'):
                A3 = 313.45 + 481266.484*t
                if abs(A3) > 360: A3 %= 360
                if A3 < 0: A3 += 360 
                # print('A3',A3)
                A3 = Angles.deg_to_rad(A3)
                
                lat = lat2*np.sin(cd2*d + cmsun2*msun + cm2*m + cf2*f)
                lat *= np.where(abs(cmsun2)==1,e,1)*np.where(abs(cmsun2)==2,e**2,1)
                lat = lat.sum()
                lat += -2235*np.sin(l) + 382*np.sin(A3) + 350*np.sin(a1)*np.cos(f)+127*np.sin(l-m)-115*np.sin(l+m)
                # print('Sb',lat)

                lat *= 1e-6
                alat = np.append(alat,lat)
        
        if ('lon' in sel) or (sel == 'all'): 
            if len(T) == 1:
                alon = alon[0]
            results += [Angles(alon,'deg')]
        if ('lat' in sel) or (sel == 'all'): 
            if len(T) == 1:
                alat = alat[0]
            results += [Angles(alat,'deg')]
        if ('dist' in sel) or (sel == 'all'): 
            if len(T) == 1:
                adist = adist[0]
            results += [adist]
        if len(results) == 1: results = results[0]
        elif len(results) == 0: raise Exception(f'Error in selection!\n`{sel}` is not allowed')

        return results 

    def equat_coor(self,date: Date, sel: str = 'all', app_val: bool = True) -> HAngles | Angles | tuple[HAngles,Angles]:
        e = mean_obliquity(date)
        lon, lat = self.get_coor(date,'lonlat')
        if app_val:
            Dlon, De = nutation_corr(date)
            lon = lon + Dlon
            e = e + De
        coor = eclipt_to_equat(Ecliptical(lon,lat),e)
        # Dra, Ddec = Angles(0,'deg'),Angles(0,'deg')
        Dra, Ddec = Sun().aberration_corr(date,coor)
        if   sel == 'ra' : return coor.alpha + Dra
        elif sel == 'dec': return coor.delta + Ddec
        elif sel == 'all': return coor.alpha + Dra, coor.delta + Ddec
        else: raise Exception(f'Error in selection!\n`{sel}` is not allowed')

    def parallax(self, date: Date) -> Angles:
        dist = self.get_coor(date,sel='dist')
        para = np.arcsin(6378.14/dist)
        return Angles(para,'rad')

    def ang_diameter(self, date: Date, alt: Angles) -> Angles:
        dist = self.get_coor(date,'dist')
        para = self.parallax(date)
        s = 716946800/dist * (1+np.sin(alt.rad)*np.sin(para.rad))
        return Angles(s/3600,'deg')*2

    def rise_set(self, date: Date, obs_pos: GeoPos, results: bool = True, iter: int = 3) -> Time | bool:
        date = date.copy()
        lon = obs_pos.lon
        lat = obs_pos.lat
        height = obs_pos.h
        # computing the GST for the meridian in date at 0h UT
        tmpdate = Date(date.date,0.)
        Dt = time_correction(tmpdate.date[0])
        GST0 = Green_ST(tmpdate,True)
        # computing r.a. and dec. on nearest days
        # (previous and consecutive one) at 0h TD     
        tmpdate = Date(date.date,0.,timetype='TD')
        jd = tmpdate.jd
        a1, d1 = self.equat_coor(tmpdate-1)
        a2, d2 = self.equat_coor(tmpdate)
        a3, d3 = self.equat_coor(tmpdate+1)
        h0 = refraction_corr(Angles(0.,'deg',lim=90),height, alt0=True) + 0.7275*self.parallax(tmpdate).deg
        del tmpdate
        # transit computation in decimals of a day
        mt = (a2 + lon - GST0).deg / 360 
        # checking the value       
        if abs(mt) > 1:
            mt -= np.sign(mt)
        elif mt < 0:
            mt += 1
        # rising and setting
        # true altitude for a zero apparent altitude  
        # computing approximated value for HA
        cosH0 = (np.sin(h0.rad) - np.sin(lat.rad)*np.sin(d2.rad)) / (np.cos(lat.rad)*np.cos(d2.rad))
        # checking the presence of rising and setting
        if abs(cosH0) <= 1:
            H0 = HAngles(np.arccos(cosH0),'rad',lim=180)
            # time of rising and setting in decimals of a day
            mr = mt - H0.deg/360
            ms = mt + H0.deg/360
            # generalizing the method 
            m = np.array([mr,ms])
            # checking the value       
            m = np.where(abs(m) > 1, m-np.sign(m), m)
            m = np.where(m < 0, m+1, m)
            # routine
            for k in range(iter):
                # computing the corresponding GST
                GST = GST0.deg + 360.985647*m
                # checking the value       
                GST = np.where(GST > 360, GST - 360*(GST/360).astype(int), GST)               
                GST = HAngles(GST,'deg')
                # converting from TD to UT
                n = m + Dt/Time.DAYSEC
                # computing the corresponging r.a.
                a = np.array([interpole_three([a1.deg,a2.deg,a3.deg],ni+jd,[jd-1,jd,jd+1]) for ni in n])
                a = HAngles(a,'deg')
                # computing the corresponging dec.
                d = np.array([interpole_three([d1.deg,d2.deg,d3.deg],ni+jd,[jd-1,jd,jd+1]) for ni in n])
                d = HAngles(d,'deg',lim=90) 
                # computing the HA
                H = GST - lon - a
                # computing corresponding altitude
                h = np.arcsin(np.sin(lat.rad)*np.sin(d.rad) + np.cos(lat.rad)*np.cos(d.rad)*np.cos(H.rad))
                h = Angles(h,'rad',lim=90)
                # correcting the previous computation
                Dm = (h-h0).deg / (360 * (np.cos(d.rad)*np.cos(lat.rad)*np.sin(H.rad)))
                m += Dm
                # checking the value       
                m = np.where(abs(m) > 1, m-np.sign(m), m)
                m = np.where(m < 0, m+1, m)
            time = date.time
            if mt < m[0]:
                m[0] -= 1
            if mt > m[1]:
                m[1] += 1    
            # checking the value       
            # m = np.where(m*24 < time.hour(), m+1, m)
            m = Time(m*Time.DAYSEC)    
            # condition to print the results
            if results:
                mr, ms = m.val
                rising  = Date(date.date,Time(mr),calendar=date.calendar,timezone=date.timezone,dl_save=date.dls)
                setting = Date(date.date,Time(ms),calendar=date.calendar,timezone=date.timezone,dl_save=date.dls)

                event = [rising,setting]
                names = ['rising ','setting']
                for i in range(2):
                    print(names[i] + ':\t' + event[i].print_date())
            return m
        else:
            m = (cosH0 < 0)
            if results:
                res_str = 'sets' if m else 'rises'
                print('Moon never ' + res_str)
            return m
    
    def ill_fract(self, date: Date, sunlon: Angles, sundist: float | np.ndarray) -> float | np.ndarray:
        lon, lat, dist = self.get_coor(date)
        Dlon, _ = nutation_corr(date)
        lon = lon + Dlon
        geo_el = np.arccos(np.cos(lat.rad)*np.cos(lon.rad-sunlon.rad))
        if type(geo_el) == np.ndarray:
            geo_el = np.where(geo_el > FLAT.rad, geo_el % FLAT.rad, geo_el)
            geo_el = np.where(geo_el < 0, geo_el + FLAT.rad, geo_el)
        else:
            if geo_el > FLAT.rad: geo_el %= FLAT.rad
            if geo_el < 0: geo_el += FLAT.rad
        sundist *= Sun().AU
        i = np.arctan(sundist*np.sin(geo_el)/(dist-sundist*np.cos(geo_el)))
        if type(i) == np.ndarray:
            i = np.where(i > FLAT.rad, i % FLAT.rad, i)
            i = np.where(i < 0, i + FLAT.rad, i)
        else:
            if i > FLAT.rad: i %= FLAT.rad
            if i < 0: i += FLAT.rad
        return (1 + np.cos(i))/2




## Target
class Target():
    """Class to store informations about a target star

    This object collects the name of the target, its equatorial
    coordinates, the standard epoch of the coordinates and 
    proper motion in right ascension and declination if one has 
    these informations

    The attributes of the class are:

    :ivar name: name of the target star
    :vartype name: str
    :ivar coor: equatorial coordinates of the target
    :vartype coor: Equatorial
    :ivar epoch: standard epoch
    :vartype epoch: str
    :ivar mua: velocity in r.a. [mas/yrs]
    :vartype mua: float
    :ivar mud: velocity in dec. [mas/yrs]
    :vartype mud: float
    
    """
    def __init__(self, name: str, coor: Equatorial | list[HAngles | float | list], epoch: str = 'J2000.0', prmt: list[float] | None = None) -> None:
        """Constructor of the class

        One can pass coordinates as `Equatorial` object or a list of `HAngles` object, 
        of list angle or of float.

        :param name: the name of the target
        :type name: str
        :param coor: equatorial coordinates
        :type coor: Equatorial | list[HAngles  |  float  |  list]
        :param epoch: standard epoch of coordinates, defaults to `'J2000.0'`
        :type epoch: str, optional
        :param prmt: informations of proper motion, defaults to `None`
        :type prmt: list[float] | None, optional
        """
        # condition for not-`Equatorial` object 
        if type(coor) == list:
            ra, dec = coor
            coor = Equatorial(ra,dec)
        
        self.name = name
        self.coor = coor.copy()
        self.epoch = epoch

        if prmt is None:
            prmt = [None, None]

        self.mua = prmt[0]*1e-3 if prmt[0] is not None else None
        self.mud = prmt[1]*1e-3 if prmt[0] is not None else None
    
    def copy(self):
        """Function to get an exact copy of the object

        :return: the exact copy of the object
        :rtype: Target
        """
        # checking if proper motion informations are present
        prmt = [self.mua, self.mud] if self.mua is not None else None
        return Target(self.name,self.coor,self.epoch,prmt)

    def obj_info(self) -> None:
        """Function to print target informations
        """
        print('Target: ' + self.name)
        print(self.coor.print_values(eph=self.epoch))
        # checking if proper motion informations are present
        if self.mua is not None:
            print('Proper motion')
            print(f'mu_alpha =\t{self.mua*1e3} mas/yr\nmu_delta =\t{self.mud*1e3} mas/yr')
        print()

    def getcoor(self, sel: str = 'all') -> HAngles | tuple[HAngles]:
        """Function to extract the stored coordinates of the star

        One can select which coordinate value to extract 
        through `sel` parameter:

          * `sel = 'ra'`: only right ascension of the star is returned 
          * `sel = 'dec'`: only declination of the star is returned 
          * `sel = 'all'`: both 

        :param sel: selection parameter, defaults to `'all'`
        :type sel: str, optional
        
        :return: right ascension or/and declination of the target star
        :rtype: HAngles | tuple[HAngles]
        
        :raises Exception: only values in docstring are allowed for `sel` parameter
        """
        Target = self.copy()
        if   sel == 'ra' : return Target.coor.alpha
        elif sel == 'dec': return Target.coor.delta
        elif sel == 'all': return Target.coor.alpha, Target.coor.delta
        else: raise Exception(f'Wrong selection!\n`{sel}` is not allowed.')

    def getprmt(self, sel: str = 'all') -> float | list[float] | None:
        """Function to extract the stored information on proper motion

        One can select which proper motion value to extract 
        through `sel` parameter:

          * `sel = 'mua'`: only proper motion in right ascension is returned 
          * `sel = 'mud'`: only proper motion in declination is returned 
          * `sel = 'all'`: both 

        :param sel: selection parameter, defaults to 'all'
        :type sel: str, optional
        
        :return: proper motion in r.a. and/or in dec. if there is
        :rtype: float | list[float] | None

        :raises Exception: only values in docstring are allowed for `sel` parameter
        """
        Target = self.copy()
        # checking the presence of proper motion informations
        if self.mua is None and self.mud is None: return None
        elif sel == 'mua': return Target.mua
        elif sel == 'mud': return Target.mud
        elif sel == 'all': return [Target.mua, Target.mud]
        else: raise Exception(f'Wrong selection!\n`{sel}` is not allowed.')

    def coor_in_date(self, date: Date, sel: str = 'all') -> HAngles | tuple[HAngles, HAngles]:
        """Function to extract equatorial coordinates of the target at a precise date

        Coordinates of `Target` class are referred to the 
        stored standard epoch. This function computes the 
        coordinates for an arbitrary epoch from `self.epoch`. 

        One can select which coordinate value to extract 
        through `sel` parameter:

          * `sel = 'ra'`: only right ascension of the star is returned 
          * `sel = 'dec'`: only declination of the star is returned 
          * `sel = 'all'`: both 

        :param date: date for which to compute coordinates
        :type date: Date
        :param sel: selection parameter, defaults to 'all'
        :type sel: str, optional

        :return: right ascension and/or declination for chosen date
        :rtype: HAngles | tuple[HAngles, HAngles]
        
        :raises Exception: only values in docstring are allowed for `sel` parameter
        """
        # getting coordinates computed for standard epoch
        ra, dec = self.getcoor()
        # getting standard epoch
        epoch = self.epoch
        # getting informations about proper motion
        prmt = self.getprmt()
        # computing precession and nutation corrections
        ra, dec  = precession_corr(date,ra,dec,epoch,prmt)
        Dlon, De = nutation_corr(date)
        e = De + mean_obliquity(date)
        
        aberr = Sun().aberration_corr(date,Equatorial(ra,dec),sel)
        # condition for right ascension only
        if sel == 'ra':
            Dra = Dlon*(np.cos(e.rad)+np.sin(e.rad)*np.tan(dec.rad)) - De*(np.cos(ra.rad)*np.tan(dec.rad))
            ra = ra + Dra + aberr
            return ra
        # condition for declination only
        elif sel == 'dec':
            Ddec = Dlon*(np.sin(e.rad)*np.cos(ra.rad))+ De*np.sin(ra.rad)
            Ddec.lim = 90
            dec = dec + Ddec + aberr
            return dec
        # condition for both results
        elif sel == 'all':
            Dra = Dlon*(np.cos(e.rad)+np.sin(e.rad)*np.tan(dec.rad)) - De*(np.cos(ra.rad)*np.tan(dec.rad))
            Ddec = Dlon*(np.sin(e.rad)*np.cos(ra.rad))+ De*np.sin(ra.rad)
            Ddec.lim = 90
            Dra1, Ddec1 = aberr
            ra = ra + Dra + Dra1
            dec = dec + Ddec + Ddec1
            return ra, dec
        else: raise Exception(f'Wrong selection!\n`{sel}` is not allowed.')
    
    def lha(self, date: Date, lon: Angles) -> HAngles:
        """Function to compute the local hour angle

        Terrestrial longitude is taken positive 
        west and negative east from Greenwich

        :param date: date for which to compute ha
        :type date: Date
        :param lon: terrestrial longitude
        :type lon: Angles
        
        :return: local hour angle
        :rtype: HAngles
        """
        ra = self.coor_in_date(date, sel='ra')
        LST = local_ST(date,lon,True,self.epoch)
        return LST - ra



def compute_alt(date: Date, obs_pos: GeoPos, obj: Target, refcor: bool = False) -> Angles:
    """Computing the altitude of a target star

    Function computes corrections in r.a. and dec. due to precession and nutation.

    Through the parameter `refcor` one can choose to return the apparent altitude 
    by correcting for refraction effect.

    :param date: date for which to compute altitude
    :type date: Date
    :param obs_pos: observatory location 
    :type obs_pos: GeoPos
    :param obj: target object
    :type obj: Target
    :param refcor: if `True` apparent altitude is returned, defaults to `False`
    :type refcor: bool, optional
    
    :return: computed altitude
    :rtype: Angles
    """
    date = date.copy()
    
    lat = obs_pos.lat
    lon = obs_pos.lon

    dec = obj.coor_in_date(date,sel='dec')     
    HA = obj.lha(date,lon).rad

    # changing name to variables
    phi = lat.rad
    delta = dec.rad
    # computing altitude
    alt = np.arcsin( np.sin(phi)*np.sin(delta) + np.cos(phi)*np.cos(delta) * np.cos(HA)  )
    alt = Angles(alt,'rad',lim=90)
    # condition to compute apparent altitude
    if refcor:
        height = obs_pos.h
        alt += refraction_corr(alt, height)
    return alt 


def trajectory(date: Date, obs_pos: GeoPos, obj: Target, numpoint: int = 1000) -> tuple[Angles, np.ndarray]:
    """Function to trace the trajectory of a star during a day

    Function returns the computed apparent altitudes and the corrispondent 
    Julian days. These variables are array and the length is set by the 
    value of `numpoint` parameter. 

    :param date: date from which to start the trajectory elaboration
    :type date: Date
    :param obs_pos: observatory location 
    :type obs_pos: GeoPos
    :param obj: target object
    :type obj: Target
    :param numpoint: length of the arrays, defaults to `1000`
    :type numpoint: int, optional

    :return: apparent altitudes and Julian days
    :rtype: tuple[Angles, np.ndarray]
    """
    # building the array of Julian days
    dayrange = date.jd + np.linspace(0,1,numpoint)
    dayrange = Date(jd=dayrange,timetype=date.time.tytime,calendar=date.calendar)
    # computing the corrisponding apparent altitudes
    alt = compute_alt(dayrange,obs_pos,obj,refcor=True)
    return alt, dayrange.jd




def tran_ris_set(date: Date, obs_pos: GeoPos, obj: Target, results: bool = False, iter: int = 3) -> Time | list[Time | bool]:
    """Computing the time of transit, rising and setting

    Firsly, function computes the time of the transit and checks the presence of
    setting and rising (that is the target star is circumpolar or not). 

    The routine to compute these times is a recursive method, then one can change
    the number of iterations through parameter `iter`.   

    :param date: date for which to compute the analisys
    :type date: Date
    :param obs_pos: observatory location
    :type obs_pos: GeoPos
    :param obj: target star
    :type obj: Target
    :param results: if `True` results are printed, defaults to `False`
    :type results: bool, optional
    :param iter: number of iterations, defaults to `3`
    :type iter: int, optional

    :return: time for transit and (if there are) setting and rising 
    :rtype: Time
    """
    date = date.copy()
    lon = obs_pos.lon
    lat = obs_pos.lat
    height = obs_pos.h
    # computing the GST for the meridian in date at 0h UT
    tmpdate = Date(date.date,0.)
    Dt = time_correction(tmpdate.date[0])
    GST0 = Green_ST(tmpdate,True,obj.epoch)
    del tmpdate
    # computing r.a. and dec. on nearest days
    # (previous and consecutive one) at 0h TD     
    tmpdate = Date(date.date,0.,timetype='TD')
    jd = tmpdate.jd
    a1, d1 = obj.coor_in_date(tmpdate-1)
    a2, d2 = obj.coor_in_date(tmpdate)
    a3, d3 = obj.coor_in_date(tmpdate+1)
    del tmpdate

    # transit computation in decimals of a day
    mt = (a2 + lon - GST0).deg / 360 
    # checking the value       
    if abs(mt) > 1:
        mt -= np.sign(mt)
    elif mt < 0:
        mt += 1
    # routine
    for k in range(iter):
        # computing the corresponding GST
        GST = GST0.deg + 360.985647*mt
        # checking the value       
        if GST > 360:
            GST -= 360*int(GST/360)
        GST = HAngles(GST,'deg')
        # converting from TD to UT
        n = mt + Dt/Time.DAYSEC
        # computing the corresponging r.a.
        a = interpole_three([a1.deg,a2.deg,a3.deg],n+jd,[jd-1,jd,jd+1])
        a = HAngles(a,'deg')
        # computing the HA
        H = GST - lon - a
        # correcting the previous computation
        Dmt = -H.deg/360
        mt += Dmt
        # checking the value       
        if abs(mt) > 1:
            mt -= np.sign(mt)
        elif mt < 0:
            mt += 1
    # saving the result in seconds 
    m = Time(mt*Time.DAYSEC)
    # condition to print the result
    if results:
        print()
        transit = Date(date.date,m,calendar=date.calendar,timezone=date.timezone,dl_save=date.dls)
        print('transit:\t' + transit.print_date())
    
    # rising and setting
    # true altitude for a zero apparent altitude  
    h0 = refraction_corr(Angles(0.,'deg',lim=90),height, alt0=True)
    # computing approximated value for HA
    cosH0 = (np.sin(h0.rad) - np.sin(lat.rad)*np.sin(d2.rad)) / (np.cos(lat.rad)*np.cos(d2.rad))
    # checking the presence of rising and setting
    if abs(cosH0) <= 1:
        H0 = HAngles(np.arccos(cosH0),'rad',lim=180)
        # time of rising and setting in decimals of a day
        mr = mt - H0.deg/360
        ms = mt + H0.deg/360
        # generalizing the method 
        m = np.array([mr,ms])
        # checking the value       
        m = np.where(abs(m) > 1, m-np.sign(m), m)
        m = np.where(m < 0, m+1, m)
        # routine
        for k in range(iter):
            # computing the corresponding GST
            GST = GST0.deg + 360.985647*m
            # checking the value       
            GST = np.where(GST > 360, GST - 360*(GST/360).astype(int), GST)               
            GST = HAngles(GST,'deg')
            # converting from TD to UT
            n = m + Dt/Time.DAYSEC
            # computing the corresponging r.a.
            a = np.array([interpole_three([a1.deg,a2.deg,a3.deg],ni+jd,[jd-1,jd,jd+1]) for ni in n])
            a = HAngles(a,'deg')
            # computing the corresponging dec.
            d = np.array([interpole_three([d1.deg,d2.deg,d3.deg],ni+jd,[jd-1,jd,jd+1]) for ni in n])
            d = HAngles(d,'deg',lim=90) 
            # computing the HA
            H = GST - lon - a
            # computing corresponding altitude
            h = np.arcsin(np.sin(lat.rad)*np.sin(d.rad) + np.cos(lat.rad)*np.cos(d.rad)*np.cos(H.rad))
            h = Angles(h,'rad',lim=90)
            # correcting the previous computation
            Dm = (h-h0).deg / (360 * (np.cos(d.rad)*np.cos(lat.rad)*np.sin(H.rad)))
            m += Dm
            # checking the value       
            m = np.where(abs(m) > 1, m-np.sign(m), m)
            m = np.where(m < 0, m+1, m)
        # if mt < m[0]:
        #     m[0] -= 1
        # if mt > m[1]:
        #     m[1] += 1
        # checking the value       
        # m = np.where(m*24 < time.hour(), m+1, m)
        m = Time(m*Time.DAYSEC)    
        # condition to print the results
        if results:
            mr, ms = m.val
            rising  = Date(date.date,Time(mr),calendar=date.calendar,timezone=date.timezone,dl_save=date.dls)
            setting = Date(date.date,Time(ms),calendar=date.calendar,timezone=date.timezone,dl_save=date.dls)

            event = [rising,setting]
            names = ['rising ','setting']
            for i in range(2):
                print(names[i] + ':\t' + event[i].print_date())
        # updating the variable of the results
        m.val = np.append(mt*Time.DAYSEC,m.val)
    else:
        rs = False if cosH0 > 1 else True
        m = [m,rs]
        if results:
            name = obj.name if obj.name != '' else 'Target'
            if rs:
                print(name + ' is circumpolar')
            else:
                print(name + ' never rises')
    return m


def visibility_plot(truedate: Date, date: Date, obj: Target, obs: GeoPos, SUN: Sun, MOON: Moon, m: float | np.ndarray, tw: np.ndarray | bool, k: float | None = None, dist: tuple[np.ndarray, np.ndarray] | None = None, altmu: float | None = None, not_vis: str | None = None, save_fig: bool = True, window: np.ndarray | None = None, win_par: bool = True) -> None:
    """Function to plot visibility data

    :param truedate: _description_
    :type truedate: Date
    :param date: _description_
    :type date: Date
    :param obj: _description_
    :type obj: Target
    :param obs: _description_
    :type obs: GeoPos
    :param SUN: _description_
    :type SUN: Sun
    :param MOON: _description_
    :type MOON: Moon
    :param m: _description_
    :type m: float | np.ndarray
    :param tw: _description_
    :type tw: np.ndarray | bool
    :param k: _description_, defaults to None
    :type k: float | None, optional
    :param dist: _description_, defaults to None
    :type dist: float | None, optional
    :param altmu: _description_, defaults to None
    :type altmu: float | None, optional
    :param not_vis: _description_, defaults to None
    :type not_vis: str | None, optional
    :param save_fig: _description_, defaults to True
    :type save_fig: bool, optional
    :param window: _description_, defaults to None
    :type window: np.ndarray | None, optional
    :param win_par: _description_, defaults to True
    :type win_par: bool, optional
    """
    ## PLOT
    # factors for the figure adjustment
    top = 0.92
    left = 0.15
    right = 0.95
    # condition to adjust the image
    if ((window is not None) and (not_vis is None)) and (window.shape[0] > 1):
        top = 0.9
    if win_par or not_vis is not None:
        left = 0.08
    plt.figure(figsize=[13,12]).subplots_adjust(left=left,top=top,right=right)
    ax = plt.axes()
    ax.set_facecolor('black')
    # True title
    plt.suptitle(obj.name + ' on ' + truedate.print_date(sel='date'),fontsize=20)
    # condition to show the visibility window of time if any
    if window is not None and not_vis is None:
        subtitle = 'visible '
        # counting the numer of windows
        cnt = 0
        for w in window:
            if cnt > 0:
                subtitle += '\nand '
            if w[0] > w[1]: w = w[::-1]
            subtitle += 'from ' + Date(jd=w[0],timetype=date.time.tytime,timezone=date.timezone,dl_save=date.dls,calendar=date.calendar).print_date() + ' to ' + Date(jd=w[1],timetype=date.time.tytime,timezone=date.timezone,dl_save=date.dls,calendar=date.calendar).print_date()
            cnt += 1
        plt.title(subtitle,fontsize=16)

    ## Target
    # computing trajectory
    alt, dayrange = trajectory(date,obs,obj)
    # computing time of transit and if any rising and setting
    event = Date(jd=m,timetype=date.time.tytime,timezone=date.timezone,dl_save=date.dls,calendar=date.calendar)
    # computing the corresponding apparent altitude
    ealt = compute_alt(event,obs,obj,True)    
    # plotting all
    plt.plot(dayrange,alt.deg,'b',label=obj.name)
    if not win_par:
        if type(m) != np.ndarray:
            plt.plot(event.jd,ealt.deg,'vg',label='transit')
        else:
            plt.plot(event.jd[0],ealt.deg[0],'vg',label='transit')
            plt.plot(event.jd[1],ealt.deg[1],'vy',label='rising')
            plt.plot(event.jd[2],ealt.deg[2],'vr',label='setting')
    
    # computing the date for each Julian Day used to sample trajectory
    dates = Date(jd=dayrange,timetype=date.time.tytime,timezone=date.timezone,dl_save=date.dls,calendar=date.calendar)    
    # getting observatory position
    lon = obs.lon
    lat = obs.lat

    ## SUN
    # computing ra and dec of the Sun
    ra, dec = SUN.corr_coor(dates)
    # computing the local HA
    lha = local_ST(dates,lon,True) - ra
    # converting to AltAz coordinates
    altaz = equat_to_altaz(Equatorial(ra,dec),lha,lat)
    # computing the angular dimension of the Sun
    sunrad = SUN.ang_diameter(dates)/2
    # adding refraction correction
    sunalt = altaz.alt + refraction_corr(altaz.alt,obs.h)
    # computing the twilight time if any
    if type(tw) != bool: 
        suntw = Date(jd=tw,timetype=date.time.tytime,timezone=date.timezone,dl_save=date.dls,calendar=date.calendar)
    # plotting all
    plt.plot(dayrange,sunalt.deg,'y',label='Sun')
    plt.plot(dayrange,sunalt.deg - sunrad.deg,'y',alpha=0.5)
    plt.plot(dayrange,sunalt.deg + sunrad.deg,'y',alpha=0.5)
    if type(tw) != bool: 
        plt.axvline(suntw.jd[0],ymin=0,ymax=1,linestyle=(0,(5,10)),alpha=0.4,color='gold',label='twilight')
        plt.axvline(suntw.jd[1],ymin=0,ymax=1,linestyle=(0,(5,10)),alpha=0.4,color='gold')

    ## MOON
    # computing ra and dec of the Sun
    ra, dec = MOON.equat_coor(dates)
    # computing the local HA
    lha = local_ST(dates,lon,True) - ra
    # converting to AltAz coordinates
    altaz = equat_to_altaz(Equatorial(ra,dec),lha,lat)
    # computing the angular dimension of the Moon
    moonalt = altaz.alt + refraction_corr(altaz.alt,obs.h)
    # plotting it
    plt.plot(dayrange,moonalt.deg,'w',label='Moon')

    ## plot stuff
    hour = date.time.hour()
    # factor to have a nice display
    N = 25 if hour <= (int(hour)+0.5) else 26
    # ticks of the plot in jds
    ticks = Date(date.date, (np.arange(0,N,1)+int(date.time.hour()))*3600,timetype=date.time.tytime,timezone=date.timezone,dl_save=date.dls,calendar=date.calendar)
    # label of ticks is in hours
    labelticks = np.round(ticks.time.hour(),0).astype(int)
    labelticks = np.where(labelticks >= 24, labelticks-24, labelticks)
    # condition for 00:00:00 UT
    midnight = ticks.jd[np.where(labelticks==0)[0]]
    for mn in midnight:
        # if midnight is present and it is not the beginning
        # or the ending of the x-axis then it is marked to
        # visualize the change of the day 
        if mn != ticks.jd[-1] and mn != ticks.jd[0]:
            plt.axvline(mn,ymin=0,ymax=1,alpha=0.4,color='white')
    # horizon
    plt.axhline(0,xmin=0,xmax=1,alpha=0.5,color='white',label='horizon')
    # condition to plot the airmass limit
    if altmu is not None:
        plt.axhline(altmu,xmin=0,xmax=1,linestyle='dashed',alpha=0.5,color='pink')
        # stuff for the text to display the value of the airmass
        textpos = 2 
        if (altmu+2 - alt.deg[0]) < 0: 
            textpos = -4
        if window is not None and win_par:
            startpos = window.min()
        else:
            startpos = dayrange[0]
        plt.annotate('$X$ = 3',(startpos,altmu),(startpos-0.03,altmu+textpos),fontsize=12,color='white')

    plt.xticks(ticks.jd,labelticks,fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel('UT [h]',fontsize=15)
    plt.ylabel('alt$_0$ [deg]',fontsize=15)
    plt.grid(linestyle='dotted',color='gray',alpha=0.5)
    # information about the place of the observatory
    location_string = obs.place_info(True)
    plt.text(0.01, 0.02, location_string, fontsize=14, transform=plt.gcf().transFigure)
    # condition to show the illuminated fraction of the Moon's disk
    if k is not None:
        info_str = f'Mean Moon ill. frac.: {k*100:.2f} %'
        # condition to show the distance from Moon
        if dist is not None:
            # extracting infotmation
            dist, dtime = dist
            # index of minimum distance
            idx = dist.argmin()
            # string for minimum distance
            info_str += f'\nMin dist. from Moon: {dist[idx]:.2f} deg'
            # condition to plot distances
            if win_par:
                # computing alt for the time of minimum distance
                dalt = compute_alt(Date(jd=dtime[idx],timezone=date.timezone,dl_save=date.dls,calendar=date.calendar),obs,obj,True).deg
                if dalt > 0:
                    plt.plot(dtime[idx],dalt,'X',color='violet',label='min dist',zorder=10)
                # rescaling arrays
                n = len(dist) // 6 
                dist  = dist[::n]
                dtime = dtime[::n]
                dalt = compute_alt(Date(jd=dtime,timezone=date.timezone,dl_save=date.dls,calendar=date.calendar),obs,obj,True).deg
                plt.plot(dtime,dalt,'|',color='white')
                # showing the distance value for each point
                for i in range(len(dtime)):
                    if dalt[i] > 0:
                        alt_pos = -3
                        if dalt[i]+alt_pos <= 0:
                            alt_pos = 2
                        plt.annotate(f'{dist[i]:.0f}$^\circ$',(dtime[i],dalt[i]),(dtime[i]-0.005,dalt[i]+alt_pos),color='white',fontsize=12)
        # displaying minimum distance
        plt.text(0.7, 0.02, info_str, fontsize=14, transform=plt.gcf().transFigure)
    # case for a not-visible target
    if not_vis is not None:
        props = dict(boxstyle='round', facecolor='red', alpha=0.8)
        plt.text(0.38, 0.5, not_vis, fontsize=23, transform=plt.gcf().transFigure, bbox=props)
    pos_leg = (-0.04, 0.9)
    # condition to plot only the visibility window
    if window is not None and win_par:
        xlims = (window.min()-1/24,window.max()+1/24)
        max_alt = ealt.deg[0] if type(ealt.deg) == np.ndarray else ealt.deg
        max_sunalt  = max(sunalt.deg)
        max_moonalt = max(moonalt.deg)
        plt.xlim(*xlims)
        plt.ylim(0,max(max_alt,max_sunalt,max_moonalt)+20)
        pos_leg = (0.85,0.8)
    if not_vis is None:
        plt.legend(numpoints=1,facecolor='black',labelcolor='w',bbox_to_anchor=pos_leg,fontsize=13).get_frame().set_alpha(None)
    # condtion to save automatically the figure
    if save_fig:
        from .stuff import RESULTS_FOLDER, ph
        namefig = obj.name.replace(' ','-') + f'_{truedate.date[0]}-{truedate.date[1]}-{truedate.date[2]:.0f}_' + obs.name.replace(' ','-')
        if not_vis is not None:
            namefig += '_not-vis'
        namefig += '.png'
        plt.savefig(ph.join(RESULTS_FOLDER,namefig), format='png')
    plt.show()



def visibility(date: Date, obj: Target, obs: GeoPos, airmass: float = 3., numpoint: int = 300, display_plots: bool = True, save_fig: bool = True, win: bool = True) -> None | tuple[list[float], float]:
    """_summary_

    :param date: _description_
    :type date: Date
    :param obj: _description_
    :type obj: Target
    :param obs: _description_
    :type obs: GeoPos
    :param airmass: _description_, defaults to 3.
    :type airmass: float, optional
    :param numpoint: _description_, defaults to 300
    :type numpoint: int, optional
    :param display_plots: _description_, defaults to True
    :type display_plots: bool, optional
    :param save_fig: _description_, defaults to True
    :type save_fig: bool, optional
    :param win: _description_, defaults to True
    :type win: bool, optional
    :return: _description_
    :rtype: None | tuple[list[float], float]
    """
    SEP = '-'*10 + '\n'

    date = date.copy()
    # initializing the visibility window
    edges = np.array([date.jd,date.jd+1])
    # setting the starting point
    start_point = edges[0]
    
    # getting the name
    name = obj.name if obj.name != '' else 'Target'
    # printing the start of the routine
    intro_str = '\n' + SEP + 'Visibility of target on ' + date.print_date()
    if date.time.tytime == 'local':
        intro_str += ' ==> ' + date.ut_to_local('str')
    intro_str += '\nat ' + obs.place_info(True) +'\n'
    print(intro_str)
    # printing target info
    obj.obj_info()
    
    # computing trans., ris. and set. of the target
    m = tran_ris_set(date,obs,obj,False)
    # initializing Sun
    SUN = Sun()
    # computing twilights
    tw, t_tw, mtw = SUN.twilight(date,obs,results=False,add_res=True)

    ## Visibility Routine
    # checking the Sun
    if type(tw) == bool:
        # Sun never sets
        if tw:
            not_vis = '\nNight never occurs'
            if display_plots:
                MOON = Moon()
                start_point = Date(jd=start_point)
                if type(m) == list:
                    m = Date(date.date,time=m[0]).jd
                else:
                    m = Date(date.date,time=m).jd 
                visibility_plot(date,start_point,obj,obs,SUN,MOON,m,tw,not_vis=not_vis,save_fig=save_fig,win_par=False)
            print(not_vis+'\nTarget is not visible!\n'+SEP)
            return None
        # Sun never rises
        else:
            sun_str = '\nDaylight never occurs'
            window = np.copy([edges])
    else:
        # getting twilight, sunset and sunrise times in jds 
        tw  = Date(date.date,time=tw).jd
        mtw = Date(date.date,time=mtw).jd
        # shifting the time to get a window during 
        # a single night
        if tw[0] < start_point:
            shift = 1
            if abs(start_point-tw[0]) < abs(start_point-tw[1]):
                shift = -1
            newtw, t_tw, newmtw = SUN.twilight(date+shift,obs,results=False,add_res=True)
            pos = (1-shift)//2
            tw[pos]  = Date((date+shift).date,newtw.val[pos]).jd    
            mtw[pos] = Date((date+shift).date,newmtw.val[pos]).jd
            del pos, shift
        else:
            if tw[1] > tw[0]:
                newtw, t_tw, newmtw = SUN.twilight(date-1,obs,results=False,add_res=True)
                tw[1]  = Date((date-1).date,newtw.val[1]).jd    
                mtw[1] = Date((date-1).date,newmtw.val[1]).jd
            if tw[0] > tw[1]+1:
                newtw, t_tw, newmtw = SUN.twilight(date-1,obs,results=False,add_res=True)
                tw[0]  = Date((date-1).date,newtw.val[0]).jd    
                mtw[0] = Date((date-1).date,newmtw.val[0]).jd
        start_point = mtw[1] 
        window = np.copy([tw])
        # storing information about Sun
        sun_str = f'\n* * SUN * *\nrising:  \t{Date(jd=mtw[0]).print_date()}\nsetting:\t{Date(jd=mtw[1]).print_date()}\n-- Twilight --\nDuration: {t_tw:.2f} h\nmorning:\t{Date(jd=tw[0]).print_date()}\nevening:\t{Date(jd=tw[1]).print_date()}'
    # checking the target
    if type(m) == list:
        # target never rises
        if not m[1]:
            not_vis = '\nTarget is always below the horizon'
            m = Date(date.date,time=m[0])
            # printing information about target and Sun
            print('transit:\t'+m.print_date()+'\n')
            print(sun_str)
            if display_plots:
                MOON = Moon()
                start_point = Date(jd=start_point)
                m = m.jd
                visibility_plot(date,start_point,obj,obs,SUN,MOON,m,tw,not_vis=not_vis,save_fig=save_fig,window=window,win_par=win)
            print(not_vis+'\nTarget is not visible!\n'+SEP)
            return None
    else:
        # extracting transit, rising and setting times in jds
        tran, rise, set = Date(date.date,time=m.val).jd
        # shifting the time to get a window during 
        # a single night
        if tran > start_point+1:
            # tran -= 1
            tran = Date((date-1).date,tran_ris_set(date-1,obs,obj).val[0]).jd
        elif tran < start_point:
            tran = Date((date+1).date,tran_ris_set(date+1,obs,obj).val[0]).jd
            # tran += 1
        if set > start_point+1:
            set = Date((date-1).date,tran_ris_set(date-1,obs,obj).val[2]).jd
            # set -= 1
        # studying the visibility
        if rise < start_point:
            rise = Date((date+1).date,tran_ris_set(date+1,obs,obj).val[1]).jd
            # rise += 1
            if set <= tw[1]:
                # target does not rise in the window
                if rise >= tw[0]:
                    not_vis = '\nTarget is below the horizon'
                    start_point = Date(jd=start_point)
                    # printing information about target and Sun
                    print(f'transit:\t{Date(jd=tran).print_date()}')
                    print(f'rising: \t{Date(jd=rise).print_date()}')
                    print(f'setting:\t{Date(jd=set).print_date()}')
                    print(sun_str)
                    if display_plots:
                        MOON = Moon()
                        m = Date(jd=np.array([tran,rise,set])).jd
                        visibility_plot(date,start_point,obj,obs,SUN,MOON,m,tw,not_vis=not_vis,save_fig=save_fig,window=window,win_par=win)
                    print(not_vis+'\nTarget is not visible!\n'+SEP)
                    return None
                else:
                    # updating the window
                    window = np.array([[max(rise,tw[1]),tw[0]]])
            else:
                if rise >= tw[0]:
                    # updating the window
                    window = np.array([[tw[1],min(set,tw[0])]])
                else:
                    # updating the window
                    window = np.array([[tw[1],set],[rise,tw[0]]])
        else:
            if rise > start_point+1:
                rise = Date((date-1).date,tran_ris_set(date-1,obs,obj).val[1]).jd
            if rise >= tw[0]:
                # target does not rise in the window
                if set > rise or set <= tw[1]:
                    not_vis = '\nTarget is below the horizon'
                    start_point = Date(jd=start_point)
                    # printing information about target and Sun
                    print(f'transit:\t{Date(jd=tran).print_date()}')
                    print(f'rising: \t{Date(jd=rise).print_date()}')
                    print(f'setting:\t{Date(jd=set).print_date()}')
                    m = Date(jd=np.array([tran,rise,set])).jd
                    print(sun_str)
                    if display_plots:
                        MOON = Moon()
                        visibility_plot(date,start_point,obj,obs,SUN,MOON,m,tw,not_vis=not_vis,save_fig=save_fig,window=window,win_par=win)
                    print(not_vis+'\nTarget is not visible!\n'+SEP)
                    return None
                else:
                    # updating the window
                    window = np.array([[tw[1],min(set,tw[0])]])
            else:
                if set <= tw[1] or set >= tw[0]:
                    # updating the window
                    window = np.array([[max(tw[1],rise),tw[0]]])
                else:
                    # updating the window
                    window = np.array([[tw[1],set],[rise,tw[0]]])                    
    # loading the coordinates of obesrvatory
    lon, lat, height = obs.coor()    
    # computing ra and dec of the target
    ra, dec = obj.coor_in_date(date)
    # initializing the Moon object
    MOON = Moon()
    # computing the altitude for the selected value of airmass
    altmu = Angles.rad_to_deg(np.arcsin(1/airmass))
    
    # initializing variables for the routine
    hvis = 0        # it stores how much time the target is visible (in hh)
    tmpcnt = 0      # it counts how many routines function does
    results = []    # it stores the results which function returns
    strresult = ''  # it stores the information about visibility
    dists = []      # it stores the array of distances from Moon
    dtime = []      # it stores the corresponding value of time in jds
    # studying the visibility in the window
    for i in range(window.shape[0]):
        # checking the order
        if window[i,0] > window[i,1]: window[i,:] = window[i,::-1]
        # renaming
        w = window[i]
        # sampling the range of time
        aw = np.linspace(w[0],w[1],numpoint)    
        dates = Date(jd=aw,calendar=date.calendar)
        # computes the illuminated fraction of the Moon's disk
        k = MOON.ill_fract(dates,SUN.app_lon(dates),SUN.distance(dates))
        # max value for illuminated fraction
        maxk  = k.max()
        dmaxk = dates.jd[k.argmax()]
        # averaging 
        k = k.sum()/len(k)
        # computing ra and dec of the Moon
        moonra, moondec = MOON.equat_coor(dates)
        lst = local_ST(dates,lon,True)
        lha = lst - moonra
        # converting them in AltAz
        moonalt = equat_to_altaz(Equatorial(moonra,moondec),lha,lat).alt
        del lst, lha
        # computing Moon angular size
        r = MOON.ang_diameter(dates,moonalt).deg/2
        mindist = r
        # computing ra and dec of the target
        ra, dec = obj.coor_in_date(dates)
        # computing the angular distance between Moon and target 
        dist = ang_sep(Equatorial(ra,dec),Equatorial(moonra,moondec)).deg
        dist = np.where(dist < 0, 180 + dist, dist)
        # dist = abs(np.where(dist > 180, 360-dist, dist))
        # checking where target is visible
        idx = np.where(dist > mindist)[0]
        # target is behind Moon then it is not visible
        if len(idx) == 0:
            # updating the counter of routines
            tmpcnt += 1
            # condition to stop the routine
            if tmpcnt == window.shape[0]:
                not_vis = '\nTarget is behind Moon'
                # printing information about target and Sun
                if type(m) == list:
                    m = Date(date.date,m[0])
                    print(f'transit:\t{m.print_date()}')
                    m = m.jd
                else:
                    print(f'transit:\t{Date(jd=tran).print_date()}')
                    print(f'rising: \t{Date(jd=rise).print_date()}')
                    print(f'setting:\t{Date(jd=set).print_date()}')
                    m = Date(jd=np.array([tran,rise,set])).jd
                print(sun_str)
                if display_plots:
                    visibility_plot(date,start_point,obj,obs,SUN,MOON,m,tw,not_vis=not_vis,save_fig=save_fig,window=window,win_par=win)
                print(not_vis+"\nIt is not visible!"+SEP)
                return None
        else:
            # storing distances information in the window of visibility
            dists += [dist[idx]]
            dtime += [dates.jd[idx]]
            # defining the start and the end of the visibility window in time
            start = dates.jd[idx[0]]
            end = dates.jd[idx[-1]]
            # storing the results
            results += [start,end]
            # updating the duration in hours
            hvis += np.diff([start,end]).sum()*Time.DAYSEC
            # print(Date(jd=start,calendar=date.calendar,timezone=date.timezone,dl_save=date.dls).print_date())
            # collecting information to print
            strresult += ('From:\t' + Date(jd=start,calendar=date.calendar,timezone=date.timezone,dl_save=date.dls).print_date()) + '\n'
            strresult += ('To:  \t' + Date(jd=end,calendar=date.calendar,timezone=date.timezone,dl_save=date.dls).print_date()) + '\n'
            # print(Date(jd=end,calendar=date.calendar,timezone=date.timezone,dl_save=date.dls).print_date())
    # getting 1-D arrays
    dists = np.array(dists).flatten()
    dtime = np.array(dtime).flatten()
    # storing the index of minimum distance
    idx = dists.argmin()
    # printing information about target, Sun and Moon
    if type(m) == list:
        m = Date(date.date,m[0])
        print(f'transit:\t{m.print_date()}')
        m = m.jd
    else:
        # print(f'transit:\t{Date(jd=tran).print_date()}')
        # print(f'rising: \t{Date(jd=rise).print_date()}')
        # print(f'setting:\t{Date(jd=set).print_date()}')
        m = Date(jd=np.array([tran,rise,set])).jd
    # print(sun_str)
    # print(f"\nMax illuminated fraction of Moon's disk: {maxk*100:.2f} % on {Date(jd=dmaxk).print_date()}")
    # print(f"Mean illuminated fraction of Moon's disk: {k*100:.2f} %")
    # print(f'\nMinimum distance from Moon: {dists[idx]:.3f} deg at {Date(jd=dtime[idx]).print_date()}')
    # printing information about the visibility
    print(f'** {maxk*100:.2f} %\t{Date(jd=dmaxk).print_date()} **')
    print('\n' + name + f" is visible for {Time.str_time(hvis,'',sep=['h ','m ','s'])[:-1]}")
    print(strresult + SEP)
    if display_plots:
        start_point = Date(jd=start_point)
        visibility_plot(date,start_point,obj,obs,SUN,MOON,m,tw,k,(dists,dtime),altmu=altmu,save_fig=save_fig,window=window,win_par=win)
    return results, hvis


def initialize_data(file_name: str, sel: int | np.ndarray | slice = slice(None)):
    obj, obs, odate = import_data(file_name, sel)
    
    # objects
    name, ra, dec, prmt, epoch = obj
    obj = [ Target(name[i],[ra[i],dec[i]],epoch[i],prmt[i]) for i in range(len(name))]
    del name, ra, dec, prmt, epoch

    # observatories
    name, lat, lon, h = obs
    obs = [GeoPos(lon[i],lat[i],h[i],name[i]) for i in range(len(name))]
    del name, lat, lon, h

    # dates
    date = odate
    odate = [Date(date=date[i],time=Time(0.,'UT')) for i in range(len(date))]
    del date

    return obj, obs, odate