from test_struc import *
from test_angle_class import *
from test_JD import julian_day, calendar_date

class Equatorial():

    def __init__(self, alpha: HAngles | float | list, delta: HAngles | float | list) -> None:
        if type(alpha) != HAngles:
            alpha = HAngles(alpha,'hms')
        if type(delta) != HAngles:
            delta = HAngles(delta,'deg',lim=90)
        self.alpha = alpha
        self.delta = delta
    
    def ha(self, ST: HAngles) -> HAngles:
        return ST - self.alpha

    def print_values(self, sel: str = 'all'):
        alpha_str = 'alpha:\t'+self.alpha.print_angle('hms')
        delta_str = 'delta:\t'+self.delta.print_angle('deg')
        if sel == 'alpha': return alpha_str
        elif sel == 'delta': return delta_str
        elif sel == 'all': return 'Equatorial Coordinates:\n' + alpha_str + delta_str
        else: raise Exception(f'!Error: `{sel}` is not allowed!\nSee the docstring of `print_values()` function')


class GeoPos():
    def __init__(self,lon: Angles | float | list, lat: Angles | float | list, h: float) -> None:
        if type(lon) != Angles:
            lon = Angles(lon,'deg')
        if type(lat) != Angles:
            lat = Angles(lat,'deg',lim=90)
        self.lon = lon
        self.lat = lat
        self.h = h


MONTHS = {
    1 : 'January',
    2 : 'February',
    3 : 'March',
    4 : 'April',
    5 : 'May',
    6 : 'June',
    7 : 'July',
    8 : 'August',
    9 : 'September',
    10 : 'October',
    11 : 'November',
    12 : 'December'
}

J2000 = 2451545.
B1950 = 2433282.4235
B1900 = 2415020.3135

# Mean ST at Greenwich
def mean_Green_HA(year: int, 
                  month: int, 
                  day: float
                  ) -> HAngles:
    JD = julian_day(year,month,day)
    T = (JD - 2451545) / 36525
    theta0 = 280.46061837 + 360.98564736629 * (JD - 2451545) + 3.87933e-4 * T**2 - T**3 / 3871e4
    theta0 = theta0 / 15
    #! Da capire e rivedere questa parte
    if np.abs(theta0) > 24:
        numbdays = np.trunc(theta0 / 24).astype(int)
        theta0 -= numbdays*24
    if theta0 < 0:
        theta0 += 24
    #!
    return HAngles(theta0,'hms')



def precession_corr(coor: Equatorial, date: list, epoch: str = 'J2000.0', prmt: list | None = None):
    ra0  = coor.alpha
    dec0 = coor.delta
    
    jd = julian_day(*date)
    if epoch[0] == 'J': 
        if epoch[:5] == 'J2000': jd0 = J2000  
        else: J2000 - (2000 - float(epoch[1:]))*365.25  
    elif epoch[0] == 'B':
        jd0 = B1900 if epoch == 'B1900' else B1950

    t = (jd - jd0) / 36525

    if prmt is not None:
        ra0  = ra0  + HAngles(prmt[0]/3600,'deg') * (t*100)
        dec0 = dec0 + HAngles(prmt[1]/3600,'deg') * (t*100)

    xi = 2306.2181 * t + 0.30188 * t**2 + 0.017998 * t**3    
    z  = 2306.2181 * t + 1.09468 * t**2 + 0.018203 * t**3    
    th = 2004.3109 * t + 0.42665 * t**2 + 0.041833 * t**3
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

    ra = np.arctan(A/B) + z.rad

    if dec0.decformat('deg') < 80:
        C = np.sin(th.rad) * np.cos(dec0.rad) * np.cos((ra0 + xi).rad) + np.cos(th.rad) * np.sin(dec0.rad)
        
        dec = np.arcsin(C)  
    else:
        dec = np.arccos(np.sqrt(A**2 + B**2))
    
    return HAngles(ra,'rad'), HAngles(dec,'rad',lim=90)
        


def trajectory(obs_pos: GeoPos, eq_pos: Equatorial, date: list, epoch: str = 'J2000.0', prmt: list | None = None) -> ArrAngle:
    lat = obs_pos.lat
    lon = obs_pos.lon

    ra, dec = precession_corr(eq_pos,date,epoch=epoch,prmt=prmt)
    
    LST = mean_Green_HA(*date) - lon
    T = np.linspace(0,24,48) + LST.decformat('hms')
    alpha = ra.decformat('hms')
    H = T - alpha
    HA = ArrAngle(H,'hms')

    phi = lat.rad
    delta = dec.rad

    alt = np.arcsin( np.sin(phi)*np.sin(delta) + np.cos(phi)*np.cos(delta) * np.cos(HA.rad)  )

    return ArrAngle(alt,'rad',lim=90)


if __name__ == '__main__':
    starting_test('TEST FOR TRAJECTORY IN THE SKY')

    try:

        lat = 43.93039980730644
        lon = 10.907521220934262 
        h = 65
        obs = GeoPos(lon,lat,h)

        ##
        print('Altair')
        ra = HAngles(['+',[19,50,46.99855]],'hms')
        dec = HAngles(['+',[8,52,05.9563]],'deg',lim=90)
        obj = Equatorial(ra,dec)
        prmt = [536.23e-3, 385.29e-3]
        print(obj.print_values())
        ##

        date = [1950,6,15]
        JDdate = julian_day(*date)
        alt = trajectory(obs,obj,date,prmt=prmt,epoch='B1900')

        hour = (date[2]-np.floor(date[2]))*24

        dayrange = np.array([calendar_date(day/24+JDdate)[-1] for day in np.linspace(0,24,48)])#+np.modf(JDdate)[0]*24

        plt.figure()
        plt.title(f'{date[0]}  {MONTHS[date[1]]}  {int(date[2])}  {hour:.1f}h at\nlat: {obs.lat.deg[0]}{obs.lat.deg[1]}, lon: {obs.lon.deg[0]}{obs.lon.deg[1]}')
        plt.grid(axis='x',linestyle='dotted',color='gray',alpha=0.7)
        plt.plot((dayrange-date[-1])*24,alt.deg)
        plt.axhline(0,xmin=0,xmax=1,color='k',linewidth=1,alpha=0.7)
        plt.xticks(np.arange(0,25,2),np.arange(0,25,2)+hour)
        plt.axvline(6.+hour,ymin=0,ymax=1,linestyle='dotted',alpha=0.5,color='y')
        plt.axvline(19.-hour,ymin=0,ymax=1,linestyle='dotted',alpha=0.5,color='y')
        plt.xlabel('t [h]')
        plt.ylabel('alt [deg]')
        
        alt0 = alt

        date = [2020,10,16]
        JDdate = julian_day(*date)
        alt = trajectory(obs,obj,date,prmt=prmt)

        hour = (date[2]-np.floor(date[2]))*24

        dayrange = np.array([calendar_date(day/24+JDdate)[-1] for day in np.linspace(0,24,48)])#+np.modf(JDdate)[0]*24

        plt.figure()
        plt.title(f'{date[0]}  {MONTHS[date[1]]}  {int(date[2])}  {hour:.1f}h at\nlat: {obs.lat.deg[0]}{obs.lat.deg[1]}, lon: {obs.lon.deg[0]}{obs.lon.deg[1]}')
        plt.grid(axis='x',linestyle='dotted',color='gray',alpha=0.7)
        plt.plot((dayrange-date[-1])*24,alt.deg)
        plt.axhline(0,xmin=0,xmax=1,color='k',linewidth=1,alpha=0.7)
        plt.xticks(np.arange(0,25,2),np.arange(0,25,2)+hour)
        plt.xlabel('t [h]')
        plt.ylabel('alt [deg]')
        plt.axvline(6.+hour,ymin=0,ymax=1,linestyle='dotted',alpha=0.5,color='y')
        plt.axvline(19.-hour,ymin=0,ymax=1,linestyle='dotted',alpha=0.5,color='y')

        plt.show()

       

        ending_test()
    except:
        test_error()