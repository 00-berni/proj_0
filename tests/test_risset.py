from test_struc import starting_test, ending_test, test_error,import_pack

import_pack()
from visibility.stuff import *
from visibility.Angles import *
from visibility.Time.Tclasses import *
from visibility.Time.dates import *
from visibility.coor import *
from visibility.sky import *


def tran_ris_set(date: Date, obs_pos: GeoPos, obj: StarObj, test: bool = False,stat: bool = False, results: bool = False) -> None | tuple:
        date = date.copy()
        lon = obs_pos.lon
        lat = obs_pos.lat
        # h0 = Angles(-0.5667,'deg',lim=90)
        h0 = refraction_corr(Angles(0.,'deg',lim=90),obs_pos.h, alt0=True)
        tmpdate = Date(date.date)
        Dt = time_correction(tmpdate.date[0])
        GHA = mean_Green_HA(tmpdate)
        tmpdate.time.tytime = 'TD'
        day = tmpdate.date[-1]
        if test:
            stat = True
            a1 = HAngles(40.68021,'deg')
            d1 = HAngles(18.04761,'deg')
            a2 = HAngles(41.73129,'deg')
            d2 = HAngles(18.44092,'deg')
            a3 = HAngles(42.78204,'deg')
            d3 = HAngles(18.82742,'deg')
        else:    
            a2, d2 = precession_corr(tmpdate,obj)
            tmpdate.date[-1] = day - 1
            a1, d1 = precession_corr(tmpdate,obj)
            tmpdate.date[-1] = day + 1
            a3, d3 = precession_corr(tmpdate,obj)
        del tmpdate

        cosH0 = (np.sin(h0.rad) - np.sin(lat.rad)*np.sin(d2.rad)) / (np.cos(lat.rad)*np.cos(d2.rad))
        
        if stat:
            results = True
            print(f'h0 = {h0.deg} deg')
            print(f'Dt = {Dt} s')
            print('GHA = ' + GHA.str_angle('hms',True) + ' =',GHA.deg,'deg')
            print(f'cosH0 = {cosH0}')

        if abs(cosH0) > 1:
            return None
        else:
            H0 = HAngles(np.arccos(cosH0),'rad',lim=180)
            m = [(a2 + lon - GHA).deg / 360]
            m += [m[0] - H0.deg/360]
            m += [m[0] + H0.deg/360]
            m = np.array(m)
            if stat:
                print(f'H0 = {H0.deg}')
            for k in range(3):
                if stat:
                    print(f'\nIteration {k}:')
                LST = []
                for i in range(3):
                    if abs(m[i]) > 1:
                        m[i] -= np.sign(m[i])
                    if m[i] < 0:
                        m[i] += 1
                    LST += [GHA.deg + 360.985647*m[i]]
                    if LST[i] > 360:
                        LST[i] -= 360*(LST[i]//360)
                    if test:
                        print(f'LST[{i+1}] =', LST[i])
                if stat:
                    print('m = ', m)
                LST = HAngles(np.array(LST),'deg')
                n = m + Dt/86400
                a = np.array([interpole_three([a1.deg,a2.deg,a3.deg],ni+day,[day-1,day,day+1]) for ni in n])
                a = HAngles(a,'deg')
                d = np.array([interpole_three([d1.deg,d2.deg,d3.deg],ni+day,[day-1,day,day+1]) for ni in n[1:]])
                d = HAngles(d,'deg',lim=90) 
                H = LST - lon - a
                h = np.arcsin(np.sin(lat.rad)*np.sin(d.rad) + np.cos(lat.rad)*np.cos(d.rad)*np.cos(H.rad[1:]))
                h = Angles(h,'rad',lim=90)
                Dm = np.array([-H.deg[0]/360])
                Dm = np.append(Dm,((h-h0).deg / 360 / (np.cos(d.rad)*np.cos(lat.rad)*np.sin(H.rad[1:]))))
                if stat:
                    print('alpha = ', a.deg) 
                    print('delta = ', d.deg) 
                    print('H = ',H.deg)
                    print('h = ',h.deg)
                    print('Dm = ', Dm)
                    print('m = ', m)
                m += Dm
                for i in range(3):
                    if abs(m[i]) > 1:
                        m[i] -= np.sign(m[i])
                    if m[i] < 0:
                        m[i] += 1
            if results:
                transit = Date(date.date,Time(m[0]*86400),epoch=date.epoch)
                rising  = Date(date.date,Time(m[1]*86400),epoch=date.epoch)
                setting = Date(date.date,Time(m[2]*86400),epoch=date.epoch)
                print()
                time = date.time.val
                event = [transit,rising,setting]
                names = ['transit','rising ','setting']
                for i in range(3):
                    if event[i].time.val > 86400:
                        days = event[i].time.val // 86400
                        event[i].date[-1] += days
                        event[i].time.val -= days * 86400
                    elif event[i].time.val < time:
                        event[i].date[-1] += 1
                    print(names[i] + ':\t' + event[i].str_date())
            return m*86400


def nullobj():
    return StarObj('',Equatorial(None,None),None)


if __name__ == '__main__':
    starting_test('TEST RISING, SETTING, TRANSIT')
    try:
        obj = nullobj()
        date = Date([1988,3,20])
        obs = GeoPos(['+',[71,5,0]],['+',[42,20,0]],0)
        m = tran_ris_set(date,obs,obj,test=True)
        # print('\ntransit:\t'+Time(m[0]).str_time())
        # print('rising :\t'+Time(m[1]).str_time())
        # print('setting:\t'+Time(m[2]).str_time())

        # - - -
        print()
        lat = 43.93039980730644
        lon = 10.907521220934262 
        h = 65
        obs = GeoPos(lon,lat,h)
        name = 'Pollux'
        ra  = HAngles(['+',[7,45,16.42]],'hms')
        dec = HAngles(['+',[28,1,34.5]],'deg',lim=90)
        # proper motion
        mu_a = -626.55e-3 # as/yr
        mu_d = -45.80e-3  # as/yr
        # ris 06:40
        # tra 14:46
        # set 22:51
        obj = StarObj(name,[ra,dec],[mu_a,mu_d])
        obj.obj_info()
        print(date.str_date())
        print()
        
        date = Date([2023,6,26],[11,5,20])

        m = tran_ris_set(date,obs,obj,stat=True)
        # print('\ntransit:\t'+Time(m[0]).str_time())
        # print('rising :\t'+Time(m[1]).str_time())
        # print('setting:\t'+Time(m[2]).str_time())
                
        ending_test()
    except:
        test_error()