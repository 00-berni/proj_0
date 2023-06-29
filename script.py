import matplotlib.pyplot as plt
from visibility.stuff import get_data
from visibility.Time.Tclasses import Date 
from visibility.sky import *

def import_data(filename, sel: int | slice = slice(None), delimiter: str =',', sep: str = ':'):
    names, ras, decs, muas, muds = get_data(filename,delimiter=delimiter)
    
    prmts = []
    for i in range(len(names)):
        alpha = ras[i]
        delta = decs[i]
        mua = muas[i]
        mud = muds[i]

        ras[i] = [alpha[0],[float(val) for val in alpha[1:].split(sep)]]
        decs[i] = [delta[0],[float(val) for val in delta[1:].split(sep)]]
        
        mua = float(mua) if mua != 'None' else None
        mud = float(mud) if mud != 'None' else None
        prmts += [[mua,mud]]
    
    return names[sel], ras[sel], decs[sel], prmts[sel]




if __name__ == '__main__':

    Npnt = 200

    lat = 43.93039980730644
    lon = 10.907521220934262 
    h = 65

    obs = GeoPos(lon,lat,h)

    name, ra, dec, prmt = import_data('targets.csv', sel=4)

    obj = StarObj(name,[ra,dec],prmt)
    obj.obj_info()

    # ---

    print()
    date = Date([1898,6,23],[20,28,0])
    print('Observation on '+date.str_date())
    alt, dayrange = trajectory(date,obs,obj,Npnt) 
    m = tran_ris_set(date,obs,obj,True)
    event = Date(date.date,m)
    ealt = compute_alt(event,obs,obj)


    hour = date.time.hour() 

    # dayrange = np.array([calendar_date(day/24+JDdate).daydec() for day in np.linspace(0,24,numpoint)])

    N = int(max(dayrange))+1 if hour <= (int(hour)+0.5) else int(max(dayrange))+2
    ticks = np.arange(int(min(dayrange)), N,1)

    plt.figure(figsize=[12,8])
    ax = plt.axes()
    ax.set_facecolor('black')
    plt.title(obj.name + ' on ' + date.str_date() + f" at\nlat: {obs.lat.str_angle('deg')}, lon: {obs.lon.str_angle('deg')}, h: {obs.h} m a.s.l.")
    
    plt.plot(dayrange,alt.deg,label='trajectory')
    plt.plot(dayrange[0],alt.deg[0],'ow',label='start')
    if type(m.val) != np.ndarray:
        plt.plot(m.hour(),ealt.deg,'vg',label='transit')
    else:
        plt.plot(m.hour()[0],ealt.deg[0],'vg',label='transit')
        plt.plot(m.hour()[1],ealt.deg[1],'vy',label='rising')
        plt.plot(m.hour()[2],ealt.deg[2],'vr',label='setting')


    plt.axhline(0,xmin=0,xmax=1,linestyle='dashed',alpha=0.5,color='white',label='horizon')
    plt.xticks(ticks,np.where(ticks < 24, ticks, ticks-24))
    plt.xlabel('UT [h]')
    plt.ylabel('alt$_0$ [deg]')
    plt.grid(axis='x',linestyle='dotted',color='gray',alpha=0.7)
    plt.legend(numpoints=1)


    print('\n - - -')            


    print()
    date = Date([2023,6,26],[11,5,20])
    print('Observation on '+date.str_date())
    alt, dayrange = trajectory(date,obs,obj,Npnt)
    m = tran_ris_set(date,obs,obj,True)
    event = Date(date.date,m)
    ealt = compute_alt(event,obs,obj)


    hour = date.time.hour() 

    # dayrange = np.array([calendar_date(day/24+JDdate).daydec() for day in np.linspace(0,24,numpoint)])

    N = int(max(dayrange))+1 if hour <= (int(hour)+0.5) else int(max(dayrange))+2
    ticks = np.arange(int(min(dayrange)), N,1)

    plt.figure(figsize=[12,8])
    ax = plt.axes()
    ax.set_facecolor('black')
    plt.title(obj.name + ' on ' + date.str_date() + f" at\nlat: {obs.lat.str_angle('deg')}, lon: {obs.lon.str_angle('deg')}, h: {obs.h} m a.s.l.")
    
    plt.plot(dayrange,alt.deg,label='trajectory')
    plt.plot(dayrange[0],alt.deg[0],'ow',label='start')
    if type(m.val) != np.ndarray:
        plt.plot(m.hour(),ealt.deg,'vg',label='transit')
    else:
        plt.plot(m.hour()[0],ealt.deg[0],'vg',label='transit')
        plt.plot(m.hour()[1],ealt.deg[1],'vy',label='rising')
        plt.plot(m.hour()[2],ealt.deg[2],'vr',label='setting')

    plt.axhline(0,xmin=0,xmax=1,linestyle='dashed',alpha=0.5,color='white',label='horizon')
    plt.xticks(ticks,np.where(ticks < 24, ticks, ticks-24))
    plt.xlabel('UT [h]')
    plt.ylabel('alt$_0$ [deg]')
    plt.grid(axis='x',linestyle='dotted',color='gray',alpha=0.7)
    plt.legend(numpoints=1)


    plt.show()
