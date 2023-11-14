import numpy as np
import matplotlib.pyplot as plt
from visibility.sky import initialize_data, visibility
import visibility.sky as sky

if __name__ == '__main__':
    data_file = 'targets.csv'
    st, end = 0,6#,None
    sel = slice(st,end)
    window = True
    save_fig = False
    displot = False
    targets, observatories, dates = initialize_data(data_file,sel=sel)    


    # for (obj,obs,date) in zip(targets,observatories,dates):
    #     _ = visibility(date,obj,obs,win=window,save_fig=save_fig,display_plots=displot)

    obs = sky.GeoPos(['-',[151,51,19.44]],['-',[27,47,51.72]],682,'Mt Kent Obs')
    MOON = sky.Moon()
    SUN = sky.Sun()
    k = np.array([])
    days = np.arange(1,32)
    k0 = np.array([70.4,79.0,86.3,92.2,96.4,99.0,99.8,99.0,96.4,92.2,86.5,79.4,71.2,61.9,51.8,41.3,30.8,20.9,12.2,5.4,1.3,0.2,2.4,7.4,14.8,23.8,33.9,44.2,54.5,64.2,73.2])
    for i in range(len(days)):
        date = sky.Date([2023,1,days[i]],0)
        kn = MOON.ill_fract(date,SUN.app_lon(date),SUN.distance(date))
        print(f'{kn*100:.2f} %\t\t{k0[i]:.2f} %\t\t'+date.print_date())
        k = np.append(k,kn*100)
    plt.figure()
    plt.subplot(2,1,1)
    plt.title('January 2023')
    plt.plot(days,k,'.-',label='$k$')
    plt.plot(days,k0,'.-',label='$k_0$')
    plt.grid(linestyle='--',alpha=0.3)
    plt.xticks(days,np.array([]*len(days)))
    plt.ylabel('$k$ [%]')
    plt.legend()

    plt.subplot(2,1,2)
    diff = abs(k-k0)
    plt.plot(days,diff,'.-',color='g')
    plt.grid(linestyle='--',alpha=0.3)
    plt.xticks(days,days)
    plt.ylabel('$|k-k_0|$')

    plt.figure()
    plt.subplot(2,1,1)
    plt.title('January 2023')
    plt.plot(days,k,'.-',label='$k$')
    plt.plot(days,k0,'.-',label='$k_0$')
    plt.grid(linestyle='--',alpha=0.3)
    plt.xticks(days,np.array([]*len(days)))
    plt.ylabel('$k$ [%]')
    plt.legend()

    plt.subplot(2,1,2)
    diff = diff/k
    plt.plot(days,diff,'.-',color='violet')
    plt.grid(linestyle='--',alpha=0.3)
    plt.xlabel('day')
    plt.ylabel('$\\frac{|k-k_0|}{k}$ [%]')
    plt.xticks(days,days)


    plt.show()
