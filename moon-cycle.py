import numpy as np
import matplotlib.pyplot as plt
import visibility.sky as sky

if __name__ == '__main__':
    month = 5
    numday = 31
    obs = sky.GeoPos(['-',[151,51,19.44]],['-',[27,47,51.72]],682,'Mt Kent Obs')
    MOON = sky.Moon()
    SUN = sky.Sun()
    k = np.array([])
    days = np.arange(1,numday+1)
#    k0 = np.array([70.4,79.0,86.3,92.2,96.4,99.0,99.8,99.0,96.4,92.2,86.5,79.4,71.2,61.9,51.8,41.3,30.8,20.9,12.2,5.4,1.3,0.2,2.4,7.4,14.8,23.8,33.9,44.2,54.5,64.2,73.2])
    k0 = np.array([78.8,85.6,91.9,96.5,99.3,99.9,98.2,94.0,87.5,79.0,68.8,57.6,46.0,34.7,24.2,15.1,8.0,3.1,0.5,0.2,2.2,6.1,11.7,18.6,26.6,35.3,44.5,54.0,63.4,72.4,80.9])
    for i in range(len(days)):
        date = sky.Date([2023,month,days[i]],0)
        kn = MOON.ill_fract(date,SUN.app_lon(date),SUN.distance(date))
        print(f'{kn*100:.2f} %\t\t{k0[i]:.2f} %\t\t'+date.print_date())
        k = np.append(k,kn*100)
    
    
    plt.figure(figsize=(10,5))
    plt.subplots_adjust(wspace=0.05)
    plt.subplot(2,1,1)
    plt.title(f'$k$ on {sky.Date.MONTHS[month]} 2023 at midnight')
    plt.plot(days,k,'.-',label='$k$')
    plt.plot(days,k0,'.-',label='$k_0$')
    plt.grid(linestyle='--',alpha=0.3)
    plt.xticks(days,np.array([]*len(days)))
    plt.ylabel('$k$ [%]')
    plt.legend()
    plt.subplots_adjust(wspace=0.05)

    plt.subplot(2,1,2)
    diff = abs(k-k0)
    plt.plot(days,diff,'.-',color='g')
    plt.grid(linestyle='--',alpha=0.3)
    plt.xticks(days,days)
    plt.xlabel('day')
    plt.ylabel('$|k-k_0|$ [%]')
    plt.subplots_adjust(wspace=0.05)

    plt.figure(figsize=(10,5))
    plt.subplot(2,1,1)
    plt.title(f'$k$ on {sky.Date.MONTHS[month]} 2023 at midnight')
    plt.plot(days,k,'.-',label='$k$')
    plt.plot(days,k0,'.-',label='$k_0$')
    plt.grid(linestyle='--',alpha=0.3)
    plt.xticks(days,np.array([]*len(days)))
    plt.ylabel('$k$ [%]')
    plt.legend()

    plt.subplot(2,1,2)
    diff = diff/k
    plt.plot(days,diff*100,'.-',color='violet')
    plt.grid(linestyle='--',alpha=0.3)
    plt.xlabel('day')
    plt.ylabel('$\\delta k$ [%]')
    plt.xticks(days,days)


    plt.show()
