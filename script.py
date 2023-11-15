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


    for (obj,obs,date) in zip(targets,observatories,dates):
        _ = visibility(date,obj,obs,win=window,save_fig=save_fig,display_plots=displot)