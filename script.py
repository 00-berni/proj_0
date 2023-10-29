import numpy as np
from visibility.sky import initialize_data, visibility

if __name__ == '__main__':
    data_file = 'targets.csv'
    st, end = 0,None
    sel = slice(st,end)
    window = True
    save_fig = False
    targets, observatories, dates = initialize_data(data_file,sel=sel)    


    for (obj,obs,date) in zip(targets,observatories,dates):
        _ = visibility(date,obj,obs,win=window,save_fig=save_fig)