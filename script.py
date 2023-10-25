import numpy as np
from visibility.sky import initialize_data, visibility

if __name__ == '__main__':
    data_file = 'targets.csv'
    sel = slice(4,None)
    targets, observatories, dates = initialize_data(data_file,sel=sel)    


    for (obj,obs,date) in zip(targets,observatories,dates):
        _ = visibility(date,obj,obs)