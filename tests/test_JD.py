import numpy as np
from test_struc import *
import_pack()
from visibility.Time.Tclasses import Date
from visibility.Time.dates import julian_day, calendar_date



if __name__ == '__main__':

    starting_test('TEST FOR JULIAN DAY')

    try:
        years, months, days, cal, JD0 = get_test_data('test_JD_data.csv')

        JD = np.array([julian_day(Date([years[i],months[i],days[i]]),cal[i]) for i in range(len(years))])
        print('> Compute the JD')
        print('JD0')
        cnt = 0
        for jd, jd0 in zip(JD,JD0):
            if jd == jd0:
                print(f'{jd0}: OK!')
            else:
                cnt += 1
                print(f'{jd0}: NO MATCH --> DIFF = {jd-jd0}\n\tJD = {jd}')
        if cnt == 0:
            print('> All results match')
        else:
            print(f'> {cnt} results do not match')
            raise Exception('No matches')
        print('> Test the inverse function')
        for i in range(len(years)):
            print(Date([years[i],months[i],days[i]]).str_date()+' ==> '+calendar_date(JD[i]).str_date())

        ending_test()
    except:
        test_error()