
import os.path as ph
import numpy as np

TEST_DIR = ph.dirname(ph.realpath(__file__))

SEP = lambda obj : '------' + obj + '------\n'

def sep_line():
    print(SEP('----------'))

def starting_test(title: str) -> None:
    print(SEP(title))
    print('> START TEST')

def ending_test() -> None:
    print('\n> TEST COMPLETE')
    sep_line()

def test_error() -> None:
    print('\n> TEST FAILD')
    sep_line()
    raise

def get_test_data(filename: str, unpack: bool = True) -> np.ndarray:
    TEST_FILE = ph.join(TEST_DIR,filename)
    print('> Open the test file: ' + TEST_FILE)
    ext = filename[-3:]
    if ext == 'csv':
        from pandas import read_csv
        data = read_csv(TEST_FILE).to_numpy().transpose()
    elif ext == 'txt':
        data = np.loadtxt(TEST_FILE,unpack=unpack)
    else: raise Exception(f'.{ext} is not allowed!\nOnly .txt or .csv files')
    return data
