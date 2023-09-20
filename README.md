
# VISIBILITY OF A STAR

**Table of Contents**<a id='toc0_'></a> 

- [**Description of the Project**](#toc1_)
- [**Project tree**](#toc2_)
    - [**Makefile**](#toc2_1_)


---

## <a id='toc1_'></a>[Description of the Project](#toc0_)

This project is an university assignment.

```
Author: Bernardo Vettori
Date:   2023/09/20
Ver:    v23
```

### <a id='toc1_1_'></a>[Task](#toc0_)

The aim of the exercise is to compute the visibility of a star from a set place and on a chosen date. The required procedure can be summed up in 5 steps:

1. compute the rise and the set of the star from the observatory
2. compute the rise and the set of the Sun
3. compute the twilight
4. compute the distance of the target from the Moon and the moon phase
5. estimate the visibility window, if any

The python library `visibility` collects all the functions to implement this tasks. For more details see the section about the library itself or the script.

### <a id='toc1_2_'></a>[Quick commands](#toc0_)

- The list of used python packages and their versions is in the file requirements.txt. Hence, one have to install them to run the script preventing possible errors of compatibility.

    A quick way to do that is to use the command: `$ make requirements`.

- The fastest and safe procedure to run the script is to use the command: `$ make script`.

For other commands see the section about makefile.

### <a id='toc1_3_'></a>[Project directory overview](#toc0_)

In addition to the library (that is the aim of the exercise), the project directory contains:

- a jupiter notebook to store and to explain in details every single step of the implementation (see the notebook section)
- a compilable file to do some operation quickly (see commands section or makefile section)
- a directory for tests (see the test section)


## <a id='toc2_'></a>[Project tree](#toc0_)

- [***data/***](data) - data directory

    - [`atm_data.csv`](data/atm_data.csv) : table of standard atmospheric parameters [[1]](#1_smith)

    - [`dT_data.csv`](data/dT_data.csv) : table of $TD - UT$ measures in year range $[1620,1998]$ [[3]](#3_meus)

    - [`moon1_data.csv`](data/moon1_data.csv) : table of periodic terms for the longitude and the distance from the Earth of the Moon [[3]](#3_meus) 

    - [`moon1_data.csv`](data/moon1_data.csv) : table of periodic terms for the latitude of the Moon [[3]](#3_meus)

    - [`targets.csv`](data/targets.csv) : list of needed parameters of the chosen targets

- [***notebook/***](notebook) - folder for implementation

    - [`implementation_notebook.ipynb`](notebook/implementation_notebook.ipynb) : jupyter notebook

- [***visibility/***](visibility) - the library

    - [***Time/***](visibility/Time) - additional package

        - [`__init__.py`](visibility/Time/__init__.py)

        - [`dates.py`](visibility/Time/dates.py) : functions to compute HA and LST

        - [`Tclasses.py`](visibility/Time/Tclasses.py) : classes for time and date

    - [`__init__.py`](visibility/__init__.py)

    - [`Angles.py`](visibility/Angles.py) : classes for angle types

    - [`coor.py`](visibility/coor.py) : classes for celestial coordinates and functions to convert them

    - [`sky.py`](visibility/sky.py) : functions to compute trajectory along the sky
    
    - [`stuff.py`](visibility/stuff.py) : functions to import files and to interpolate data

- [***tests/***](tests) - tests directory

- [`.gitignore`](.gitignore)

- [`LICENCE.md`](LICENCE.md) : the free licence

- [`Makefile`](Makefile) : compilable file for useful commands

- [`README.md`](README.md) : informations about the project     

- [`README.txt`](README.txt) : same file as this one in `.txt`

- [`requirements.txt`](requirements.txt) : require packages abd versions

- [`script.py`](script.py)
        

### <a id='toc2_1_'></a>[visibility](#toc0_)

The library is made up 1 package and 4 modules:

- **_Time/_** package
    - `dates.py`
    
        This module collects the functions to compute the Greenwich hour angle and the local sidereal time.

    - `Tclasses.py`

        This module implements two classes:

        - The `Time` class, which stores the value of the time in seconds and the type (TD or UT). It allows also some manipulations such as converting TD in UT and vice versa.

        - The `Date` class, which is very useful due to the fact that it stores a date with time, its corresponding julian day and if the date is a civil local one, the class converts it in UT, storing informations about the time zone and the daylight saving. It is possible to pass directly a julian day, too.

- `Angles.py`

    This module implements two classes for the manipulation of angles. The class `Angles` collecs indeed the value of any angle in both degrees and radiants and it simplifies the mathematical operation between angles. The class `HAngles` is merely the same class, but it stores the value in hours, too. 

- `coor.py` 
- `sky.py` 
- `stuff.py` 

### <a id='toc2_2_'></a>[Makefile](#toc0_)


### <a id='toc2_4_'></a>[script.py](#toc0_)

#### <a id='toc1_1_'></a>[Inputs](#toc0_)

The script takes **3 main inputs**: 

1. the **coordinates** and the **proper motion** informations of the star
2. the **coordinates** and the **height** of the observatory location on the Earth
3. the **date** of the observation

#### <a id='toc1_2_'></a>[Outputs](#toc0_)

The script returns 

## <a id='toc4_'></a>[References](#toc0_)

1. <a id='1_smith'></a> P. Duffett-Smith and J. Zwart., _Practical Astronomy with your Calculator or Spreadsheet_, Cambridge University Press, 2011.
2. <a id='2_icao'></a> _Manual of the ICAO Standard Atmosphere_, ICAO, 3rd edition, 1993.
3. <a id='3_meus'></a> J. Meeus, _Astronomical Algorithms_, Willmann-Bell Inc., 1998.
4. <a id='4_notes'></a> _Notes of Astrofisica Osservativa course_, 2021-2022.
5. <a id='5_almanac'></a> P. K. Seidelmann, _Explanatory supplement to the Astronomical almanac_,University Science Books, 2006
