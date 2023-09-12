
# VISIBILITY OF A STAR

**Table of Contents**<a id='toc0_'></a> 

- [**Description of the Project**](#toc1_)
- [**Project tree**](#toc2_)
    - [**Makefile**](#toc2_1_)


---

## <a id='toc1_'></a>[Description of the Project](#toc0_)

This project is an university exercise. 

The `visibility` directory is the python package of functions for the visibility of a star: set a place from which to observe and the date of the observation, it is possible to compute if a target object is visible or not. To have more details see the section about [<u>the script</u>](#toc2_4_).

In addition to this package (that is the aim of the exercise), the project directory contains:

- a jupiter notebook to store and explain in details every single step of the implementation 
- 



## <a id='toc2_'></a>[Project tree](#toc0_)

- [`README.md`](README.md) : informations about the project 

- [`README.txt`](README.txt) : same file as this one in `.txt`

- [`LICENCE.md`](LICENCE.md) : the free licence

- [`Makefile`](Makefile) : compilable file

- [`requirements.txt`](requirements.txt) : require packages

- [`script.py`](script.py) : the script

- [***visibility/***](visibility) - the package

    - [`__init__.py`](visibility/__init__.py) : script

    - [`Angles.py`](visibility/Angles.py) : function script

    - [`coor.py`](visibility/coor.py) : function script

    - [`sky.py`](visibility/sky.py) : function script
    
    - [`stuff.py`](visibility/stuff.py) : function script

    - [***Time/***](visibility/Time) - 
    
        - [`__init__.py`](visibility/Time/__init__.py) : script

        - [`dates.py`](visibility/Time/dates.py) : 

        - [`Tclasses.py`](visibility/Time/Tclasses.py) :

- [***data/***](data) - folder for data

    - [`atm_data.csv`](data/atm_data.csv) : 

    - [`dT_data.csv`](data/dT_data.csv) :

    - [`moon1_data.csv`](data/moon1_data.csv) : 

    - [`moon1_data.csv`](data/moon1_data.csv) : 

    - [`targets.csv`](data/targets.csv) :

- [***notebook/***](notebook) - folder for implementation

    - [`implementation_notebook.ipynb`](notebook/implementation_notebook.ipynb) : jupyter notebook
    
    
- [***tests/***](tests) - tests directory



### <a id='toc2_1_'></a>[Makefile](#toc0_)

Use `proj_2$ make requirements` to install all the necessary packages.


### <a id='toc2_4_'></a>[script.py](#toc0_)

#### <a id='toc1_1_'></a>[Inputs](#toc0_)

The script takes **3 main inputs**: 

1. the **coordinates** and the **proper motion** informations of the star
2. the **coordinates** and the **height** of the observatory location on the Earth
3. the **date** of the observation

#### <a id='toc1_2_'></a>[Outputs](#toc0_)

The script returns 
