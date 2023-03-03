# Visibility

## Table of Contents

1. [Structure of Project Folder](#structure-of-project-folder)
2. [Task](#task)

---

## Structure of Project Folder

- **_./_**
  - **`README.md`**
    File to summerise the aim of the script and the structure of the project
  - **`requirements.txt`**
    List of required packages of Python 3
  - **`Makefile`**
  - **_notebook/_**
    - **`notebook.ipynb`**
      Notebook used to write the code
  - **_script/_**
    - **`script.py`**
      The script for the visibility project
    - **`functions.py`**
      File with my own functions, implemented for the code
  - **_tests/_**
    - **`test01.py`**

---

## Task

**Input**

- Coordinates (alpha, delta) of the object
- Coordinates (long, lat, alt, zone time) of obs position on the Earth 
- The time of the obs [YY, MM, DD, hh, mm, ss.ss]

**Output**

- Coordinates (HA, zenit) of the object
- Coordinates (HA, zenit) of the Sun
- Coordinates (HA, zenit) of the Moon
- Airmass

The visibility plot

## References

 1. _Explanatory supplement to the Astronomical almanac_
 2. _Textbook on Spherical Astronomy_