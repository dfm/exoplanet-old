#!/bin/tcsh

#make ersatz planets of all e and m combinations

rm *madeplanets

python -u ersatz_exoplanet.py 1000  100 30 dum            0 1 0
python -u ersatz_exoplanet.py 1000  100 30 dum            0 1 1
python -u ersatz_exoplanet.py 1000  100 30 dum            0 1 2

python -u ersatz_exoplanet.py 1000  100 30 dum            0 0 0
python -u ersatz_exoplanet.py 1000  100 30 dum            0 0 1
python -u ersatz_exoplanet.py 1000  100 30 dum            0 0 2
