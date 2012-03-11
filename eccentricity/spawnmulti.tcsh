#!/bin/tcsh

rm e*ml.*log
rm e*pd.*log

#here's how to make new planets
#
#makeplanmulti.tcsh

( python -u ersatz_exoplanet.py 400 100000 30 e0m1ml 1 0 1 > e0m1ml.log ) >& e0m1ml.errlog &
( python -u ersatz_exoplanet.py 400 100000 30 e0m1pd 0 0 1 > e0m1pd.log ) >& e0m1pd.errlog &
( python -u ersatz_exoplanet.py 400 100000 30 e0m2ml 1 0 2 > e0m2ml.log ) >& e0m2ml.errlog &
( python -u ersatz_exoplanet.py 400 100000 30 e0m2pd 0 0 2 > e0m2pd.log ) >& e0m2pd.errlog &
