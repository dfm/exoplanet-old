# Plot distributions of parameters (truth)
# over 10000 planets.
# plot eccentricity for both realistic and
# strawmen distributions
#
# intellectual property:
#   Copyright 2010 David W. Hogg and Adam D. Myers.
#   NYU and UIUC
#   All rights reserved.         Watch this space
#   for an open-source license to be inserted later.
#
# V1.0 ADM August 10th, 2010

# this rc block must be before the matplotlib import?
from matplotlib import rc
rc('font',**{'family':'serif','serif':'Computer Modern Roman','size':18})
rc('text', usetex=True)
# now import matplotlib
import matplotlib
matplotlib.use('Agg')
import os
from ersatz_exoplanet import ersatz_planets

if __name__ == '__main__':
 
    planets = 1000

    print 'plotting distributions'
    fp = ersatz_planets(planets,1,crazy=False)
    print 'initialized normal ersatz_planet class'
    #parameter distributions
    fp.plot_edist(hardcopy=1)
    fp.plot_Mpdist(hardcopy=1)
    fp.plot_Tdist(hardcopy=1)
    fp.plot_idist(hardcopy=1)
    fp.plot_Kdist(hardcopy=1)

    print 'initializing crazy ersatz_planet class'
    fp = ersatz_planets(planets,1,crazy=True)
    print 'initialized  crazy ersatz_planet class'
    #parameter distributions
    fp.plot_edist(hardcopy=1)
