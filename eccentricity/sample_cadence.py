import numpy
_DAY= 86400.
_MONTH=27.321661*_DAY
def sample_cadence(pbadweather=0.01):
    """
    NAME:
       sample_cadence
    PURPOSE:
       same a MARVELS like CADENCE
    INPUT:
       pbadweather - chance that a target cannot be observed
    OUTPUT:
    HISTORY:
       2010-07-08 - Written - Bovy (NYU)
    """
    #MARVELS CADENCE:
    # 15 visits during first two of first seven months
    # 1 visit during other of first seven months
    # 2 visits / month the next year for 6.5 months
    # ----------------------------------------------
    # 33
    #
    # How are the 15 visits arranged? assume rather uniform

    firstN= 15
    firstStep= 2.*_MONTH/firstN
    cadence= [ii*firstStep+round((2.*numpy.random.uniform()-1))*_DAY
              for ii in range(firstN)]
    for ii in range(firstN):
        if numpy.random.uniform() < pbadweather:
            cadence.pop(ii)

    secondN= 5
    for ii in range(secondN):
        if numpy.random.uniform() >= pbadweather:
            cadence.append(cadence[-1]+_MONTH+
                            round((8.*numpy.random.uniform()-4)))
    #One year later ...
    cadence.append(cadence[-1]+6.*_MONTH+
                   round((6.*numpy.random.uniform()-3)))
    thirdN= 12
    for ii in range(thirdN):
        if numpy.random.uniform() >= pbadweather:
            cadence.append(cadence[-1]+_MONTH/2.+
                           round((6.*numpy.random.uniform()-3)))
            
    return cadence
    
