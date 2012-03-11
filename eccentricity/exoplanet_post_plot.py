# To plot things after running the code, rather than during code execution
# targdir is where to write the plot html pages to
# if clobber is true then overwrite existing directory
# always makes hardcopies and follows standard html page structuring
# of exoplanet_distribution.py
#
# intellectual property:
#   Copyright 2010 David W. Hogg and Adam D. Myers.
#   NYU and UIUC
#   All rights reserved.         Watch this space
#   for an open-source license to be inserted later.
#
# V1.0 ADM August 4th, 2010

from glob import glob
import cPickle
from time import time
import os
from numpy import *
#import matplotlib       #uncomment these two lines for better
#matplotlib.use('TkAgg') #plotting if you don't have a default matplotlib.rc
import matplotlib.pyplot as plt

def post_plot(targdir,clobber=False):

    plotstart = time()

    if clobber:
        remfile = targdir +'planet*/*png'
        os.system('/bin/rm -rvf %s' % remfile)

    os.chdir(targdir)

    filelist = glob('*pickled*')
    for listedfile in filelist:

        print 'Doing planet',listedfile,' t....', time()-plotstart,'secs'

        outplandir = "".join(listedfile.split('.')[0:2])
        
        file = open(listedfile, 'r')
        dic = cPickle.load(file)
        #radial velocity plots

        os.chdir(outplandir)

        K = sqrt(dic["truthKsin"]**2.+dic["truthKcos"]**2.)
        e = sqrt(dic["truthesin"]**2.+dic["truthecos"]**2.)
        Mp = dic["truthMp"]
        Ms = dic["truthMs"]
        i = dic["truthi"]
        T = exp(dic["truthlnT"])
        pomega = arctan2(dic["truthesin"],dic["truthecos"])
        phi = arctan2(dic["truthKsin"],dic["truthKcos"]) - pomega
        phi = arctan2(sin(phi), cos(phi))
        ptime = dic["time"]
        Vdata = dic["Vdata"]
        Verr = dic["Verr"]
        Vmodel = dic["Vmodel"]
        timegod = dic["timegod"]
        Vmodelgod = dic["Vmodelgod"]

        stringy = "Mp="+str(round(Mp,1))+"Mjup, "
        stringy += "Ms="+str(round(Ms,1))+"Mjup, "
        stringy += "i="+str(round(i*180./pi,1))+"deg, "
        stringy += "T="+str(round(T,1))+"days, "
        stringy += "phi="+str(round(phi,2))+", "
        stringy += "e="+str(round(e,2))+", "
        stringy += "pomega = "+str(round(pomega,2))

        #plot of true radial velocity for planet
        plt.clf()
        plt.suptitle(stringy)
        if len(ptime) > 0:
            plt.errorbar(ptime,Vdata,Verr,fmt='kx')
        plt.xlabel("Time (days)")
        plt.ylabel("Radial Velocity (km/s)")
        plt.plot(timegod,Vmodelgod,color='k')
        plt.savefig('god_plot_Vt.png')

        #plot of data radial velocity for planet
        plt.clf()
        plt.suptitle(stringy)
        if len(ptime) > 0:
            plt.errorbar(ptime,Vdata,Verr,fmt='kx')
        plt.xlabel("Time (days)")
        plt.ylabel("Radial Velocity (km/s)")
        plt.savefig('plot_Vt.png')


        #plots of parameter samplings and histograms
        truth = [dic["truthlnT"], dic["truthKcos"],dic["truthKsin"],dic["truthecos"],dic["truthesin"],dic["truthS"],dic["truthV0"],e,K,T,pomega,phi,log(K)]
        pardict = {0:'lnT',1:'Kcos',2:'Ksin',3:'ecos',4:'esin',5:'S',6:'V0',7:'e',8:'K',9:'T',10:'pomega',11:'phi',12:'lnK'}

        K = sqrt(dic["Ksin"]**2.+dic["Kcos"]**2.)
        e = sqrt(dic["esin"]**2.+dic["ecos"]**2.)
        T = exp(dic["lnT"])
        pomega = arctan2(dic["esin"],dic["ecos"])
        phi = arctan2(dic["Ksin"],dic["Kcos"]) - pomega
        phi = arctan2(sin(phi), cos(phi))

        likechain = array(dic["lnLikelihood"])
        plotchain = column_stack([dic['lnT'],dic['Kcos'],dic['Ksin'],dic['ecos'],dic['esin'],dic['S'],dic['V0'],e,K,T,pomega,phi,log(K)])
        #thin to ~10000 for plotting
        thinstep = ceil(len(plotchain)/10000)
        if thinstep > 1:
            print "thinning chain for plotting by a factor of %d" % thinstep
            print shape(plotchain)
            plotchain = plotchain[::thinstep,:]
            print shape(plotchain)
            print "thinning likechain for plotting by a factor of %d" % thinstep
            print shape(likechain)
            likechain = likechain[::thinstep]
            print shape(likechain)
            #the index of the maximum likelihood
        ML = argmax(likechain)

        for parameter in range(13):
            parst = str(pardict[parameter])
            plt.clf()
            #sampling
            plt.plot(plotchain[:,parameter],'.',ms=1.2,color='k')
            #maximum likelihood
            plt.axhline(plotchain[ML,parameter],color='k')
            #truth
            plt.axhline([truth[parameter]],color='0.75')
            plt.xlabel('Link')
            plt.ylabel(parst)
            plt.savefig(parst+'_plot_sampling.png')
            
            plt.clf()
            plt.hist(plotchain[:,parameter],20,facecolor='0.75')
            #maximum likelihood
            plt.axvline(plotchain[ML,parameter],color='k',ls='-.')
            #truth
            plt.axvline([truth[parameter]],color='0.3',ls='-')
            plt.xlabel(parst)
            plt.ylabel("N("+parst+")")
            plt.savefig(parst+'_histo_sampling.png')

            print 'done plot',parameter,'of 12'

        #plot of likelihood sampling
        likeensemble = sum(likechain,axis=1)
        plt.clf()
        stringy = 'Likelihood Progression'
        plt.suptitle(stringy)
        #The maximum likelihood value
        plt.plot([max(likeensemble)]*len(likechain),'k--')
        #The likelihood sampling
        plt.plot(likeensemble,'.',ms=1.2,color='k')
        plt.axis([0,len(likeensemble),median(likeensemble)-30,median(likeensemble)+30])
        plt.xlabel('Link')
        plt.ylabel('Ensemble ln(Likelihood)')
        plt.savefig('plot_like.png')

        os.chdir('..')
    os.chdir('..')

    return

if __name__ == '__main__':
    from sys import argv
    
    if len(argv) < 2:
        print '***************************************************************************'
        print 'USAGE -'
        print 'python exoplanet_post_plot.py targetdirectory'
        print '     +  (clobberflag (1/0))'
        print '***************************************************************************'
        raise IOError

    targdir = argv[1]
    
    clobberflag = 0

    if len(argv) > 2:
        clobberflag = int(argv[2])

    print 'targdir', targdir, 'clobberflag', clobberflag


    post_plot(targdir,clobberflag)
