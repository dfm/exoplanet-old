
# To make plots for a specific planet
# if clobber is true then overwrite existing plots
# writes hardcopy plots to "."
# idea is that one can use the webpages to find a nice
# planet and then plot it here
# must have an established directory structure from a
# run of ersatz_planets.py
#
# intellectual property:
#   Copyright 2010 David W. Hogg and Adam D. Myers.
#   NYU and UIUC
#   All rights reserved.         Watch this space
#   for an open-source license to be inserted later.
#
# V1.0 ADM August 10th, 2010

from glob import glob
import cPickle
from time import time
import os
from numpy import *

# this rc block must be before the matplotlib import?
from matplotlib import rc
rc('font',**{'family':'serif','serif':'Computer Modern Roman','size':18})
rc('text', usetex=True)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ersatz_exoplanet import *

suffix = '.pdf'

#targdir is the directoy in which the planets are stored
#planet is the exact string representation of the planet we want to plot
def fig_plot_planet(normaldir,likedir,planet,clobber=False):

    plotstart = time()

    if clobber:
        remfile = '*planet'+planet+'*'+suffix
        os.system('/bin/rm -rvf %s' % remfile)

    file1 = open(likedir+'/planet.'+planet+'.pickled', 'r')
    print 'Finding ML index for planet',planet,' t....', time()-plotstart,'secs'
    dic = cPickle.load(file1)
    K = sqrt(dic["Ksin"]**2.+dic["Kcos"]**2.)
    e = sqrt(dic["esin"]**2.+dic["ecos"]**2.)
    T = exp(dic["lnT"])
    pomega = arctan2(dic["esin"],dic["ecos"])
    phi = arctan2(dic["Ksin"],dic["Kcos"]) - pomega
    phi = arctan2(sin(phi), cos(phi))
    ML = argmax(array(dic["lnLikelihood"]))
    plotchainML = column_stack([dic['lnT'],dic['Kcos'],dic['Ksin'],dic['ecos'],dic['esin'],dic['S'],dic['V0'],e,K,T,pomega,phi,log(K)])
    V0ML = (dic['V0'])[ML]
    KML = K[ML]
    eML = e[ML]
    TML = T[ML]
    pomegaML = pomega[ML]
    phiML = phi[ML]
    plotchainML = plotchainML[ML,:]

    file2 = open(normaldir+'/planet.'+planet+'.pickled', 'r')
    print 'Doing planet',planet,' t....', time()-plotstart,'secs'
    dic = cPickle.load(file2)

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

    MML = mean_anomaly_from_time(timegod,TML,phiML)
    EML = eccentric_anomaly_from_mean_anomaly(MML,eML)
    fML = true_anomaly_from_eccentric_anomaly(EML,eML)
    VmodelML = V0ML+rad_vel(KML,fML,pomegaML,eML)

    #plot of true radial velocity for planet
    plt.clf()
    godcolor = 'k'
    godlw = 1.05
    godls = '-'
    godalpha = 0.5
    MLcolor = 'k'
    MLls = '--'
    plt.plot(timegod,VmodelML,color=MLcolor,ls=MLls,lw=godlw)
    plt.plot(timegod,Vmodelgod,color=godcolor,ls=godls,lw=godlw,alpha=godalpha)
    if len(ptime) > 0:
        plt.errorbar(ptime,Vdata,Verr,fmt='kx')
    plt.xlabel("time (d)")
    plt.ylabel("radial velocity (km\,s$^{-1}$)")
    plt.savefig('god_plot_Vt_planet'+planet+suffix)

    #plots of parameter samplings and histograms
    truth = [dic["truthlnT"], dic["truthKcos"],dic["truthKsin"],dic["truthecos"],dic["truthesin"],dic["truthS"],dic["truthV0"],e,K,T,pomega,phi,log(K)]
    pardict = {0:'lnT',1:'Kcos',2:'Ksin',3:'ecos',4:'esin',5:'S',6:'V0',7:'e',8:'K',9:'T',10:'pomega',11:'phi',12:'lnK'}

    K = sqrt(dic["Ksin"]**2.+dic["Kcos"]**2.)
    e = sqrt(dic["esin"]**2.+dic["ecos"]**2.)
    T = exp(dic["lnT"])
    pomega = arctan2(dic["esin"],dic["ecos"])
    phi = arctan2(dic["Ksin"],dic["Kcos"]) - pomega
    phi = arctan2(sin(phi), cos(phi))
    
    plotchain = column_stack([dic['lnT'],dic['Kcos'],dic['Ksin'],dic['ecos'],dic['esin'],dic['S'],dic['V0'],e,K,T,pomega,phi,log(K)])
    #thin to ~10000 for plotting
    thinstep = ceil(len(plotchain)/10000)
    if thinstep > 1:
        print "thinning chain for plotting by a factor of %d" % thinstep
        print shape(plotchain)
        plotchain = plotchain[::thinstep,:]
        print shape(plotchain)

    #plot lnT versus eccentricity
    parx = 7
    pary = 0
    plt.clf()
    #sampling
    plt.plot(plotchain[:,parx],plotchain[:,pary],'k.',ms=1.5,alpha='0.5')
    plt.axvline(plotchainML[parx],color=MLcolor,ls=MLls,lw=godlw)    
    plt.axhline(plotchainML[pary],color=MLcolor,ls=MLls,lw=godlw)
    plt.axvline([truth[parx]],color=godcolor,ls=godls,lw=godlw,alpha=godalpha)
    plt.axhline([truth[pary]],color=godcolor,ls=godls,lw=godlw,alpha=godalpha)
    my= median(plotchain[:,pary])
    plt.axis([0.,1.,my-0.2,my+0.2])
    plt.xlabel('eccentricity $e$')
    plt.ylabel('log period $\ln(T / [1~\mathrm{d}])$')
    plt.savefig('plot_lnT_e_planet'+planet+suffix)

    #this parameter is the eccentricity
    parameter = 7
    parst = str(pardict[parameter])
    plt.clf()
    n, bins, patches = plt.hist(plotchain[:,parameter],20,facecolor='0.5',edgecolor='0.5')
    plt.axvline(plotchainML[parameter],color=MLcolor,ls=MLls,lw=godlw)
    plt.axvline([truth[parameter]],color=godcolor,ls=godls,lw=godlw,alpha=godalpha)
    plt.axis([0.,1.,0.,1.1*max(n)])
    plt.xlabel('eccentricity $e$')
    plt.ylabel("number of samples per bin")
    plt.savefig(parst+'_histo_sampling_planet'+planet+suffix)
    return

if __name__ == '__main__':
    fig_plot_planet('normal','likemaxim','03',0)
    fig_plot_planet('normal','likemaxim','58',0)
