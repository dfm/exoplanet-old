# name:
#   eccentricity-distribution.py
# purpose:
#   Read in eccentricity samplings, return eccentricity distribution.
# author:
#   David W. Hogg (NYU) & Adam Myers (UIUC)
# license:
#   Copyright 2010 Hogg & Myers; all rights reserved (for now).

import numpy as np
# this rc block must be before the matplotlib import?
from matplotlib import rc
rc('font',**{'family':'serif','serif':'Computer Modern Roman','size':18})
rc('text', usetex=True)
# now import matplotlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import cPickle as pickle
import os
import glob
from deconvolve_sampling import *

def pickle_to_file(data, fn):
    print 'writing %s' % fn
    f = open(fn, 'wb')
    # MAGIC -1: highest pickle protocol
    pickle.dump(data, f, -1)
    f.close()

def unpickle_from_file(fn):
    print 'reading %s' % fn
    f = open(fn, 'rb')
    data = pickle.load(f)
    f.close()
    return data

def weighted_mean(x,w):
	return (np.sum(x * w) / np.sum(w))

def concatenate_pickle_files(keyword, path):
    outfn = path + 'e.pickled'
    if not os.path.exists(outfn):
        baseurl = 'http://geordie.astro.illinois.edu/forhogg2/exoplanets/' + keyword + '/'
        cmd = 'mkdir -p %s' % path
        print cmd
        os.system(cmd)
        filelist = []
        nPlanet = 300
        eML = np.zeros(nPlanet)
        eTrue = np.zeros(nPlanet)
        eChainSet = []
        j = 0
        k = 0
        while j < nPlanet:
            print j
            fn = 'planet.%2.2d.pickled' % k
            cmd = 'curl -o %s%s %s%s' % (path, fn, baseurl, fn)
            print cmd
            os.system(cmd)
            filelist.append(path+fn)
            exop = unpickle_from_file(path+fn)
            cmd = '/bin/rm -vf %s' % (path+fn)
            print cmd
            os.system(cmd)
            K = np.sqrt(exop['Kcos']**2 + exop['Ksin']**2)
            KMean = np.mean(K)
            KSNR = KMean / np.std(K)
            print 'K mean and s/n:', KMean, KSNR
            e = np.sqrt(exop['ecos']**2 + exop['esin']**2)
            eMean = np.mean(e)
            print 'e mean and s/n:', eMean, eMean / np.std(e)
            if True: # KSNR > 1.:
                eChainSet.append(e)
                eML[j] = e[np.argmax(exop['lnLikelihood'])]
                eTrue[j] = np.sqrt(exop['truthecos']**2 + exop['truthesin']**2)
                j += 1
            else:
                print 'K s/n too low'
            k += 1
        pickle_to_file((eTrue, eML, eChainSet), outfn)
    return outfn

def eHistogram(e):
    nBins = 20
    eIndx = np.floor(e * nBins).astype(int)
    alpha = np.zeros(nBins)
    for i in range(nBins):
        alpha[i] = len((np.where(eIndx == i))[0])
    alpha[0] += len((np.where(eIndx < 0))[0])
    alpha[nBins-1] += len((np.where(eIndx >= nBins))[0])
    return alpha.astype(float) * nBins / len(e)

def download_infer_and_plot(keyword, Gaussian=False, Turner=False, fit=True, nObj=300, nSampleMax=None, beta=False):
    path = 'deconvolve/' + keyword + '/'
    parstr = 'alpha'
    if beta:
        parstr += '-beta'
    outfn = path + '%s.%d.pickled' % (parstr, nObj)
    suffix = '_%s_%d.pdf' % (keyword, nObj)
    if beta:
        suffix = '_beta' + suffix
    if not nSampleMax is None:
        outfn = path + 'alpha.%d.%d.pickled' % (nSampleMax, nObj)
        suffix = '_%d%s' % (nSampleMax, suffix)
    fn = concatenate_pickle_files(keyword, path)
    (eTrue, eML, eChainSet) = unpickle_from_file(fn)
    eTrue = eTrue[:nObj]
    eML = eML[:nObj]
    eChainSet = eChainSet[:nObj]
    if fit:
        if not os.path.exists(outfn):
            if not nSampleMax is None:
                eNewChainSet = []
                for e in eChainSet:
                    nlink = len(e)
                    thinstep = np.floor(float(len(e)+1)/float(nSampleMax)).astype('int')
                    thinstep = max(thinstep, 1)
                    eNewChainSet.append(e[::thinstep])
                eChainSet = eNewChainSet
            if beta:
                nLink = 2000
                proposalWidth = 0.15
                if nObj < 100:
                    proposalWidth = 0.30
            else:
                nLink = 10000
                proposalWidth = 0.15
            (alphaBest, alphaChain) = sampleAlpha(eChainSet, beta=beta, nLink=nLink, proposalWidth=proposalWidth)
            pickle_to_file((alphaBest, alphaChain), outfn)
        (alphaBest, alphaChain) = unpickle_from_file(outfn)
    plt.clf()
    trueLinestyle = 'k:'
    trueAlpha = 0.75
    trueLinewidth = 1.5
    if Gaussian:
        xG = np.arange(1001).astype(float) / 1000.
        mean = 0.3
        variance = 0.05**2
        yG = 1./np.sqrt(2. * np.pi * variance) * np.exp(-0.5 * (xG - mean)**2 / variance)
        plt.plot(xG, yG, trueLinestyle, alpha=trueAlpha, lw=trueLinewidth)
    if Turner:
        nT = 1000
        xT = np.arange(nT+1).astype(float) / float(nT)
        yT = (1./((1.+xT)**4.))-(xT/(2.**4.))
        yT = yT / np.sum(yT / nT)
        plt.plot(xT, yT, trueLinestyle, alpha=trueAlpha, lw=trueLinewidth)
    plotStep(eHistogram(eTrue))
    plt.axvline(np.sum(eTrue) / len(eTrue), color='k', alpha=0.5, lw=0.5)
    plt.xlim(0.,1.)
    if Gaussian:
        plt.ylim(0.,10.)
    if Turner:
        plt.ylim(0.,5.)
    plt.xlabel('eccentricity $e$')
    plt.ylabel('frequency $p(e)$')
    titlePrefix = '%d stars' % nObj
    plt.title('%s / truth' % titlePrefix)
    plt.savefig(path + 'eTruth' + suffix)
    plt.clf()
    if Gaussian:
        plt.plot(xG, yG, trueLinestyle, alpha=trueAlpha, lw=trueLinewidth)
    if Turner:
        plt.plot(xT, yT, trueLinestyle, alpha=trueAlpha, lw=trueLinewidth)
    plotStep(eHistogram(eML))
    plt.axvline(np.sum(eML) / len(eML), color='k', alpha=0.5, lw=0.5)
    plt.xlim(0.,1.)
    if Gaussian:
        plt.ylim(0.,10.)
    if Turner:
        plt.ylim(0.,5.)
    plt.xlabel('eccentricity $e$')
    plt.ylabel('frequency $p(e)$')
    plt.title('%s / ML estimates' % titlePrefix)
    plt.savefig(path + 'eML' + suffix)
    plt.clf()
    plt.plot(eTrue, eML, 'ko', alpha=0.5)
    plt.plot([0,1], [0,1], 'k-', alpha=0.5)
    plt.xlim(0.,1.)
    plt.ylim(0.,1.)
    plt.xlabel('true eccentricity $e$')
    plt.ylabel('maximum-likelihood eccentricity $e$')
    plt.title(titlePrefix)
    plt.savefig(path + 'eCompareML' + suffix)
    if fit:
        plt.clf()
        if Gaussian:
            plt.plot(xG, yG, trueLinestyle, alpha=trueAlpha, lw=trueLinewidth)
        if Turner:
            plt.plot(xT, yT, trueLinestyle, alpha=trueAlpha, lw=trueLinewidth)
        thinstep = np.floor(len(alphaChain)/8.).astype(int)
        for (lnP, alpha) in alphaChain[thinstep-1::thinstep]:
            if beta:
                plotBeta(alpha, alpha=0.5)
            else:
                plotStep(alpha, dots=False, alpha=0., edgecolor='0.5')
        if beta:
            dx = 0.001
            egrid = np.arange(0.5*dx,1.,dx)
            numerator = 0. * egrid
        else:
            numerator = 0.*alpha
        denominator = 0.
        for (lnP, alpha) in alphaChain:
            if beta:
                numerator += rawBeta(egrid, alpha)
            else:
                numerator += alpha
            denominator += 1.
        if beta:
            ygrid = numerator / denominator
            plt.plot(egrid, ygrid, 'k-', lw=2)
            plt.axvline(np.sum(egrid * ygrid) / np.sum(ygrid), color='k', alpha=0.5, lw=0.5)
        else:
            marginalizedAlpha = numerator / denominator
            plotStep(marginalizedAlpha)
            nAlpha = len(marginalizedAlpha)
            # next line must be synchronized with deconvolve_sampling.py
            xAlpha = (np.arange(nAlpha) + 0.5) / float(nAlpha)
            plt.axvline(np.sum(xAlpha * marginalizedAlpha) / np.sum(marginalizedAlpha), color='k', alpha=0.5, lw=0.5)
        plt.xlim(0.,1.)
        if Gaussian:
            plt.ylim(0.,10.)
        if Turner:
            plt.ylim(0.,5.)
        plt.xlabel('eccentricity $e$')
        plt.ylabel('frequency $p(e)$')
        if not nSampleMax is None:
            titlePrefix = '%s / %s samples per star' % (titlePrefix, nSampleMax)
        plt.title('%s / inferred distribution' % titlePrefix)
        plt.savefig(path + 'eInferred' + suffix)

        if not beta:
            # NOW: go hierarchical
            nBins = len(marginalizedAlpha)
            # woah, check this out:
            eMOS = [weighted_mean(e,np.exp(marginalizedAlpha[np.floor(e*nBins).astype(int)])) for e in eChainSet]
            plt.clf()
            plt.plot(eTrue, eMOS, 'ko', alpha=0.5)
            plt.plot([0,1], [0,1], 'k-', alpha=0.5)
            plt.xlim(0.,1.)
            plt.ylim(0.,1.)
            plt.xlabel('true eccentricity $e$')
            plt.ylabel('mean-of-sampling eccentricity $e$')
            plt.title(titlePrefix)
            plt.savefig(path + 'eCompareMOS' + suffix)
    return

if __name__ == '__main__':
    for nObj in [30, 300]:
        download_infer_and_plot('likemaxim', Turner=True, fit=False, nObj=nObj)
        download_infer_and_plot('crazylikemaxim', Gaussian=True, fit=False, nObj=nObj)
        download_infer_and_plot('normal', Turner=True, nObj=nObj)
        download_infer_and_plot('normal', Turner=True, nObj=nObj, beta=True)
        download_infer_and_plot('normal', Turner=True, nObj=nObj, nSampleMax=50)
        download_infer_and_plot('crazy', Gaussian=True, nObj=nObj)
        download_infer_and_plot('crazy', Gaussian=True, nObj=nObj, nSampleMax=50)
