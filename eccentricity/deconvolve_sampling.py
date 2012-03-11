# name:
#   deconvolve-sampling.py
# purpose:
#   Take a set of posterior samplings of x, return true x distribution.
# author:
#   David W. Hogg (NYU)
# license:
#   Copyright 2010 Hogg & Myers; all rights reserved (for now).

import numpy as np
from scipy.special import gammaln
# this rc block must be before the matplotlib import?
from matplotlib import rc
rc('font',**{'family':'serif','serif':'Computer Modern Roman','size':18})
rc('text', usetex=True)
# now import matplotlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt

# input: sampling from one object, step-function model parameters
# output: likelihood of the model for that one object
# assumes x ranges from 0 to 1
def oneObjectLikelihoodStep(xList, alpha):
    nStep = len(alpha)
    return np.mean(alpha[np.floor(xList*nStep).astype('int')])

def lnLikelihoodStep(xListSet, alpha):
    likeList = [oneObjectLikelihoodStep(xList, alpha) for xList in xListSet]
    return np.sum(np.log(likeList))

def rawBeta(x, alpha):
    return np.exp( gammaln(alpha[0]+alpha[1])
                   - gammaln(alpha[0]) - gammaln(alpha[1])
                   + np.log(x) * (alpha[0]-1.)
                   + np.log(1.-x) * (alpha[1]-1.) )

def lnLikelihoodBeta(xListSet, alpha):
    likeList = [np.mean(rawBeta(xList, alpha)) for xList in xListSet]
    return np.sum(np.log(likeList))

# this implements a Gaussian-process-like prior
# - epsilon should be a passed-in parameter
# - technically improper
def lnPriorStep(alpha):
    epsilon = 2.0
    lnA = np.log(alpha)
    return -0.5 * epsilon * np.sum((lnA[:-1]-lnA[1:])**2)

# improper; hack
def lnPriorBeta(alpha):
    if alpha[0] < 0.:
        return -1e9
    if alpha[1] < 0.:
        return -1e9
    return 0.

def lnPosterior(xListSet, alpha, beta=False):
    if beta:
        return lnPriorBeta(alpha) + lnLikelihoodBeta(xListSet, alpha)
    return lnPriorStep(alpha) + lnLikelihoodStep(xListSet, alpha)

# input: set of amplitudes alpha of a histogram model
# output: renormalized model
# assumes x ranges from 0 to 1.
def normalizeStep(alpha):
    dx = 1. / len(alpha)
    norm = np.sum(alpha * dx)
    return alpha / norm

def proposeNewStep(alpha, dLnAlpha):
    nStep = len(alpha)
    delta = dLnAlpha * np.random.normal(0.,1.,nStep)
    newalpha = np.exp(np.log(alpha) + delta)
    return normalizeStep(newalpha)

def proposeNewBeta(alpha, dalpha):
    return alpha + dalpha * np.random.normal(0.,1.,2)

# input: current amplitudes alpha, proposal width
# output: new amplitudes alpha
def proposeNew(alpha, proposalWidth, beta=False):
    if beta:
        return proposeNewBeta(alpha, proposalWidth)
    return proposeNewStep(alpha, proposalWidth)

# input: data, current parameters and log probability, step size, temperature
# output: next parameters and log probability
def MetropolisHastingsStep(xListSet, alpha, lnProb, proposalWidth,
                           temperature=1., beta=False):
    newAlpha = proposeNew(alpha, proposalWidth, beta=beta)
    newLnProb = lnPosterior(xListSet, newAlpha, beta=beta)
    accept = 1
    if np.log(np.random.uniform()) > ((newLnProb - lnProb) / temperature):
        newAlpha = alpha.copy()
        newLnProb = lnProb.copy()
        accept = 0
    return (newAlpha, newLnProb, accept)

def sampleAlpha(xListSet, nLink=10000, thinstep=1, proposalWidth=None,
                beta=False):
    if beta:
        alpha = np.zeros(2) + 1.
        if proposalWidth is None:
            proposalWidth = 0.15
    else:
        alpha = np.zeros(20) + 1.
        alpha = normalizeStep(alpha)
        if proposalWidth is None:
            proposalWidth = 0.15
    lnP = lnPosterior(xListSet, alpha, beta=beta)
    lnPBest = lnP
    alphaBest = alpha
    for i in range(2): # run twice, first time is a burn-in
        naccept = 0.
        alphaChain = []
        for link in range(nLink):
            (alpha, lnP,
             accept) = MetropolisHastingsStep(xListSet, alpha, lnP,
                                              proposalWidth, beta=beta)
            naccept += accept
            if lnP > lnPBest:
                lnPBest = lnP
                alphaBest = alpha
            if beta:
                print link, alpha, lnP, naccept/(link+1)
            else:
                print link, lnP, naccept/(link+1)
            if i:
                if ((link+1) % thinstep) == 0:
                    alphaChain.append((lnP, alpha))
    return (alphaBest, alphaChain)

def plotStep(amp, dots=True, alpha=0.2, edgecolor='k'):
    nStep = len(amp)
    xStep = (np.arange(nStep) + 0.5) / float(nStep)
    xplot = np.zeros([2 * nStep + 2, 2])
    # make staircase function
    xplot[0:2*nStep:2,0] = xStep - 0.5 / float(nStep)
    xplot[1:2*nStep:2,0] = xStep + 0.5 / float(nStep)
    xplot[0:2*nStep:2,1] = amp
    xplot[1:2*nStep:2,1] = amp
    # close the polygon
    xplot[-2,0] = 1.
    xplot[-1,0] = 0.
    xplot[-2,1] = 0.
    xplot[-1,1] = 0.
    # plot polygon
    ax = plt.gca()
    p1 = plt.Polygon(xplot, ec='none', fc='k', alpha=alpha)
    p2 = plt.Polygon(xplot, ec=edgecolor, fc='none', alpha=1., lw=0.5)
    ax.add_artist(p1)
    ax.add_artist(p2)
    # plot dots
    if dots:
        plt.plot(xStep, amp, 'ko', alpha=1.)
    plt.xlim(0.,1.)
    plt.ylim(0.,1.05*max(amp))

def plotBeta(params, alpha=1.):
    dx = 0.0001
    xplot = np.arange(0.5*dx,1.,dx)
    yplot = rawBeta(xplot, params)
    plt.plot(xplot,yplot,'k-', alpha=alpha)
    plt.xlim(0.,1.)
    plt.ylim(0.,1.05*max(yplot))

# functional test
if __name__ == '__main__':
    plt.clf()
    plotBeta([-3, 1.69339392])
    plt.xlabel('$x$')
    plt.savefig('xTest.png')
    os.exit(0)
    nTrue = 300
    xTrue = np.random.beta(1.5,5.0,size=nTrue)
    plt.clf()
    plt.hist(xTrue,20)
    plt.xlim(0.0,1.)
    plt.xlabel('true $x$')
    plt.savefig('xTrue.png')
    nSample = 1000
    xListSet = []
    for xT in xTrue:
        uncertainty = np.random.uniform(0.0,0.1)
        xList = xT + uncertainty * np.random.normal(0.,1.,nSample)
        bad = np.where(np.logical_or(xList < 0.,xList > 1.))
        nBad = len(bad[0])
        while nBad > 0:
            xList[bad] = xT + uncertainty * np.random.normal(0.,1.,nBad)
            bad = np.where(np.logical_or(xList < 0.,xList > 1.))
            nBad = len(bad[0])
        xListSet.append(xList)
    print len(xListSet)
    print len(xListSet[0])
    for i in range(10):
        plt.clf()
        plt.hist(xListSet[i],20)
        plt.xlim(0.0,1.0)
        plt.xlabel('observed $x$')
        plt.savefig('xObs%d.png' % i)
    plt.clf()
    (alphaBest, foo) = sampleAlpha(xListSet, nLink=1000, beta=True)
    plotBeta(alphaBest)
    plt.xlabel('inferred $x$')
    plt.savefig('xInf.png')
    (alphaBest, foo) = sampleAlpha(xListSet, nLink=1000)
    plotStep(alphaBest)
