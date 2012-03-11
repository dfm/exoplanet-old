# name:
#   monotonic_function.py
# purpose:
#   Fit a monotonically decreasing function to an iid data sample.
# author:
#   David W. Hogg (NYU)
# license:
#   Copyright 2010 Hogg; all rights reserved (for now).

import numpy as np
# this rc block must be before the matplotlib import?
from matplotlib import rc
rc('font',**{'family':'serif','serif':'Computer Modern Roman','size':18})
rc('text', usetex=True)
# now import matplotlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import markovpy as mcmc

BAD_VALUE = -1.0

# input: sampling, bin heights
# output: likelihood of the model
# assumes x ranges from 0 to 1
def lnLikelihood(pars, x):
    if len(x) < 1:
        return 0.
    alpha = unpack(pars)
    nStep = len(alpha)
    return np.sum(np.log(alpha[np.floor(x*nStep).astype('int')]))

# enforce non-negativity and monotonicity here
def unpack(pars):
    alpha = np.cumsum(np.exp(pars[::-1]))[::-1]
    return alpha * len(alpha) / np.sum(alpha)

def lnPrior(pars):
    epsilon = 1.0
    return -1. * (epsilon / len(pars)) * np.sum(pars**2)

def lnPosterior(pars, x):
    return lnPrior(pars) + lnLikelihood(pars, x)

def sampleAlpha(x):
    nStep = 30
    pars = -2. * np.zeros(nStep)
    if len(x) > 0:
        imax = np.ceil(np.max(x) * nStep).astype('int')
        pars[imax:] = pars[imax:] - 2.
    p0 = []
    for i in range(len(pars)):
        p0.append([pars[i]-0.001,pars[i]+0.001])
    parsChain,lnPChain,afrac = mcmc.mcfit(lnPosterior,p0,burnin=0,N=3000,
                                           args=(x,),outfile='foo.mcmc')
    return parsChain

def plotStep(amp, dots=True, fill=True, fillalpha=0.2, edgealpha=1.0, lw=0.5):
    nStep = len(amp)
    xStep = (np.arange(nStep) + 0.5) / float(nStep)
    xplot = np.zeros([2 * nStep + 2, 2])
    # make staircase function
    xplot[0:2*nStep:2,0] = xStep - 0.5 / float(nStep)
    xplot[1:2*nStep:2,0] = xStep + 0.5 / float(nStep)
    xplot[0:2*nStep:2,1] = amp
    xplot[1:2*nStep:2,1] = amp
    plt.plot(xplot[:-2,0], xplot[:-2,1], 'k-', alpha=edgealpha, lw=lw)
    # close the polygon
    xplot[-2,0] = 1.
    xplot[-1,0] = 0.
    xplot[-2,1] = 0.
    xplot[-1,1] = 0.
    # plot polygon
    ax = plt.gca()
    if fill:
        p1 = plt.Polygon(xplot, ec='none', fc='k', alpha=fillalpha)
        ax.add_artist(p1)
    # plot dots
    if dots:
        plt.plot(xStep, amp, 'ko', alpha=edgealpha)

def main(x):
    np.random.seed(42)
    pChain = sampleAlpha(x)
    (nLink,nPar) = pChain.shape
    nStep = nPar
    alphaChain = np.zeros((nLink/2,nStep))
    for i in range(nLink / 2):
        alphaChain[i,:] = unpack(pChain[nLink/2+i,:])
    plotLink = np.random.randint(nLink/2, size=(8))
    print plotLink
    plt.clf()
    plt.plot(x, 0.1+np.zeros(len(x)), 'ko', alpha=0.5)
    for I in plotLink:
        plotStep(alphaChain[I,:], dots=False, fill=False, edgealpha=0.25, lw=2.)
    meanAlpha = np.mean(alphaChain, axis=0)
    plotStep(meanAlpha, dots=False, edgealpha=1., lw=0.5)
    plt.xlabel('x')
    plt.xlim(0.,1.)
    plt.ylim(0.,6.)

# functional test
if __name__ == '__main__':
    x = np.array([])
    main(x)
    plt.savefig('monotonic_function_0.png')
    nData = 20
    x = 0.2 * np.random.uniform(size=(nData))
    main(x)
    plt.savefig('monotonic_function_1.png')
    x[0:2] = x[0:2]+0.7
    main(x)
    plt.savefig('monotonic_function_2.png')
