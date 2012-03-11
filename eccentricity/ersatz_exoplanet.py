# ersatz_exoplanet.py
#
#   class radial_velocity
#   Generate ersatz V(t) radial velocity data for an exoplanet of given input orbital parameters
#   Make a likelihood function to test orbital parameters against that planet
#
#   class ersatz_planets
#   Generate an ensemble of ersatz planets from sensible distributions of orbital parameters
#
#   class mcmc_planets
#   Generate a Markov Chain of lenght Nlinks for a set of N ersatz planets over N ersatz times
#
# intellectual property:
#   Copyright 2010 David W. Hogg and Adam D. Myers.
#   NYU and UIUC
#   All rights reserved.     Watch this space
#   for an open-source license to be inserted later.
#
# COMMENTS:
#
# v1.0 ADM:
#
# we start the chain at close to truth
# across all planets
#
# 1 million chain links will take about 80 hours for 400 planets!
#
# mcmc_planets can be tested by sending zero observations (ntimes=0) and
# checking the sampling plots are as they should be
#
# 
# KNOWN BUGS:
#
# TO DO LIST:
#     Multiple planets
#     Maybe discard arrays (multiple planet functionality over the top for this module)
#     Make order of sampling of individual parameters random not uniform
#     Make plot functionality separate
#     hardcode binning of histograms to 0 left edge, over 20 bins


from numpy import *
from glob import glob
import cPickle
import matplotlib.pyplot as plt

from time import time

# make a Markov Chain for ersatz planets
# based on the product of the likelihood
# across all planets
# Nplanets: number of planets
# Nlinks: number of links in Chain
# Ntimes: number of times sampled for each planet
# verbose: if this is 1 write out info to screen
# if god=1 then add a "godplanets" instance to plot truth
# if passplanets==1 then use the passed planet rather than generating
# passedplanets is the passed planet
# a new ersatz one
# in general this slows things down for large numbers of planets
# so should only be used to plot truth
# parameters are sampled one-at-a-time
#
# e.g.
# from ersatz_exoplanet import mcmc_planets
# 
# a = mcmc_planets(2,10) #create a chain of 10 for a distribution of 2 planets
#make sure numpy is imported in your shell for slicing
# print a.chain[0] # link 0 in the chain for all planets
# print a.chain[9] # link 9 in the chain for all planets
# print a.chain[10] # fails because there are only 10 links in the chain
# print a.likechain # the values of the likelihood function
#
# print a.chain[:,0] #all parameters in the chain for planet 0
# print a.chain[:,1] #all parameters in the chain for planet 1
# print a.chain[:,2] #fails because there is no third planet
#
# print a.chain[:,0][:,0] #parameter 0 (the period) in the chain for planet 0
# print a.chain[:,1][:,0] #parameter 0 (the period) in the chain for planet 1
# print a.planets[0].T #the true period of planet 0
# print a.planets[1].T #the true period of planet 1
# print a.chain[:,0][:,3] #parameter 3 (e*cos(pomega)) in the chain for planet 0
# print a.chain[:,1][:,3] #parameter 3 (e*cos(pomega)) in the chain for planet 1
# print a.planets[0].e*cos(a.planets[0].pomega) #the true e*cos(pomega) for planet 0
#
# what would the highest likelihood be for an individual planet?
# p = a.planets[0] #p is the instance of planet 0
# print p.lnlikelihood(p.K, p.T, p.phi, p.e, p.pomega, 0., 0.)
# p = a.planets[1] #p is the instance of planet 1
# print p.lnlikelihood(p.K, p.T, p.phi, p.e, p.pomega, 0., 0.)
#
# plot the progression of the Ensemble Likelihood with each link and
# make a hardcopy
# a.plot_like(hardcopy=1)
# 

class mcmc_planets:
    def __init__(self,Nlinks,Nplanets=0,Ntimes=0,passedplanets=0,verbose=True,god=True,test=False,eCrazy=False,likemaxim=False,passplanets=True,anneal=False):

        start = time()

        self.Nlinks = Nlinks

        #set up N planets
        if passplanets:
            fp = passedplanets
            self.Nplanets = len(fp.planets)
        else:
            fp = ersatz_planets(Nplanets,Ntimes,god=god,test=test,eCrazy=eCrazy)
            self.Nplanets = Nplanets

        self.planets = fp.planets

        if god:
            self.godplanets = fp.godplanets

        plan = self.planets[0]
        print "Truth for this planet"
        print 'T:', plan.T, 'K:', plan.K , 'phi:', plan.phi, 'pomega:', plan.pomega, 'e:', plan.e, 'V0:', 0., 'S:', 0., 'And observed',plan.ntimes, 'times'
    
        #note we're working the chain in the following parameters
        #T
        #Kcos(phi+pomega)
        #Ksin(phi+pomega)
        #ecos(pomega)
        #esin(pomega)
        #Vo
        #S
        #we'll cheat and start things at truth (or very close)

        self.Tstart = log(fp.Tdist)
        self.Kcosstart = fp.Kdist*cos(fp.phidist+fp.pomegadist)
        self.Ksinstart = fp.Kdist*sin(fp.phidist+fp.pomegadist)
        self.ecosstart = fp.edist*cos(fp.pomegadist)
        self.esinstart = fp.edist*sin(fp.pomegadist)
        self.Sstart = array(self.Nplanets*[0.])
        self.V0start = array(self.Nplanets*[0.])

        #the initial widths of the step functions to consider for each parameter

        self.stepwidths = zeros([self.Nplanets,7]) + [0.01, 0.01, 0.01, 0.1, 0.1, 0.01, 0.01]

        #iterate to set the step widths of the Gaussian sampler for each planet
        
        #dummy conditional
        toolo = [[1,1],[1,1]]
        toohi = [[1,1],[1,1]]
        count = 0
        while len(toolo[0]) + len(toohi[0]) > 0:
            if count == 30:
                print 'kicking'
                #give it a weird push
                toolo = where(accfrac < 0.3)
                self.stepwidths[toolo[0],toolo[1]] *= 1.3
                toohi = where(accfrac > 0.6)
                self.stepwidths[toohi[0],toohi[1]] /= 1.3
            if count > 60:
                raise StopIteration
            count+=1
            print count
            accfrac = self.setstepwidths(likemaxim=likemaxim)
            if verbose:
                print '........................................................'
                print 'Setting step widths...iteration',count
                print 'Step widths', str(around(self.stepwidths,5))
                print 'Accepted fractions', str(around(accfrac,3))
            waytoolo = where(accfrac < 0.05)
            self.stepwidths[waytoolo[0],waytoolo[1]] /= 10.
            toolo = where(accfrac < 0.3)
            self.stepwidths[toolo[0],toolo[1]] /= 2.
            toohi = where(accfrac > 0.6)
            self.stepwidths[toohi[0],toohi[1]] *= 1.5
            waytoohi = where(accfrac > 0.95)
            self.stepwidths[waytoohi[0],waytoohi[1]] *= (20/1.5)

        #self.stepwidths = zeros([self.Nplanets,7]) + [1., 0.003, 0.003, 0.2, 0.2, 0.003, 0.003]

        #generate array of starting parameters
        self.oldpars = column_stack([self.Tstart,self.Kcosstart,self.Ksinstart,self.ecosstart,self.esinstart,self.Sstart,self.V0start])

        #just for now call this newpars purely to evaluate the likelihood
        self.newpars = self.oldpars
        (oldlnp, oldlike) = self.posterior_evaluate(likemaxim=likemaxim)

        chain = []
        likechain = []
        lnpchain = []

        #to keep track of which parameter we're working with
        parcount = 0
        #counters for how many successful steps the sampler takes
        self.accepted = zeros([self.Nplanets,7])
        self.rejected = zeros([self.Nplanets,7])

        #set the seed of the random number generator so the chain is reproduceable
        random.seed(1413)

        temperature = 1.0
        for i in range(self.Nlinks):
            if anneal:
                temperature = float(self.Nlinks - i) / float(self.Nlinks)
            #get new parameters. Note that likelihood_step always works on self.oldpars
            self.newpars = self.likelihood_step(parcount)

            #and their posteriors. Note that posterior_evaluate always works on self.newpars
            (newlnp, newlike) = self.posterior_evaluate(likemaxim=likemaxim)
            #probabilities are log() so we take difference not ratio and ln the random number

            #remember that likelihoods are all arrays across N planets
            lnpratio = (newlnp-oldlnp)/temperature
            randnum = log(random.uniform(size=self.Nplanets))

            # M-H step
            wreject = where(lnpratio < randnum)
            waccept = where(logical_not(lnpratio < randnum))
            self.newpars[wreject] = self.oldpars[wreject]
            newlnp[wreject] = oldlnp[wreject]
            newlike[wreject] = oldlike[wreject]
            self.rejected[wreject,parcount] += 1
            self.accepted[waccept,parcount] += 1
            #append the newpars to the chain. Also append eccentricity
            chain.append(column_stack([self.newpars]))
            #keep the likelihoods for debugging
            likechain.append(newlike)
            lnpchain.append(newlnp)

            self.oldpars = self.newpars
            oldlnp = newlnp
            oldlike = newlike

            parcount +=1
            if parcount == 7:
                parcount = 0

            if verbose:
                if (20.*i/Nlinks)-(20*i/Nlinks) < 1e-10:
                    print "mcmc_planets", 100.*i/Nlinks,"percent complete"
                    print "t....", time()-start,"secs"


        self.chain = array(chain)
        self.likechain = array(likechain)
        self.lnpchain = array(lnpchain)

        if verbose: print "Generated", Nlinks, "links for", Nplanets, "planets"

    #Determine reasonable step widths for each planet to set the
    #size of the Gaussian sampler for each parameter
    def setstepwidths(self,likemaxim=False):

        #generate array of starting parameters
        self.oldpars = column_stack([self.Tstart,self.Kcosstart,self.Ksinstart,self.ecosstart,self.esinstart,self.Sstart,self.V0start])
        #just for now call this newpars purely to evaluate the likelihood
        self.newpars = self.oldpars
        (oldlnp, oldlike) = self.posterior_evaluate(likemaxim=likemaxim)

        #to keep track of which parameter we're currently working on
        parcount = 0
       #counters to record how many successful steps the sampler takes
        accepted = zeros([self.Nplanets,7])
        rejected = zeros([self.Nplanets,7])

        #set the seed of the random number generator so the chain is reproduceable
        random.seed(1413)

        for i in range(10000):
            #get new parameters. Note that likelihood_step always works on self.oldpars
            self.newpars = self.likelihood_step(parcount)

            #and their likelihood. Note that likelihood_evaluate always works on self.newpars
            (newlnp, newlike) = self.posterior_evaluate(likemaxim=likemaxim)
            #likelihoods are lnL so we take difference not ratio and ln the random number

            lnpratio = newlnp-oldlnp
            randnum = log(random.uniform(size=self.Nplanets))

            #don't take a step if the newlike is less than the oldlike
            #unless their ratio is less than a random number
            #wreject = where((logical_and((newlike < oldlike),(likeratio < randnum))))
            #waccept = where(logical_not(logical_and((newlike < oldlike),(likeratio < randnum))))
            wreject = where(lnpratio < randnum)
            waccept = where(logical_not(lnpratio < randnum))

            self.newpars[wreject] = self.oldpars[wreject]
            newlnp[wreject] = oldlnp[wreject]
            rejected[wreject,parcount] += 1
            #we did take a step where newlnp was greater than oldlnp etc.
            accepted[waccept,parcount] += 1

            self.oldpars = self.newpars
            oldlnp = newlnp
                        
            parcount +=1
            if parcount == 7:
                parcount = 0

        accfrac = 1.*accepted/(accepted+rejected)
        return accfrac


    #Plot the progression of the ensemble likelihood with each iteration
    #if hardcopy=1 a png is made
    def plot_like(self,hardcopy=0):

        #to subsample if the plots would get busy
        #me = self.Nlinks/10000

        #make an ensemble likelihood across all planets
        likeensemble = sum(self.likechain,axis=1)
        plt.clf()
        stringy = 'Likelihood Progression ('
        stringy += str(len(self.chain[0,:]))+' planets)'
        plt.suptitle(stringy)
        #The maximum likelihood value
        plt.plot([max(likeensemble)]*self.Nlinks,'r--')
        #The likelihood sampling
        plt.plot(likeensemble,'.',ms=1.2)
        plt.axis([0,len(likeensemble),median(likeensemble)-30,median(likeensemble)+30])
        plt.xlabel('Link')
        plt.ylabel('Ensemble ln(Likelihood)')
        if hardcopy:
            plt.savefig('plot_like.png')
        else:
            plt.show()

    # plot the parameter sampling for any planet and any parameter
    # two plots are created..parameter as a function of chain link
    # and a histogram of the parameter sampling
    # parameter order is as in pardict
    # if hardcopy=1 then make a png
    def plot_sampling(self,parameter,hardcopy=0):

        #to subsample if the plots would get busy
        #me = self.Nlinks/10000

        plan = self.planets[0]
        truth = [log(plan.T), plan.K*cos(plan.phi+plan.pomega),plan.K*sin(plan.phi+plan.pomega),plan.e*cos(plan.pomega),plan.e*sin(plan.pomega),0.,0.,plan.e,plan.K,plan.T,plan.pomega,plan.phi,log(plan.K)]
        pardict = {0:'lnT',1:'Kcos',2:'Ksin',3:'ecos',4:'esin',5:'S',6:'V0',7:'e',8:'K',9:'T',10:'pomega',11:'phi',12:'lnK'}

        ecc = array(sqrt(self.chain[:,0,3]**2+self.chain[:,0,4]**2))
        K =  array(sqrt(self.chain[:,0,1]**2+self.chain[:,0,2]**2))
        T = exp(self.chain[:,0,0])
        pomega = arctan2(self.chain[:,0,4],self.chain[:,0,3])
        phi = arctan2(self.chain[:,0,2],self.chain[:,0,1]) - pomega
        phi = arctan2(sin(phi), cos(phi))

        plotchain = column_stack([self.chain[:,0],ecc,K,T,pomega,phi,log(K)])
        likechain = self.likechain*1.
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

        parst = str(pardict[parameter])

        plt.clf()

        #sampling
        plt.plot(plotchain[:,parameter],'.',ms=1.2)
        #maximum likelihood
        plt.axhline(plotchain[ML,parameter],color='k')
        #truth
        plt.axhline([truth[parameter]],color='r')
        plt.xlabel('Link')
        plt.ylabel(parst)
        if hardcopy:
            plt.savefig(parst+'_plot_sampling.png')
        else:
            plt.show()

        plt.clf()
        plt.hist(plotchain[:,parameter],20)
        #maximum likelihood
        plt.axvline(plotchain[ML,parameter],color='k')
        #truth
        plt.axvline([truth[parameter]],color='r')
        plt.xlabel(parst)
        plt.ylabel("N("+parst+")")
        if hardcopy:
            plt.savefig(parst+'_histo_sampling.png')
        else:
            plt.show()

    #Step the parameters to a new place to sample the likelihood function
    #Like likelihood_step but steps in the parameters one-at-a-time
    #parcount is the parameter to step in (0 for Tstep, 1 for Kcosstep etc.)
    def likelihood_step(self,parcount):

        stepwidth = self.stepwidths[:,parcount]

        #lambda function to represent the step to take
        step = lambda x: x + random.normal(0,stepwidth,self.Nplanets)

        #the *1. insists this is a new array, not a repeat instance of oldpars
        newpars = self.oldpars*1.

        #move the correct parameter on according to its step width
        newpars[:,parcount] = step(newpars[:,parcount])

        return newpars

    #Take a list of
    #[T, Kcos(phi+pomega),Ksin(phi+pomega),ecos(pomega),esin(pomega),S,V0]
    #convert to K,T,phi,e,pomega,S,V0 and call likelihood function using these parameters
    #then to sum the likelihood across all of the planets
    #!!evaluates based on whatever is currently called newpars!!
    #returns ln(Likelihood) totalled across the ensemble of all passed planets
    def posterior_evaluate(self,likemaxim=False):

        parsarr = self.newpars

        lnP = zeros([self.Nplanets])
        lnL = zeros([self.Nplanets])

        T = exp(parsarr[:,0])
        K = sqrt(parsarr[:,1]**2.+parsarr[:,2]**2.)
        e = sqrt(parsarr[:,3]**2.+parsarr[:,4]**2.)
        pomega = arctan2(parsarr[:,4],parsarr[:,3])
        phi = arctan2(parsarr[:,2],parsarr[:,1]) - pomega
        S = parsarr[:,5]
        V0 = parsarr[:,6]

        for i in range(self.Nplanets):
            (foo, bar) = self.planets[i].lnposterior(K[i],T[i],phi[i],e[i],pomega[i],S[i],V0[i],likemaxim=likemaxim)
            lnP[i] = foo
            lnL[i] = bar

        return (lnP, lnL)

# draw N ersatz planet observations from a sensible set
# of distribution functions of the orbital parameters
# N is the input number of planets and n is the time sampling of each planet
# n is currently uniform
# try:
# import ersatz_exoplanet
# a = ersatz_exoplanet.ersatz_planets(10000)
# a.planets[1].plot_Vt()
# a.planets[199].plot_Vt()
# print a.planets[1].T #The period of planet 1 (the second planet)
# print a.planets[1].e #The eccentricity of planet 1
# a.plot_edist()
# a.plot_likedist()
# a = ersatz_exoplanet.ersatz_planets(100,god=True)
# a.godplanets[1].plot_Vt(god=True)
# if test = 1 then we make everything a reasonably "simple" well-behaved planet

class ersatz_planets:
    def __init__(self,N,ntimes,god=True,test=False,eCrazy=False,mflag=0):

        godplan = []
        plan = []
        Mp = generate_m_planets(N, mflag=mflag)
        Ms = generate_m_stars(N)
        T = generate_periods(N)
        phi = generate_phases(N)
        e = generate_eccentricities(N, eCrazy=eCrazy)
        i = generate_inclinations(N)
        pomega = generate_pomegas(N)
        self.ntimes = zeros(N)+ntimes
        self.eCrazy = eCrazy

        edist = []
        Kdist = []
        Mpdist = []
        Tdist = []
        idist = []
        phidist = []
        pomegadist = []
        likedist = []

        testy = test

        for j in range(N):
            if god:
                planet = radial_velocity(Mp[j],Ms[j],i[j],T[j],phi[j],e[j],pomega[j],self.ntimes[j],god=True,test=testy)
                godplan.append(planet)
            planet = radial_velocity(Mp[j],Ms[j],i[j],T[j],phi[j],e[j],pomega[j],self.ntimes[j],test=testy)
            plan.append(planet)
            edist.append(planet.e)
            Kdist.append(planet.K)
            Mpdist.append(planet.Mp)
            Tdist.append(planet.T)
            idist.append(planet.i)
            phidist.append(planet.phi)
            pomegadist.append(planet.pomega)
            (lnP, lnL) = planet.lnposterior(planet.K, planet.T, planet.phi, planet.e, planet.pomega, 0., 0.)
            likedist.append(lnL)

        self.planets = plan
        self.godplanets = godplan

        self.Kdist = array(Kdist)
        self.edist = array(edist)
        self.Mpdist = array(Mpdist)
        self.Tdist = array(Tdist)
        self.idist = array(idist)
        self.phidist = array(phidist)
        self.pomegadist = array(pomegadist)
        self.likedist = array(likedist)

    def plot_edist(self,hardcopy=0):
        plt.clf()
        plt.xlabel("Eccentricity, $e$")
        plt.ylabel("N($e$)")
        if self.eCrazy:
            n, bins, patches = plt.hist(self.edist,20,facecolor='0.75',edgecolor='0.75',normed=True)
            plt.axis([0,1,0,1.1*max(n)])
            mean = 0.3
            sdev = 0.05
            var = sdev*sdev
            gaussian = lambda x: (1./sqrt(2*pi*var))*exp(-((x-mean)**2./2./var))
            x = array(range(1001))*0.001
            plt.plot(x, gaussian(x),color='0.35')
            if hardcopy:
                plt.savefig('plot_edist_crazy.png')
            else:
                plt.show()
        else:
            n, bins, patches = plt.hist(self.edist,20,facecolor='0.5',edgecolor='0.5',normed=True)
            plt.axis([0,1,0,1.1*max(n)])
            efunc = lambda e: (1./((1.+e)**4.))-(e/(2.**4.))
            x = array(range(1001))*0.001
            #the 3.84 factor is the reciprocal of the integral under efunc (to normalize)
            plt.plot(x, 3.840043008481695*efunc(x),color='0.35')                
            if hardcopy:
                plt.savefig('plot_edist.png')       
            else:
                plt.show()

    def plot_Mpdist(self,hardcopy=0):
        plt.clf()
        plt.hist(log10(self.Mpdist),20,facecolor='0.5',edgecolor='0.5',normed=True)
        plt.xlabel("Planet Mass, log($m / {\mathrm{M_{jup}}}$)")
        plt.ylabel("N($m$)")
        #plot a uniform log distribution from -1 to 1
        x = (array(range(1001))*0.002)-1
        #0.5 here is because of normalization (integral under -1 to 1)
        plt.plot(x, 0.5+(0*x),color='0.35')                
        if hardcopy:
            plt.savefig('plot_Mpdist.png')
        else:
            plt.show()

    def plot_Tdist(self,hardcopy=0):
        plt.clf()
        n, bins, patches = plt.hist(log10(self.Tdist),20,facecolor='0.5',edgecolor='0.5',normed=True)
        plt.xlabel("Period, log($T$ / 1 d)")
        plt.ylabel("N($T$)")
        #plot a uniform log distribution from 0.3 to 3.3
        x = (array(range(1001))*0.003)+0.3
        #0.3333333 here is because of normalization (integral under -1 to 1)
        plt.plot(x, 0.3333333+(0*x),color='0.35')                
        plt.axis([0.3,3.3,0,1.1*max(n)])        
        if hardcopy:
            plt.savefig('plot_Tdist.png')
        else:
            plt.show()

    def plot_Kdist(self,hardcopy=0):
        plt.clf()
        plt.hist(log10(self.Kdist),20,facecolor='0.5',edgecolor='0.5',normed=True)
        plt.xlabel("Velocity Amplitude, log($\kappa / {\mathrm{m s}}^{-1})$")
        plt.ylabel("N($\kappa$)")
        if hardcopy:
            plt.savefig('plot_Kdist.png')
        else:
            plt.show()

    def plot_idist(self,hardcopy=0):
        plt.clf()
        plt.hist(self.idist*180./pi,20,facecolor='0.5',edgecolor='0.5',normed=True)
        plt.xlabel("Inclination, $i$ (degrees)")
        plt.ylabel("N($i$)")
        #plot a uniform log distribution in cos
        x = (array(range(1001))*0.001)*90
        #remember to normzlize from degrees to radians (1/180/pi)
        #sinx is same thing as cos(90-x)
        plt.plot(x, sin(x*pi/180.)*pi/180.,color='0.35')                 
        if hardcopy:
            plt.savefig('plot_idist.png')
        else:
            plt.show()

    def plot_likedist(self,hardcopy=0):
        plt.clf()
        plt.hist(self.likedist,20,facecolor='0.5',edgecolor='0.5',normed=True)
        plt.suptitle("ln Likelihood Distribution")
        plt.xlabel("ln(Likelihood) ($L$)")
        plt.ylabel("N($L$)")
        if hardcopy:
            plt.savefig('plot_likedist.png')
        else:
            plt.show()


# class to generate radial velocity data for a planet of
# given orbital parameters
# retains all of the inputted information for the planet
# generates a random time array to sample
# the orbit and random velocities and errors at those times
# generates a function to evaluate the likelihood of different
# orbital parameters for the planet
# Phase angles are in radians.
# Mp - Mass of planet (Jupiter masses)
# Ms - Mass of star (Jupiter masses)
# i - inclination (degrees)
# T - period (days)
# phi - phase
# e - eccentricity
# pomega  - eccentric longitude (radians)
# retrievable when the class is instantiated
# tday - time of observation (days)
# rvel - radial velocity (km/s)
# Try:
# import ersatz_exoplanet
# Mp,Ms,i,T,phi,e,pomega = 1.,1000.,87.6,100.,0.,0.05,0.
# a=ersatz_exoplanet.radial_velocity(Mp,Ms,i,T,phi,e,pomega)
# a.plot_Vt()
# a=ersatz_exoplanet.radial_velocity(Mp,Ms,i,T,phi,e,pomega,god=True)
# a.plot_Vt(god=True)
#if test is True we set the (hopefully simple) default test planet

class radial_velocity:
    def __init__(self,Mp,Ms,i,T,phi,e,pomega,ntimes,god=False,test=False):

        if test:
            self.Mp = 1.
            self.Ms = 1000.
            self.i = 87.6
            self.T = 100.
            self.phi = 0.
            self.e = 0.05
            self.pomega = 0.
            self.ntimes = 30
        else:
            self.Mp = float(Mp)
            self.Ms = float(Ms)
            self.i = float(i*pi/180.)
            self.T = float(T)
            self.phi = float(phi)
            self.e = float(e)
            self.pomega = float(pomega)
            self.ntimes = float(ntimes)

        if god:
            #generate time series
            self.timegod = arange(5000)/5.
            #generate intermediaries
            M = mean_anomaly_from_time(self.timegod,self.T,self.phi)
            E = eccentric_anomaly_from_mean_anomaly(M,self.e)
            f = true_anomaly_from_eccentric_anomaly(E,self.e)
            K = rad_vel_amp_from_mass_one_planet(self.Mp,self.Ms,self.i,self.e,self.T)

            #final radial velocity distribution
            self.Vmodelgod = rad_vel(K,f,self.pomega,self.e)


        #generate time series
        if test:
            self.time = array([1.,20.,43.,67.,81.,90.,112.,137.,160.,165.,185.,201.,233.,248.,271.,290.,301.,322.,340.,367.,385.,400.,411.,439.,455.,465.,475.,483.,490.,499])
        else:
            self.time = generate_times(ntimes)
        #generate intermediaries
        self.M = mean_anomaly_from_time(self.time,self.T,self.phi)
        self.E = eccentric_anomaly_from_mean_anomaly(self.M,self.e)
        self.f = true_anomaly_from_eccentric_anomaly(self.E,self.e)
        self.K = rad_vel_amp_from_mass_one_planet(self.Mp,self.Ms,self.i,self.e,self.T)
        self.M2 = mean_anomaly_from_eccentric_anomaly(self.E,self.e)

        #final radial velocity distribution
        self.Vmodel = rad_vel(self.K,self.f,self.pomega,self.e)

        #generate errors on velocity
        if test:
            self.Verr = array([0.005])
        else:
            self.Verr = generate_verrs(ntimes)

        #generate a series of offsets based
        #on Gaussian sampling of the errors
        if test:
            self.Voffset = array(ntimes*[0.])
        else:
            self.Voffset = random.normal(0,1,len(self.Verr))
        #then normalize each offset by its errors
        self.Voffset*=self.Verr

        #Now make the data
        self.Vdata = self.Vmodel+self.Voffset

    def plot_Vt(self,god=False,hardcopy=0):
        stringy = "Mp="+str(round(self.Mp,1))+"Mjup, "
        stringy += "Ms="+str(round(self.Ms,1))+"Mjup, "
        stringy += "i="+str(round(self.i*180./pi,1))+"deg, "
        stringy += "T="+str(round(self.T,1))+"days, "
        stringy += "phi="+str(round(self.phi,2))+", "
        stringy += "e="+str(round(self.e,2))+", "
        stringy += "pomega = "+str(round(self.pomega,2))

        if god:
            plt.clf()
            plt.suptitle(stringy)
            if len(self.time) > 0:
                plt.errorbar(self.time,self.Vdata,self.Verr,fmt='bx')
            plt.xlabel("Time (days)")
            plt.ylabel("Radial Velocity (km/s)")
            plt.plot(self.timegod,self.Vmodelgod)
            if hardcopy:
                plt.savefig('god_plot_Vt.png')
            else:
                plt.show()

        else:
            plt.clf()
            plt.suptitle(stringy)
            if len(self.time) > 0:
                plt.errorbar(self.time,self.Vdata,self.Verr,fmt='bx')
            plt.xlabel("Time (days)")
            plt.ylabel("Radial Velocity (km/s)")
            if hardcopy:
                plt.savefig('plot_Vt.png')
            else:
                plt.show()


# for a given planet instantiated in the class
# construct a function that defines the likelihood
# of a set of parameters against that planet
# The things that are passed here are *not*
# instantiated in the class...they are dummy variables
# to call the function with at a later point
#  K       - radial velocity amplitude guessed by observational astronomer
#  T       - period (days) guessed by observational astronomer
#  phi     - phase (radians) guessed by observational astronomer
#  e       - eccentricity guessed by observational astronomer
#  pomega  - eccentric longitude (radians) guessed by observational astronomer
#  S       - Jitter parameter guessed by observational astronomer
#  V0      - Velocity offset guessed by observational astronomer
#
# NOTE THAT WE ALWAYS RETURN lnLn instead of L
# AS IT WILL ALWAYS BE MORE PRECISE
# TO SUM ln(L) THAN TO MULTIPLY THE LIKELIHOOD
# THUS THE USER SHOULD BE AWARE THAT THE
# OUTPUT QUANTITIES SHOULD BE SUMMED NOT MULTIPLIED
# NOTE ALSO THAT WE TAKE THE NUMBER OF DOF OF
# -2ln(L) FOR THE SAME REASON

    def lnprior(self,K,T,phi,e,pomega,S,V0,likemaxim=False):
        sigmaV0 = 10 # km / s
        Kmin = 1e-5 # km / s
        lnp = 0.0
        if not likemaxim:
            # convert from "natural" $p \propto e$ prior to $p \propto e^0$
            lnp += log(1./e)
            # convert from "natural" $p \propto K$ to $p \propto 1/K$
            lnp += log(1./K**2)
        # Gaussian prior in V0
        lnp += -0.5 * log(2 * pi * sigmaV0**2) - 0.5 * (V0 / sigmaV0)**2
        if e > 1 or e < 0:
            #print 'etrip', e
            lnp = -1e10
        if T < 0.1 or T > 3.6525e5: # days
            #print 'Ttrip'
            lnp = -1e10
        if K < Kmin or K > 10.0: # km / s
            #print 'Ktrip'
            lnp = -1e10
        if S**2 > 1.0:
            #print 'Strip'
            lnp = -1e10
        if lnp != lnp:
            print K,T,phi,e,pomega,S,V0
            raise IOError

        return lnp

    def lnlikelihood(self,K,T,phi,e,pomega,S,V0):
    #Dear God: what is truth for this planet?
    #Dear User: self.Vdata and self.Verr and self.time
    #are the velocity and error at a given time
        if len(self.time) < 1:
            return 0.

        #based on passed parameters, determine
        #output of radial velocity equation at same times as data (self.time)

        M = mean_anomaly_from_time(self.time,T,phi)
        E = eccentric_anomaly_from_mean_anomaly(M,e)
        f = true_anomaly_from_eccentric_anomaly(E,e)
        V = rad_vel(K,f,pomega,e)

        #now we have the radial velocity based on the sample parameters
        #generated at the times for the planet represented by this class
        #we can calculate the likelihood based on the planetary data generated
        #in this class. The generated data and errors are self.Vdata and self.Verr
        #note that these (and the V calculated above) are arrays at times self.time

        errjit = (self.Verr**2.+S**2)

        # The likelihood function
        minus2lnL = sum(log(errjit)) + sum((V0 + V - self.Vdata)**2./errjit)

        #this is the mean error on our sampled exoplanet distribution in
        #the function generate_verrs(). This informs the normalization
        #of the likelihood function. verrmax and verrmin are in km/s
        #        verrmax = 0.01
        #        verrmin = 0.001
        #        merr = sqrt((verrmax**2+verrmin**2.)/2.)

        #normalization constant
        #        Q = len(self.time)*(log(merr**2.)+1.)
        if minus2lnL != minus2lnL:
            print "fuckup",K,T,phi,e,pomega,S,V0
            minus2lnL = -1.0e10
        return -0.5*(minus2lnL)

    def lnposterior(self,K,T,phi,e,pomega,S,V0,likemaxim=False):
        lnp= self.lnprior(K,T,phi,e,pomega,S,V0,likemaxim=likemaxim)
        if lnp < -1.0e8:
            return (lnp, lnp)
        lnL = self.lnlikelihood(K,T,phi,e,pomega,S,V0)
        return (lnp+lnL, lnL)

# convert day of observation to mean anomaly
# assumed that time starts on day 0
# time = time of observation since t0 (days)
# T = period (days)
# phi = phase (radians)
def mean_anomaly_from_time(time,T,phi):
    return (2.*pi*time/T)+phi

# convert eccentric anomaly to mean anomaly
#  E       - eccentric anomaly (radians)
#  eccentricity is instantiated in class
#  e       - eccentricity
#  return  - mean anomaly (radians)
def mean_anomaly_from_eccentric_anomaly(E,e):
    return (E - e * sin(E))

# convert mean anomaly to eccentric anomaly
#  M       - [array of] mean anomaly (radians)
#  e       - eccentricity
# These are instantiated in the class
#  return  - eccentric anomaly (radians)
def eccentric_anomaly_from_mean_anomaly(M,e,verbose=False):
    E = M + e * sin(M)
    iteration = 0
    deltaM = 100.0
    mi = 1000.
    #should never hit this limit
    tol = 1e-15
    #maybe incorporate tolerance in the class if runtime grows?
    while (iteration < mi) and all(abs(deltaM) > tol):
        deltaM = (M - mean_anomaly_from_eccentric_anomaly(E,e))
        E = E + deltaM / (1. - e * cos(E))
        iteration += 1
        if verbose: print 'eccentric anomaly iterations:',iteration
    return E

# convert eccentric anomaly to true anomaly
#  E       - eccentric anomaly (radians)
#  e       - eccentricity
#  These are instantiated in class
#  return  - true anomaly (radians)
def true_anomaly_from_eccentric_anomaly(E,e):
    f = arccos((cos(E) - e) / (1. - e * cos(E)))
    f *= (sign(sin(f)) * sign(sin(E)))
    return f

# compute radial velocity amplitude (K)
# only strictly true in the one-planet sense
#  Mp      - the mass of the planet (jupiter masses)
#  Ms      - the mass of the star (jupiter masses)
#  T       - the orbital period (days)
#  i       - the inclination (radians)
#  e       - eccentricity
#  These are instantiated in class
#  return  - radial velocity amplitude in km/s
def rad_vel_amp_from_mass_one_planet(Mp,Ms,i,e,T):
    GinSIunits = 6.67428e-11
    secperday = 86400.
    Mjup = 1.8986e27
    # G in km/s cubed per Jupiter mass per year
    G = GinSIunits*1e-9*Mjup/secperday
    # K is a product of three terms
    K = (2*pi*G/T)**(1./3)
    K *= Mp*sin(i)/((Mp+Ms)**(2/3.))
    K *= 1./sqrt(1.-(e*e))
    return K

# compute radial velocity
#  K       - radial velocity amplitude
#  f       - true anomaly (radians)
#  e       - eccentricity
#  pomega  - eccentric longitude (radians)
#  These are instantiated in class
#  return  - radial velocity (same units as K)
def rad_vel(K,f,pomega,e):
    return K * (sin(f + pomega) + e * sin(pomega))

# currently a sampling of 0 to 30 random days over 1000 days
def generate_times(n):
    return random.uniform(0,1000,n)

# currently a sampling of 30 at 10 to 100 m/s in the square
def generate_verrs(n):
    squareverr = random.uniform(10,100,n)
    #remember to convert to km/s for main code 
    return sqrt(squareverr)/1000.

#generate a sample of planet masses
# from an unrealistically narrow distribution
def generate_crazy_m_planets(N):
    return 10.**(0.0 + 0.2*random.normal(0.,1.,N))

#generate a sample of planet masses
# currently a Schechter function
# note acceptance-rejection logic
def generate_schechter_m_planets(N):
    alpha = -1.125
    m0 = 5.0
    m = zeros(N)
    mmax = 100.0
    mmin = 0.01
    bad = (where(m < 1.0))[0]
    nbad = len(bad)
    while nbad > 0:
        m[bad] = random.uniform(mmin**(alpha+1.),mmax**(alpha+1.),nbad)**(1./(alpha+1.))
        bad = bad[(where(logical_or(random.uniform(0.,1.,nbad) > exp(-m[bad]/m0),m[bad] < mmin))[0])]
        nbad = len(bad)
    return m

#generate a sample of planet masses
# currently a sampling of N drawn in the log from 0.1 to 10Mjup
def generate_m_planets(N, mflag=0):
    if (mflag == 1):
        return generate_schechter_m_planets(N)
    if (mflag == 2):
        return generate_crazy_m_planets(N)
    return 10**random.uniform(-1.,1.,N)

#generate a sample of star masses
#currently a sampling of N all of 1 solar mass
def generate_m_stars(N):
    sol_to_jup_mass = 1047.56
    return random.uniform(1.*sol_to_jup_mass,1.*sol_to_jup_mass,N)

#generate a sample of orbital periods
#currently a sampling of N drawn in the log from 2 to 2000days
def generate_periods(N):
    logperiod = random.uniform(0.3,3.3,N)
    return 10**logperiod

def generate_phases(N):
    theta = random.uniform(0.,2.*pi,N)
    return arctan2(sin(theta),cos(theta))

def generate_pomegas(N):
    theta = random.uniform(0.,2.*pi,N)
    return arctan2(sin(theta),cos(theta))

#generate a sample of eccentricities
# from an unrealistically narrow distribution
# note logic to deal with e<0 and 1<e.
def generate_crazy_eccentricities(N):
    e = zeros(N)-0.1
    bad = where(logical_or((e < 0.),(e > 1.)))
    while len(bad[0]) > 0:
        e[bad] = 0.3 + 0.05*random.normal(0.,1.,N)
        bad = where(logical_or((e < 0.),(e > 1.)))
    return e

#generate a sample of eccentricities
# currently a sampling of N drawn by a bizarre derivative function
# function is from Shen \& Turner
# acceptance / rejection method used is INSANE
def generate_eccentricities(N, eCrazy=False):
    if eCrazy:
        return generate_crazy_eccentricities(N)

    #The bizarre derivative function
    efunc = lambda e: (1./((1.+e)**4.))-(e/(2.**4.))

    sampsize = N*2
    mi = 1000
    iteration = 0
    while (iteration < mi):
        #draw sampsize eccentricities in 0,1
        e = random.uniform(0.,1.,sampsize)
        #evaluate the derivative function
        Q = efunc(e)
        #draw sampsize test points in 0,1
        test = random.uniform(0.,1.,sampsize)
        #this is e where the efunc beats test, zero otherwise
        tested = (Q > test)*e
        #remove the zeroes
        tested = tested[where(tested != 0)]
        if len(tested) > N:
            return tested[0:N]
        else:
            #if we didn't generate enough random points sampsize was
            #too small
            sampsize*=2
            iteration+=1
    print '**********************************************'
    print 'Too many iterations in generate_eccentricities'
    print '**********************************************'
    raise IOError

#generate a sample of inclinations
#currently a sampling of N drawn in cos i from 0 to 1
def generate_inclinations(N):
    cosi = random.uniform(0.,1.,N)
    return arccos(cosi)*180./pi


#from a premade list of ersatz_planets, call MCMC chain and make webpages
#planetlist is a list of premade planets
#chainlinks is the length of the chain
#targdir is where to write the webpage
#if likemaxim is True use the likelihood-maximizing method rather than the posterior-sampling
#if plot = False don't make the plots to save time
#if clobber = True overwrite any existing directory structure
def do_everything(planetlist, chainlinks, targdir, likemaxim=False, plot=False, clobber=True, anneal=False):
    mainstart = time()

    if clobber:
        if os.system('/bin/rm -rvf %s' % targdir):
            print 'I hope you wanted to do that!'

    print 'will make diagnostics.html in directory',targdir, 'to open in any browser'

    #make and change to the target directory
    os.mkdir(targdir)
    os.chdir(targdir)

    today = str(date.today())
    mainfile = open('diagnostics.html','w')
    mainfile.write("""<!doctype html public "-//w3c//dtd html 4.0 transitional//en">\n<html>\n<h1>Created by ersatz_exoplanet.py on """+today)
    mainfile.write("""<body>\n\n""")
    mainfile.write("""<h2>Planets</h2>\n\n""")
    
    distfile = open('distributions.html','w')

    distfile.write("""<!doctype html public "-//w3c//dtd html 4.0 transitional//en">\n<html>\n<h1>Created by ersatz_exoplanet.py on """+today)
    distfile.write("""<body>\n\n""")
    distfile.write("""<h2>Distributions of the Initial Parameters over 5000 planets</h2>\n\n""")

    distfile.write("""<img SRC="plot_edist.png" height=380 width=570>\n""")
    distfile.write("""<img SRC="plot_Mpdist.png" height=380 width=570>\n""")
    distfile.write("""<br><br>\n""")
    distfile.write("""<img SRC="plot_Tdist.png" height=380 width=570>\n""")
    distfile.write("""<img SRC="plot_idist.png" height=380 width=570>\n""")
    distfile.write("""<br><br>\n""")
    distfile.write("""This is the initial likelihoods of the planets evaluated at truth<br><br>\n""")
    for nbsp in range(4):
        distfile.write('''&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\n''')
    distfile.write("""<img SRC="plot_likedist.png" height=570 width=770><br><br>\n""")
    
    distfile.write("""</body>\n</html>\n""")

    distfile.close()

    #the number of passed planets in the planet list
    totalplanets = len(planetlist)
    #loop counter through N planets
    i = 0
    #count of planets that fail to converge
    badcount = 0
    while i < totalplanets:
        print 'making planet', i,' t....', time()-mainstart,'secs'
        try:
            mainfile.write('''<a href = "planet%2.2d''' %i +'''/diagnostics.html">planet%2.2d''' %i +'''</a><br>\n''')
            outfile = 'planet.%2.2d.pickled' % i
            outf = open(outfile, 'w')

            #make mcmcchain

            a = mcmc_planets(chainlinks,passedplanets=planetlist[i],passplanets=True,likemaxim=likemaxim,anneal=anneal)
            plan = a.planets[0]
            godplan = a.godplanets[0]

            # thin the chain -- maybe?
            thinstep = ceil(chainlinks/100000)
            if thinstep > 1:
                print "thinning chain by a factor of %d" % thinstep
                print shape(a.chain)
                a.chain = a.chain[::thinstep,:,:]
                print shape(a.chain)
                print "thinning likechain by a factor of %d" % thinstep
                print shape(a.likechain)
                a.likechain = a.likechain[::thinstep,:]
                print shape(a.likechain)


            #write out chain and truth to file
            print 'writing file t....', time()-mainstart,'secs'
            truth = [log(plan.T), plan.K*cos(plan.phi+plan.pomega),plan.K*sin(plan.phi+plan.pomega),plan.e*cos(plan.pomega),plan.e*sin(plan.pomega),0.,0.]
            outdict = {'lnT':a.chain[:,:,0],'Kcos':a.chain[:,:,1],'Ksin':a.chain[:,:,2],'ecos':a.chain[:,:,3],'esin':a.chain[:,:,4],'S':a.chain[:,:,5],'V0':a.chain[:,:,6],'lnPosterior':a.lnpchain,'lnLikelihood':a.likechain,'truthlnT':truth[0],'truthKcos':truth[1],'truthKsin':truth[2],'truthecos':truth[3],'truthesin':truth[4],'truthS':truth[5],'truthV0':truth[6],'truthMp':plan.Mp,'truthMs':plan.Ms,'truthi':plan.i,'timegod':godplan.timegod,'Vmodelgod':godplan.Vmodelgod,'time':plan.time,'Vmodel':plan.Vmodel,'Vdata':plan.Vdata,'Verr':plan.Verr}
            cPickle.dump(outdict,outf)

            outf.close()

            #make a directory for each planet to store webpage
            pdir = 'planet%2.2d' % i

            os.mkdir(pdir)

            #write the page. This is filthy now, could generalize later
            file = open(pdir+'/diagnostics.html','w')

            file.write("""<!doctype html public "-//w3c//dtd html 4.0 transitional//en">\n<html>\n<h1>Created by ersatz_exoplanet.py on """+today+"""</h1>""")

            file.write("""<body>\n\n""")

            if i > 0:
                file.write('''<a href="../planet%2.2d''' %(i-1) + '''/diagnostics.html">previous planet</a><br>''')
            if i+1 < totalplanets:
                file.write('''<a href="../planet%2.2d''' %(i+1) + '''/diagnostics.html">next planet</a><br>''')

            file.write("""<h2>Planet's Radial Velocity Data and Truth</h2>\n\n""")

            file.write('''<img SRC="plot_Vt.png" height=380 width=570>\n''')
            file.write('''<img SRC="god_plot_Vt.png" height=380 width=570>\n''')
            file.write('''<br><br>\n''')

            file.write('''<h2>Distributions of the Initial Parameters over 5000 planets</h2>\n\n''')
            file.write('''<a href="../distributions.html">are here...</a>''')


            file.write('''<h2>Change in Ensemble Likelihood and the 7 Parameters with Chain Link</h2>\n\n''')

            file.write('''The first plot is the progression of the likelihood sampling<br><br>\n''')
            for nbsp in range(4):
                file.write('''&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\n''')
            file.write('''<img SRC="plot_like.png" height=570 width=770><br><br>\n''')
            file.write('''Now the progression of the parameter samplings<br><br>\n''')

            pardict = {0:'lnT',1:'Kcos',2:'Ksin',3:'ecos',4:'esin',5:'S',6:'V0',7:'e',8:'K',9:'T',10:'pomega',11:'phi',12:'lnK'}

            for value in pardict.values():
                file.write('''<img SRC="'''+value+'''_plot_sampling.png" height=380 width=570>\n''')
                file.write('''<img SRC="'''+value+'''_histo_sampling.png" height=380 width=570>\n''')
                file.write('''<br><br>\n''')
           
            file.write('''<br><br>\n''')

            if plot:
                #the likelihood as a function of link
                a.plot_like(hardcopy=1)
                #radial velocity data plots
                a.planets[0].plot_Vt(hardcopy=1)
                #radial velocity truth plots
                a.godplanets[0].plot_Vt(god=1,hardcopy=1)
                for j in range(13):
                    print 'plot', j, 'of', 12,'t....', time()-mainstart,'secs'
                    a.plot_sampling(j,hardcopy=1)
                os.system('mv *sampling*png planet%2.2d' %i)
                os.system('mv *Vt*png planet%2.2d' %i)
                os.system( 'mv plot_like.png planet%2.2d' %i)
                #the class that retains the distributions
                # Ntimes = 1 here, we're just recovering the distributions

            pardict = {0:'lnT',1:'Kcos',2:'Ksin',3:'ecos',4:'esin',5:'S',6:'V0'}

            file.write("<PRE>\n")
            file.write("     PARAMETERS:   ["+" ".join(pardict.values())+"]\n")
            file.write("     ACCEPTED:     "+str(a.accepted)+"\n")
            file.write("     REJECTED:     "+str(a.rejected)+"\n")
            file.write("     FRACTION:     "+str(around((1.*a.accepted/(a.accepted+a.rejected)),3))+"\n")
            file.write("<PRE>\n\n")

            if i > 0:
                file.write('''<a href="../planet%2.2d'''%(i-1) +'''/diagnostics.html">previous planet</a><br>''')
            if i+1 < totalplanets:
                file.write('''<a href="../planet%2.2d'''%(i+1) +'''/diagnostics.html">next planet</a><br>''')

            file.write("""</body>\n</html>\n""")

            file.close()

            print str(around((1.*a.accepted/(a.accepted+a.rejected)),3))
            print str(around(a.stepwidths,4))

            #loop to the next planet
            print 'done planet', i+1, 'of',totalplanets,'t....', time()-mainstart,'secs'
            i = i +1
        except StopIteration:
            badcount += 1
            print 'Failed to converge starting step widths on',badcount,'occasions so far'
            pass

    mainfile.close()
    os.chdir('..')
    print 'Failed to converge starting step widths on',badcount,'occasions'
    return

# If makeplanets = 1, write some planets to cache, if not, read in premade planets from cache and apply to those.
# If makeplanets = 1 and eflag = 1, we write out to crazycache
# If makeplanets = 1 and mflag = 1, we write out to mass1cache
# etc
if __name__ == '__main__':
    from datetime import date
    import os
    from sys import argv

    if len(argv) < 5:
        print '***************************************************************************'
        print 'USAGE -'
        print 'python ersatz_exoplanet.py Nplanets Nlinks Nobs outputdirectory'
        print ' + likemaximflag (1/0) eflag (1/0) mflag (2/1/0) plotflag (1/0) clobberflag (1/0)'
        print '***************************************************************************'
        raise IOError

    #default flags
    likemaximflag = 0
    eflag = 0
    mflag = 0
    plotflag = 0
    clobberflag = 1
   
    Nplanets = int(argv[1])
    Nlinks = int(argv[2])
    Nobs = int(argv[3])
    outputdirectory = argv[4]

    if len(argv) > 5:
        likemaximflag = int(argv[5])
    if len(argv) > 6:
        eflag = int(argv[6])
    if len(argv) > 7:
        mflag = int(argv[7])
    if len(argv) > 8:
        plotflag = int(argv[8])
    if len(argv) > 9:
        clobberflag = int(argv[9])

    outputdirectory += '/'

    print 'Nplanets',Nplanets, 'Nlinks',Nlinks, 'Nobs',Nobs, 'outputdirectory',outputdirectory,'likemaximflag',likemaximflag,'eflag', eflag, 'mflag', mflag, 'plotflag',plotflag, 'clobberflag',clobberflag

    predir = 'e%dm%dmadeplanets/' % (eflag, mflag)
    prefilename = predir+'premadeplanets.pickled'

    #is the made planets directory cache exists, then don't make new planets
    makeplanets = 1
    if os.path.exists(predir):
        makeplanets = 0

    if makeplanets:
        os.mkdir(predir)
        print 'making',Nplanets,'planets to put in',predir
        planetlist = []
        for i in range(Nplanets):
            print 'Planet', i
            planetlist.append(ersatz_planets(1,Nobs,eCrazy=eflag,mflag=mflag,god=True))

        prefile = open(prefilename,'w')

        cPickle.dump(planetlist,prefile)
        prefile.close()
    else:
        prefile = open(prefilename,'r')
        planetlist = cPickle.load(prefile)

        if Nplanets > len(planetlist):
            print '******************************************************************'
            print 'Number of requested planets',Nplanets,'but number stored only',len(planetlist)
            print 'need to make more observations in',predir,'by removing the madeplanets cache directory'
            print '******************************************************************'
            raise IndexError
        planetlist = planetlist[:Nplanets]

        for plan in planetlist:
            if Nobs > len(plan.planets[0].time):
                print '*********************************************************************'
                print 'Number of needed observations',Nobs,'but number stored only',len(plan.planets[0].time)
                print 'need to make more observations in',predir,'by removing the madeplanets cache directory'
                print '*********************************************************************'
                raise IndexError

            plan.planets[0].ntimes = Nobs
            plan.planets[0].time = plan.planets[0].time[0:Nobs]
            plan.planets[0].Verr = plan.planets[0].Verr[0:Nobs]
            plan.planets[0].Vmodel = plan.planets[0].Vmodel[0:Nobs]
            plan.planets[0].Vdata = plan.planets[0].Vdata[0:Nobs]
            plan.planets[0].Voffset = plan.planets[0].Voffset[0:Nobs]

        do_everything(planetlist,Nlinks,outputdirectory,likemaxim=likemaximflag,plot=plotflag,clobber=clobberflag,anneal=likemaximflag)

        print 'plotting distributions'
        os.chdir(outputdirectory)
        fp = ersatz_planets(5000,1,eCrazy=eflag,mflag=mflag)
        #parameter distributions
        fp.plot_edist(hardcopy=1)
        fp.plot_Mpdist(hardcopy=1)
        fp.plot_Tdist(hardcopy=1)
        fp.plot_idist(hardcopy=1)
        fp.plot_likedist(hardcopy=1)
        os.chdir('..')

        print 'DONE'
