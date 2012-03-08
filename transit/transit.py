# encoding: utf-8
"""


"""

from __future__ import division

__all__ = ['Transit']

import numpy as np

class Transit(object):
    def __init__(self, R, Rs, Rp, incl):
        self.R  = R
        self.Rs = Rs
        self.p  = Rp/Rs
        self.incl = incl
        self.mu = np.cos(incl)

    def evaluate(self, phase):
        phase = np.atleast_1d(phase)
        m = np.cos(self.incl)*np.cos(phase) > 0
        z = np.sqrt(np.sin(phase)**2 + np.sin(self.incl)**2*np.cos(phase)**2) \
                * self.R/self.Rs
        if len(phase) == 1:
            if m:
                return float(self.occultation(z))
            return 1
        ret = np.ones_like(phase)
        ret[m] = self.occultation(z[m])
        return ret

    def occultation(self, z):
        p = self.p
        z = np.atleast_1d(np.abs(z))

        m1 = z >= 1+p
        m2 = (np.abs(1-p) < z) * (z < 1+p)
        m3 = z <= 1-p
        m4 = z <= p-1

        p2 = p**2
        z2 = z[m2]**2
        tmp = z2-p2+1
        k1 = np.arccos(0.5*tmp/z[m2])
        k2 = np.arccos(0.5*(z2+p2-1)/z[m2]/p)
        k3 = np.sqrt(z2-0.25*tmp**2)

        z[m1] = 0
        z[m2] = (k1 + p**2 * k2 - k3)/np.pi
        z[m3] = p**2
        z[m4] = 1

        return 1-z

    def animate(self, dphase=0.1):
        import time
        import matplotlib.pyplot as pl
        import matplotlib.patches as mpatches

        phase = 0

        pl.ion()
        pl.figure(figsize=(8,8))
        ax = pl.axes([0,0.3,1,0.7], aspect='equal', frameon=False,
                xticks=[], yticks=[])
        ax2 = pl.axes([0.1,0.1,0.8,0.2])

        ph = np.linspace(-np.pi, np.pi, 2000)
        lc = self.evaluate(ph)
        ax2.plot(ph/np.pi, lc, 'k')
        ax2.set_ylim(0.9*np.min(lc), 1.1)
        pt, = ax2.plot(0, self.evaluate(0), "or")

        star = mpatches.Circle([0, 0], 1, color="w", ec="k")
        ax.add_patch(star)

        planet = mpatches.Circle([0, self.R/self.Rs*np.sqrt(1-self.mu**2)],
                self.p, color="gray", ec="k")
        ax.add_patch(planet)
        ax.set_xlim((self.R/self.Rs+self.p)*np.array([-1,1]))
        ax.set_ylim([-2,2])

        for i in range(100):
            phase += dphase
            x = np.sin(phase) * self.R/self.Rs
            y = -np.sin(self.incl) * np.cos(phase) * self.R/self.Rs
            planet.center = [x,y]

            z = np.cos(self.incl)*np.cos(phase)
            if z < 0:
                star2 = mpatches.Circle([0, 0], 1, color="w", ec="k")
                ax.add_patch(star2)

            pt.set_xdata((phase/np.pi+1)%2-1)
            pt.set_ydata(self.evaluate(phase))

            pl.draw()
            if z < 0:
                star2.remove()
            time.sleep(0.05)

if __name__ == '__main__':
    transit = Transit(10, 1, 0.8, np.radians(2))
    transit.animate()

