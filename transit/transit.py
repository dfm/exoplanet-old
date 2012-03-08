# encoding: utf-8
"""


"""

from __future__ import division

__all__ = ['Transit']

import numpy as np

class Transit(object):
    def __init__(self, Rs, Rp, incl):
        self.Rs = Rs
        self.p  = Rp/Rs
        self.incl = incl

    def occultation(self, z):
        p = self.p
        z = np.atleast_1d(np.abs(z))

        m1 = z > 1+p
        m2 = (np.abs(1-p) <= z) * (z <= 1+p)
        m3 = z < 1-p
        m4 = z < p-1

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

if __name__ == '__main__':
    import matplotlib.pyplot as pl

    transit = Transit(1, 0.5, 10)
    z = np.linspace(-2, 2, 500)

    pl.plot(z, transit.occultation(z), "k")

    pl.show()

