# encoding: utf-8
"""


"""


__all__ = ['Transit']

import numpy as np
import scipy.optimize as op

# Planet parameter priors.

Rmin = 0.1    # M_earth
Rmax = 50.0   # M_earth
nR = -1.5


def occultation(p, z):
    z = np.atleast_1d(np.abs(z))

    m1 = z >= 1 + p
    m2 = (np.abs(1 - p) < z) * (z < 1 + p)
    m3 = z <= 1 - p
    m4 = z <= p - 1

    p2 = p ** 2
    z2 = z[m2] ** 2
    tmp = z2 - p2 + 1
    k1 = np.arccos(0.5 * tmp / z[m2])
    k2 = np.arccos(0.5 * (z2 + p2 - 1) / z[m2] / p)
    k3 = np.sqrt(z2 - 0.25 * tmp ** 2)

    z[m1] = 0
    z[m2] = (k1 + p ** 2 * k2 - k3) / np.pi
    z[m3] = p ** 2
    z[m4] = 1

    return 1 - z


class Planet(object):

    def __init__(self, radius, period, incl, phi, e, a):
        self.radius = radius
        self.T = period
        self.incl = incl
        self.phi = phi
        self.e = e
        self.a = a

    def eccentric_anomaly(self, t):
        """
        Numerically solve for the eccentric anomaly as a function of time.

        """
        e = self.e
        wt = 2 * np.pi * np.atleast_1d(t) / self.T + self.phi
        psi0s = wt + e * np.sin(wt)
        f = lambda psi, wt0: wt0 - psi + e * np.sin(psi)
        return np.array(map(lambda (psi0, wt0): op.newton(f, psi0,
                                                args=(wt0,)), zip(psi0s, wt)))

    def coords(self, t, i0=0.0):
        psi = self.eccentric_anomaly(t)
        cpsi = np.cos(psi)

        e = self.e
        d = 1 - e * cpsi

        cth = (cpsi - e) / d
        r = self.a * d

        # In the plane of the orbit.
        x, y = r * cth, np.sign(np.sin(psi)) * r * np.sqrt(1 - cth * cth)

        # Rotate into observable coordinates.
        i = i0 + self.incl
        return x * np.cos(i), y, x * np.sin(i)


class PlanetarySystem(object):

    def __init__(self, radius, incl):
        self.radius = radius
        self.incl = incl
        self.planets = []

    @classmethod
    def sample_prior(cls, Np=1):
        r0 = 150 + 20 * np.random.randn()
        i0 = 0.0 * np.radians(5 * np.random.randn())

        ps = cls(r0, i0)

        radius = ((Rmax ** (nR + 1) - Rmin ** (nR + 1)) * np.random.rand(Np)
                + Rmin ** (nR + 1)) ** (1.0 / (1 + nR))
        period = 0.5 + 500 * np.random.rand(Np)
        incl = 0.0 * np.radians(0.05 * np.random.randn(Np))
        phase = 2 * np.pi * np.random.rand(Np)
        e = 0.3 * np.random.rand(Np)
        a = 500 + 100 * np.random.randn(Np)

        for n in range(Np):
            ps.add_planet(radius[n], period[n], incl[n], phase[n], e[n], a[n])

        return ps

    def add_planet(self, r, T, i, phi, e, a):
        r0 = self.radius
        self.planets.append(Planet(r / r0, T, i, phi, e, a / r0))

    def lightcurve(self, t):
        lc = np.ones_like(t)
        for p in self.planets:
            x, y, z = p.coords(t, i0=self.incl)
            m = x > 0
            D = np.sqrt(y[m] ** 2 + z[m] ** 2)
            lc[m] *= occultation(p.radius, D)
        return lc


if __name__ == '__main__':
    import matplotlib.pyplot as pl

    np.random.seed(42)
    ps = PlanetarySystem(10.0, -1.3)

    # Hot Jupiter.
    T = 0.05
    ps.add_planet(1.2, T, 1.3, np.pi, 0.01, 108)

    # Jupiter.
    # ps.add_planet(1.0, 12.0, 1.3, 0.0, 0.05, 11000)

    # Saturn.
    # ps.add_planet(0.85, 30.0, 2.5, 0.0, 0.05, 21000)

    t = np.linspace(0, 15.8, 8491)
    lc = ps.lightcurve(t)

    pl.plot(t % T, lc, ".k")
    pl.savefig("lightcurve.png")
