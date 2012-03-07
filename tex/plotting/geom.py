#!/usr/bin/env python
# encoding: utf-8
"""
Plot the transit geometry.

"""

import numpy as np

import matplotlib.pyplot as pl
import matplotlib.patches as mpatches

p = 0.5
z = 1.2

x = (z**2+1-p**2)/(2.0*z)
y = np.sqrt(1 - x**2)

fig = pl.figure(figsize=(10,10))
ax = fig.add_subplot(111, aspect='equal') #, frameon=False, xticks=[], yticks=[])

# Add the star
ax.plot(0, 0, "ok")
star = mpatches.Circle([0, 0], 1, color="none", ec="k")
ax.add_patch(star)

# Add the planet
ax.plot(z, 0, "ok")
planet = mpatches.Circle([z, 0], p, color="gray", ec="k", alpha=0.5)
ax.add_patch(planet)

# Add some lines
ax.plot([x,0,x], [-y,0,y], "k")
ax.plot([x,z,x], [-y,0,y], "k")
ax.plot([x,x], [-y,y], "k", lw=2)
ax.text(x, 0, "$a$", va="center", ha="right", size=18)
ax.text(x/2, y/2, "$R_*$", va="bottom", ha="right", size=18)

ax.set_xlim(-1.5,2)
ax.set_ylim(-1.5,1.5)

pl.show()

