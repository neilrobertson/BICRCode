#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

data = [-0.29640676173133207, 0.043827985784286941, 0.0057229950591030383, 0.089159498207098495, 0.14429545132536603, 0.22271779350660062, 0.062620059923404733, 0.095328514966641298, 0.030274115183646717, 0.00045913687181425725, 0.13025104671503734, 0.051121066556439186, 0.12583756779882799, -0.0045724491565433669, -0.094735977861178988, -0.15993521839209429, -0.11331752839177754, -0.063472751077248132, -0.076083703548705375, -0.037618492958842231]

width = 0.5
ind = np.arange(len(data))


plt.bar(ind, data, width)

plt.xticks(ind + width / 2., ('G1', 'G2', 'G3', 'G4') )
plt.show()