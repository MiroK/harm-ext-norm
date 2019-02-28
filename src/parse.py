import matplotlib.pyplot as plt
import numpy as np
import glob
import os


dx_value = lambda name: (name, float(name[name.index('dx_')+3:]))
dx_pairs = map(dx_value, glob.glob('results/steklov_circle_dx_*'))

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

dxs, mins, maxs = [], [], []
for path, dx in sorted(dx_pairs, key=lambda p: p[1], reverse=True):
    data = np.loadtxt(path)
    
    if not len(data): continue
    
    min_, max_ = data[:, 3].T, data[:, 4].T
    print '\t', dx
    print min_
    print max_

    ax1.semilogy(dx*np.ones_like(min_), 1./min_, marker='o') 
    ax2.semilogy(dx*np.ones_like(max_), max_, marker='x') 
    
    dxs.append(dx)
    mins.append(min_[-1])
    maxs.append(max_[-1])

ax1.semilogy(dxs, 1./np.array(mins))
ax2.semilogy(dxs, maxs)

plt.show()
