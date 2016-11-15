# Makefile program for running Project4_d.cpp


# Read data from 'Ising2D_d.txt' and plot energy histogram,
# i.e. energy probability distribution.

from math import *
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

def read_E(filename):
    infile = open(filename, 'r')
    # Elements to be read in file:
    E = [];
    # Read lines except for the first one:
    lines = infile.readlines()
    first_line = lines[0]
    words_in_first_line = first_line.split()
    T = float(words_in_first_line[-1])
    for line in lines[2:-1]:
    	words = line.split()
        E.append(float(words[0]))
    final_line = lines[-1]
    E_variance = float(final_line.split()[-1])
    infile.close()
    return T, E, E_variance

# Fetching data by a call on read_E_and_M:
T, E, E_variance = read_E('4d_omg')

print "Energy per spin std.: ", E_variance

# Plot results:
fig, ax = plt.subplots(1)
plt.rcParams.update({'font.size': 14})
n, bins, patches = plt.hist(E, 5, normed=1, label='$\sigma_E^2 = %9f$' % E_variance)
(mu, sigma) = norm.fit(E)
y = mlab.normpdf( bins, mu, sigma)
variance_fit = sigma*sigma
l = plt.plot(bins, y, 'k--', linewidth=2, label='Best tilpasset Gausskurve, $\sigma^2 = %.9f$' % variance_fit)

# Checking the returned variables:
print "n: ", n
print "Bins: ", bins
print "Sum of bin values: ", np.sum(n)
print "Sum of bin areas: ", np.sum(n*np.diff(bins))

plt.xlabel('$E/n_{\mathrm{spinn}}$')
plt.ylabel('$P(E)$')
plt.legend(loc='upper right',fancybox='True')
plt.title('Energi histogram ved T = %.1f [kT/J]' % T)
#ax.set_ylim([0.0,4.0])
plt.show()