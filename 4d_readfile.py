import numpy as np
from scipy.stats import norm
import matplotlib.mlab as mlab

def read(filename):
	#function that reads the x and y coordinate
	infile = open(filename, 'r')
	T = []
	E_mean_NN = []  #Energy expectation value/Nspins/Nspins   <E>
	E_var_TT = []        #Energy variance/T/T                      (sigma_E)^2
	M_mean_NN =[]   #Magnetization expectation value/Nspins/Nspins              <M>
	M_var_T = []    #Magnetization variamce/T				 (sigma_M)^2
	M_abs_mean_NN = [] #Magnetization absolute exp value/Nspins/Nspins  <|M|>
	MCcycle = [] #number of Monte Carlo cycles
	Accepted = [] #number of accetped configurations

	relevant_lines = infile.readlines()[1:] #skips the irrelevant lines
	for line in relevant_lines:
	    data = line.split()
	    T.append(float(data[0]))
	    E_mean_NN.append(float(data[1]))
	    E_var_TT.append(float(data[2]))
	    M_mean_NN.append(float(data[3]))
	    M_var_T.append(float(data[4]))
	    M_abs_mean_NN.append(float(data[5]))
	    MCcycle.append(float(data[6]))
	    #Accepted.append(float(data[7]))
	infile.close()
	T = np.array(T)
	E_mean_NN = np.array(E_mean_NN)
	E_var_TT = np.array(E_var_TT)
	M_mean_NN = np.array(M_mean_NN)
	M_var_T = np.array(M_var_T)
	M_abs_mean_NN = np.array(M_abs_mean_NN)
	MCcycle = np.array(MCcycle)
	Accepted = np.array(Accepted)

	return T, E_mean_NN, E_var_TT, M_mean_NN, M_var_T, M_abs_mean_NN, MCcycle, Accepted

T, E_mean_NN, E_var_TT, M_mean_NN, M_var_T, M_abs_mean_NN, MCcycle, Accepted = read('helpme4dT1')
T2, E_mean_NN2, E_var_TT2, M_mean_NN2, M_var_T2, M_abs_mean_NN2, MCcycle2, Accepted2 = read('4dT2420')


E_relevant = E_mean_NN[19:] #after the steady state situation has been reached
#E_relevant = E_mean_NN2[38:] #after the steady state situation has been reached

E_variance = 0.023237777;
E_variance2 = 1.4091419;


import matplotlib.pyplot as plt


#4d: Histogram when T=1 and T=2.4
weights = np.ones_like(E_relevant)/float(len(E_relevant))
fig, ax = plt.subplots(1)
plt.rcParams.update({'font.size': 14})
n, bins, patches = plt.hist(E_relevant, weights = weights, bins= 20, normed=1, facecolor='red', label='$\sigma_E^2 = %f$' % E_variance2)
(mu, sigma) = norm.fit(E_relevant)
y = mlab.normpdf( bins, mu, sigma)
l = plt.plot(bins, y, 'k--', linewidth=2, label='Best tilpassede Gausskurve')

# Checking the returned variables:
print "n: ", n
print "Bins: ", bins
print "Sum of bin values: ", np.sum(n)
print "Sum of bin areas: ", np.sum(n*np.diff(bins))


ax = plt.gca()
ax.get_xaxis().get_major_formatter().set_useOffset(False)
plt.draw()

plt.xlabel('Energi per spinn, $E/n_{\mathrm{spinn}}$')
plt.ylabel('Frekvens')
plt.legend(loc='upper right',fancybox='True')
plt.title('Energi histogram ved T = %.1f [kT/J]' % T2[1])
plt.show()