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
	    Accepted.append(float(data[7]))
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

T, E_mean_NN, E_var_TT, M_mean_NN, M_var_T, M_abs_mean_NN, MCcycle, Accepted = read('4crandomAcceptvsT20')


import matplotlib.pyplot as plt

#4c
plt.plot((MCcycle), E_mean_NN) 
plt.rcParams.update({'font.size': 14})
ax = plt.gca()
ax.get_yaxis().get_major_formatter().set_useOffset(False)
ax.set_xscale('log')
plt.draw()
plt.xlabel('Monte Carlo cycles')
plt.ylabel('<E>')
plt.show()

plt.plot((MCcycle), M_abs_mean_NN)
plt.rcParams.update({'font.size': 14})
ax = plt.gca()
ax.get_yaxis().get_major_formatter().set_useOffset(False)
ax.set_xscale('log')
plt.draw()
plt.xlabel('Monte Carlo cycles')
plt.ylabel('<|M|>')
plt.show()

#4c accepted moves as funtion of MCC  (file 1lattice20 (when T=1))
plt.plot((MCcycle), Accepted)
plt.rcParams.update({'font.size': 14})
ax = plt.gca()
ax.get_yaxis().get_major_formatter().set_useOffset(False)
ax.set_xscale('log')
plt.draw()
plt.xlabel('Monte Carlo cycles')
plt.ylabel('Accepted moves')
plt.show()

#4c accepted moves as funtion of T (file latticeT20 (when T=1))
plt.plot(T, Accepted)
plt.xlabel('Temperature, [kT/J]')
plt.ylabel('Accepted moves')
plt.show()

