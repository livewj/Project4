import numpy as np
from scipy.stats import norm
import matplotlib.mlab as mlab

def read(filename):
	#function that reads the x and y coordinate
	infile = open(filename, 'r')
	T = []
	E_mean_NN = []     # <E>
	Cv = []            # Cv
	M_mean_NN =[]      # <M>
	X = []             # X
	M_abs_mean_NN = [] #<|M|>

	relevant_lines = infile.readlines()[1:] #skips the irrelevant lines
	for line in relevant_lines:
	    data = line.split()
	    T.append(float(data[0]))
	    E_mean_NN.append(float(data[1]))
	    Cv.append(float(data[2]))
	    M_mean_NN.append(float(data[3]))
	    X.append(float(data[4]))
	    M_abs_mean_NN.append(float(data[5]))
	infile.close()
	T = np.array(T)
	E_mean_NN = np.array(E_mean_NN)
	Cv = np.array(Cv)
	M_mean_NN = np.array(M_mean_NN)
	X = np.array(X)
	M_abs_mean_NN = np.array(M_abs_mean_NN)

	return T, E_mean_NN, Cv, M_mean_NN, X, M_abs_mean_NN

T40, E_mean_NN40, Cv40, M_mean_NN40, X40, M_abs_mean_NN40 = read('4eup40')
T60, E_mean_NN60, Cv60, M_mean_NN60, X60, M_abs_mean_NN60 = read('4eup60')
T100, E_mean_NN100, Cv100, M_mean_NN100, X100, M_abs_mean_NN100 = read('4eup40')
T140, E_mean_NN140, Cv140, M_mean_NN140, X140, M_abs_mean_NN140 = read('4eup40')

import matplotlib.pyplot as plt

#<E>
plt.plot(T40, E_mean_NN40, T60, E_mean_NN60) 
plt.rcParams.update({'font.size': 14})
plt.legend(['L=40','L=60'])
plt.xlabel('Temperature [kT/J]')
plt.ylabel('<E>')
plt.show()

#<|M|>
plt.plot(T40, M_abs_mean_NN40, T60, M_abs_mean_NN60)
plt.rcParams.update({'font.size': 14})
plt.legend(['L=40', 'L=60'])
plt.xlabel('Temperature [kT/J]')
plt.ylabel('<|M|>')
plt.show()

#Cv
plt.plot(T40, Cv40, T60, Cv60)
plt.rcParams.update({'font.size': 14})
plt.legend(['L=40', 'L=60'])
plt.xlabel('Temperature [kT/J]')
plt.ylabel('$C_V$')
plt.show()

#X
plt.plot(T40, X40, T60, X60)
plt.rcParams.update({'font.size': 14})
plt.legend(['L=40', 'L=60'])
plt.xlabel('Temperature [kT/J]')
plt.ylabel('$\chi(T)$')
plt.show()


