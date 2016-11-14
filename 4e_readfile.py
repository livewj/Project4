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

T40, E_mean_NN40, Cv40, M_mean_NN40, X40, M_abs_mean_NN40 = read('helpme40')
T60, E_mean_NN60, Cv60, M_mean_NN60, X60, M_abs_mean_NN60 = read('helpme60')
T100, E_mean_NN100, Cv100, M_mean_NN100, X100, M_abs_mean_NN100 = read('helpme100vol2')
T140, E_mean_NN140, Cv140, M_mean_NN140, X140, M_abs_mean_NN140 = read('helpme140')

import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

#<E>
f1 = UnivariateSpline(T40, E_mean_NN40)
f2 = UnivariateSpline(T60, E_mean_NN60)
f3 = UnivariateSpline(T100, E_mean_NN100)
f4 = UnivariateSpline(T140, E_mean_NN140)
plt.plot(T40, f1(T40)/max(f1(T40))*max(E_mean_NN40), T60, f2(T60)/max(f2(T60))*max(E_mean_NN60), T100, f3(T100)/max(f3(T100))*max(E_mean_NN100), T140, f4(T140)/max(f4(T140))*max(E_mean_NN140))

#plt.plot(T40, E_mean_NN40, T60, E_mean_NN60, T100, E_mean_NN100, T140, E_mean_NN140) 
plt.rcParams.update({'font.size': 14})
plt.legend(['L=40','L=60', 'L=100','L=140'])
plt.xlabel('Temperature [kT/J]')
plt.ylabel('<E>')
plt.show()

#<|M|>
f1 = UnivariateSpline(T40, M_abs_mean_NN40)
f2 = UnivariateSpline(T60, M_abs_mean_NN60)
f3 = UnivariateSpline(T100, M_abs_mean_NN100)
f4 = UnivariateSpline(T140, M_abs_mean_NN140)
plt.plot(T40, f1(T40)/max(f1(T40))*max(M_abs_mean_NN40), T60, f2(T60)/max(f2(T60))*max(M_abs_mean_NN60), T100, f3(T100)/max(f3(T100))*max(M_abs_mean_NN100), T140, f4(T140)/max(f4(T140))*max(M_abs_mean_NN140))

#plt.plot(T40, M_abs_mean_NN40, T60, M_abs_mean_NN60, T100, M_abs_mean_NN100, T140, M_abs_mean_NN140)
plt.rcParams.update({'font.size': 14})
plt.legend(['L=40', 'L=60', 'L=100', 'L=140'])
plt.xlabel('Temperature [kT/J]')
plt.ylabel('<|M|>')
plt.show()

#Cv
f1 = UnivariateSpline(T40, Cv40)
f2 = UnivariateSpline(T60, Cv60)
f3 = UnivariateSpline(T100, Cv100)
f4 = UnivariateSpline(T140, Cv140)
plt.plot(T40, f1(T40)/max(f1(T40))*max(Cv40), T60, f2(T60)/max(f2(T60))*max(Cv60), T100, f3(T100)/max(f3(T100))*max(Cv100),T140, f4(T140)/max(f4(T140))*max(Cv140))

#plt.plot(T40, Cv40, T60, Cv60, T100, Cv100, T140, Cv140)
plt.rcParams.update({'font.size': 14})
plt.legend(['L=40', 'L=60','L=100','L=140'])
plt.xlabel('Temperature [kT/J]')
plt.ylabel('$C_V$')
plt.show()

#X
f1 = UnivariateSpline(T40, X40)
f2 = UnivariateSpline(T60, X60)
f3 = UnivariateSpline(T100, X100)
f4 = UnivariateSpline(T140, X140)
f4.set_smoothing_factor(5)
plt.plot(T40, f1(T40)/max(f1(T40))*max(X40),  T60, f2(T60)/max(f2(T60))*max(X60), T100, f3(T100)/max(f3(T100))*max(X100), T140, f4(T140)/max(f4(T140))*max(X140))

#plt.plot(T40, X40, T60, X60, T100, X100, T140, X140)
plt.rcParams.update({'font.size': 14})
plt.legend(['L=40', 'L=60','L=100','L=140'])
plt.xlabel('Temperature [kT/J]')
plt.ylabel('$\chi(T)$')
plt.show()


