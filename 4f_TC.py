# Read data from all the different txt-files and make a least square
# fit to find the best estimate of T_c.
# We use the maxima of C_V for the different values
# of L to estimate the critical temperature.

from math import *
import numpy as np
import matplotlib.pyplot as plt

def read_quantities(filename):
    infile = open(filename, 'r')
    # Elements to be read in file:
    T = []; E = []; Cv = []; Chi = []; absM = [] 
    # Read lines except for the first one:
    lines = infile.readlines()
    for line in lines:
    	words = line.split()
        T.append(float(words[0]))
        E.append(float(words[1]))
        Cv.append(float(words[2]))
        # Skip the storage of M here.
        Chi.append(float(words[4]))
        absM.append(float(words[5]))
    infile.close()
    return T, E, Cv, Chi, absM

def find_maximum(T,CV):
    maxCV = max(CV)
    i = CV.index(maxCV)
    return T[i]

def curve(a,Tc,reciprocal_L):  #reciprocal = inverse = 1/L
    return Tc+a*reciprocal_L

# Simple brute force approach to find the lest square fit:
def best_fit(a_min, a_max, Tc_min, Tc_max, L_list, T_list, N=101):
    A = np.linspace(a_min, a_max, N)
    TC = np.linspace(Tc_min, Tc_max, N)
    TC_best = 0.0; A_best = 0.0; current_min_S = 10**6
    for k in range(len(A)):
        A_test = A[k]
        for l in range(len(TC)):
            TC_test = TC[l]
            S = 0.0
            for i in range(len(T_list)):
                T_i = T_list[i]
                L_i = L_list[i]
                S += ( curve(A_test,TC_test,L_i) - T_i )**2
            if S < current_min_S:
                A_best = A_test
                TC_best = TC_test
                current_min_S = S

    return A_best, TC_best

# Some exact results provided by Onsager (1944):
Tc_exact = 2.0/log(1.0+sqrt(2.0))

# Use txt files with refined resolution in T:
T1, E1, CV1, Chi1, absM1 = read_quantities('helpme40')
T2, E2, CV2, Chi2, absM2 = read_quantities('helpme60')
T3, E3, CV3, Chi3, absM3 = read_quantities('helpme100vol2')
T4, E4, CV4, Chi4, absM4 = read_quantities('helpme140')

T1c = find_maximum(T1,CV1)
T2c = find_maximum(T2,CV2)
T3c = find_maximum(T3,CV3)
T4c = find_maximum(T4,CV4)

one_over_L_list = [1.0/20,1.0/40,1.0/60,1.0/80]  #1/L
T_list = [T1c,T2c,T3c,T4c]

# Finding the best values of the slope and the constant value:
A_b, TC_b = best_fit(0, 3.0, 2.0, 2.5, one_over_L_list, T_list, N=401)
print "Best fit slope: ", A_b
print "Best fit critical temperature: ", TC_b

y = []; reciprocal_refined = one_over_L_list[:]
for val in reciprocal_refined:
    y.append(curve(A_b, TC_b, val))

# Adding the final point manually for visual display:
reciprocal_refined.append(0)
y.append(TC_b)


plt.rcParams.update({'font.size': 14})
fig, ax = plt.subplots(1)
ax.plot(one_over_L_list, T_list,'go',label='Datapunkter fra $C_V$')
ax.plot(reciprocal_refined, y,'b-',label='Minste kvadraters metode, $T_c(\infty) = %.9f$' % TC_b)
ax.legend(loc='best',fancybox='True',shadow='True')
ax.set_ylim([TC_b-0.01,max(T_list)+0.01])
ax.set_xlabel('$1/L$')
ax.set_ylabel('Maxumim av varmekapasitet $\max{(C_V(L))}$')
ax.grid()
plt.show()