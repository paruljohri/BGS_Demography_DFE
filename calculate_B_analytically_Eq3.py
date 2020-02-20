#This script is to get the slope of recovery of pi and the number of  base pairs that would lead to 50%, 75% and 90% recovery for a given recombination rate.

import sys
import math
import numpy as np

#Define variables and constants:
g = 0.0 #rate of gene conversion
r = 1.0*1e-8 #rate of recombination
l = 6000.0 #(*Length of genomic element*)
u = 1.25*1e-8 #(*Mutation rate*)
U = l*u
Ne =20000.0 #(*Effective population size, only required to calculate expected nucleotide diversity under neutrality*)
pi = 4*Ne*u #(*Expected nucleotide diversity under neutrality*)
f0 = 0.22 #(*Proportion of effectively neutral mutations with 0 <= |2Nes| < 1 *)
f1 = 0.27 #(*Proportion of weakly deleterious mutations with 1 <= |2Nes| < 10 *)
f2 = 0.13 #(*Proportion of moderately deleterious mutations with 10 <= |2Nes| < 100 *)
f3 = 0.38 #(*Proportion of strongly deleterious mutations with |2Nes| >= 100 *)
#(*Note that the number of classes can easily be increased to whatever is required to approximate the continuous DFE *)
h = 0.5 #(* dominance coefficient *)
#(*Now we define the boundaries of the fixed intervals over which we will integrate *)
t0 = 0.0
t1 = h*(1/(2*Ne))
t1half = h*(5/(2*Ne)) #(* This is the cut-off value of 2Nes=5. This derivation assumes that all mutations with 2Nes<5 will not contribute to BGS *)
t2 = h*(10/(2*Ne))
t3 = h*(100/(2*Ne))
t4 = h*1.0
posn_specific = 1 #Some specific position that you want to get B and pi at

def calculate_exponent(t_start, t_end, posn):
    a = g + r*posn
    b = g + r*(posn + l)
    E1 = (U/(r*l*(1-a))) * (1 + (a/((1-a)*(t_end-t_start))) * math.log((a+(t_start*(1-a)))/(a + (t_end*(1-a)))))
    E2 = -1.0*(U/(r*l*(1-b)))*(1 + (b/((1-b)*(t_end-t_start)))*math.log((b + ((1-b)*t_start))/(b + ((1-b)*t_end))))
    E = E1 + E2
    return (E)

def calculate_pi(posn):
    E_f1 = calculate_exponent(t1half, t2, posn) 
    E_f2 = calculate_exponent(t2, t3, posn)
    E_f3 = calculate_exponent(t3, t4, posn)
    pi_posn = f0*pi + f1*0.5*pi + f1*0.5*pi*math.exp(-1.0*E_f1) + f2*pi*math.exp(-1.0*E_f2) + f3*pi*math.exp(-1.0*E_f3)
    return (pi_posn)

def calculate_B(posn):
    B_posn = calculate_pi(posn)/pi
    return (B_posn)

def calculate_pi_window(win_start, win_end):
    i=win_start
    pi_sum = 0.0
    while i <= win_end:
        pi_sum = pi_sum + float(calculate_pi(i))
        i = i + 1
    return(pi_sum/(win_end - win_start+1))

#print (calculate_pi(10000))

def fit_log_curve(l_posns, l_pi):
    x_data = np.array(l_posns)
    y_data = np.array(l_pi)
    log_x_data = np.log(x_data)
    curve_fit = np.polyfit(log_x_data, y_data, 1)
    return (curve_fit[0], curve_fit[1])

#Printing pi and B at a specific position:                                      
print ("pi at position " + str(posn_specific) + ":" + '\t' + str(calculate_pi(posn_specific)))
print ("B at position " + str(posn_specific) + ":" + '\t' + str(calculate_B(posn_specific)))
                                                                                
#Getting an average pi over a window:

#Fit a log curve to from position=1 to 100kb
posn_first = 1
posn_last = 10000 #100000
l_posns, l_pi = [], []
j=posn_first
while j <= posn_last:
    l_posns.append(j)
    l_pi.append(calculate_pi(j))
    j += 1
log_fit = fit_log_curve(l_posns, l_pi)
slope = log_fit[0]
intercept = log_fit[1]

#Calculate numbp50, 75 and 90
pi_50 = float(calculate_pi(1)) + ((pi - float(calculate_pi(1)))*0.5)
pi_75 = float(calculate_pi(1)) + ((pi - float(calculate_pi(1)))*0.75)
pi_90 = float(calculate_pi(1)) + ((pi - float(calculate_pi(1)))*0.9)

numbp50 = math.exp((pi_50-intercept)/slope)
numbp75 = math.exp((pi_75-intercept)/slope)
numbp90 = math.exp((pi_90-intercept)/slope)

print ("Number of bases for 50% recovery of pi:" + '\t' + str(numbp50))
print ("Number of bases for 75% recovery of pi:" + '\t' + str(numbp75))
print ("Number of bases for 90% recovery of pi:" + '\t' + str(numbp90))
print("done")


