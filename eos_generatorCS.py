#EoS generator for cs parametrization

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import quad

output_file= 'eos_aleCS'

p0= 2.6791429e14
c= 2.9979e10
mn= 1.6749e-24  #cgs
mng=939.565  #MeV
n0g=0.16 #fm^-3
finalPNV= 8.81700e+32
P= np.zeros(6+75)

EFTpressure = [0.4470*1.6022e33, 0.7162*1.6022e33, 0.9094*1.6022e33,1.154*1.6022e33, 1.464*1.6022e33, 1.851*1.6022e33]

#EFTenergy

egcrust = [87.9*1.7827e12, 108.2*1.7827e12, 119.5*1.7827e12, 131.5*1.7827e12, 144.4*1.7827e12, 158.0*1.7827e12]

eg= np.linspace(2.99137e+14, 1.88148e+15, 75)

egfinal= np.concatenate((egcrust, eg))

cscsqr = np.zeros(len(egfinal))

#eg was calculated using polytropic equations for the limits and using interpolation for the rest of the equation

def constant(P1, density1):
    ''' generates a constant according to the dimension of the
        pressure at the transition density, a6 was found
        by isolating the constant from the linear integral'''

    def A(e):
        x = e / (mng * n0g)
        return a1 * np.exp((-0.5 * np.power((x - a2), 2)) / np.power(a3, 2))

    def B(e):
        x = e / (mng * n0g)
        return (1 / 3) / (1 + np.exp(-a5 * (x - a4)))

    def C(e):
        return 1

    def D(e):
        x = e / (mng * n0g)
        return 1 / (1 + np.exp(-a5 * (x - a4)))

    Aint1 = quad(A, 0, density1)[0]
    Bint1 = quad(B, 0, density1)[0]
    Cint1 = quad(C, 0, density1)[0]
    Dint1 = quad(D, 0, density1)[0]


    a6= (P1-Aint1-Bint1)/(Cint1-Dint1)

    return a6

def pressure(density):

    def csc(e):
        x = e / (mng * n0g)
        return a1*np.exp((-0.5*np.power((x-a2),2))/np.power(a3, 2)) + a6 + ((1/3)-a6)/(1+np.exp(-a5*(x-a4)))

    P= quad(csc, 0, density)[0]
    return P

while True:

    '''Check if the values are correct according to equation 4'''
    a1= np.random.uniform(0.1, 1.5)
    a2= np.random.uniform(1.5, 12)
    a3= np.random.uniform(0.05, 2) * a2
    a4= np.random.uniform(1.5, 37)
    a5= np.random.uniform(0.1, 1)

    for i in range (len(egcrust)):

        P[i]= EFTpressure[i]

    a6= constant(2.163, 167.8)

    CSCsqr= a1*np.exp((-0.5*np.power((((229.8) / (mng * n0g))-a2),2))/np.power(a3, 2)) + a6 + ((1/3)-a6)/(1+np.exp(-a5*(((229.8) / (mng * n0g))-a4)))

    print(a6)

    valid = True

    if CSCsqr > 0.163 or CSCsqr<0:

        print("Values discarded")

        continue

    for i in range (len(egfinal)):

        cscsqr[i]= (a1*np.exp((-0.5*np.power((((egfinal[i]/1.7827e12) / (mng * n0g))-a2),2))/np.power(a3, 2)) + a6 + ((1/3)-a6)/(1+np.exp(-a5*(((egfinal[i]/1.7827e12) / (mng * n0g))-a4))))

        if cscsqr[i] < 0:

            print("Values discarded")

            valid= False

            break #exit for loop

    if not valid:
        continue

    break

for i in range (len(eg)):

    P[i+6]= pressure(eg[i]/1.7827e12)*1.60218e33

# Read the data from data.txt
with open('pVNV.txt', 'r') as file:
    data_lines = file.readlines()

# Write the data to the output file
with open(output_file, 'w') as f:

    for line in data_lines:
        f.write(line)

    for j in range(len(P)):
        f.write('{:.5e} {:.5e}\n'.format(egfinal[j], P[j]))


for i in range (len(egfinal)):

    cscsqr[i]= (a1*np.exp((-0.5*np.power((((egfinal[i]/1.7827e12) / (mng * n0g))-a2),2))/np.power(a3, 2)) + a6 + ((1/3)-a6)/(1+np.exp(-a5*(((egfinal[i]/1.7827e12) / (mng * n0g))-a4))))

plt.plot(np.log10(egfinal), np.log10(P))
plt.xlabel('Log10(Energy Density)')
plt.ylabel('Log10(Pressure)')
plt.title('Equation of State for Piecewise Polytropes')
plt.show()

#the shape is consistent with figure 1 but not the value of a6

plt.plot(egfinal, cscsqr)
plt.xlabel('energy_density')
plt.ylabel('(cs/c)^2')
plt.xscale('log')
plt.title('Equation of State for Piecewise Polytropes')
plt.show()
