#EoS generator for piecewise polytropes

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import quad

output_file= 'eos_alePP'
output_file2= 'eos_units'


p0= 2.6791429e14 #transition mass density in cgs #2.9468e14
c= 2.9979e10 #speed of light in cgs
finalPNV= 3.92900e+32
final_e= 1.00000e+14

#creating piecewise structure of a neutron star using different intervals of mass density

gamma1=4
gamma2=3
gamma3=2.5
p12=3*p0
p23=4.5*p0

crust= np.linspace(p0*0.5792, p0*1.0371, 6)

core1= np.linspace(p0*1.1, p12, 30)

core2= np.linspace(3.1*p0, p23, 30)

core3= np.linspace(4.6*p0, 5.4*p0, 30)

#initializing arrays
P0= np.zeros(len(crust))
e0= np.zeros(len(crust))


P1= np.zeros(len(core1))
e1= np.zeros(len(core1))


P2= np.zeros(len(core2))
e2= np.zeros(len(core2))


P3= np.zeros(len(core3))
e3= np.zeros(len(core3))


#Calculates Pressure value between the last value of pressure in the NV EoS and the first value at p0

EFTpressure = [0.4470*1.6022e33, 0.7162*1.6022e33, 0.9094*1.6022e33,1.154*1.6022e33, 1.464*1.6022e33, 1.851*1.6022e33]

EFTenergy = [87.9*1.7827e12, 108.2*1.7827e12, 119.5*1.7827e12, 131.5*1.7827e12, 144.4*1.7827e12, 158.0*1.7827e12]

#fill the pressure array with the interpolated values

for i in range (len(crust)):

    P0[i] = EFTpressure[i]

    e0[i] = EFTenergy[i]



P1[0]= 2.163*1.6022e33
e1[0]= 167.8*1.7827e12


k1= P1[0]/np.power(core1[0], gamma1)

a1= (e1[0]/core1[0])-1-((k1/np.power(c, 2))/(gamma1-1))*np.power(core1[0], gamma1-1)



for i in range (len(core1)-1):

    P1[i+1] = k1 * np.power(core1[i+1], gamma1)

    e1[i+1] = (1+a1)*core1[i+1] + (1/(gamma1-1)) * (P1[i+1]/np.power(c, 2))



    if np.sqrt((gamma1*(P1[i+1])/e1[i+1]+(P1[i+1]/np.power(c, 2))))>c:
        print("Hold your horses! You are faster than light!")

k2= P1[len(core1)-1]/np.power(p12, gamma2)

a2= (e1[len(core1)-1]/p12)-1-((k2/np.power(c, 2))/(gamma2-1))*np.power(p12, gamma2-1)



for i in range (len(core2)):

    P2[i] = k2 * np.power(core2[i], gamma2)

    e2[i] = (1+a2)*core2[i] + (1/(gamma2-1)) * (P2[i]/np.power(c, 2))



    if np.sqrt((gamma2*(P2[i])/(e2[i]+(P2[i]/np.power(c, 2)))))>c:
        print("Hold your horses! You are faster than light!")

k3= P2[len(core2)-1]/np.power(p23, gamma3)

a3= (e2[len(core2)-1]/p23)-1-((k3/np.power(c, 2))/(gamma3-1))*np.power(p23, gamma3-1)

for i in range (len(core3)):

    P3[i]= k3 * np.power(core3[i], gamma3)

    e3[i]= (1+a3)*core3[i] + (1/(gamma3-1)) * (P3[i]/np.power(c, 2))

    if np.sqrt((gamma3*(P3[i])/(e3[i]+(P3[i]/np.power(c, 2)))))>c:
        print("Hold your horses! You are faster than light!")

P= np.concatenate((P0, P1, P2, P3))
e= np.concatenate((e0, e1, e2, e3))


# Read the data from the NVeos file
with open('PVNVPP.txt', 'r') as file:
    data_lines = file.readlines()

with open(output_file, 'w') as f:

    for line in data_lines:
        f.write(line)

    for j in range(len(P)):
        f.write('{:.5e} {:.5e}\n'.format(e[j], P[j]))

mass_density = np.concatenate((crust, core1, core2, core3))

plt.plot(np.log10(mass_density), np.log10(P))
plt.xlabel('Log10(Energy Density)')
plt.ylabel('Log10(Pressure)')
plt.title('Equation of State for Piecewise Polytropes')
plt.show()
