# First let's set up our packages
import numpy as np
from matplotlib import pyplot as plt
from scipy import integrate

H0kmsmpc = 70.  # Hubble constant in km/s/Mpc
H0s = H0kmsmpc * 3.2408e-20 # H0 in inverse seconds is H0 in km/s/Mpc * (3.2408e-20 Mpc/km)
H0y = H0s * 3.154e7 * 1.e9

om0 = 0.3
ol0 = 0.7
o_r0 = 8.24e-4

astart = 0.0001
astop = 100
astep = 0.0001

bstart = 0.0
bstop = 100
bstep = 1.0

b_arr = np.arange(bstart,bstop,bstep)
a_arr = np.arange(astart,astop,astep)
t_Gyr = np.zeros(len(a_arr))


    
def adotinv_flatmatter(a):
    return np.sqrt(a) 
for i,a_end in enumerate(a_arr): # enumerate adds an index to each value
    t_Hubble,uncert = integrate.quad(adotinv_flatmatter,0,a_end)
    t_Gyr[i] = t_Hubble/H0y

a_arr = np.arange(astart,astop,astep)
a = np.zeros(len(a_arr))
om = np.zeros(len(a_arr))
ol = np.zeros(len(a_arr))
o_r = np.zeros(len(a_arr))
tot = np.zeros(len(a_arr))

for i ,a in enumerate(a_arr): # enumerate adds an index to each value
    om[i] = om0*np.power(a,-3);
    ol[i] = ol0;
    o_r[i] = o_r0*np.power(a,-4);
    tot[i] = om[i]+ol[i]+o_r[i]
    
omp = om/tot
olp = ol/tot
orp = o_r/tot

x=t_Gyr
plt.figure(dpi=1200) 
plt.plot( x,omp,label='omega mass') 
plt.plot( x,olp,label='omega lambda') 
plt.plot( x,orp,label='omega radiation') 
#plt.plot( a_arr,omp) 
#plt.plot( a_arr,olp) 
#plt.plot( a_arr,orp) 
plt.xscale("log")
plt.xlabel('Gyr')
plt.ylabel('factor')
plt.show()

plt.figure(dpi=1200) 
plt.plot( om,x)
plt.xscale("log")
plt.xlabel('omega mass')
plt.ylabel('Gyr')
plt.show()  
 
plt.figure(dpi=1200) 
plt.plot( o_r,x) 
plt.xscale("log")
plt.xlabel('omega r')
plt.ylabel('Gyr')
plt.show()  

plt.figure(dpi=1200) 
plt.plot( omp,x)
plt.xscale("log")
plt.xlabel('omega mass over total')
plt.ylabel('Gyr')
plt.show() 
 
plt.figure(dpi=1200) 
plt.plot( orp,x) 
plt.xscale("log")
plt.xlabel('omega r over total')
plt.ylabel('Gyr')
plt.show()  
