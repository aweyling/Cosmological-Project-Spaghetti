# First let's set up our packages
import numpy as np
from matplotlib import pyplot as plt
from scipy import integrate

# And set some constants
c = 299792.458 # km/s (speed of light)

H0kmsmpc = 70.  # Hubble constant in km/s/Mpc
H0s = H0kmsmpc * 3.2408e-20 # H0 in inverse seconds is H0 in km/s/Mpc * (3.2408e-20 Mpc/km)
H0y = H0s * 3.154e7 * 1.e9 # H0 in inverse Giga years is H0 in inverse seconds * (3.154e7 seconds/year) * (1e9 years / Giga year)
print('Hubble time = ',1/H0y,'Gyr')

# Let's format that to a more appropriate number of significant figures.  
# The first % means put the variable here.  The .2f means make it a float with 2 decimal places.
print('Hubble time = %.2f Gyr'%(1/H0y))

# Write a function for the integrand, i.e. $1/\dot{a}$ (adot inverse)
# I'm going to keep the constant H_0 outside the integral and add it later.
def adotinv_flatmatter(a):
    return np.sqrt(a) 

# First note some weird python notation, if you say "x, y = 1, 2" it's the same as saying "x=1" and "y=2".

# Now integrate that equation from scalefactor 0 to 1 (as an example) which will give the age of the universe today (because we define a so that a=1 today):
age_Hubble, uncert = integrate.quad(adotinv_flatmatter,0,1) 
# This returns age in Hubble times, and the uncertainty in the numerical integral

# Convert that to billions of years
age_Gyr = age_Hubble/H0y # Age in Gyr

print('Age of the universe in Hubble times = %.3f which is %.2f Gyr'%(age_Hubble,age_Gyr))
print('Does this match the analytic prediction of 2/3 Hubble times?  Yes!')
print('')
print('Uncertainty on the numerical integral is ', uncert, 'so should be negligible.') # It's good to check this to check that the integral has succeeded.

# Note, to get just the age without the uncertainty you can use [0] to get the first element that's returned,
# and then you can also divide by H0y immediately to get the age in one step
age_Gyr_easier = integrate.quad(adotinv_flatmatter,0,1)[0] / H0y 

print('Just checking that the one-step age calculation (%.2f Gyr) matches the two step one (%.2f Gyr).'%(age_Gyr_easier,age_Gyr))

age_half_Gyr = integrate.quad(adotinv_flatmatter,0,0.5)[0] / H0y 
age_double_Gyr = integrate.quad(adotinv_flatmatter,0,2)[0] / H0y 
print(age_half_Gyr)
print(age_double_Gyr)

#Scalefactor vs time
# Start by making an array of scalefactors
astart = 0.0
astop = 2.1
astep = 0.001 # Make this finer to make the plot smoother
a_arr = np.arange(astart,astop,astep)
print('a_arr=',a_arr)

# First set up an array of times (initially all set to zero) into which we'll put our calculated times
t_Gyr = np.zeros(len(a_arr))  # len(a_arr) gives the length of the a_arr 

# Make a loop and do that integral for every final value of a (starting at a=0 every time)
for i,a_end in enumerate(a_arr): # enumerate adds an index to each value
    t_Hubble,uncert = integrate.quad(adotinv_flatmatter,0,a_end)
    t_Gyr[i] = t_Hubble/H0y

# It's good practise to compare it to the analytic calculation when available
t_analytic_Gyr = 'Add the analytic equation from above here (and uncomment the plotting line for it below)'
    
plt.plot(t_Gyr,a_arr)
plt.plot(age_Gyr, 1.0,'o') # Put a dot at the current time
#plt.plot(t_analytic_Gyr, a_arr,':',color='red')
plt.xlabel('Time (Gyr)')
plt.ylabel('Scalefactor')
plt.show()

astart = 0.0
astop = 2.1
astep = 0.001
a_arr = np.arange(astart,astop,astep)
t_lookback_Gyr = np.zeros(len(a_arr))
for i,a_end in enumerate(a_arr): # enumerate adds an index to each value
    t_lookback_Hubble,uncert = integrate.quad(adotinv_flatmatter,1,a_end)
    t_lookback_Gyr[i] = t_lookback_Hubble/H0y

# Uncomment the commands below to plot your result
plt.plot(t_lookback_Gyr,a_arr) 
plt.axvline(x=0,linestyle=':') # Plot some crosshairs 
plt.axhline(y=1,linestyle=':')
plt.plot(0.0, 1.0,'o') # Put a dot at the current time
plt.xlabel('Lookback time (Gyr)')
plt.ylabel('Scalefactor')
plt.show()

# First calculate the index corresponding to a=1.0.  (Find when |a-1.0| is minimum.  You could also do this by just redoing the integral from 0<a<1, but the way I've set it up above we know we have an a=1 in the array, so we can just find what we've already calculated.)
index_today = np.argmin(np.abs(a_arr - 1.0))

# Then find the time corresponding to that index  (note you do NOT have to do the integral again)
age_today = t_Gyr.item(index_today)

t2_lookback_Gyr = t_Gyr-age_today

# Uncomment the commands below to plot your result
plt.plot(t2_lookback_Gyr,a_arr) 
plt.axvline(x=0,linestyle=':') # Plot some crosshairs 
plt.axhline(y=1,linestyle=':')
plt.plot(0.0, 1.0,'o') # Put a dot at the current time
plt.xlabel('Lookback time (Gyr)')
plt.ylabel('Scalefactor')
plt.show()

# First write a function that takes as input a, Omega_M (om), and Omega_Lambda (ol) and outputs 1/adot
def adotinv(a,om,ol):
    adot=1/(np.sqrt(om*a**-3+(1-om-ol)*a**-2+ol)*a)
    return 1.0/adot

# Calculate for the universe we think we live in, with approximately matter density 0.3 and cosmological constant 0.7
om = 0.1
ol = 1.3

# Note that when you integrate something with more than one argument you pass it with args=(arg1,arg2) in the integrate function
# e.g. "integrate.quad(adotinv, lower_limit, uper_limit, args=(om,ol))""
t_lookback_Gyr,uncert = np.arrange([integrate.quad(adotinv, 1.0, a_end, args=(om,ol))])/H0y

# Uncomment these commands to plot this new model (note I've added a label that can be used in the legend)
plt.plot(t_lookback_Gyr,a_arr,label='$(\Omega_M,\Omega_\Lambda)$=(%.2f,%.2f)'%(om,ol)) 
plt.axvline(x=0,linestyle=':') # Plot some crosshairs 
plt.axhline(y=1,linestyle=':')
plt.xlabel('Lookback time (Gyr)')
plt.ylabel('Scalefactor')
plt.legend(loc='lower right',frameon=False)
plt.show()