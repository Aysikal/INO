import numpy as np
from functions import airmass_function , get_fli , calculate_sky_magnitude

#date and time
print("ATTENTION: Time and date entries MUST be UTC")
year = int(input("Enter the year: "))
month = int(input("Enter the month (1-12): "))
day = int(input("Enter the day (1-31): "))
hour = int(input("Enter the hour (0-23): "))
minute = int(input("Enter the minute (0-59): "))

#object inputs:
print("RA MUST be in the form of HH:MM:SS")
RA = input("RA: ")
print("DEC MUST be in the form of DD:MM:SS")
DEC = input("DEC: ")
print("enter magnitude in the chosen filter. note that the magnitude should be in the AB system.")
m = float(input("magnitude: "))

#system 
print("ATTENTION: exposure time should be in seconds")
t = int(input("Enter exposure time: "))
print("Binning is either 1x1 and 2x2, enter either 1 or 2 for each respectively.")
binning = int(input("enter binning: "))
filter = ("filter choose from u , g, r, i, z")

#fixed values:
extiction = 0.5
seeing = 1
pixel_scale = 0.047
dc = 0.08
h = 6.62620 * 10**(-34)
c = 2.9979250 * 10**(8)
readnoise = 3.7
D = 3.4 #m
d = 0.6 #m
S = np.pi*(D/2)**2 - np.pi*(d/2)**2 
#effective wavelengths (Central Wavelength)
CW = {'u': 3560 *10**(-10),
     'g': 4825 * 10**(-10),
     'r': 6261 * 10**(-10),
     'i': 7672 * 10**(-10),
     'z': 9097 * 10**(-10)}

# Fukugita et al. (1996)
band_width = {'u': 463 *10**(-10),
                'g': 988 * 10**(-10),
                'r': 1340* 10**(-10),
                'i': 1064 * 10**(-10),
                'z': 1248 * 10**(-10)}
# Fukugita et al. (1996)
f0 = {'u': 859.5*(10**(-11)), #ergs/s/cm2/A
      'g': 466.9*(10**(-11)),
      'r': 278.0*(10**(-11)),
      'i': 185.2*(10**(-11)),
     'z': 131.5*(10**(-11))}

# Fukugita et al. (1996)
f0 = {'u': 8590.5*(10**(-11)), #watt/m2/um
      'g': 4660.9*(10**(-11)),
      'r': 2780.0*(10**(-11)),
      'i': 1850.2*(10**(-11)),
     'z': 1310.5*(10**(-11))}

#problem values : ----------------------------------------------------------------------------------------------------------------------------------------------------------------#

E =   {'u': 0.1,
      'g': 0.7,
      'r': 0.75,
      'i': 0.55,
      'z': 0.2 }
 
#signal -----------------------------------------------------------------------------------------------------------------------------------------------------#
airmass = airmass_function(year, month, day, hour, minute, RA, DEC)
npix = (np.pi*((seeing/pixel_scale)**2)) / (binning)**2
P = ((h*c)/CW[filter])
m_corrected = m + (airmass * extiction)
f_nu = 10 ** (-0.4 *(m_corrected + 48.6))
f_lambda = ((f_nu *c)/(CW[filter]**2)) * (10 **(-10)) #erg/cm2/A
signal = (f_lambda*10**(-7) * (band_width[filter]*(10**(10)))* t * E[filter] * S *(10**4) ) / P

#sky ------------------------------------------------------------------------------------------------------------------------------------------------------------#
fli = get_fli(year, month, day, hour, minute)
offset = {'u': 22,
        'g': 22,
        'r': 22,
        'i': 22,
        'z': 22}
sky_mag = calculate_sky_magnitude(offset[filter], fli)
f_nu_s = 10 ** (-0.4 *(sky_mag + 48.6))
f_lambda_s = ((f_nu_s *c)/(CW[filter]**2)) * (10 **(-10)) #erg/cm2/A
N_sky = (f_lambda_s*10**(-7) * (band_width[filter]*(10**(10))) * t * E[filter] * S *(10**4) * (pixel_scale**2)) / P

#noise---------------------------------------------------------------------------------------------------------------------------------------
noise = np.sqrt(signal + npix *(N_sky + readnoise**2) + t*dc*npix)
print("SNR is" , (signal/noise))
