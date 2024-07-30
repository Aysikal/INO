import numpy as np
import astropy 
import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt
from scipy import ndimage
import pandas as pd
from scipy.optimize import curve_fit 
from sklearn.linear_model import LinearRegression
from photutils.aperture import CircularAperture, aperture_photometry
from photutils.aperture import CircularAnnulus
from photutils.aperture import ApertureStats
import functions as func
from scipy import ndimage, misc 
from matplotlib import pyplot as plt
from astropy.visualization import LogStretch
from matplotlib.patches import Circle
#input
gain = 16.5
t = 20
binning = 1
seeing = 0.5
airmass = 1.2
extiction = 0.5
m = 15
background_mean = 0
filter = 'r'
#fixed input:
pixel_scale = 0.047
dc = 0.08
h = 6.62620 * 10**(-34)
c = 2.9979250 * 10**(8)
readnoise = 3.7
sigma_f = 0.289
pix_size = 81 * (10**(-12)) #m
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
'''
zero_flux_mag = {'u': 24.63,
                 'g': 25.11,
                 'r': 24.80,
                 'i':24.36,
                 'z':22.83}

b = {'u': 1.4*10**(-10),
     'g':0.9 * 10**(-10),
     'r':1.2 * 10**(-10),
     'i':1.8 * 10**(-10),
     'z':7.8 * 10**(-10)} 

m_10b = {'u': 22.122,
        'g': 22.60,
        'r': 22.29,
        'i':21.85,
        'z':20.32}
'''
#----------------------------------------------------------------------------------------------------------------------------------------------------------------#
#problem values : 
QE = {'u': 0.1,
      'g': 0.7,
      'r': 0.75,
      'i': 0.55,
      'z': 0.2 }
 
delta_i = 0.135*10**(-3) 
#----------------------------------------------------------------------------------------------------------------------------------------------------------------#
npix = 2 *(seeing/pixel_scale)

P = ((h*c)/CW[filter])
m_corrected = m 
f_nu = 10 ** (-0.4 *(m_corrected + 48.6))
print(f_nu)
#f_lambda = f_nu * (1/3.34)*10**(-4)*(10**(-10)/CW[filter])**2
f_lambda = (f_nu *c )/(CW[filter]**2)  #erg/cm2/A
#flux = (10**(-0.4*m_corrected) * f0[filter])  # erg/s/cm2/A
signal = (f_lambda*10**(-7) * band_width[filter]* t * 0.3 * S*(10**4) ) / P
noise = np.sqrt(signal + binning *(background_mean + readnoise**2) +  t*dc*npix)
print(signal/noise) 

#sky ------------------------------------------------------------------------------------------------------------------------------------------------------------#
t = 60
m_sky = 21.8
f_nu_s = 10 ** (-0.4 *(m_sky + 48.6))
f_lambda_s = (f_nu_s *c)/(CW[filter]**2)  #erg/cm2/A
N_sky = (f_lambda_s*10**(-7) * band_width[filter]* t * 0.3 * S*(10**4)*(pixel_scale**2)) / P
print(N_sky)