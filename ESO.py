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
extiction = 0
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
                'r': 955 * 10**(-10),
                'i': 1064 * 10**(-10),
                'z': 1248 * 10**(-10)}
# Fukugita et al. (1996)
f0 = {'u': 859.5*10**(-11), #ergs/s/cm2/A
      'g': 466.9*10**(-11),
      'r': 278.0*10**(-11),
      'i': 185.2*10**(-11),
     'z': 131.5*10**(-11)}

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
#----------------------------------------------------------------------------------------------------------------------------------------------------------------#
#problem values : 
QE = {'u': 0.1,
      'g': 0.7,
      'r': 0.75,
      'i': 0.55,
      'z': 0.2 }
 
delta_i = 0.135*10**(-3) 
#----------------------------------------------------------------------------------------------------------------------------------------------------------------#
#mag = -(2.5/np.ln(10))*[np.asinh((f/f0)/2*b[filter])+np.ln(b[filter])]
npix = 2 *(seeing/pixel_scale)

P = ((h*c)/CW[filter])
m_corrected = m - 0.2 * 1.2


flux = (10**(-0.4*m_corrected) * f0[filter])  # erg/s/cm2/A

'''
m_ab = 15 - (1.2*( 0.11))
#Janskey method
f_jy = 10**(-0.4*(m_ab+8.90)) 
jy_to_w_m2_hz = 10 **(-26)
jy_to_w_m2_um = jy_to_w_m2_hz * ((lambda_g)**2/c)
flux = f_jy * jy_to_w_m2_um
'''

'''
f_lambda = f_jy * (3.34)**(-1) * 10**(-4) * ((10**(-10)/lambda_a)**2) * 10**(12)

f_lambda = (f_jy * c * (550 * 10**(-9))**(-2))/10**6

print(f_lambda)
'''
signal = (flux*10**(-7) * band_width[filter]*(10**10) * t * QE[filter] * S*(10**4) ) / P
noise = np.sqrt(signal + binning *(background_mean + readnoise**2) +  t*dc*npix)

print(signal/noise)
