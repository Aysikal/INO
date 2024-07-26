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
t = 60
binning = 1
seeing = 0.5
airmass = 1.2
g_extiction = 0
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
ix_size = 81 * (10**(-12)) #m
D = 3.4 #m
d = 0.6 #m
S = np.pi*(D/2)**2 - np.pi*(d/2)**2
lambda_g = 0.470033 * (10**-6)
b = {'u': 1.4*10**(-10),
     'g':0.9 * 10**(-10),
     'r':1.2 * 10**(-10),
     'i':1.8 * 10**(-10),
     'z':7.8 * 10**(-10)}
zero_flux_mag = {'u': 24.63,
                 'g': 25.11,
                 'r': 24.80,
                 'i':24.36,
                 'z':22.83}
m_10b = {'u': 22.122,
        'g': 22.60,
        'r': 22.29,
        'i':21.85,
        'z':20.32}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------#
#problem values : 
QE_G = 	0.365
delta_i = 0.135*10**(-3) 
#----------------------------------------------------------------------------------------------------------------------------------------------------------------#
#formula:
f0_r = 10**(-0.4**zero_flux_mag['r'])

#mag = -(2.5/np.ln(10))*[np.asinh((f/f0)/2*b[filter])+np.ln(b[filter])]
npix = 2 *(seeing/pixel_scale)

P = ((h*c)/550 * 10**(-9)) * 10**11
m_corrected = m - g_extiction * airmass
m_ab = 15 - (1.2*( 0.11))
f_jy = 10**(-0.4*(m_ab+8.90)) 
'''
jy_to_w_m2_hz = 10 **(-26)
jy_to_w_m2_um = jy_to_w_m2_hz * ((lambda_g)**2/c)
flux = f_jy * jy_to_w_m2_um
'''

'''
f_lambda = f_jy * (3.34)**(-1) * 10**(-4) * ((10**(-10)/lambda_a)**2) * 10**(12)
'''
f_lambda = (f_jy * c * (550 * 10**(-9))**(-2))/10**6

print(f_lambda)
signal = (f_lambda * delta_i * t * QE_G * S) / P
noise = np.sqrt(signal + binning *(background_mean + readnoise**2) +  t*dc*npix)

print(signal/noise)
