############################################################################################################################################################
#This code was written by Aysan Hemmati. on the summer of 2024. 
#you can contact me for any additional questions or information via Email 
#email address :aysanhemmatiortakand@gmail.com
#github = https://github.com/Aysikal
############################################################################################################################################################
import numpy as np
from functions import airmass_function , get_fli , calculate_sky_magnitude
import math

#date and time
mode = input("choose calculator mode. enter either (snr) for snr calculator or (exp) for exposure time calculator.")
print(mode)
if mode == 'snr': 
      print("ATTENTION: exposure time should be in seconds")
      t = float(input("Enter exposure time: "))
if mode == 'exp':
      snr = int(input("snr: "))


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
filter = input("filter choose from u , g, r, i, z: ")
print("enter magnitude in the chosen filter. note that the magnitude should be in the AB system.")
m = float(input("magnitude: "))
#system 
full_well = 70000
print("Binning is either 1x1 and 2x2, enter either 1 or 2 for each respectively.")
binning = int(input("enter binning: "))


#fixed values:
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

#problem values : ----------------------------------------------------------------------------------------------------------------------------------------------------------------#
#should be updated based on INO
extiction = {'u': 0.2,
      'g': 0.2,
      'r': 0.2,
      'i': 0.2,
      'z': 0.2 }

reflectivity = 0.6
filter_reflection = 0.8

E =   {'u': 0.1*(reflectivity)**2 * filter_reflection,
      'g': 0.7*(reflectivity)**2 * filter_reflection,
      'r': 0.75*(reflectivity)**2 * filter_reflection,
      'i': 0.55*(reflectivity)**2 * filter_reflection,
      'z': 0.2*(reflectivity)**2 * filter_reflection }

#SKY BACKGROUND MAGNITUDE WHEN NO MOON IS PRESENT.
offset = {'u': 22,
        'g': 22,
        'r': 22,
        'i': 22,
        'z': 22}

#---------------------------------------------------------------------------------------------------------------------------------------------------------
def calculate_snr(year, month, day, hour, minute, RA, DEC, seeing, pixel_scale, binning, h, c, CW, filter, m, extiction, band_width, t, E, S, get_fli, offset, calculate_sky_magnitude, readnoise):
    # Signal calculation
    airmass = airmass_function(year, month, day, hour, minute, RA, DEC)
    npix = (np.pi * ((seeing / pixel_scale) ** 2)) / (binning ** 2)
    P = ((h * c) / CW[filter])
    m_corrected = m + (airmass * extiction[filter])
    f_nu = 10 ** (-0.4 * (m_corrected + 48.6))
    f_lambda = ((f_nu * c) / (CW[filter] ** 2)) * (10 ** (-10))  # erg/cm2/A
    '''
    b = f_lambda * (band_width*(10**10))*(10**(-7)) * (10**4) #watt/m2 * S
    center_pix = 2 * b * math.erf(pixel_scale/(2*(2**(-1/2)) * seeing))
    if center_pix * t >= full_well:
         print("CCD SATURATED!")
    '''
    A = (f_lambda * 10 ** (-7) * (band_width[filter] * (10 ** (10))) *  E[filter] * S * (10 ** 4)) / P
    signal = A * t

    # Sky calculation
    fli = get_fli(year, month, day, hour, minute)
    sky_mag = calculate_sky_magnitude(offset[filter], fli)
    f_nu_s = 10 ** (-0.4 * (sky_mag + 48.6))
    f_lambda_s = ((f_nu_s * c) / (CW[filter] ** 2)) * (10 ** (-10))  # erg/cm2/A
    C = (f_lambda_s * 10 ** (-7) * (band_width[filter] * (10 ** (10))) * E[filter] * S * (10 ** 4) * (pixel_scale ** 2)) / P
    N_sky = C * t

    # Noise calculation
    B = npix * (N_sky + readnoise ** 2)
    noise = np.sqrt(A * t  + B)

    signal_to_noise = signal / noise
    return signal_to_noise


def solve_for_t(A, npix, C, readnoise, s):
    a = A**2
    b = -s**2 * (A + npix * C)
    c = -s**2 * npix * readnoise**2

    discriminant = b**2 - 4 * a * c

    if discriminant < 0:
        return None  # No real solution

    t1 = (-b + math.sqrt(discriminant)) / (2 * a)
    t2 = (-b - math.sqrt(discriminant)) / (2 * a)
    if t1 >= 0 and t2 >= 0:
        return min(t1, t2)
    elif t1 >= 0:
        return t1
    elif t2 >= 0:
        return t2


def calculate_exposure_time(snr, year, month, day, hour, minute, RA, DEC, seeing, pixel_scale, binning, h, c, CW, filter, m, extiction, band_width, E, S, get_fli, offset, calculate_sky_magnitude, readnoise):
    airmass = airmass_function(year, month, day, hour, minute, RA, DEC)
    npix = (np.pi * ((seeing / pixel_scale) ** 2)) / (binning ** 2)
    P = ((h * c) / CW[filter])
    m_corrected = m + (airmass * extiction[filter])
    f_nu = 10 ** (-0.4 * (m_corrected + 48.6))
    f_lambda = ((f_nu * c) / (CW[filter] ** 2)) * (10 ** (-10))  # erg/cm2/A
    A = (f_lambda * 10 ** (-7) * (band_width[filter] * (10 ** (10))) *  E[filter] * S * (10 ** 4)) / P
    fli = get_fli(year, month, day, hour, minute)
    sky_mag = calculate_sky_magnitude(offset[filter], fli)
    f_nu_s = 10 ** (-0.4 * (sky_mag + 48.6))
    f_lambda_s = ((f_nu_s * c) / (CW[filter] ** 2)) * (10 ** (-10))  # erg/cm2/A
    C = (f_lambda_s * 10 ** (-7) * (band_width[filter] * (10 ** (10))) * E[filter] * S * (10 ** 4) * (pixel_scale ** 2)) / P
    exposure_time = solve_for_t(A,npix,C,readnoise,snr)
    return exposure_time


if mode == 'snr':
      SNR = calculate_snr(year, month, day, hour, minute, RA, DEC, seeing, pixel_scale, binning, h, c, CW, filter, m, extiction, band_width, t, E, S, get_fli, offset, calculate_sky_magnitude, readnoise)
      print("SNR is", SNR)

if mode == 'exp':
      t = calculate_exposure_time(snr, year, month, day, hour, minute, RA, DEC, seeing, pixel_scale, binning, h, c, CW, filter, m, extiction, band_width, E, S, get_fli, offset, calculate_sky_magnitude, readnoise)
      print("Exposure time t is:", t)
