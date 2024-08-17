############################################################################################################################################################
#This the source code for INO exposure time calculator. 
#This code was written by Aysan Hemmati. on the summer of 2024. 
#you can contact me for any additional questions or information via Email 
#email address :aysanhemmatiortakand@gmail.com
#github = https://github.com/Aysikal
#documentation : https://docs.google.com/document/d/1D06Oxp0n3AuIE-N1idVbhhsnLQNO6vtncn_TaSSiS_Q/edit
############################################################################################################################################################
import numpy as np
import math

#date and time
mode = input("choose calculator mode. enter either (snr) for snr calculator or (exp) for exposure time calculator.")
print(mode)
if mode == 'snr': 
      print("ATTENTION: exposure time should be in seconds")
      t = float(input("Enter exposure time: "))
if mode == 'exp':
      snr = int(input("snr: "))

object_type = input("choose object type. choose between (point source) or (extended): ")
if object_type == "extended":
       omega_i = float(input("solid angle: "))

airmass = float(input("Enter airmass: "))
fli = float(input("Enter FLI: "))
print("choose from the following options for turbulence: optimal turbulence (0.6 - 0.8) , minimal turbulence (0.8 - 1) , moderate turbulence (1 - 1.3), high turbulence (1.3 - 1.5), very high turbulence (1.5 - 2)")
seeing_conditions = input("Only type the name not the range. Example : minimal ")
seeing_dict = {"optimal": 8,
          "minimal" : 1,
          "moderate": 1.3,
          "high": 1.5,
          "very high":2}

seeing = seeing_dict[seeing_conditions]
filter = input("filter choose from u , g, r, i, z: ")
print("enter magnitude in the chosen filter. note that the magnitude should be in the AB system.")
m = float(input("magnitude: "))
#system 
full_well = 70000
print("Binning is either 1x1 and 2x2, enter either 1 or 2 for each respectively.")
binning = int(input("enter binning: "))


#fixed values:
pixel_scale = 0.047
dc = 0.08
h = 6.62620 * 10**(-34)
c = 2.9979250 * 10**(8)
readnoise = 3.7
D = 3.4 #m
d = 0.6 #m
S = np.pi*(D/2)**2 - np.pi*(d/2)**2 
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
#filter data
#effective wavelengths (Central Wavelength)
CW = {'u': 3530 *10**(-10),
     'g': 4860 * 10**(-10),
     'r': 6260 * 10**(-10),
     'i': 7670* 10**(-10),
     'z': 9097 * 10**(-10)}

# Fukugita et al. (1996) / ??
band_width = {'u': 630 *10**(-10),
                'g': 1530 * 10**(-10),
                'r': 1340* 10**(-10),
                'i': 1330 * 10**(-10),
                'z': 1248 * 10**(-10)}

#problem values : ----------------------------------------------------------------------------------------------------------------------------------------------------------------#
#should be updated based on INO

extiction = {'u': 0.2,
      'g': 0.2,
      'r': 0.2,
      'i': 0.2,
      'z': 0.2 }

reflectivity = 0.6 #mirror reflectivity


#how much of the light goes though the filter. most probabley depends on the filter at use:
filter_reflection = {'u': 0.8,
                    'g': 0.8,
                    'r': 0.8,
                    'i': 0.8,
                    'z': 0.8 }   

filter_reflection = 0.8 #or it may be a constant for all 


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

def calculate_sky_magnitude(offset , fli):
    # FLI ranges from 0 to 1 (0 = New Moon, 1 = Full Moon)
    # Assume a constant offset for sky brightness (adjust as needed)
    offset = 21.0  # Example value (mag/arcsec^2)
    if fli ==0: 
          sky_magnitude = offset
    # Calculate sky brightness (in magnitudes per square arcsecond)
    else: sky_magnitude = offset + 2.5 * math.log10(1-fli)

    return sky_magnitude
#---------------------------------------------------------------------------------------------------------------------------------------------------------
def calculate_snr(airmass,fli, seeing, pixel_scale, binning, h, c, CW, filter, m, extiction, band_width, t, E, S, offset, calculate_sky_magnitude, readnoise):
    # Signal calculation
    npix = (np.pi * ((seeing / pixel_scale) ** 2)) / (binning ** 2)
    P = ((h * c) / CW[filter])
    m_corrected = m + (airmass * extiction[filter])
    f_nu = 10 ** (-0.4 * (m_corrected + 48.6))
    f_lambda = ((f_nu * c) / (CW[filter] ** 2)) * (10 ** (-10))  # erg/cm2/A
    A = (f_lambda * 10 ** (-7) * (band_width[filter] * (10 ** (10))) *  E[filter] * S * (10 ** 4)) / P
    signal = A * t
    # Sky calculation
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

def calculate_snr_extended(omega, airmass,fli, seeing, pixel_scale, binning, h, c, CW, filter, m, extiction, band_width, t, E, S, offset, calculate_sky_magnitude, readnoise):
    # Signal calculation
    npix = (np.pi * ((seeing / pixel_scale) ** 2)) / (binning ** 2)
    P = ((h * c) / CW[filter])
    m_corrected = m + (airmass * extiction[filter])
    f_nu = 10 ** (-0.4 * (m_corrected + 48.6))
    f_lambda = ((f_nu * c) / (CW[filter] ** 2)) * (10 ** (-10))  # erg/cm2/A
    A = (f_lambda * 10 ** (-7) * (band_width[filter] * (10 ** (10))) *  E[filter] * S * (10 ** 4) * omega) / P
    signal = A * t
    # Sky calculation
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

def calculate_exposure_time(snr, airmass,fli, seeing, pixel_scale, binning, h, c, CW, filter, m, extiction, band_width, E, S, offset, calculate_sky_magnitude, readnoise):
    npix = (np.pi * ((seeing / pixel_scale) ** 2)) / (binning ** 2)
    P = ((h * c) / CW[filter])
    m_corrected = m + (airmass * extiction[filter])
    f_nu = 10 ** (-0.4 * (m_corrected + 48.6))
    f_lambda = ((f_nu * c) / (CW[filter] ** 2)) * (10 ** (-10))  # erg/cm2/A
    A = (f_lambda * 10 ** (-7) * (band_width[filter] * (10 ** (10))) *  E[filter] * S * (10 ** 4)) / P
    sky_mag = calculate_sky_magnitude(offset[filter], fli)
    f_nu_s = 10 ** (-0.4 * (sky_mag + 48.6))
    f_lambda_s = ((f_nu_s * c) / (CW[filter] ** 2)) * (10 ** (-10))  # erg/cm2/A
    C = (f_lambda_s * 10 ** (-7) * (band_width[filter] * (10 ** (10))) * E[filter] * S * (10 ** 4) * (pixel_scale ** 2)) / P
    exposure_time = solve_for_t(A,npix,C,readnoise,snr)
    return exposure_time
def calculate_exposure_time_extended(omega, snr, airmass,fli, seeing, pixel_scale, binning, h, c, CW, filter, m, extiction, band_width, E, S, offset, calculate_sky_magnitude, readnoise):
    npix = (np.pi * ((seeing / pixel_scale) ** 2)) / (binning ** 2)
    P = ((h * c) / CW[filter])
    m_corrected = m + (airmass * extiction[filter])
    f_nu = 10 ** (-0.4 * (m_corrected + 48.6))
    f_lambda = ((f_nu * c) / (CW[filter] ** 2)) * (10 ** (-10))  # erg/cm2/A
    A = (f_lambda * 10 ** (-7) * (band_width[filter] * (10 ** (10))) *  E[filter] * S * (10 ** 4)) * omega/ P
    sky_mag = calculate_sky_magnitude(offset[filter], fli)
    f_nu_s = 10 ** (-0.4 * (sky_mag + 48.6))
    f_lambda_s = ((f_nu_s * c) / (CW[filter] ** 2)) * (10 ** (-10))  # erg/cm2/A
    C = (f_lambda_s * 10 ** (-7) * (band_width[filter] * (10 ** (10))) * E[filter] * S * (10 ** 4) * (pixel_scale ** 2)) / P
    exposure_time = solve_for_t(A,npix,C,readnoise,snr)
    return exposure_time



if object_type == "point source":
    if mode == 'snr':
      SNR = calculate_snr(airmass,fli, seeing, pixel_scale, binning, h, c, CW, filter, m, extiction, band_width, t, E, S, offset, calculate_sky_magnitude, readnoise)
      print("SNR is", SNR)

    if mode == 'exp':
      t = calculate_exposure_time(snr,airmass,fli,seeing, pixel_scale, binning, h, c, CW, filter, m, extiction, band_width, E, S, offset, calculate_sky_magnitude, readnoise)
      print("Exposure time t is:", t)

if object_type == "extended": 
    if mode == "snr": 
        SNR = calculate_snr_extended(omega_i , airmass,fli, seeing, pixel_scale, binning, h, c, CW, filter, m, extiction, band_width, t, E, S, offset, calculate_sky_magnitude, readnoise)
        print("SNR is", SNR)
    if mode == "exp":
        t = calculate_exposure_time_extended(omega_i, snr,airmass,fli,seeing, pixel_scale, binning, h, c, CW, filter, m, extiction, band_width, E, S, offset, calculate_sky_magnitude, readnoise)
        print("Exposure time t is:", t)