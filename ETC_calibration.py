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
from photutils.aperture import CircularAperture
from photutils.aperture import CircularAnnulus
from photutils.aperture import ApertureStats
import functions as func
from scipy import ndimage, misc 
from matplotlib import pyplot as plt
from astropy.visualization import LogStretch

gain = 16.5
readnoise = 1
box_size = 40

sample_file = func.open_fits(r"C:/Users\AYSAN\Desktop/project/INO_airmass\data\BW/G/1/light-g-2023_10_10-exp00.02.00.000-1x1_High_1.fit")
center = [1590, 1600]

box = sample_file[int(center[1]) - box_size : int(center[1]) + box_size, int(center[0]) - box_size : int(center[0]) + box_size]
bkg_snr = 1/(np.std(box))**(1/2)
print(bkg_snr)
plt.imshow(box, origin="lower")
plt.show()
star_center = ndimage.center_of_mass(box)
box_size = box_size - 2

centered_box = box[int(star_center[1]) - box_size : int(star_center[1]) + box_size, int(star_center[0]) - box_size : int(star_center[0]) + box_size]
def gaussian(x, a, x0, sigma, c):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2)) + c 

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
import numpy as np

def radial_profile(data, center):
    data = data/np.max(data)
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(int)
    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    radialprofile = radialprofile[:box_size]
    popt, _ = curve_fit(gaussian, np.arange(len(radialprofile)), radialprofile)
    # Get the FWHM from the fitted Gaussian
    fwhm = 2.355 * popt[2]
    return radialprofile , fwhm/2 , popt

def get_radius(image , center , HWHM , radius_step = 0.5 , inner_radius=2.5, outer_radius=3.0):
        radius_min = HWHM *1/2
        radius_max = HWHM*2

        y, x = np.indices((image.shape))
        dist = np.sqrt((x - center[1])**2 + (y - center[0])**2)
        image = image.astype(np.float64)

        max_snr = 0
        snrs = []
        radii = []
   
        for radius in np.arange(radius_min, radius_max + radius_step, radius_step):
        
         star_data = image[dist <= radius]
         sum_brightness = np.sum(star_data)

         background_data = (image[(dist > inner_radius * radius)
                                & (dist <= outer_radius * radius)])
         background_brightness = np.mean(background_data) * len(star_data)
         
         corrected_brightness = sum_brightness - background_brightness

         noise = np.sqrt(sum_brightness + len(star_data) * readnoise**2 + len(star_data)*gain*np.mean(background_data))
         snr = corrected_brightness / noise
         snrs.append(snr)
         radii.append(radius)
         if snr > max_snr:
             max_snr = snr
             best_radius = radius
             best_brightness = corrected_brightness

        return best_radius , max_snr , snrs , radii

'''
# Example usage:
radial_profile_values, hwhm, popt = radial_profile(centered_box, star_center)
print(hwhm)
print(f"HWHM distance: {hwhm} pixels")
# Plot the radial profile
plt.plot(radial_profile_values, label='Radial Profile')
plt.plot(gaussian(np.arange(len(radial_profile_values)), *popt), label='Gaussian Fit')
plt.xlabel('Distance from center (pixels)')
plt.ylabel('Intensity')
plt.legend()
plt.show()

best_radius , max_snr, snrs, radii = get_radius(centered_box, star_center, hwhm)
plt.plot(radii , snrs)
plt.xlabel("radius from the center of the star")
plt.ylabel("snr")
plt.show()
'''

