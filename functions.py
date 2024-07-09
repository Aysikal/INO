from astropy.io import fits
import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt
from astropy.stats import SigmaClip
from photutils.background import Background2D,SExtractorBackground
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from scipy import ndimage
from scipy.ndimage import gaussian_filter
import cv2
import astroalign as aa
from scipy.optimize import curve_fit 

def open_fits(path):
    fitsfile = fits.open(path)
    file = fitsfile[0].data
    return file 
def get_boxes(images,center,box_size):
    # Slice the array
    box_size = int(box_size/2)
    image_boxes = []
    for i in range(0 , len(images)):
        box = images[i][int(center[1]) - box_size : int(center[1]) + box_size, int(center[0]) - box_size : int(center[0]) + box_size]
        image_boxes.append(box)
    return image_boxes

def gaussian(x, a, x0, sigma, c):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2)) + c 

def radial_profile(data, center, box_size):
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

def get_radius(image , center , HWHM , readnoise , gain, radius_step = 0.5 , inner_radius=2.5, outer_radius=3.0):
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