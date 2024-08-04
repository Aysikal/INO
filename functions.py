from astropy.io import fits
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.stats import SigmaClip
from photutils.background import Background2D, SExtractorBackground
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from scipy import ndimage
from scipy.ndimage import gaussian_filter
import cv2
import astroalign as aa
from scipy.optimize import curve_fit
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import datetime, time
import pytz
import ephem
from astropy.coordinates import get_body
import math


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
    return radialprofile , fwhm , popt

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

def airmass_function(year, month, day , hour, minute, RA, DEC):
    # INO location
    elevation = 3600
    observer_location = EarthLocation(lat=33.674*u.deg, lon=51.3188*u.deg , height = elevation*u.meter)
    ra = RA.split(':')
    ra_hours = float(ra[0])
    ra_minutes = float(ra[1])
    ra_seconds = float(ra[2])
    ra_deg = (ra_hours * 15) + (ra_minutes * 0.25) + (ra_seconds * 0.00417)

    dec = DEC.split(':')
    dec_degrees = float(dec[0])
    dec_minutes = float(dec[1])
    dec_seconds = float(dec[2])
    dec_deg = dec_degrees + dec_minutes/60 + dec_seconds/3600
    star_position = SkyCoord(ra=ra_deg*u.degree, dec=dec_deg*u.degree, frame='icrs')
    date_time_string = f"{year}-{month:02d}-{day:02d} {hour:02d}:{minute:02d}:00"
    observing_time = Time(date_time_string, scale='utc', location=observer_location)
    # Convert to local horizontal coordinates (AltAz frame)
    star_altaz = star_position.transform_to(AltAz(obstime=observing_time, location=observer_location))
    print("altitude of the star in degrees", star_altaz.alt.degree)
    zenith_angle = 90*u.deg - star_altaz.alt
    zenith_angle_deg = zenith_angle.to_value(u.deg)
    z= np.radians(zenith_angle)
    X = 1 / (np.cos(z) + 0.50572 * (6.07995 + 90 - zenith_angle_deg) ** (-1.6364))
    print("airmass" , X)
    return X


def get_fli(year, month, day, hour, minute):
    # INO location
    latitude = '33:40:28'
    longitude = '51:19:08'
    elevation = 3600
    location = EarthLocation(lat=33.674*u.deg, lon=51.3188*u.deg , height = elevation*u.meter)


    # Create the date_time string
    date_time_string = f"{year}-{month:02d}-{day:02d} {hour:02d}:{minute:02d}:00"
    time = observing_time = Time(date_time_string, scale='utc', location=location)

    #moon data
    moon_position = get_body("moon", time, location)
    ra_hours = moon_position.ra.hour
    ra_degrees = ra_hours*15
    dec_degrees = float(moon_position.dec.deg)
    moon_position = SkyCoord(ra=(ra_degrees)*u.degree, dec=dec_degrees*u.degree, frame='icrs')
    moon_altaz = moon_position.transform_to(AltAz(obstime=observing_time, location=location))
    print("moon altitude", moon_altaz.alt.deg)
    print("moon azimuth", moon_altaz.az.deg )

    #sun data
    sun_position = get_body("sun",time, location)
    sun_ra_hours = sun_position.ra.hour
    sun_ra_degrees = sun_ra_hours*15
    sun_dec_degrees = float(sun_position.dec.deg) 
    sun_coords = SkyCoord(ra=sun_ra_degrees*u.degree, dec=sun_dec_degrees*u.degree)


    # Calculate ecliptic longitudes
    sun_lon = sun_coords.barycentrictrueecliptic.lon.deg
    moon_lon = moon_position.barycentrictrueecliptic.lon.deg

    # Phase angle (0°–360°)
    phase_angle = (moon_lon - sun_lon) % 360.0

    # Convert phase angle to radians
    phase_angle_radians = np.radians(phase_angle)

    # Calculate fraction illuminated (fli)
    fli = 1 - (1 + np.cos(phase_angle_radians)) / 2

    # Convert to percentage
    percentage_illuminated = fli * 100

    # moon phase names:
    if phase_angle < 1.0 or phase_angle >= 359.0:
        phase_name = "New Moon"
    elif 1.0 <= phase_angle < 89.0:
        phase_name = "Waxing Crescent"
    elif 89.0 <= phase_angle <= 91.0:
        phase_name = "First Quarter"
    elif 91.0 < phase_angle < 179.5:
        phase_name = "Waxing Gibbous"
    elif 179.5 <= phase_angle <= 180.5:
        phase_name = "Full Moon"
    elif 180.5 < phase_angle < 269.5:
        phase_name = "Waning Gibbous"
    elif 269.5 <= phase_angle <= 270.5:
        phase_name = "Last Quarter"
    else:
        phase_name = "Waning Crescent"

    return fli

def calculate_sky_magnitude(offset , fli):
    # FLI ranges from 0 to 1 (0 = New Moon, 1 = Full Moon)
    # Assume a constant offset for sky brightness (adjust as needed)
    offset = 21.0  # Example value (mag/arcsec^2)
    if fli ==0: 
          sky_magnitude = offset
    # Calculate sky brightness (in magnitudes per square arcsecond)
    else: sky_magnitude = offset + 2.5 * math.log10(1-fli)

    return sky_magnitude

'''
def airmass(RA, DEC, lon, lat, utc_date_time):
    observer_location = EarthLocation(lat=-29.2563*u.deg, lon=-70.738*u.deg)

    RA = "13:26:47.28"
    DEC = "-47:28:46.092"

    ra = RA.split(':')
    ra_hours = float(ra[0])
    if ra_hours < 0:
        sgn = -1
    else:
        sgn = 1

    ra_minutes = float(ra[1])
    ra_seconds = float(ra[2])
    ra_deg = (ra_hours * 15) + sgn*(ra_minutes * 0.25) + sgn*(ra_seconds * 0.00417)

    if dec_deg < 0: 
        sgn2 = -1
    else:
        sgn2 = 1

    dec = DEC.split(':')
    dec_degrees = float(dec[0])
    dec_minutes = float(dec[1])
    dec_seconds = float(dec[2])
    dec_deg = dec_degrees + sgn2*(dec_minutes/60) + sgn2*(dec_seconds/3600)

    # Example: RA = 10.625 degrees, Dec = 41.2 degrees
    star_position = SkyCoord(ra=ra_deg*u.degree, dec=dec_deg*u.degree, frame='icrs')
    observing_time = Time(utc_date_time, scale='utc', location=observer_location)

    # Convert to local horizontal coordinates (AltAz frame)
    star_altaz = star_position.transform_to(AltAz(obstime=observing_time, location=observer_location))
    zenith_angle = 90*u.deg - star_altaz.alt
    zenith_angle_deg = zenith_angle.to_value(u.deg)

    z= np.radians(zenith_angle)
    X = 1 / (np.cos(z) + 0.50572 * (6.07995 + 90 - zenith_angle_deg) ** (-1.6364))
    return zenith_angle, X
'''