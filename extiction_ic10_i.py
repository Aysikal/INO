import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt
import reza
from scipy import ndimage
import pandas as pd
from scipy.optimize import curve_fit 
from sklearn.linear_model import LinearRegression
from photutils.aperture import CircularAperture
from photutils.aperture import CircularAnnulus
from photutils.aperture import ApertureStats
import functions
from datetime import datetime
import pytz

# HWHM 
pixel_number = 65
gain = 16.5
readnoise = 3.7
exp2 = 120
exp3 = 180
'''
#import darks: 
dark_file = fits.open(r"C:/Users\AYSAN\Desktop/project\background\blackwidow/2min\background_2min_r.fits")
dark = dark_file[0].data

dark_file_3 = fits.open(r"C:/Users\AYSAN\Desktop/project\background\blackwidow/3min\background_3min_R.fits")
dark3 = dark_file_3[0].data

#import gaintable : 
gain_file = fits.open(r"C:/Users\AYSAN\Desktop/project/masterflats/r-gaintable.fits")
gaintable = gain_file[0].data
'''
directory = r"C:/Users\AYSAN\Desktop/project/INO/new-data/ic-10 field/Infrared/30s"
directory_3 = r"C:/Users\AYSAN\Desktop/project/INO/new-data/ic-10 field/Infrared/200s"

corrected_lights = []
for file in os.listdir(directory):
   filepath = os.path.join(directory, file)
        # Open the fits file
   light_file = fits.open(filepath)
        # Get the data from the first extension
   light_data = light_file[0].data
   corrected_light = (light_data)
   corrected_lights.append(corrected_light)

corrected_lights_3 = []
for file in os.listdir(directory_3):
   filepath = os.path.join(directory_3, file)
        # Open the fits file
   light_file = fits.open(filepath)
        # Get the data from the first extension
   light_data = light_file[0].data
   corrected_light_3 = (light_data)
   corrected_lights_3.append(corrected_light_3)


# functions: 
def gaussian(x, a, x0, sigma): 
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) 

def radial_profile(data, center):
    normal_data = data/np.max(data)
    y, x = np.indices((normal_data.shape))
    r = np.sqrt((x - center[1])**2 + (y - center[0])**2)
    r = r.astype(int)
    tbin = np.bincount(r.ravel(), normal_data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    peak = np.max(normal_data[r <= 2*pixel_number])
    x_ax = np.linspace( 0 , 2*pixel_number - center[1] , num= len(radialprofile))
    
    return radialprofile , tbin, np.min(r[normal_data <= (peak / 2)])*2 , nr , x_ax

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def HWHM(starimages , starcenters):
 HWHM = []
 for i in range( 0 , len(starimages)):
   star_radial = radial_profile(starimages[i],starcenters[i])
   rad = star_radial[0]
   x_ax = star_radial[-1]
   # Executing curve_fit on noisy data 
   popt, pcov = curve_fit(gaussian, x_ax, rad) 
   #popt returns the best fit values for parameters of the given model
   hwhm = 2* popt[-1]
   HWHM.append(hwhm)
 return HWHM

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

        return best_radius , snr , snrs , radii

def get_best_radius_for_star(images , centers , HWHM):
    list_of_radii = []
    for i in range( 0 , len(images)):
        image = images[i]
        center = centers[i]
        best_radius , snr , snrs , radii = get_radius(image , center , HWHM[i])
        list_of_radii.append(best_radius)
    best_radius_for_star = np.max(list_of_radii)
    return best_radius_for_star , list_of_radii

def get_magnitudes(images , centers , radius , exp ,inner_radius=2.5, outer_radius=3.0):
     magnitudes = []
     snrs = []
     errors = []
     for i in range ( 0 , len(images)):
        image = images[i]
        center = centers[i]
        y, x = np.indices((image.shape))
        image = image.astype(np.float64)
        circlestat = ApertureStats(image , CircularAperture(center , radius) )
        sum_brightness = circlestat.sum
        aperstats = ApertureStats(image, CircularAnnulus(center , inner_radius*radius , outer_radius*radius))
        background_brightness_density = aperstats.median
        background_brightness = (background_brightness_density)*(np.pi*(radius**2))
        corrected_brightness = sum_brightness - background_brightness
        shot_noise = (sum_brightness/gain)**(1/2)
        noise = (np.sqrt(((shot_noise)**2) *gain) + (np.pi *(radius**2)) * readnoise**2)
        snr = corrected_brightness / noise
        snrs.append(snr)
        error = 1.08 / snr
        errors.append(error)
        magnitude = (-2.5 * np.log10((corrected_brightness/exp)))
        magnitudes.append(magnitude)
     return  magnitudes , snrs , errors

def plots(starimage , starcenter  ,HWHM):
   radial_data = radial_profile(starimage , starcenter)
   radial_x = radial_data[-1] 
   radial_y = radial_data[0]
   snr_data = get_radius(starimage , starcenter , HWHM)
   snr_ax = snr_data[2]
   radius_ax = snr_data[3]
   return radial_x , radial_y , radius_ax , snr_ax 
def line(x, a, b):
    return (a * x) + b

def dms_to_decimal(dms): #in the format "45°30'15.5\""
  
    try:
        # Split the input string into components
        degrees, rest = dms.split("°")
        minutes, seconds = rest.split("'")

        # Convert to integers or floats
        degrees = int(degrees)
        minutes = int(minutes)
        seconds = float(seconds[:-1])  # Remove the trailing double quote

        # Calculate decimal degrees
        decimal_degrees = degrees + minutes / 60 + seconds / 3600
        return decimal_degrees
    except ValueError:
        return None  

def convert_iran_to_utc(iran_time_str):
    # Define the Tehran timezone
    tehran_tz = pytz.timezone('Asia/Tehran')

    # Parse the datetime string to a datetime object
    tehran_time = datetime.strptime(iran_time_str, "%Y-%m-%dT%H:%M:%S.%f")

    # Localize the datetime object to Tehran timezone
    tehran_time = tehran_tz.localize(tehran_time)

    # Convert the datetime object to UTC
    utc_time = tehran_time.astimezone(pytz.utc)

    # Extract the year, month, day, hour, and minute
    year = utc_time.year
    month = utc_time.month
    day = utc_time.day
    hour = utc_time.hour
    minute = utc_time.minute

    return year, month, day, hour, minute


def extinction(magnitude_data,RA, DEC, local_image_time):
    star_magnitudes = magnitude_data[0]
    star_errors = magnitude_data[2]
    year , month, day, hour1, minute1 = convert_iran_to_utc(local_image_time[0])
    year , month, day, hour2, minute2 = convert_iran_to_utc(local_image_time[1])
    year , month, day, hour3, minute3 = convert_iran_to_utc(local_image_time[2])
    year , month, day, hour4, minute4 = convert_iran_to_utc(local_image_time[3])
    year , month, day, hour5, minute5 = convert_iran_to_utc(local_image_time[4])
    year , month, day, hour6, minute6 = convert_iran_to_utc(local_image_time[5])
    year , month, day, hour7, minute7 = convert_iran_to_utc(local_image_time[6])
    year , month, day, hour8, minute8 = convert_iran_to_utc(local_image_time[7])

    x1 = functions.airmass_function(year, month, day, hour1, minute1, RA, DEC)
    x2 = functions.airmass_function(year, month, day, hour2, minute2, RA, DEC)
    x3 = functions.airmass_function(year, month, day, hour3, minute3, RA, DEC)
    x4 = functions.airmass_function(year, month, day, hour4, minute4, RA, DEC)
    x5 = functions.airmass_function(year, month, day, hour5, minute5, RA, DEC)
    x6 = functions.airmass_function(year, month, day, hour6, minute6, RA, DEC)
    x7 = functions.airmass_function(year, month, day, hour7, minute7, RA, DEC)
    x8 = functions.airmass_function(year, month, day, hour8, minute8, RA, DEC)


    x = [x1 ,x2 , x3, x4, x5, x6, x7, x8]
    x = np.asanyarray(x)
    # plotting differential magnitudes in respect to secz
    # gives a and b values and their respective errors:
    popt, pcov = curve_fit(line, x, star_magnitudes, sigma=star_errors, absolute_sigma=True)
    line_fit = line(x, *popt)
    plt.xlabel("airmass")
    plt.ylabel("instrumental magnitude")
    
    plt.errorbar(x, star_magnitudes, yerr = star_errors ,fmt='none',ecolor = 'blue',color='yellow') 

    # plotting a weighted fitted line:
    plt.plot(x, line_fit, label='Weighted fit (WLS)')
    # plotting the data:
    plt.plot(x, star_magnitudes, '.')
    plt.legend(loc='lower center')
    plt.show()
    
    # calculate regression:
    x = np.array(x).reshape((-1, 1))
    y = np.array(star_magnitudes)

    model = LinearRegression().fit(x, y)
    r_sq = model.score(x, y)
    print(f"coefficient of determination (Regression): {r_sq}")
    print(f"intercept: {popt[1]}")
    print(f"slope: {popt[0]}")
    return {'Weighted fit parameters:':popt,
            "covariants error" : np.sqrt(np.diag(pcov))} , line_fit 


#star locations: 
loc_2 = [[[3096,2435],[3101,2432],[3096,2425]]]
number_of_stars = len(loc_2)
centers = []
stars = []
for j in range(0,number_of_stars):
     for i in range(0,len(corrected_lights)):
        image = corrected_lights[i]
        star = image[loc_2[j][i][1]-pixel_number : loc_2[j][i][1]+pixel_number , loc_2[j][i][0]-pixel_number : loc_2[j][i][0]+pixel_number]
        #plt.imshow(star)
        #plt.show()
        stars.append(star)
        center = ndimage.center_of_mass(star)
        centers.append(center)
        
loc_3 = [[[3095,2426],[3094,2828],[3098,3231],[2678,3225],[2270,3220]]]
stars_3 = []
centers_3 = []    
for j in range(0,number_of_stars):
     for i in range(0,len(corrected_lights_3)):
        image = corrected_lights_3[i]
        star = image[loc_3[j][i][1]-pixel_number : loc_3[j][i][1]+pixel_number , loc_3[j][i][0]-pixel_number : loc_3[j][i][0]+pixel_number]
        #plt.imshow(star)
        #plt.show()
        stars_3.append(star)
        center = ndimage.center_of_mass(star)
        centers_3.append(center)

#30s exposure:
star1images = stars[0:3]
star1centers = centers[0:3]
'''
star2images = stars[3:6]
star2centers = centers[3:6]
star3images = stars[6:9]
star3centers = centers[6:9]
star4images = stars[9:12]
star4centers = centers[9:12]
star5images = stars[12:15]
star5centers = centers[12:15]
star6images = stars[15:18]
star6centers = centers[15:18]
star7images = stars[18:21]
star7centers = centers[18:21]
star8images = stars[21:24]
star8centers = centers[21:24]
'''
HWHM1 = HWHM(star1images,star1centers)
'''
HWHM2 = HWHM(star2images,star2centers)
HWHM3 = HWHM(star3images,star3centers)
HWHM4 = HWHM(star4images,star4centers)
HWHM5 = HWHM(star5images,star5centers)
HWHM6 = HWHM(star6images,star6centers)
HWHM7 = HWHM(star7images,star7centers)
HWHM8 = HWHM(star8images,star8centers)
'''
HWHMs = HWHM1# + HWHM2 + HWHM3 + HWHM4 + HWHM5 + HWHM6 + HWHM7 + HWHM8

star1radius = get_best_radius_for_star(star1images , star1centers , HWHM1)[0]
'''
star2radius = get_best_radius_for_star(star2images , star2centers , HWHM2)[0]
star3radius = get_best_radius_for_star(star3images , star3centers , HWHM3)[0]
star4radius = get_best_radius_for_star(star4images , star4centers , HWHM4)[0]
star5radius = get_best_radius_for_star(star5images , star5centers , HWHM5)[0]
star6radius = get_best_radius_for_star(star6images , star6centers , HWHM6)[0]
star7radius = get_best_radius_for_star(star7images , star7centers , HWHM7)[0]
star8radius = get_best_radius_for_star(star8images , star8centers , HWHM8)[0]
'''
star1magnitues = get_magnitudes(star1images , star1centers , star1radius , exp2)
'''
star2magnitues = get_magnitudes(star2images , star2centers , star2radius , exp2)
star3magnitues = get_magnitudes(star3images , star3centers , star3radius , exp2)
star4magnitues = get_magnitudes(star4images , star4centers , star4radius , exp2)
star5magnitues = get_magnitudes(star5images , star5centers , star5radius , exp2)
star6magnitues = get_magnitudes(star6images , star6centers , star6radius , exp2)
star7magnitues = get_magnitudes(star7images , star7centers , star7radius , exp2)
star8magnitues = get_magnitudes(star8images , star8centers , star8radius , exp2)
'''
#3:20 min exp :
star1images_3 = stars_3[0:5]
star1centers_3 = centers_3[0:5]
'''
star2images_3 = stars_3 [4:8]
star2centers_3 = centers_3[4:8]
star3images_3 = stars_3[8:12]
star3centers_3 = centers_3[8:12]
star4images_3 = stars_3[12:16]
star4centers_3 = centers_3[12:16]
star5images_3 = stars_3[16:20]
star5centers_3 = centers_3[16:20]
star6images_3 = stars_3[20:24]
star6centers_3 = centers_3[20:24]
star7images_3 = stars_3[24:28]
star7centers_3 = centers_3[24:28]
star8images_3 = stars_3[28:32]
star8centers_3 = centers_3[28:32]
'''
HWHM1_3 = HWHM(star1images_3,star1centers_3)
'''
HWHM2_3 = HWHM(star2images_3,star2centers_3)
HWHM3_3 = HWHM(star3images_3,star3centers_3)
HWHM4_3 = HWHM(star4images_3,star4centers_3)
HWHM5_3 = HWHM(star5images_3,star5centers_3)
HWHM6_3 = HWHM(star6images_3,star6centers_3)
HWHM7_3 = HWHM(star7images_3,star7centers_3)
HWHM8_3 = HWHM(star8images_3,star8centers_3)
'''
HWHMs_3 = HWHM1_3# + HWHM2_3 + HWHM3_3 + HWHM4_3 + HWHM5_3 + HWHM6_3 + HWHM7_3 + HWHM8_3

star1radius_3 = get_best_radius_for_star(star1images_3 , star1centers_3 , HWHM1_3)[0]
'''
star2radius_3 = get_best_radius_for_star(star2images_3 , star2centers_3 , HWHM2_3)[0]
star3radius_3 = get_best_radius_for_star(star3images_3 , star3centers_3 , HWHM3_3)[0]
star4radius_3 = get_best_radius_for_star(star4images_3 , star4centers_3 , HWHM4_3)[0]
star5radius_3 = get_best_radius_for_star(star5images_3 , star5centers_3 , HWHM5_3)[0]
star6radius_3 = get_best_radius_for_star(star6images_3 , star6centers_3 , HWHM6_3)[0]
star7radius_3 = get_best_radius_for_star(star7images_3 , star7centers_3 , HWHM7_3)[0]
star8radius_3 = get_best_radius_for_star(star8images_3 , star8centers_3 , HWHM8_3)[0]
'''
star1magnitues_3 = get_magnitudes(star1images_3 , star1centers_3 , star1radius_3 , exp3)
'''
star2magnitues_3 = get_magnitudes(star2images_3 , star2centers_3 , star2radius_3 , exp3)
star3magnitues_3 = get_magnitudes(star3images_3 , star3centers_3 , star3radius_3 , exp3)
star4magnitues_3 = get_magnitudes(star4images_3 , star4centers_3 , star4radius_3 , exp3)
star5magnitues_3 = get_magnitudes(star5images_3 , star5centers_3 , star5radius_3 , exp3)
star6magnitues_3 = get_magnitudes(star6images_3 , star6centers_3 , star6radius_3 , exp3)
star7magnitues_3 = get_magnitudes(star7images_3 , star7centers_3 , star7radius_3 , exp3)
star8magnitues_3 = get_magnitudes(star8images_3 , star8centers_3 , star8radius_3 , exp3)
'''

local_image_time = ["2023-10-10T00:57:37.575","2023-10-10T00:58:13.382" ,"2023-10-10T00:59:04.649","2023-10-10T01:02:36.497", "2023-10-10T01:06:07.580", "2023-10-10T01:09:46.340","2023-10-10T01:13:38.056","2023-10-10T01:17:16.486"]

star1magnitues_list= star1magnitues + star1magnitues_3
print(star1magnitues_list)
extinction(star1magnitues_list, "00:20:25.50","59:17:17.0", local_image_time) 
