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
pixel_number = 65
gain = 16.5
readnoise = 3.7
exp =   96.9


directory = r"C:/Users\AYSAN\Desktop/project\data/green"
number_of_lights = len(os.listdir(directory))

corrected_lights = []
for file in os.listdir(directory):
   filepath = os.path.join(directory, file)
        # Open the fits file
   light_file = fits.open(filepath)
        # Get the data from the first extension
   light_data = light_file[0].data
   corrected_light = light_data
   corrected_lights.append(corrected_light)



#star locations: 
locations = [[[658, 1142],[665, 1148],[1149, 1588],[1246, 2238],[1052, 2327],[665, 2503]]
             ,[[2490, 753],[2498, 748],[2981, 1198],[3081, 1848],[2885,1939],[2499, 2115]]] #star1
             

number_of_stars = len(locations)
centers = []
stars = []
for j in range(0,number_of_stars):
     for i in range(0,len(corrected_lights)):
        image = corrected_lights[i]
        star = image[locations[j][i][1]-pixel_number : locations[j][i][1]+pixel_number , locations[j][i][0]-pixel_number : locations[j][i][0]+pixel_number]
        stars.append(star)
        center = ndimage.center_of_mass(star)
        centers.append(center)

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

def HWHM(starimages , starcenters ):
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
    
         corrected_brightness = (sum_brightness - background_brightness)
         noise = np.sqrt(sum_brightness+ len(star_data) * readnoise**2 + len(star_data)*gain*np.mean(background_data))
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

        noise = np.sqrt(corrected_brightness + (np.pi *(radius**2)) * readnoise**2 + (np.pi *(radius**2))*gain*background_brightness_density)
        
        #shot_noise = (sum_brightness/gain)**(1/2)
        #noise = np.sqrt(((shot_noise)**2 * gain) + (np.pi*(radius**2)) * readnoise**2)
        
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

from datetime import datetime
import pytz

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

local_image_time = ["2023-10-09T23:08:32.276","2023-10-09T23:11:33.522" ,"2023-10-09T23:14:07.129","2023-10-09T23:16:40.234", "2023-10-09T23:19:25.813", "2023-10-09T23:21:39.124"]

def extinction(magnitude_data,RA, DEC, local_image_time):
    star_magnitudes = magnitude_data[0]
    star_errors = magnitude_data[2]
    year , month, day, hour1, minute1 = convert_iran_to_utc(local_image_time[0])
    year , month, day, hour2, minute2 = convert_iran_to_utc(local_image_time[1])
    year , month, day, hour3, minute3 = convert_iran_to_utc(local_image_time[2])
    year , month, day, hour4, minute4 = convert_iran_to_utc(local_image_time[3])
    year , month, day, hour5, minute5 = convert_iran_to_utc(local_image_time[4])
    year , month, day, hour6, minute6 = convert_iran_to_utc(local_image_time[5])

    x1 = functions.airmass_function(year, month, day, hour1, minute1, RA, DEC)
    x2 = functions.airmass_function(year, month, day, hour2, minute2, RA, DEC)
    x3 = functions.airmass_function(year, month, day, hour3, minute3, RA, DEC)
    x4 = functions.airmass_function(year, month, day, hour4, minute4, RA, DEC)
    x5 = functions.airmass_function(year, month, day, hour5, minute5, RA, DEC)
    x6 = functions.airmass_function(year, month, day, hour6, minute6, RA, DEC)

    x = [x1 ,x2 , x3, x4, x5, x6]
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

#2 minutes exposure:
star1images = stars[0:6]
star1centers = centers[0:6]
star2images = stars[6:12]
star2centers = centers[6:12]
'''
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
HWHM2 = HWHM(star2images,star2centers)
'''
HWHM3 = HWHM(star3images,star3centers)
HWHM4 = HWHM(star4images,star4centers)
HWHM5 = HWHM(star5images,star5centers)
HWHM6 = HWHM(star6images,star6centers)
HWHM7 = HWHM(star7images,star7centers)
HWHM8 = HWHM(star8images,star8centers)
'''
HWHMs = HWHM1 + HWHM2 #+ HWHM3 + HWHM4 + HWHM5 + HWHM6 + HWHM7 + HWHM8

star1radius = get_best_radius_for_star(star1images , star1centers , HWHM1)[0]
star2radius = get_best_radius_for_star(star2images , star2centers , HWHM2)[0]
'''
star3radius = get_best_radius_for_star(star3images , star3centers , HWHM3)[0]
star4radius = get_best_radius_for_star(star4images , star4centers , HWHM4)[0]
star5radius = get_best_radius_for_star(star5images , star5centers , HWHM5)[0]
star6radius = get_best_radius_for_star(star6images , star6centers , HWHM6)[0]
star7radius = get_best_radius_for_star(star7images , star7centers , HWHM7)[0]
star8radius = get_best_radius_for_star(star8images , star8centers , HWHM8)[0]
'''
star1magnitues = get_magnitudes(star1images , star1centers , star1radius , exp)
star2magnitues = get_magnitudes(star2images , star2centers , star2radius , exp)
'''
star3magnitues = get_magnitudes(star3images , star3centers , star3radius , exp)
star4magnitues = get_magnitudes(star4images , star4centers , star4radius , exp)
star5magnitues = get_magnitudes(star5images , star5centers , star5radius , exp)
star6magnitues = get_magnitudes(star6images , star6centers , star6radius , exp)
star7magnitues = get_magnitudes(star7images , star7centers , star7radius , exp)
star8magnitues = get_magnitudes(star8images , star8centers , star8radius , exp)
'''
'''
for i in range(0 , len(stars)):
     x1 , y1 , x2 , y2 = plots(stars[i], centers[i] , HWHMs[i])

     #radial profile:
     fig = plt.figure() 
     ax = fig.add_subplot(111) 
     ax.plot(x1, y1, c='k') 

     popt, pcov = curve_fit(gaussian, x1, y1) 

     fit = gaussian(x1, popt[0], popt[1], popt[2]) 
     ax.plot(x1, fit, c='r', label='Best fit for image') 
     ax.legend() 
     
     if i < 32:
         fig.savefig('radialprofiles_plots/exp = 120, radialprofile %i.png'%(i+1))
     else: 
         fig.savefig('radialprofiles_plots/exp = 180, radialprofile %i.png'%(i+1))
    

     plt.title("radial profile for image %i"%(i))
     plt.xlabel('distance from the center of the star')
     plt.ylabel('brightness')
     plt.show()
     #snr profile
     fig2 = plt.figure() 
     plt.plot(x2,y2)
     plt.title("snr with respect to aperture radius for image %i"%(i))
     plt.xlabel('aperture radius')
     plt.ylabel('SNR')
    
     if i < 32:
         fig2.savefig('snr_plots/exp = 120, SNR vs. r %i.png'%(i+1))
     else: 
         fig2.savefig('snr_plots/exp = 180, SNR vs. r %i.png'%(i+1))
     
     plt.show()
'''

extinction(star1magnitues,"22:15:30.59","51:34:53.727", local_image_time)
extinction(star2magnitues,"22:15:30.59","51:34:53.727", local_image_time)

