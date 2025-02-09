############################################################################################################################################################
#This code was written by Aysan Hemmati. on the summer of 2024. 
#you can contact me for any additional questions or information via Email 
#email address :aysanhemmatiortakand@gmail.com
#github = https://github.com/Aysikal
############################################################################################################################################################

from astropy.io import fits
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from astropy.coordinates import get_body
import math

def open_fits(path):
    fitsfile = fits.open(path)
    file = fitsfile[0].data
    return file 

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
    #print("altitude of the star in degrees", star_altaz.alt.degree)
    zenith_angle = 90*u.deg - star_altaz.alt
    zenith_angle_deg = zenith_angle.to_value(u.deg)
    z= np.radians(zenith_angle)
    X = 1 / (np.cos(z) + 0.50572 * (6.07995 + 90 - zenith_angle_deg) ** (-1.6364))
    #print("airmass" , X)
    return X


def get_fli(year, month, day, hour, minute):
    # INO location
    elevation = 3600
    #iran
    #location = EarthLocation(lat=33.674*u.deg, lon=51.3188*u.deg , height = elevation*u.meter)
    #eso
    location = EarthLocation(lat=-24.6274*u.deg, lon=-70.4039*u.deg , height = elevation*u.meter)

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




