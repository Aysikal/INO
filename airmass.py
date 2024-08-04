from astropy import units as u
from astropy.coordinates import SkyCoord , EarthLocation, AltAz
from astropy.time import Time
import numpy as np

#date and time inputs:
print("ATTENTION: Time and date entries MUST be UTC")
year = int(input("Enter the year: "))
month = int(input("Enter the month (1-12): "))
day = int(input("Enter the day (1-31): "))
hour = int(input("Enter the hour (0-23): "))
minute = int(input("Enter the minute (0-59): "))

#ra and dec inputs:
print("RA MUST be in the form of HH:MM:SS")
RA = input("RA: ")
print("DEC MUST be in the form of DD:MM:SS")
DEC = input("DEC: ")

def airmass_function(year, month, day , hour, minute, RA, DEC):
    # INO location
    latitude = '33:40:28'
    longitude = '51:19:08'
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

