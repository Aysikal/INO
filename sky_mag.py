import datetime, time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pytz
import ephem
from astropy.coordinates import get_body
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
import math 

#date and time inputs:
print("ATTENTION: Time and date entries MUST be UTC")
year = int(input("Enter the year: "))
month = int(input("Enter the month (1-12): "))
day = int(input("Enter the day (1-31): "))
hour = int(input("Enter the hour (0-23): "))
minute = int(input("Enter the minute (0-59): "))

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
print(f"Moon RA (hours): {ra_hours:.6f}")
print(f"Moon DEC (degrees): {dec_degrees:.6f}")

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

print(f"Phase angle: {phase_angle:.1f}°")

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

import math

def calculate_sky_brightness(fli):
    # FLI ranges from 0 to 1 (0 = New Moon, 1 = Full Moon)
    # Assume a constant offset for sky brightness (adjust as needed)
    offset = 22.0  # Example value (mag/arcsec^2)
    if fli ==0: 
          sky_brightness = offset
    # Calculate sky brightness (in magnitudes per square arcsecond)
    else: sky_brightness = offset + 2.5 * math.log10(1-fli)

    return sky_brightness

sky_mag = calculate_sky_brightness(fli)
print("phase: ", phase_name)
print("FLI: ", fli)
print(f"Estimated sky brightness: {sky_mag:.2f} mag/arcsec^2")






