from astropy.coordinates import EarthLocation, AltAz, get_sun, get_moon, solar_system_ephemeris
from astropy.time import Time
import astropy.units as u
import math 
from astropy.coordinates import Angle
import math

observation_time = '2024-06-20 00:00:00'
RA_hours= 18.0
dec_deg = 30.0 

lat_deg = 51.19
lon_deg = 33.40
# Define your observer's location (INO)
observer_location = EarthLocation(lat_deg,lon_deg)
# Get the current UTC time
utc_time = Time(observation_time)

# Calculate the local sidereal time
lst = utc_time.sidereal_time('apparent', 'greenwich')


# Assuming you have the right ascension of an object (in hours), compute the hour angle
hour_angle = lst - RA_hours * u.hourangle # format : H M S 

def hour_angle_to_degrees(hour_angle):
    hour_angle_str = hour_angle.to_string()
    # Extract hours, minutes, and seconds
    components = hour_angle_str.split('h')[1].split('m')
    hours, minutes, seconds = float(hour_angle_str.split('h')[0]), float(hour_angle_str.split('h')[1].split('m')[0]),float(hour_angle_str.split('h')[1].split('m')[1].split('s')[0])
    
    # Convert to degrees
    total_degrees = (hours * 15) + (minutes * 0.25) + (seconds * 0.00417)
    
    return total_degrees

hour_angle_degrees = hour_angle_to_degrees(hour_angle)

def calculate_zenith_angle(HA_deg, longitude_deg, latitude_deg, declination_deg):
    
    # Calculate local hour angle (LHA)
    LHA = HA_deg + longitude_deg
    
    # Convert latitude and declination to radians
    lat_rad = math.radians(latitude_deg)
    declination_rad = math.radians(declination_deg)
    
    # Compute zenith angle
    cos_Z = math.sin(lat_rad) * math.sin(declination_rad) + math.cos(lat_rad) * math.cos(declination_rad) * math.cos(math.radians(LHA))
    Z = math.degrees(math.acos(cos_Z))
    
    return Z

zenith_angle = calculate_zenith_angle(hour_angle_degrees, lon_deg, lat_deg, dec_deg)
print(f"Zenith angle: {zenith_angle:.2f} degrees")






