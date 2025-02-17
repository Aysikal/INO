import numpy as np
import os
from astropy.io import fits
import matplotlib.pyplot as plt

# Variables
mode = "low"
dark_folder_path_1min_bin1_low = r"C:\Users\AYSAN\Desktop\project\INO\darks\dark_1min_bin1\Low"
dark_folder_path_1s_bin1_low = r"C:\Users\AYSAN\Desktop\project\INO\darks\dark_1s\Low"
dark_folder_path_2s_bin1_low = r"C:\Users\AYSAN\Desktop\project\INO\darks\dark_2s\Low"
dark_folder_path_4s_bin1_low = r"C:\Users\AYSAN\Desktop\project\INO\darks\dark_4s\Low"
dark_folder_path_5min_bin1_low = r"C:\Users\AYSAN\Desktop\project\INO\darks\dark_5min\Low"
dark_folder_path_5s_bin1_low = r"C:\Users\AYSAN\Desktop\project\INO\darks\dark_5s\Low"
dark_folder_path_10ms_bin1_low = r"C:\Users\AYSAN\Desktop\project\INO\darks\dark_10ms\Low"
dark_folder_path_20ms_bin1_low = r"C:\Users\AYSAN\Desktop\project\INO\darks\dark_20ms\Low"
dark_folder_path_30s_bin1_low = r"C:\Users\AYSAN\Desktop\project\INO\darks\dark_30s\Low"
dark_folder_path_30s_bin1_low_2 = r"C:\Users\AYSAN\Desktop\project\INO\darks\dark_30s_bin1\Low"
dark_folder_path_50ms_bin1_low = r"C:\Users\AYSAN\Desktop\project\INO\darks\dark_50ms\Low"
dark_folder_path_200s_bin1_low = r"C:\Users\AYSAN\Desktop\project\INO\darks\dark_200s\Low"
dark_folder_path_200s_bin1_low_2 = r"C:\Users\AYSAN\Desktop\project\INO\darks\dark_200s_bin1\Low"
dark_folder_path_250ms_bin1_low = r"C:\Users\AYSAN\Desktop\project\INO\darks\dark_250ms\Low"
dark_folder_path_300ms_bin1_low = r"C:\Users\AYSAN\Desktop\project\INO\darks\dark_300ms\Low"

bins = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
dark_exposures = [60, 1, 2, 4, 300, 5, 10 * (10**-3), 20 * (10**-3), 30, 30, 50 * (10**-3), 200, 200, 250 * (10**-3), 300 * (10**-3)]
dark_paths = [dark_folder_path_1min_bin1_low, dark_folder_path_1s_bin1_low, dark_folder_path_2s_bin1_low, dark_folder_path_4s_bin1_low,
              dark_folder_path_5min_bin1_low, dark_folder_path_5s_bin1_low, dark_folder_path_10ms_bin1_low, dark_folder_path_20ms_bin1_low,
              dark_folder_path_30s_bin1_low, dark_folder_path_30s_bin1_low_2, dark_folder_path_50ms_bin1_low, dark_folder_path_200s_bin1_low,
              dark_folder_path_200s_bin1_low_2, dark_folder_path_250ms_bin1_low, dark_folder_path_300ms_bin1_low]
bias_folder_path = r"C:\Users\AYSAN\Desktop\project\INO\bias\Low"

def calculate_rms(arr):
    arr_flat = arr.flatten()
    rms = np.sqrt(np.mean(arr_flat**2))
    return rms

def get_masterdark_and_masterbias(dark_path, bin):
    def open_fits(file_path):
        with fits.open(file_path, ignore_missing_simple=True) as fitsfile:
            file = fitsfile[0].data
            return file


    files = []
    for idx, filename in enumerate(os.listdir(dark_path)):
        if filename.endswith('.fit') or filename.endswith('.fits'):  # Ensure we're only reading FIT/FITS files
            file_path = os.path.join(dark_path, filename)
            data = open_fits(file_path)
            data = np.array(data)
            files.append(data)

    # Stack the 2D arrays along a new axis and calculate the mean across that axis
    stack_files = np.stack(files, axis=0)
    mean_sum_files = np.mean(stack_files, axis=0)

    box = mean_sum_files[int(1000/bin):int(3000/bin), int(1000/bin):int(3000/bin)]
    box_rms = calculate_rms(box)
    box_mean = np.mean(box)

    return files, box, box_rms, box_mean

biases, bias_box, bias_box_rms, bias_box_mean = get_masterdark_and_masterbias(bias_folder_path, 1)

means = []
means_minus_bias = []
RMSs = []
RMSs_minus_bias = []
for i in range(0, len(dark_paths)):
    darks, dark_box, dark_box_mean, dark_box_rms = get_masterdark_and_masterbias(dark_paths[i], bins[i])
    dark_box_minus_bias = dark_box - bias_box
    means.append(dark_box_mean)
    means_minus_bias.append(np.mean(dark_box_minus_bias))
    RMSs.append(dark_box_rms)
    RMSs_minus_bias.append(calculate_rms(dark_box_minus_bias))

# Fit line to mean values
mean_fit = np.polyfit(dark_exposures, means, 1)
mean_fit_line = np.poly1d(mean_fit)

# Fit line to RMS values
rms_fit = np.polyfit(dark_exposures, RMSs, 1)
rms_fit_line = np.poly1d(rms_fit)

# Plot Mean Value vs. Exposure time with fit line
plt.scatter(dark_exposures, means, label='Data')
plt.plot(dark_exposures, mean_fit_line(dark_exposures), color='red', label=f'Fit Line: y = {mean_fit[0]:.2f}x + {mean_fit[1]:.2f}')
plt.title(f"Mean Value vs. Exposure time for {mode} mode")
plt.xlabel("Exposure Time (s)")
plt.ylabel("Mean Value (counts)")
plt.legend()
plt.show()

# Plot RMS vs. Exposure time with fit line
plt.scatter(dark_exposures, RMSs, label='Data')
plt.plot(dark_exposures, rms_fit_line(dark_exposures), color='red', label=f'Fit Line: y = {rms_fit[0]:.2f}x + {rms_fit[1]:.2f}')
plt.title(f"RMS vs. Exposure time for {mode} mode")
plt.xlabel("Exposure Time (s)")
plt.ylabel("RMS")
plt.legend()
plt.show()

# Fit line to mean values (Dark minus Bias)
mean_fit_bias = np.polyfit(dark_exposures, means_minus_bias, 1)
mean_fit_line_bias = np.poly1d(mean_fit_bias)

# Fit line to RMS values (Dark minus Bias)
rms_fit_bias = np.polyfit(dark_exposures, RMSs_minus_bias, 1)
rms_fit_line_bias = np.poly1d(rms_fit_bias)

# Plot Mean Value vs. Exposure time (Dark minus Bias) with fit line
plt.scatter(dark_exposures, means_minus_bias, label='Data')
plt.plot(dark_exposures, mean_fit_line_bias(dark_exposures), color='red', label=f'Fit Line: y = {mean_fit_bias[0]:.2f}x + {mean_fit_bias[1]:.2f}')
plt.title(f"Mean Value vs. Exposure time for {mode} mode (Dark minus Bias)")
plt.xlabel("Exposure Time (s)")
plt.ylabel("Mean Value (counts)")
plt.legend()
plt.show()

# Plot RMS vs. Exposure time (Dark minus Bias) with fit line
plt.scatter(dark_exposures, RMSs_minus_bias, label='Data')
plt.plot(dark_exposures, rms_fit_line_bias(dark_exposures), color='red', label=f'Fit Line: y = {rms_fit_bias[0]:.2f}x + {rms_fit_bias[1]:.2f}')
plt.title(f"RMS vs. Exposure time for {mode} mode (Dark minus Bias)")
plt.xlabel("Exposure Time (s)")
plt.ylabel("RMS")
plt.legend()
plt.show()
