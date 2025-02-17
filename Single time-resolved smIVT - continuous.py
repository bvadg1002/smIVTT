# TIF stacks recorded on F-Trap can be opened and saved as numpy

import os
import cv2
import numpy as np
import matplotlib.pyplot as plt
import tifffile as tif


def fetch_od_and_exposure(file_name):
    # Fetch OD value from the file name
    od_position_start = file_name.find("OD") + 2
    od_position_end = od_position_start + file_name[od_position_start:].find("-")
    od_value_string = file_name[od_position_start:od_position_end]
    od_value_float = float(od_value_string.replace(",", "."))

    # Fetch exposure value [s] from the file name
    exp_position_start = od_position_end + 1
    exp_position_end = exp_position_start + file_name[exp_position_start:].find("_")
    exp_value_string = file_name[exp_position_start:exp_position_end]
    exp_value_float = float(exp_value_string.replace(",", "."))

    return od_value_float, exp_value_float


def load_tif_stack(path):
    # Load the TIF file in the stack and transpose to match the openCV indexing (height, width, number)
    im_stack = tif.imread(path)
    im_stack = np.transpose(im_stack, (1, 2, 0))

    # Plot the first TIF image in the stack
    plt.imshow(im_stack[:, :, 0], cmap="gist_heat")
    plt.title(path)
    plt.show()
    return im_stack


def measure_background(im_stack):
    bg_values = np.empty(im_stack.shape[2])
    for i in range(im_stack.shape[2]):
        image = im_stack[:, :, i]
        bg_values[i] = image.mean()
    background = bg_values.mean()
    return background


def measure_mean_signal_over_time(im_stack):
    # Get dimensions of first frame in TIF stack
S    frame_height, frame_width, num_channels = im_stack.shape[:]

    # Calculate local brightness of each pixel in the first frame
    border_size = 25
    kernel_size = 25
    cropped_im_stack = im_stack[border_size:-border_size, border_size:-border_size, 1]
    cropped_height = frame_height - 2 * border_size
    cropped_width = frame_width - 2 * border_size
    local_brightness = cv2.blur(cropped_im_stack, (kernel_size, kernel_size))

    # Find the pixel with the maximum local brightness in the first frame
    max_loc = np.unravel_index(
        np.argmax(local_brightness), (cropped_height, cropped_width)
    )
    roi_center = (max_loc[0] + border_size, max_loc[1] + border_size)

    # Calculate mean value in the rectangular ROI around the illumination spot for all images in stack
    mean_signal = np.empty(im_stack.shape[2])
    std_signal = np.empty(im_stack.shape[2])
    for i in range(im_stack.shape[2]):
        image = im_stack[:, :, i]
        roi_im = image[
            roi_center[0] - 25 : roi_center[0] + 25,
            roi_center[1] - 25 : roi_center[1] + 25,
        ]

        mean_signal[i] = roi_im.mean()
        std_signal[i] = roi_im.std()
    return mean_signal, std_signal


def adjust_signal(mean_signal, od_value, exp_value):
    # Scale to the default exposure value of 0.5 s
    mean_signal *= 0.5 / exp_value

    # Scale to the default OD value of 1.5
    mean_signal *= 10 ** (od_value - 1.5)

    return mean_signal


def fetch_timepoints_min(im_stack, exp_value, delay, preparation):
    step = delay + exp_value
    start = preparation
    stop = preparation + im_stack.shape[2] * step
    time = np.arange(start, stop, step)
    time /= 60
    return time


def fetch_timepoints_h(im_stack, exp_value, delay, preparation):
    step = delay + exp_value
    start = preparation
    stop = preparation + im_stack.shape[2] * step
    time = np.arange(start, stop, step)
    time /= 3600
    return time


def shaded_error_plot(x, y, y_err):
    print(mean_values)
    plt.plot(x, y)
    plt.fill_between(x, y - y_err, y + y_err, alpha=0.2)
    plt.xlabel("Time, h")
    plt.ylabel("Fluorescent signal, a.u.")
    plt.show()
    return


# Define the time [s] elapsed since the sample mixing till the measurement was started
prep_time = 60

# Define the time delay [s] between two frames in a time-resolved data measurement
delay_time = 300

# Define the directory containing the TIF files, the data file name and the background file name
tif_folder = r"E:\MyFiles\Experimental data\smIVTT\20231201 - 4th RT HiT7 measurement"
tif_file_data = (
    "500nM_MB-ATTO647N_10pM_BTx32DNA_T7RNAPolMix_HiScribe_0-69,42h_OD1,5-0,2_gain_50_MMStack_Pos0.ome"
)
tif_file_data += ".tif"
tif_file_bg = "MQ_OD1,5-0,2_gain_50_MMStack_Pos0.ome"
tif_file_bg += ".tif"

full_path_data = os.path.join(tif_folder, tif_file_data)
full_path_bg = os.path.join(tif_folder, tif_file_bg)

# Fetch the measurement metadata
od_filter, exposure = fetch_od_and_exposure(full_path_data)

# Load data and background measurement TIF stacks
data_stack = load_tif_stack(full_path_data)
bg_stack = load_tif_stack(full_path_bg)

# Calculate mean background
mean_bg = measure_background(bg_stack)

# Calculate the fluorescent signal over time with the subtracted background
mean_values, std_values = measure_mean_signal_over_time(data_stack)
mean_values -= mean_bg

# Adjust the fluorescent signal values to the default settings: OD = 1.5, exposure = 200 ms
mean_values = adjust_signal(mean_values, od_filter, exposure)
std_values = adjust_signal(std_values, od_filter, exposure)

# Calculate the "real" time [min or h] of each frame from the onset of the reaction
time_min = fetch_timepoints_min(data_stack, exposure, delay_time, prep_time)
time_h = fetch_timepoints_h(data_stack, exposure, delay_time, prep_time)

# Plot the Fluorescent signal over time
shaded_error_plot(time_h, mean_values, std_values)
