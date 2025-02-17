# TIF stacks recorded on F-Trap can be opened and saved as numpy

import os
import cv2
import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tifffile as tif


matplotlib.use("Qt5Agg")

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

    # # Plot the first TIF image in the stack
    # plt.imshow(im_stack[:, :, 0], cmap="gist_heat")
    # plt.title(path)
    # plt.show()
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
    frame_height, frame_width, num_channels = im_stack.shape[:]

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


def shaded_error_plot_combined(data, path, norm):
    if norm == True:
        norm_name = ' normalized'
        y_name = 'Mean normalized signal, a.u.'
        y_err_name = 'STD normalized signal, a.u.'
    else:
        norm_name = ''
        y_name = 'Mean signal, a.u.'
        y_err_name = 'STD signal, a.u.'
    #
    # for i in data.index:
    #     x = np.array(data.at[i, "Time, h"])
    #     y = np.array(data.at[i, y_name])
    #     y_err = np.array(data.at[i, y_err_name])
    #     color = "orangered"
    #     plt.plot(x, y, color=color)
    #     # plt.fill_between(x, y - y_err, y + y_err, alpha=0.1, color=color)

    # Prepare dataframe for plotting the mean of all measurements
    combined_mean_values = data[y_name].apply(pd.Series).transpose()
    data_out = pd.DataFrame(
        {
            "Time, h": np.array(data.at[1, "Time, h"]),
            "1st measurement, a.u.": np.array(combined_mean_values[0]),
            "2nd measurement, a.u.": np.array(combined_mean_values[1]),
            "3rd measurement, a.u.": np.array(combined_mean_values[2]),
            str(y_name): combined_mean_values.mean(axis=1),
            str(y_err_name): combined_mean_values.std(axis=1),
        }
    )

    # Plot the mean of all measurements
    x = data_out["Time, h"]
    y = data_out[y_name]
    y_err = data_out[y_err_name]
    label = "Mean" + norm_name + " signal"
    color = "midnightblue"
    plt.plot(
        x,
        y,
        label=label,
        color=color,
        linewidth=4,
        linestyle="solid",
    )
    plt.fill_between(x, y - y_err, y + y_err, alpha=0.1, color=color)

    plt.title("Transcription signal over time:  HiScribe T7 10 pM BTx32 DNA")
    plt.xlabel("Time, h")
    plt.ylabel("Fluorescent signal, a.u.")
    plt.legend(loc="upper left")
    plt.xlim(-1, 69)

    plt.show()
    path_out = path + "\\HiScribe T7 10 pM BTx32 DNA" + norm_name + ".xlsx"
    data_out.to_excel(path_out, index=False)
    return


# Define the time [s] elapsed since the sample mixing till the measurement was started
prep_time = 60

# Define the time delay [s] between two frames in a time-resolved data measurement
delay_time = 300

# Create a pandas DataFrame to store all data about measurements
data_combined = pd.DataFrame(
    columns=[
        "Mean normalized signal, a.u.",
        "STD normalized signal, a.u.",
        "Mean signal, a.u.",
        "STD signal, a.u.",
        "Time, min",
        "Time, h",
        "Mean background, a.u.",
        "OD",
        "Exposure, s",
        "Interval, s",
        "Prep time, s",
        "Temperature, °C",
        "Molecular Beacon concentration, nM",
        "DNA concentration",
        "Data file path",
        "Background file path",
        "Folder",
        "Measurement type",
    ]
)

root_folders = [r"E:\MyFiles\Experimental data\smIVTT\20231026 - 1st RT HiT7 measurement",
               r"E:\MyFiles\Experimental data\smIVTT\20231103 - 2nd RT HiT7 measurement",
               r"E:\MyFiles\Experimental data\smIVTT\20231117 - 3rd RT HiT7 measurement",
               r"E:\MyFiles\Experimental data\smIVTT\20231201 - 4th RT HiT7 measurement"]

norm = True

for i in range(len(root_folders)):
    data_combined.at[i, "Folder"] = root_folders[i]
    data_combined.at[i, "Temperature, °C"] = 20
    data_combined.at[i, "DNA concentration"] = 10
    data_combined.at[i, "Measurement type"] = "Continuous"
    data_combined.at[i, "Mean background, a.u."] = 100
    data_combined.at[i, "Mean normalized signal, a.u."] = []
    data_combined.at[i, "STD normalized signal, a.u."] = []
    data_combined.at[i, "Mean signal, a.u."] = []
    data_combined.at[i, "STD signal, a.u."] = []
    data_combined.at[i, "Data file path"] = []
    data_combined.at[i, "Time, min"] = []
    data_combined.at[i, "Time, h"] = []

    files = os.listdir(root_folders[i])
    sorted_files = sorted(files, reverse=True)
    for file_name in sorted_files:
        if file_name.startswith("MQ"):
            tif_file_bg = file_name
            full_path_bg = os.path.join(root_folders[i], file_name)
            data_combined.at[i, "Background file path"] = full_path_bg
            bg_stack = load_tif_stack(full_path_bg)
            mean_bg = measure_background(bg_stack)
            data_combined.at[i, "Mean background, a.u."] = mean_bg

        elif file_name.startswith("500nM"):
            tif_file_data = file_name
            full_path_data = os.path.join(root_folders[i], file_name)
            data_stack = load_tif_stack(full_path_data)

            # Fetch the measurement metadata
            od_filter, exposure = fetch_od_and_exposure(full_path_data)

            # Calculate the fluorescent signal over time with the subtracted background
            mean_values, std_values = measure_mean_signal_over_time(data_stack)
            mean_values -= data_combined.at[i, "Mean background, a.u."]

            # Adjust the fluorescent signal values to the default settings: OD = 1.5, exposure = 200 ms
            mean_values = adjust_signal(mean_values, od_filter, exposure)
            std_values = adjust_signal(std_values, od_filter, exposure)

            data_combined.at[i, "Mean normalized signal, a.u."] = mean_values / mean_values[0]
            data_combined.at[i, "STD normalized signal, a.u."] = std_values / mean_values[0]
            data_combined.at[i, "Mean signal, a.u."] = mean_values
            data_combined.at[i, "STD signal, a.u."] = std_values

            # Calculate the "real" time [min or h] of each frame from the onset of the reaction
            time_min = fetch_timepoints_min(data_stack, exposure, delay_time, prep_time)
            time_h = fetch_timepoints_h(data_stack, exposure, delay_time, prep_time)

            data_combined.at[i, "Time, min"] = time_min
            data_combined.at[i, "Time, h"] = time_h

# Plot the Fluorescent signal over time
shaded_error_plot_combined(data_combined, root_folders[2], norm)