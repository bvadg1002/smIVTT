# TIF stacks recorded on F-Trap can be opened and saved as numpy

import os
import cv2
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import tifffile as tif

matplotlib.use('Qt5Agg')

def fetch_od_and_exposure(file):
    # Fetch OD value from the file name
    od_position_start = file.find("OD") + 2
    od_position_end = od_position_start + file[od_position_start:].find("-")
    od_value_string = file[od_position_start:od_position_end]
    od_value_float = float(od_value_string.replace(",", "."))

    # Fetch exposure value [s] from the file name
    exp_position_start = od_position_end + 1
    exp_position_end = exp_position_start + file[exp_position_start:].find("_")
    exp_value_string = file[exp_position_start:exp_position_end]
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
    for k in range(im_stack.shape[2]):
        image = im_stack[:, :, k]
        bg_values[k] = image.mean()
    background = bg_values.mean()
    return background


def measure_mean_signal_over_time(im_stack):
    # Get dimensions of first frame in TIF stack
    frame_height, frame_width, num_channels = im_stack.shape[:]

    # Calculate local brightness of each pixel in the first frame
    kernel_size = 25
    local_brightness = cv2.blur(im_stack[:, :, 0], (kernel_size, kernel_size))

    # Find the pixel with the maximum local brightness in the first frame
    max_loc = np.unravel_index(np.argmax(local_brightness), (frame_height, frame_width))
    roi_center = max_loc
    roi_radius = 25
    mask = np.zeros((frame_height, frame_width), dtype=np.uint8)
    cv2.circle(mask, roi_center, roi_radius, 255, -1)
    # Calculate mean value in the rectangular ROI around the illumination spot for all images in stack
    mean_signal = np.empty(im_stack.shape[2])
    std_signal = np.empty(im_stack.shape[2])
    for k in range(im_stack.shape[2]):
        image = im_stack[:, :, k]
        roi_im = image[
            roi_center[0] - 10 : roi_center[0] + 10,
            roi_center[1] - 10 : roi_center[1] + 10,
        ]

        mean_signal[k] = roi_im.mean()
        std_signal[k] = roi_im.std()
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


def fetch_timepoint_aliquoted(path):
    time_position_end = path.rfind("min")
    time_position_start = path[:time_position_end].rfind("_")
    time = float(path[time_position_start + 1 : time_position_end])
    return time, time / 60


def fetch_temperature(path):
    temp_position_end = path.rfind(" degrees")
    temp = path[temp_position_end-2:temp_position_end]
    return temp


def fetch_dna_concentration(path):
    conc_position_start = path.rfind("\\") + 1
    conc_position_end = path[conc_position_start:].find('M') + conc_position_start + 1
    conc = path[conc_position_start:conc_position_end]
    return conc


# Define the time [s] elapsed since the sample mixing till the measurement was started
prep_time = 60

# Define the time delay [s] between two frames in a time-resolved data measurement
delay_time = 60

# Create a pandas DataFrame to store all data about measurements
data_combined = pd.DataFrame(
    columns=[
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

# Define the directory containing several measurement folders (TIF files)
root_folder = (
    r"E:\MyFiles\Scripts\smIVTT\Time-resolved transcription tracking - mixed"
)

# Iterate through the subfolders, analyze each of the data files and save the results into a dataframe
for i, (folder_path, folder_name, file_names) in enumerate(os.walk(root_folder)):
    # Skip parental directory
    if folder_path == root_folder:
        continue  # Skip the parent directory

    # For aliquoted measurements proceed with the following reading and analysis
    if "aliquots" in folder_path:
        data_combined.at[i, "Folder"] = folder_path
        data_combined.at[i, "Temperature, °C"] = fetch_temperature(folder_path)
        data_combined.at[i, "DNA concentration"] = fetch_dna_concentration(folder_path)
        data_combined.at[i, "Measurement type"] = "Aliquots"
        data_combined.at[i, "Mean background, a.u."] = 100
        data_combined.at[i, "Mean signal, a.u."] = []
        data_combined.at[i, "STD signal, a.u."] = []
        data_combined.at[i, "Data file path"] = []
        data_combined.at[i, "Time, min"] = []
        data_combined.at[i, "Time, h"] = []

        for file_name in file_names:
            if file_name.startswith("MQ"):
                tif_file_bg = file_name
                full_path_bg = os.path.join(folder_path, tif_file_bg)
                data_combined.at[i, "Background file path"] = full_path_bg
                bg_stack = load_tif_stack(full_path_bg)
                mean_bg = measure_background(bg_stack)
                data_combined.at[i, "Mean background, a.u."] = mean_bg

            if not file_name.startswith("MQ"):
                tif_file_data = file_name
                full_path_data = os.path.join(folder_path, tif_file_data)
                data_combined.at[i, "Data file path"].append(full_path_data)

        for j, file_name in enumerate(data_combined.at[i, "Data file path"]):
            full_path_data = data_combined.at[i, "Data file path"][j]
            data_stack = load_tif_stack(full_path_data)

            od_filter, exposure = fetch_od_and_exposure(full_path_data)
            data_combined.at[i, "OD"] = od_filter
            data_combined.at[i, "Exposure, s"] = exposure
            mean_values, std_values = measure_mean_signal_over_time(data_stack)
            mean_values -= data_combined.at[i, "Mean background, a.u."]
            data_combined.at[i, "Mean signal, a.u."].append(
                adjust_signal(mean_values, od_filter, exposure).mean()
            )
            data_combined.at[i, "STD signal, a.u."].append(
                adjust_signal(std_values, od_filter, exposure).mean()
            )

            time_min, time_h = fetch_timepoint_aliquoted(full_path_data)

            data_combined.at[i, "Time, min"].append(time_min)
            data_combined.at[i, "Time, h"].append(time_h)

        sorted_indices = sorted(
            range(len(data_combined.at[i, "Time, min"])),
            key=lambda k: data_combined.at[i, "Time, min"][k],
        )
        data_combined.at[i, "Mean signal, a.u."] = [
            data_combined.at[i, "Mean signal, a.u."][k] for k in sorted_indices
        ]
        data_combined.at[i, "STD signal, a.u."] = [
            data_combined.at[i, "STD signal, a.u."][k] for k in sorted_indices
        ]
        data_combined.at[i, "Time, min"] = [
            data_combined.at[i, "Time, min"][k] for k in sorted_indices
        ]
        data_combined.at[i, "Time, h"] = [
            data_combined.at[i, "Time, h"][k] for k in sorted_indices
        ]

    # For continuous measurements apply a different processing
    if "aliquots" not in folder_path or "cont" in folder_path:
        for file_name in file_names:
            data_combined.at[i, "Temperature, °C"] = fetch_temperature(folder_path)
            data_combined.at[i, "DNA concentration"] = fetch_dna_concentration(folder_path)
            data_combined.at[i, "Folder"] = folder_path
            if file_name.startswith("MQ"):
                tif_file_bg = file_name
                full_path_bg = os.path.join(folder_path, tif_file_bg)
                data_combined.at[i, "Background file path"] = full_path_bg
            if not file_name.startswith("MQ"):
                tif_file_data = file_name
                full_path_data = os.path.join(folder_path, tif_file_data)
                data_combined.at[i, "Data file path"] = full_path_data


        folder_path = data_combined.at[i, "Folder"]
        full_path_bg = data_combined.at[i, "Background file path"]
        full_path_data = data_combined.at[i, "Data file path"]
        bg_stack = load_tif_stack(full_path_bg)
        data_stack = load_tif_stack(full_path_data)

        mean_bg = measure_background(bg_stack)
        data_combined.at[i, "Mean background, a.u."] = mean_bg
        od_filter, exposure = fetch_od_and_exposure(full_path_data)
        data_combined.at[i, "OD"] = od_filter
        data_combined.at[i, "Exposure, s"] = exposure

        mean_values, std_values = measure_mean_signal_over_time(data_stack)

        mean_values -= mean_bg
        data_combined.at[i, "Mean signal, a.u."] = adjust_signal(
            mean_values, od_filter, exposure
        )
        data_combined.at[i, "STD signal, a.u."] = adjust_signal(
            std_values, od_filter, exposure
        )

        time_min = fetch_timepoints_min(data_stack, exposure, delay_time, prep_time)
        time_h = fetch_timepoints_h(data_stack, exposure, delay_time, prep_time)
        data_combined.at[i, "Time, min"] = time_min
        data_combined.at[i, "Time, h"] = time_h

        last_slash_position = fetch_dna_concentration(full_path_data)
        data_combined.at[i, "Temperature, °C"] = fetch_temperature(full_path_data)

print(data_combined.head())

for i in data_combined.index:
    x = np.array(data_combined.at[i, "Time, h"])
    y = np.array(data_combined.at[i, "Mean signal, a.u."])
    y_err = np.array(data_combined.at[i, "STD signal, a.u."])
    label_temp = data_combined.at[i, "Temperature, °C"]
    label_conc = data_combined.at[i, "DNA concentration"]
    plt.plot(x, y, label=label_conc + " BTx32 DNA at " + label_temp + " °C")
    plt.fill_between(x, y - y_err, y + y_err, alpha=0.2)


plt.title("Transcription over time")
plt.xlabel("Time, h")
plt.ylabel("Fluorescent signal, a.u.")
plt.legend()

plt.show()
