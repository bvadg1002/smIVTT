# TIF stacks recorded on F-Trap can be opened and saved as numpy

import os
import cv2
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import tifffile as tif

matplotlib.use("Qt5Agg")


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
    for k in range(im_stack.shape[2]):
        image = im_stack[:, :, k]
        roi_im = image[
            roi_center[0] - 25 : roi_center[0] + 25,
            roi_center[1] - 25 : roi_center[1] + 25,
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


def fetch_temperature_aliquoted(path):
    temp_position_end = path.rfind(" degrees")
    temp = float(path[temp_position_end - 2 : temp_position_end])
    return temp


def fetch_dna_concentration(path):
    conc_position_start = path.rfind("\\") + 1
    conc_position_end = path[conc_position_start:].find("M") + conc_position_start + 1
    conc = path[conc_position_start:conc_position_end]
    return conc


def plot_mean_signal(data, path, norm):
    # Plot individual measurements

    parts = data.at[1, "DNA concentration"].split(' ')
    supplier = parts[0]
    concentration = parts[1] + ' ' + parts[2]
    if norm == True:
        norm_name = ' normalized'
        y_name = 'Mean normalized signal, a.u.'
        y_err_name = 'STD normalized signal, a.u.'
    else:
        norm_name = ''
        y_name = 'Mean signal, a.u.'
        y_err_name = 'STD signal, a.u.'

    for i in data.index:
        x = np.array(data.at[i, "Time, min"])
        y = np.array(data.at[i, y_name])
        y_err = np.array(data.at[i, y_err_name])
        color = "orangered"
        plt.plot(x, y, color=color)
        plt.fill_between(x, y - y_err, y + y_err, alpha=0.1, color=color)

    # Prepare dataframe for plotting the mean of all measurements
    combined_mean_values = data[y_name].apply(pd.Series).transpose()

    data_out = pd.DataFrame(
        {
            "Time, min": np.array(data.at[1, "Time, min"]),
            "1st measurement, a.u.": np.array(combined_mean_values[1]),
            "2nd measurement, a.u.": np.array(combined_mean_values[2]),
            "3rd measurement, a.u.": np.array(combined_mean_values[3]),
            str(y_name): combined_mean_values.mean(axis=1),
            str(y_err_name): combined_mean_values.std(axis=1),
        }
    )

    # Plot the mean of all measurements
    x = data_out["Time, min"]
    y = data_out[y_name]
    y_err = data_out[y_err_name]
    label = "Mean" + norm_name + " signal"
    plt.errorbar(
        x,
        y,
        yerr=y_err,
        label=label,
        fmt="o",
        capsize=5,
        color="midnightblue",
        linewidth=1,
        linestyle="None",
    )

    plt.title("Transcription signal over time:  " + supplier + " T7 " + concentration + " BTx32 DNA")
    plt.xlabel("Time, min")
    plt.ylabel("Fluorescent signal, a.u.")
    plt.legend(loc="upper left")

    plt.show()

    path_out = path + "\\" + supplier + " T7 " + concentration + " BTx32 DNA" + norm_name + ".xlsx"
    return data_out, path_out


# Define the time [s] elapsed since the sample mixing till the measurement was started
prep_time = 60

# Define the time delay [s] between two frames in a time-resolved data measurement
delay_time = 60

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

# Define the directory containing several measurement folders (TIF files)
root_folder = r"E:\MyFiles\Experimental data\smIVTT\20231207 - smIVT Takara 10 pM PURE 37 degrees"

# Iterate through the subfolders, analyze each of the data files and save the results into a dataframe
for i, (folder_path, folder_name, file_names) in enumerate(os.walk(root_folder)):
    # Skip parental directory
    if folder_path == root_folder:
        continue  # Skip the parent directory

    # For aliquoted measurements proceed with the following reading and analysis
    if "aliquots" in folder_path:
        data_combined.at[i, "Folder"] = folder_path
        data_combined.at[i, "Temperature, °C"] = fetch_temperature_aliquoted(
            folder_path
        )
        data_combined.at[i, "DNA concentration"] = fetch_dna_concentration(folder_path)
        data_combined.at[i, "Measurement type"] = "Aliquots"
        data_combined.at[i, "Mean background, a.u."] = 100
        data_combined.at[i, "Mean normalized signal, a.u."] = []
        data_combined.at[i, "STD normalized signal, a.u."] = []
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

            data_combined.at[i, "Mean normalized signal, a.u."].append(
                data_combined.at[i, "Mean signal, a.u."][j]
            )

            data_combined.at[i, "STD normalized signal, a.u."].append(
                data_combined.at[i, "STD signal, a.u."][j]
            )

            time_min, time_h = fetch_timepoint_aliquoted(full_path_data)

            data_combined.at[i, "Time, min"].append(time_min)
            data_combined.at[i, "Time, h"].append(time_h)
        sorted_indices = sorted(
            range(len(data_combined.at[i, "Time, min"])),
            key=lambda k: data_combined.at[i, "Time, min"][k],
        )

        time_is_zero_index = int(np.where(np.array(data_combined.at[i, "Time, min"]) == 0.0)[0])
        norm_factor = data_combined.at[i, "Mean signal, a.u."][time_is_zero_index]

        data_combined.at[i, "Mean normalized signal, a.u."] = [
            data_combined.at[i, "Mean normalized signal, a.u."][k] / norm_factor for k in sorted_indices
        ]
        data_combined.at[i, "STD normalized signal, a.u."] = [
            data_combined.at[i, "STD normalized signal, a.u."][k] / norm_factor for k in sorted_indices
        ]
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
saved_data, saved_path = plot_mean_signal(data_combined, root_folder, norm=True)
saved_data.to_excel(saved_path, index=False)

saved_data, saved_path = plot_mean_signal(data_combined, root_folder, norm=False)
saved_data.to_excel(saved_path, index=False)

