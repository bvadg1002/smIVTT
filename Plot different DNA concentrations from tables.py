import os
import re
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.use("Qt5Agg")

root_folders = [r"E:\MyFiles\Experimental data\smIVTT\20231103 - smIVT HiScribe 1 pM 37 degrees",
        r"E:\MyFiles\Experimental data\smIVTT\20231101 - smIVT HiScribe 10 pM 37 degrees",
        r"E:\MyFiles\Experimental data\smIVTT\20231103 - smIVT HiScribe 100 pM 37 degrees"]
norm = True

pattern = r"\d+\s?pM"
concentrations = [re.search(pattern, path).group() for path in root_folders if re.search(pattern, path)]

if norm == True:
    norm_string = " normalized"
    y_name = "Mean normalized signal, a.u."
    y_err_name = "STD normalized signal, a.u."
else:
    norm_string = ""
    y_name = "Mean signal, a.u."
    y_err_name = "STD signal, a.u."

for i, folder in enumerate(root_folders):
    file = "HiScribe T7 " + concentrations[i] + " BTx32 DNA" + norm_string + ".xlsx"
    print(file)
    data = pd.read_excel(os.path.join(folder, file))
    x = data["Time, min"]
    y = data[y_name]
    y_err = data[y_err_name]
    plt.errorbar(x, y, yerr=y_err, label=concentrations[i] + " BTx32 DNA")
    plt.fill_between(x, y - y_err, y + y_err, alpha=0.2)

plt.xticks(np.arange(0, 660, 60))
plt.xlim(-30, 660)
plt.ylim(0.9, 5.2)
plt.title("Transcription over time" + norm_string)
plt.xlabel("Time, min")
plt.ylabel("Fluorescent signal, a.u.")
plt.legend()

plt.show()

for i, folder in enumerate(root_folders):
    file = "HiScribe T7 " + concentrations[i] + " BTx32 DNA" + norm_string + ".xlsx"
    data = pd.read_excel(os.path.join(folder, file))
    x = data["Time, min"]
    y = data[y_name]
    y_err = data[y_err_name]
    #plt.errorbar(x, y, yerr=y_err, label=concentrations[i] + " BTx32 DNA")
    plt.plot(x, y, label=concentrations[i] + " BTx32 DNA")
    plt.fill_between(x, y - y_err, y + y_err, alpha=0.2)

plt.title("Transcription over time" + norm_string)
plt.xlabel("Time, min")
plt.ylabel("Fluorescent signal, a.u.")
plt.legend()
plt.yscale('log')

plt.show()