import os
import re
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

matplotlib.use("Qt5Agg")

root_folders = [r"E:\MyFiles\Experimental data\smIVTT\20231101 - smIVT HiScribe 10 pM 37 degrees",
        r"E:\MyFiles\Experimental data\smIVTT\20231108 - smIVT NEB 10 pM 37 degrees",
        r"E:\MyFiles\Experimental data\smIVTT\20231108 - smIVT Invitrogen 10 pM 37 degrees",
        r"E:\MyFiles\Experimental data\smIVTT\20231108 - smIVT Roche 10 pM 37 degrees",
        r"E:\MyFiles\Experimental data\smIVTT\20231108 - smIVT ThermoFischer 10 pM 37 degrees",
        r"E:\MyFiles\Experimental data\smIVTT\20231124 - smIVT TranscriptAid 10 pM 37 degrees",
        r"E:\MyFiles\Experimental data\smIVTT\20231124 - smIVT Takara 10 pM 37 degrees",
        r"E:\MyFiles\Experimental data\smIVTT\20231124 - smIVT Promega 10 pM 37 degrees"]

#Plot normalized data (adjusted by the signal at time point 0) or unchanged fluorescent signal
signal_norm = True

#Plot data for the adjusted units of polymerase activity or unchanged fluorescent signal
unit_norm = False

pattern = r"smIVT\s+(\w+)"
suppliers = [re.search(pattern, path).group(1) for path in root_folders if re.search(pattern, path)]

if signal_norm == True:
    signal_norm_string = " normalized"
    y_name = "Mean normalized signal, a.u."
    y_err_name = "STD normalized signal, a.u."
else:
    signal_norm_string = ""
    y_name = "Mean signal, a.u."
    y_err_name = "STD signal, a.u."

if unit_norm == True:
    unit_norm_string = " per T7 activity unit"
    unit_suppliers =  {
        "HiScribe": 15,
        "NEB": 5,
        "Invitrogen": 10,
        "Roche": 2,
        "ThermoFischer": 4,
        "TranscriptAid": 20,
        "Takara": 10,
        "Promega": 4
    }
else:
    unit_norm_string = ""
    unit_suppliers = {
        "HiScribe": 1,
        "NEB": 1,
        "Invitrogen": 1,
        "Roche": 1,
        "ThermoFischer": 1,
        "TranscriptAid": 1,
        "Takara": 1,
        "Promega": 1
    }

color_suppliers =  {
        "HiScribe": "darkorange",
        "NEB": "orange",
        "Invitrogen": "darkorchid",
        "Roche": "mediumblue",
        "ThermoFischer": "red",
        "TranscriptAid": "firebrick",
        "Takara": "dodgerblue",
        "Promega": "gold"
    }

for i, folder in enumerate(root_folders):
    file = suppliers[i] + " T7 10 pM BTx32 DNA" + signal_norm_string + ".xlsx"
    data = pd.read_excel(os.path.join(folder, file))
    x = data["Time, min"]
    y = (data[y_name] - data.at[0,y_name]) / unit_suppliers[suppliers[i]] + data.at[0,y_name]
    y_err = data[y_err_name] / unit_suppliers[suppliers[i]]
    plt.plot(x, y, label=suppliers[i], color=color_suppliers[suppliers[i]])
    plt.fill_between(x, y - y_err, y + y_err, alpha=0.2, color=color_suppliers[suppliers[i]])


plt.title("Transcription over time" + unit_norm_string + signal_norm_string)
plt.xlabel("Time, min")
plt.ylabel("Fluorescent signal, a.u.")
plt.legend()

plt.show()
