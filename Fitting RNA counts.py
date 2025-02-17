import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define the linear function with intercept at (0,0)
def linear_func(x, slope):
    return slope * x

# Provided data
data = np.array([
    [9.9,    8,    3,    5],
    [21.8,    7,    1,    6],
    [6.1,    7,    4,    3],
    [57.9,    11,    7,    4],
    [7.9,    9,    6,    3],
    [90.0,    23,    7,    16],
    [7.2,    9,    7,    2]

])

# Extracting columns
volumes = data[:, 0]
rnatotal = data[:, 1]
rnacount = data[:, 3]

# Perform curve fitting for # RNAtotal
popt_total, pcov_total = curve_fit(linear_func, volumes, rnatotal)
slope_total, slope_error_total = popt_total[0], np.sqrt(np.diag(pcov_total))[0]

# Perform curve fitting for # RNAcount
popt_count, pcov_count = curve_fit(linear_func, volumes, rnacount)
slope_count, slope_error_count = popt_count[0], np.sqrt(np.diag(pcov_count))[0]

# Plotting
plt.figure(figsize=(10, 6))

# Plotting # RNAtotal data and fitted line
plt.scatter(volumes, rnatotal, color='blue', label='# RNAtotal')
plt.plot(volumes, linear_func(volumes, slope_total), color='blue', linestyle='--', label='Fitted Line (# RNAtotal)')

# Plotting # RNAcount data and fitted line
plt.scatter(volumes, rnacount, color='red', label='# RNAcount')
plt.plot(volumes, linear_func(volumes, slope_count), color='red', linestyle='--', label='Fitted Line (# RNAcount)')

plt.xlabel('Volume (um^3)')
plt.ylabel('Number of mRNA Molecules')
plt.title('Volume vs. Number of mRNA Molecules')
plt.legend()
plt.grid(True)
plt.show()

# Printing fitting equation with slope error
print("Fitting Equation for # RNAtotal: y = x * ({:.3f} ± {:.3f})".format(slope_total, slope_error_total))
print("Fitting Equation for # RNAcount: y = x * ({:.3f} ± {:.3f})".format(slope_count, slope_error_count))
