import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define the linear function with intercept at (0,0)
def linear_func(x, slope):
    return slope * x

# Provided data
data = np.array([
    [5.1, 2, 0, 2, 0, 2],
    [62.9, 5, 3, 4, 3, 2],
    [30.3, 2, 1, 1, 1, 1],
    [12.3, 2, 1, 3, 1, 1],
    [7.2, 1, 0, 1, 0, 1],
    [25.8, 3, 1, 4, 1, 2],
    [10.7, 2, 1, 3, 1, 1],
    [82.4, 5, 1, 10, 0, 4],
    [9.2, 2, 1, 4, 1, 1],
    [18.3, 3, 2, 4, 1, 1],
    [13.2, 3, 0, 4, 0, 3],
    [60.4, 9, 8, 12, 7, 1],
    [76.5, 6, 4, 8, 3, 2],
    [19.4, 1, 1, 1, 1, 0],
    [20.6, 3, 2, 7, 3, 1],
    [37.0, 4, 2, 7, 2, 2],
    [57.9, 7, 6, 6, 2, 1]
])

# Extracting columns
volumes = data[:, 0]
dnatotal = data[:, 1]
dnacount = data[:, 5]

# Perform curve fitting for # DNAtotal
popt_total, pcov_total = curve_fit(linear_func, volumes, dnatotal)
slope_total, slope_error_total = popt_total[0], np.sqrt(np.diag(pcov_total))[0]

# Perform curve fitting for # DNAcount
popt_count, pcov_count = curve_fit(linear_func, volumes, dnacount)
slope_count, slope_error_count = popt_count[0], np.sqrt(np.diag(pcov_count))[0]

# Plotting
plt.figure(figsize=(10, 6))

# Plotting # DNAtotal data and fitted line
plt.scatter(volumes, dnatotal, color='blue', label='# DNAtotal')
plt.plot(volumes, linear_func(volumes, slope_total), color='blue', linestyle='--', label='Fitted Line (# DNAtotal)')

# Plotting # DNAcount data and fitted line
plt.scatter(volumes, dnacount, color='red', label='# DNAcount')
plt.plot(volumes, linear_func(volumes, slope_count), color='red', linestyle='--', label='Fitted Line (# DNAcount)')

plt.xlabel('Volume (um^3)')
plt.ylabel('Number of DNA Molecules')
plt.title('Volume vs. Number of DNA Molecules')
plt.legend()
plt.grid(True)
plt.show()

# Printing fitting equation with slope error
print("Fitting Equation for # DNAtotal: y = x * ({:.3f} ± {:.3f})".format(slope_total, slope_error_total))
print("Fitting Equation for # DNAcount: y = x * ({:.3f} ± {:.3f})".format(slope_count, slope_error_count))
