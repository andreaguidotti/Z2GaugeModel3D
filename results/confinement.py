"""
Fit routine to extract the string tension in the confinement phase.

STAGE 1:
    For each spatial size (Xs), the potential is fitted as a function of Xt
    using either a linear (a/x + b) or exponential (a*exp(-c*x) + b) model.
    The parameter b represents the asymptotic value of the potential
    for large Xt.

STAGE 2:
    The extracted asymptotic values b(Xs) are then fitted as a function of Xs
    using both a Cornell potential model:

        V(Xs) = a * Xs + b + c / Xs

    and a simple linear model:

        V(Xs) = a * Xs + b

    The string tension is given by the slope parameter 'a'.

Rows containing NaN values are automatically ignored.
The first fit type can be choose (linear or exponential)
from the command line.

Input file format:
    Column 0: beta
    Column 1: Xs (spatial size)
    Column 2: Xt (temporal size)
    Column 3: Wilson loop average
    Column 4: Wilson loop standard deviation
    Column 5: String potential average
    Column 6: String potential standard deviation
"""
import argparse
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Command line interface setup using argparse

parser = argparse.ArgumentParser(
    description=
        ("""
         Fit routine to extract the string tension in the confinement phase.
         Input file format:
           0: beta
           1: Xs (spatial size)
           2: Xt (temporal size)
           3: Wilson loop avg
           4: Wilson loop std
           5: String potential avg
           6: String potential std
         Rows with NaN values are ignored.
         Choose 'linear' or 'exponential' as the first-fit model.
         """),

    formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument(
    "inputFile",
    help="The path to the text file containing the data to be analyzed"
)
parser.add_argument(
    "skipFirst",
    type=int,
    help="Number of initial points (in increasing Xt) to exclude" \
    " from the first fit used to determine the asymptotic value."
)
parser.add_argument(
    "skipSecond",
    type=int,
    help="Number of initial points (in increasing Xs) to exclude" \
    " from the final fit used to extract string tension value."
)
parser.add_argument(
    "--fitType",
    choices=["linear", "exponential"],
    default="linear",
    help="Type of fit to use for the first analysis: linear (a/x + b)"\
    " or exponential (a*exp(-c*Xt)+b)"
)
args = parser.parse_args()

#extract command line arguments
input_file = args.inputFile
skipXt = args.skipFirst

skipXs = args.skipSecond
fitType = args.fitType

#load data and remove row with NaN
data = np.loadtxt(input_file)
data = data[~np.isnan(data).any(axis=1)]

#extract rows from file
allXs = data[:, 1].astype(int)
unique_Xs = np.unique(allXs)  
Xt = data[:, 2]                  
uString = data[:, 5]                  
uStringerr = data[:, 6]               

#first fit functions
def extractLinear(x, a, b):
    return a / x + b

def extractExp(x, a, b, c):
    return a * np.exp(-c * x) + b

# -----First fit to get asymptotic values of string Potential for each Xs-----
colors = cm.tab10(np.linspace(0, 1, len(unique_Xs)))

b_values = []
b_errors = []
valid_Xs = []

plt.figure(figsize=(8,6))
for i, Xs in enumerate(unique_Xs):
    #filter values with desired Xs
    belongToXs = (allXs == Xs)
    #save interested data in new arrays
    xi = Xt[belongToXs]
    yi = uString[belongToXs]
    yerri = uStringerr[belongToXs]

    # filter to exclude first points from fit 
    valid = (xi > skipXt)

    xi = xi[valid]
    yi = yi[valid]
    yerri = yerri[valid]

    if len(xi) < 2:
        print(f"Skipping XS = {Xs}: not enough data" 
              "points excluding first {skipXt} points")
        continue

    #select fitting function from command line argument
    if fitType == "linear":
        fitFunction = extractLinear

    elif fitType == "exponential":
        fitFunction = extractExp
    else:
        raise ValueError("Invalid fit type. Choose 'linear' or 'exponential'.")

    #fit
    popt, pcov = curve_fit(fitFunction, xi, yi, sigma=yerri, absolute_sigma=True)
    perr = np.sqrt(np.diag(pcov))

    a = popt[0]
    b = popt[1]

    b_values.append(b)
    b_errors.append(perr[1])
    valid_Xs.append(Xs)

    # plot 
    plt.errorbar(xi, yi, yerr=yerri, fmt='.', color=colors[i % len(colors)],
                 label=f'XS = {Xs}', capsize=3)
    
    xfit = np.linspace(min(xi)-0.1, max(xi)+0.1, 100)
    plt.plot(xfit, fitFunction(xfit, *popt), '-', color=colors[i % len(colors)])

plt.xlabel('Xt')
plt.ylabel('String Potential')
plt.title(f'Fits y = a/x + b for each XS value (x ≥ {skipXt})')
plt.legend()

plt.show()

#-----Second Fit to extract string tension-----

x2 = np.array(valid_Xs, dtype=int)
b_values = np.array(b_values)
b_errors = np.array(b_errors)

# points to skip in final fit
x2_fit = x2[skipXs:]
b_values_fit = b_values[skipXs:]
b_errors_fit = b_errors[skipXs:]

# defining Cornel potential as fitting function
def Cornell_potential(x, a, b, c):
    return a * x + b + c / x

#fit 
popt_cornell, pcov_cornell = curve_fit(Cornell_potential, x2_fit, b_values_fit,
                               sigma=b_errors_fit, absolute_sigma=True)

a_ext, b_ext, c_ext = popt_cornell
a_ext_err, b_ext_err, c_ext_err = np.sqrt(np.diag(pcov_cornell))

# defining linear function as alternative
def linear_fit(x, a, b):
    return a * x + b

#fit
popt_linear, pcov_linear = curve_fit(linear_fit, x2_fit, b_values_fit,
                               sigma=b_errors_fit, absolute_sigma=True)
a_lin, b_lin = popt_linear
a_lin_err, b_lin_err = np.sqrt(np.diag(pcov_linear))

# reduced Chi squared
def chi2_reduced(y, yfit, yerr, n_params):
    chi2 = np.sum(((y - yfit) / yerr) ** 2)
    dof = len(y) - n_params
    return chi2, chi2 / dof

chi2_cornell, redchi2_cornell = chi2_reduced(b_values_fit, Cornell_potential(x2_fit, *popt_cornell), b_errors_fit, 3)
chi2_linear, redchi2_linear = chi2_reduced(b_values_fit, linear_fit(x2_fit, *popt_linear), b_errors_fit, 2)

# print results
print("\n--- Fit results ---")
print(f"Cornell:  a = {a_ext:.6f} ± {a_ext_err:.6f}, "
      f"b = {b_ext:.6f} ± {b_ext_err:.6f}, "
      f"c = {c_ext:.6f} ± {c_ext_err:.6f}")
print(f"  χ² = {chi2_cornell:.3f},  χ²_red = {redchi2_cornell:.3f}")

print(f"\nLinear: a = {a_lin:.6f} ± {a_lin_err:.6f}, "
      f"b = {b_lin:.6f} ± {b_lin_err:.6f}")
print(f"  χ² = {chi2_linear:.3f},  χ²_red = {redchi2_linear:.3f}")

# plot 
plt.figure(figsize=(8,6))
plt.errorbar(x2, b_values, yerr=b_errors, fmt='.', color='k', capsize=4, label='Data')

xfit = np.linspace(min(x2)-0.5, max(x2)+0.5, 200)
plt.plot(xfit, Cornell_potential(xfit, *popt_cornell), '-', label=f'Cornell\nχ²_red = {redchi2_cornell:.3f}', color='red')
plt.plot(xfit, linear_fit(xfit, *popt_linear), '--', label=f'Linear\nχ²_red = {redchi2_linear:.3f}', color='blue')

plt.xlabel(r'$XS$')
plt.ylabel('String potential')
plt.title(rf'Potential for asymptotic Xt vs $X_S$')
plt.legend()


plt.show()