import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
import time 
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score



'''
Usage instructions:
requirements:
Must have SPARC data downloaded.
Must have SPARC summary downloaded. Please remove everything except the data at the bottom. 
when ran in terminal: 
arg 1: this file's name (with path)
arg 2: the summary file's name (with path)
arg 3: file path to where your SPARC data is (right click a data file and click get info) 
(Ex: /Users/michaeltellis/Downloads/SPARCData/)
'''
####Constants
start_time = time.time()
SUMMARY_FILE = sys.argv[1]
FILE_PATH = sys.argv[2]
# prune galaxies with R < 8
COLUMN_NAMES = ["Galaxy Name","Hubble Type (1)", "Distance Mpc", "Mean error on D Mpc", 
"Distance Method (2)", "Inclination deg", "Mean error on Inc deg", "Total Luminosity at [3.6]_10+9solLum",
"Mean error on L[3.6]_10+9solLum", "Effective Radius at [3.6] kpc", "Effective Surface Brightness at [3.6]_solLum/pc2", 
"Disk Scale Length at [3.6]_kpc", "Disk Central Surface Brightness at [3.6]_solLum/pc2", "Total HI mass_10+9solMass", 
"HI radius at 1 Msun/pc2_kpc", "Asymptotically Flat Rotation Velocity km/s", "Mean error on Vflat km/s", "Quality Flag (3)",
"References for HI and Ha data (4)"]
HUBBLE_TYPE = {"0":"S0","1":"Sa","2":"Sab","3":"Sb", "4":"Sbc", "5":"Sc",
 "6":"Scd", "7":"Sd", "8":"Sdm", "9":"Sm", "10":"Im", "11":"BCD"}
FUNCTION_TYPE = ["negative_exponential_func", "logarithmic_func", "exponential_func", "polynomial_func", "power_func", "radical_func"]
COLSPECS = [(0,11),(12,14),(14,20),(20,25),(25,27),(27,31),(31,35),(35,42),(42,49),(49,54),(54,62),(62,67),(67,75),(75,82),(82,87),(87,92),(92,97),(97,100),(100,114)]
summary_df = pd.read_fwf(SUMMARY_FILE, header=None,names=COLUMN_NAMES, colspecs=COLSPECS)
galaxy_names = summary_df[COLUMN_NAMES[0]].tolist()



def negative_exponential_func(x, a, b, c):
	return -a * np.exp(-b * x) + c

def exponential_func(x, a, b, c):
	return a * np.exp(b*x) + c


def logarithmic_func(x, a, b):
    return a * np.log(x) + b

def polynomial_func(x, a, b, c):
    return a * x**2 + b * x + c

def power_func(x,a, b, c):
	return a * x ** b + c

def radical_func(x,a, b):
	return a + b*np.sqrt(x)
def find_best_func(R, V_obs):
	max_r_squared = 0
	best_fit = ""
	finala = 0
	finalb = 0
	finalc = 0
	for function in FUNCTION_TYPE:
		r_squared= 0
		a = 0
		b= 0
		c = 0
		#please find a better way to do this 
		if function == "negative_exponential_func":
			try:
				params, covariance = curve_fit(negative_exponential_func, R, V_obs)
				a, b , c= params
				y_pred = negative_exponential_func(R, a, b, c)
				r_squared = r2_score(V_obs, y_pred)
			except RuntimeError as e:
				pass
		elif function == "logarithmic_func":
			try:
				params, covariance = curve_fit(logarithmic_func, R, V_obs)
				a, b = params
				y_pred = logarithmic_func(R, a, b)
				r_squared = r2_score(V_obs, y_pred)
			except RuntimeError as e:
				pass
		elif function == "exponential_func":
			try:
				params, covariance = curve_fit(exponential_func, R, V_obs)
				a, b, c = params
				y_pred = exponential_func(R, a, b,c)
				r_squared = r2_score(V_obs, y_pred)
			except RuntimeError as e:
				pass
		elif function == "polynomial_func":
			try:
				params, covariance = curve_fit(polynomial_func, R, V_obs)
				a, b, c = params
				y_pred = polynomial_func(R, a, b,c)
				r_squared = r2_score(V_obs, y_pred)
			except RuntimeError as e:
				pass
		elif function == "power_func":
			try:
				params, covariance = curve_fit(power_func, R, V_obs)
				a, b, c = params
				y_pred = power_func(R, a, b,c)
				r_squared = r2_score(V_obs, y_pred)
			except RuntimeError as e:
				pass
		elif function == "radical_func":
			try:
				params, covariance = curve_fit(radical_func, R, V_obs)
				a, b = params
				y_pred = radical_func(R, a, b)
				r_squared = r2_score(V_obs, y_pred)
			except RuntimeError as e:
				pass
		else:
			print("IDK how you got here man")
		if r_squared > max_r_squared:
			max_r_squared = r_squared
			best_fit = function
			finala = a
			finalb = b
			finalc = c
	return best_fit, max_r_squared, finala, finalb, finalc


fits = {"best fit": [], "R^2":[], "a":[],"b":[],"c":[]}
for i in range(len(galaxy_names)):
	galaxy_name = galaxy_names[i]
	galaxy_file = FILE_PATH + galaxy_name + "_rotmod.dat"
	R, V_obs, error_V_obs, V_gas, V_disk, V_bulge, SB_disk, SB_bulge = np.loadtxt(galaxy_file, unpack=True)


	last_point = V_obs[len(V_obs)-1]
	best_fit, r_squared, a, b, c = find_best_func(R, V_obs)	
	if best_fit == "negative_exponential_func":
		print("best fit: " + best_fit)
		print("R^2: " + str(r_squared))
		print("a: " + str(-a))
		print("b: " + str(-b))
		print("c: " + str(c))
		fits["best fit"].append(best_fit)
		fits["R^2"].append(r_squared)
		fits["a"].append(-a)
		fits["b"].append(-b)
		fits["c"].append(c)
	else:
		print("best fit: " + best_fit)
		print("R^2: " + str(r_squared))
		print("a: " + str(a))
		print("b: " + str(b))
		print("c: " + str(c))
		fits["best fit"].append(best_fit)
		fits["R^2"].append(r_squared)
		fits["a"].append(a)
		fits["b"].append(b)
		fits["c"].append(c)
fits = pd.DataFrame(fits)
fits.to_csv("fit_functions.tsv", sep = "\t", index = False)
print("Wrote data to tsv file: fit_functions.tsv")
print(fits)
runtime = time.time() - start_time
print("runtime: " + str(runtime))

'''
	x_fit = np.linspace(min(R), max(R), 100)
	y_fit = None

	if best_fit == "negative_exponential_func":
		y_fit = negative_exponential_func(x_fit, a,b,c)
	elif best_fit == "logarithmic_func":
		y_fit = logarithmic_func(x_fit, a, b)
	elif best_fit == "exponential_func":
		y_fit = exponential_func(x_fit, a, b, c)
	elif best_fit == "polynomial_func":
		y_fit = polynomial_func(x_fit, a, b, c)
	elif best_fit == "power_func":
		y_fit = power_func(x_fit, a, b, c)
	elif best_fit == "radical_func":
		y_fit = radical_func(x_fit, a, b)
	else:
		print("I ALSO HAVE NO IDEA HOW YOU GOT HERE")

	plt.scatter(R, V_obs, label='Data points')

	plt.plot(x_fit, y_fit, label=best_fit, color='chartreuse')
	plt.legend()
	plt.show()
	break 
'''

