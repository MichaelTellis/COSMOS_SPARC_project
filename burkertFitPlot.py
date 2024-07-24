import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
import time 
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from scipy.stats import chi2_contingency
from scipy.constants import G


'''
Usage instructions:
requirements:
Must have SPARC data downloaded.
Must have SPARC summary downloaded. Please remove everything except the data at the bottom. 
when ran in terminal: 
arg 1: this file's name (with path)
arg 2: the summary file's name (with path)
arg 3: file path to where your SPARC data is (right click a data file and click get info) 
arg 4: file path to the burkert (type tsv, with path)
(Ex: /Users/michaeltellis/Downloads/SPARCData/)
'''
####Constants
start_time = time.time()

SUMMARY_FILE = sys.argv[1]
FILE_PATH = sys.argv[2]

SPARC_BURKERT_PROFILE_PATH = sys.argv[3]
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
sparc_df = pd.read_csv(SPARC_BURKERT_PROFILE_PATH, sep="\t", index_col = ["Name"])
print(sparc_df)
galaxy_names = summary_df[COLUMN_NAMES[0]].tolist()
hubble_types = summary_df[COLUMN_NAMES[1]].tolist()
for i in range(len(hubble_types)):
	hubble_types[i] = HUBBLE_TYPE[str(hubble_types[i])]

def burkert(r, r_s, C_200, V_200):
	x = r / r_s
	return np.sqrt(C_200/x)*np.sqrt((0.5*np.log(1 + x**2) + np.log(1 + x) - np.arctan(x))/(0.5*np.log(1 + C_200**2) + np.log(1 + C_200) - np.arctan(C_200)))*(V_200)


def negative_exponential_func(x, a, b, c):
	return -a * np.exp(-b * x) + c

def exponential_func(x, a, b, c):
	return a * np.exp(b*x) + c


def logarithmic_func(x, a, b):
	return a * np.log(x) + b

def polynomial_func(x, a, b, c):
	return a*(x**2)+(b*x)+c

def power_func(x,a, b, c):
	return a * x ** b + c

def radical_func(x,a, b):
	return a + b*np.sqrt(x)

def find_best_func(R, V_obs):
	least_chi2 = 0
	max_r_squared = 0
	best_fit = ""
	finala = 0
	finalb = 0
	finalc = 0
	for function in FUNCTION_TYPE:
		chi2 = 0
		a = 0
		b= 0
		c = 0
		r2_score = 0
		#please find a better way to do this 
		if function == "negative_exponential_func":
			try:
				params, covariance = curve_fit(negative_exponential_func, R, V_obs)
				a, b , c= params
				y_pred = negative_exponential_func(R, a, b, c)
				chi2, r2_score = get_chi_2(R, V_obs,a,b,c,"negative_exponential_func")
			except RuntimeError as e:
				pass
		elif function == "logarithmic_func":
			try:
				params, covariance = curve_fit(logarithmic_func, R, V_obs)
				a, b = params
				y_pred = logarithmic_func(R, a, b)
				chi2, r2_score = get_chi_2(R, V_obs,a,b,None,"logarithmic_func")
				c = None
			except RuntimeError as e:
				pass
		elif function == "exponential_func":
			try:
				params, covariance = curve_fit(exponential_func, R, V_obs)
				a, b, c = params
				y_pred = exponential_func(R, a, b,c)
				chi2, r2_score = get_chi_2(R, V_obs,a,b,c,"exponential_func")
			except RuntimeError as e:
				pass
		elif function == "polynomial_func":
			try:
				params, covariance = curve_fit(polynomial_func, R, V_obs)
				a, b, c = params
				y_pred = polynomial_func(R, a, b,c)
				chi2, r2_score = get_chi_2(R, V_obs,a,b,c,"polynomial_func")
			except RuntimeError as e:
				pass
		elif function == "power_func":
			try:
				params, covariance = curve_fit(power_func, R, V_obs)
				a, b, c = params
				y_pred = power_func(R, a, b,c)
				chi2, r2_score = get_chi_2(R, V_obs,a,b,c,"power_func")
			except RuntimeError as e:
				pass
		elif function == "radical_func":
			try:
				params, covariance = curve_fit(radical_func, R, V_obs)
				a, b = params
				y_pred = radical_func(R, a, b)
				chi2, r2_score = get_chi_2(R, V_obs,a,b,None,"radical_func")
				c = None
			except RuntimeError as e:
				pass
		else:
			print("IDK how you got here man")
		if (chi2 != None and chi2 < least_chi2) or (r2_score != None and r2_score > max_r_squared):
			least_chi2 = chi2
			max_r_squared = r2_score
			best_fit = function
			finala = a
			finalb = b
			finalc = c
	return best_fit, least_chi2, finala, finalb, finalc, max_r_squared

def get_chi_2(x, V_obs,a, b,c, function):
	chi2PerDf = None
	r2 = None
	if function == "negative_exponential_func":
		try:
			chi_squared_fit = negative_exponential_func(x, a, b, c)
			observed_independence = np.array([list(V_obs),list(chi_squared_fit)])
			chi2, p_value, dof, expected = chi2_contingency(observed_independence)
			chi2PerDf = chi2/dof
		except RuntimeError as e:
			pass
		except ValueError as e:
			y_pred = negative_exponential_func(R, a, b, c)
			r2 = r2_score(V_obs, y_pred)
	elif function == "logarithmic_func":
		try:
			chi_squared_fit = logarithmic_func(x, a, b)
			print(V_obs)
			print(chi_squared_fit)
			observed_independence = np.array([list(V_obs),list(chi_squared_fit)])
			print(observed_independence)
			chi2, p_value, dof, expected = chi2_contingency(observed_independence)
			chi2PerDf = chi2/dof
		except RuntimeError as e:
			pass
		except ValueError as e:
			y_pred = negative_exponential_func(R, a, b, c)
			r2 = r2_score(V_obs, y_pred)
	elif function == "exponential_func":
		try:
			chi_squared_fit = exponential_func(x, a, b, c)
			observed_independence = np.array([list(V_obs),list(chi_squared_fit)])
			chi2, p_value, dof, expected = chi2_contingency(observed_independence)
			chi2PerDf = chi2/dof
		except RuntimeError as e:
			pass
		except ValueError as e:
			y_pred = negative_exponential_func(R, a, b, c)
			r2 = r2_score(V_obs, y_pred)
	elif function == "polynomial_func":
		try:
			chi_squared_fit = polynomial_func(x, a, b, c)
			observed_independence = np.array([list(V_obs),list(chi_squared_fit)])
			chi2, p_value, dof, expected = chi2_contingency(observed_independence)
			chi2PerDf = chi2/dof
		except RuntimeError as e:
			pass
		except ValueError as e:
			y_pred = negative_exponential_func(R, a, b, c)
			r2 = r2_score(V_obs, y_pred)
	elif function == "power_func":
		try:
			chi_squared_fit = power_func(x, a, b, c)
			observed_independence = np.array([list(V_obs),list(chi_squared_fit)])
			chi2, p_value, dof, expected = chi2_contingency(observed_independence)
			chi2PerDf = chi2/dof
		except RuntimeError as e:
			pass
		except ValueError as e:
			y_pred = negative_exponential_func(R, a, b, c)
			r2 = r2_score(V_obs, y_pred)
	elif function == "radical_func":
		try:
			chi_squared_fit = radical_func(x, a, b)
			observed_independence = np.array([list(V_obs),list(chi_squared_fit)])
			chi2, p_value, dof, expected = chi2_contingency(observed_independence)
			chi2PerDf = chi2/dof
		except RuntimeError as e:
			pass
		except ValueError as e:
			y_pred = negative_exponential_func(R, a, b, c)
			r2 = r2_score(V_obs, y_pred)
	elif function == "burkert":
		try:
			chi_squared_fit = burkert(x, a, b, c)
			observed_independence = np.array([list(V_obs),list(chi_squared_fit)])
			chi2, p_value, dof, expected = chi2_contingency(observed_independence)
			chi2PerDf = chi2/dof
		except RuntimeError as e:
			pass
		except ValueError as e:
			y_pred = negative_exponential_func(R, a, b, c)
			r2 = r2_score(V_obs, y_pred)
	return chi2PerDf, r2




fits = {"galaxy name":[], "Chi^2/df":[], "rs":[],"c200":[],"v200":[], "data points": [], "galaxy_type": []}
for i in range(len(galaxy_names)):
	galaxy_name = galaxy_names[i]
	#galaxy_name = "NGC5055"
	galaxy_type = hubble_types[i]
	galaxy_file = FILE_PATH + galaxy_name + "_rotmod.dat"
	R, V_obs, error_V_obs, V_gas, V_disk, V_bulge, SB_disk, SB_bulge = np.loadtxt(galaxy_file, unpack=True)
	print(type(R))

	last_point = V_obs[len(V_obs)-1]
	'''
	best_fit, chi2, a, b, c, r2 = find_best_func(R, V_obs)	
	if best_fit == "negative_exponential_func":
		print("galaxy name: " + str(galaxy_name))
		print("galaxy type: " + galaxy_type)
		print("Data points: " + str(len(R)))
		print("best fit: " + best_fit)
		print("Chi^2/df: " + str(chi2))
		print("r2: " + str(r2))
		print("a: " + str(-a))
		print("b: " + str(-b))
		print("c: " + str(c))
		fits["galaxy name"].append(galaxy_name)
		fits["best fit"].append(best_fit)
		fits["Chi^2/df"].append(chi2)
		fits["a"].append(-a)
		fits["b"].append(-b)
		fits["c"].append(c)
		fits["data points"].append(len(R))
		fits["galaxy_type"].append(galaxy_type)
		fits["r2"].append(r2)

	else:
		print("galaxy name: " + str(galaxy_name))
		print("galaxy type: " + galaxy_type)
		print("Data points: " + str(len(R)))
		print("best fit: " + best_fit)
		print("Chi^2/df: " + str(chi2))
		print("a: " + str(a))
		print("b: " + str(b))
		print("c: " + str(c))
		print("r2: " + str(r2))
		fits["galaxy name"].append(galaxy_name)
		fits["best fit"].append(best_fit)
		fits["Chi^2/df"].append(chi2)
		fits["a"].append(a)
		fits["b"].append(b)
		fits["c"].append(c)
		fits["data points"].append(len(R))
		fits["galaxy_type"].append(galaxy_type)
		fits["r2"].append(r2)
	'''
	fitted_rho_0 = None
	fitted_r_0 = None
	fitted_c200 = None 
	fitted_v200 = None
	burkert_chi2 = None
	fitted_velocities = None
	fitted_r = None
	x_fit = np.linspace(min(R), max(R), 100)
	try: 
		popt, pcov = curve_fit(burkert, R, V_obs)
		fitted_r,fitted_c200, fitted_v200 = popt
		fitted_velocities = burkert(x_fit,fitted_r,fitted_c200 , fitted_v200)
		chi_squared_fit = burkert(R, fitted_r, fitted_c200, fitted_v200)
		observed_independence = np.array([list(V_obs),list(chi_squared_fit)])
		burkert_chi2, p_value, dof, expected = chi2_contingency(observed_independence)
		
		burkert_chi2 = burkert_chi2/dof
	except RuntimeError as e:
		pass
	

	

	#plt.plot(x_fit, fitted_velocities, label='Fitted Burkert Profile', color= "purple")

	print(fitted_velocities)
	print("burkert C200: " + str(fitted_c200))
	print("burkert V200: " +  str(fitted_v200))
	print("burkert r: " + str(fitted_r))
	print("burkert chi squared/df: " + str(burkert_chi2))
	print("SPARC rs: " + str(sparc_df.loc[galaxy_name,"rs - Halo scale radius - kpc"]))
	print("SPARC c200: " + str(sparc_df.loc[galaxy_name,"c200 - Halo concentration"]))
	print("SPARC v200: " + str(sparc_df.loc[galaxy_name,"v200 - Rotation velocity at r200 - km/s"]))
	fits["galaxy name"].append(galaxy_name)
	fits["Chi^2/df"].append(burkert_chi2)
	fits["rs"].append(fitted_r)
	fits["c200"].append(fitted_c200)
	fits["v200"].append(fitted_v200)
	fits["data points"].append(len(R))
	fits["galaxy_type"].append(galaxy_type)

	'''
	y_fit = None

	
	
	y_fit = burkert(x_fit, sparc_df.loc[galaxy_name,"rs - Halo scale radius - kpc"], sparc_df.loc[galaxy_name,"c200 - Halo concentration"],sparc_df.loc[galaxy_name,"v200 - Rotation velocity at r200 - km/s"])
	print(y_fit)
	plt.plot(x_fit, y_fit, label= "Sparc fitted burkert profile", color = "red")
	y_fit = burkert(x_fit, 19.55 ,12.37,176.53)
	plt.plot(x_fit, y_fit, color = "orange")
	



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

	

	plt.plot(x_fit, y_fit, label=best_fit, color='chartreuse')
	
	plt.scatter(R, V_obs, label='Data points')
	plt.title(galaxy_name)
	plt.legend()
	plt.show()
	'''
	



fits = pd.DataFrame(fits)
fits.to_csv("fit_functions.tsv", sep = "\t", index = False)
print("Wrote data to tsv file: fit_functions.tsv")
print(fits)
runtime = time.time() - start_time
print("runtime: " + str(runtime))