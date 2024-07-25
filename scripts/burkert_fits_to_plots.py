import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
import time 
from scipy.optimize import curve_fit
start_time = time.time()

# use full path for burkert table
# sparc data's path
'''
Usage Instructions:
arg 1: This file's path and name
arg 2: TSV file which includes the fitted burkert velocity parameters. 
arg 3: txt SPARC galaxy summary file.
optional arg 4: galaxy name for plotting an individual galaxy 
'''
BURKERT_FITS = sys.argv[1]
SPARC_DATA_PATH = sys.argv[2]
SUMMARY_FILE = sys.argv[3]
one_galaxy = False
if len(sys.argv) > 4:
	GALAXY_NAME = sys.argv[4]
	one_galaxy = True

#CONSTANTS
COLUMN_NAMES = ["Galaxy Name","Hubble Type (1)", "Distance Mpc", "Mean error on D Mpc", 
"Distance Method (2)", "Inclination deg", "Mean error on Inc deg", "Total Luminosity at [3.6]_10+9solLum",
"Mean error on L[3.6]_10+9solLum", "Effective Radius at [3.6] kpc", "Effective Surface Brightness at [3.6]_solLum/pc2", 
"Disk Scale Length at [3.6]_kpc", "Disk Central Surface Brightness at [3.6]_solLum/pc2", "Total HI mass_10+9solMass", 
"HI radius at 1 Msun/pc2_kpc", "Asymptotically Flat Rotation Velocity km/s", "Mean error on Vflat km/s", "Quality Flag (3)",
"References for HI and Ha data (4)"]
HUBBLE_TYPE = {"0":"S0","1":"Sa","2":"Sab","3":"Sb", "4":"Sbc", "5":"Sc",
 "6":"Scd", "7":"Sd", "8":"Sdm", "9":"Sm", "10":"Im", "11":"BCD"}
COLSPECS = [(0,11),(12,14),(14,20),(20,25),(25,27),(27,31),(31,35),(35,42),(42,49),(58,62),(58,62),(62,67),(67,75),(75,82),(82,87),(87,92),(92,97),(97,100),(100,114)]
# turns summary txt file to dataframe. 
summary_df = pd.read_fwf(SUMMARY_FILE, header=None,names=COLUMN_NAMES, colspecs=COLSPECS)
#creates list of hubble types, galaxy names, and outer radii of galaxies. 
galaxy_names = summary_df[COLUMN_NAMES[0]].tolist()
hubble_types = summary_df[COLUMN_NAMES[1]].tolist()
galaxy_radii = summary_df[COLUMN_NAMES[9]].tolist()
galaxy_radii = np.array(galaxy_radii, dtype= float)
print(galaxy_radii)
burkert_df = pd.read_csv(BURKERT_FITS, sep="\t",index_col = "galaxy name")
print(burkert_df)


#burkert velocity equation.
def burkert(r, r_s, C_200, V_200):
	x = r / r_s
	return np.sqrt(C_200/x)*np.sqrt((0.5*np.log(1 + x**2) + np.log(1 + x) - np.arctan(x))/(0.5*np.log(1 + C_200**2) + np.log(1 + C_200) - np.arctan(C_200)))*(V_200)



if one_galaxy:
	galaxy_file = SPARC_DATA_PATH + GALAXY_NAME + "_rotmod.dat"
	R, V_obs, error_V_obs, V_gas, V_disk, V_bulge, SB_disk, SB_bulge = np.loadtxt(galaxy_file, unpack=True)

	x_fit = np.linspace(min(R), max(R), 100)
	y_fit = burkert(x_fit,burkert_df.loc[GALAXY_NAME,"rs"],burkert_df.loc[GALAXY_NAME,"c200"],burkert_df.loc[GALAXY_NAME,"v200"])	
	plt.figure(figsize=(10, 6))
	plt.plot(x_fit, y_fit, label='Fitted Burkert Profile', color= "purple")
	plt.scatter(R, V_obs, label='Data points')
	plt.title("Fitted Burkert profile of " + GALAXY_NAME + " compared to observed velocities" )
	plt.xlabel("Radius: kpc")
	plt.ylabel("Rotational velocities: km/s")
	plt.legend()
	plt.show()
else:
	galaxy_names = burkert_df.index.tolist()
	
	print(galaxy_names)

	fig, axs = plt.subplots(2, 2, figsize=(10, 8))
	sm_yfit = []
	sdm_yfit = []
	sd_yfit = []
	BCD_yfit = []
	sm_count = 0
	sdm_count = 0
	sd_count = 0
	BCD_count = 0


	sm_maxr = 0
	sd_maxr = 0
	sdm_maxr = 0
	BCD_maxr = 0

	for i in range(len(galaxy_names)):
		galaxy_name = galaxy_names[i]
		galaxy_type = hubble_types[i]
		galaxy_file = SPARC_DATA_PATH + galaxy_name + "_rotmod.dat"
		R, V_obs, error_V_obs, V_gas, V_disk, V_bulge, SB_disk, SB_bulge = np.loadtxt(galaxy_file, unpack=True)

		x_fit = np.linspace(min(R), max(R), 100)
		y_fit = burkert(x_fit,burkert_df.loc[galaxy_name,"rs"],burkert_df.loc[galaxy_name,"c200"],burkert_df.loc[galaxy_name,"v200"])

		x_fit = np.divide(x_fit, np.max(R))



		if galaxy_type == 7:
			axs[0,0].plot(x_fit, y_fit, label = galaxy_name)
			sm_yfit.append(y_fit)
			sm_count +=1
			if np.max(R) > sm_maxr:
				sm_maxr = np.max(R)
		if galaxy_type == 8:
			axs[0,1].plot(x_fit, y_fit, label = galaxy_name)
			sdm_yfit.append(y_fit)
			sdm_count +=1
			if np.max(R) > sdm_maxr:
				sdm_maxr = np.max(R)
		if galaxy_type == 9:
			axs[1,0].plot(x_fit, y_fit, label = galaxy_name)
			sd_yfit.append(y_fit)
			sd_count +=1
			if np.max(R) > sd_maxr:
				sd_maxr = np.max(R)
		if galaxy_type == 11:
			axs[1,1].plot(x_fit, y_fit, label = galaxy_name)
			BCD_yfit.append(y_fit)
			BCD_count +=1
			if np.max(R) > BCD_maxr:
				BCD_maxr = np.max(R)
	sm_yfit = np.stack(sm_yfit)
	sm_avg = np.mean(sm_yfit, axis=0)
	params, cov = curve_fit(burkert, x_fit, sm_avg)
	sm_r_s, sm_c200, sm_v200 = params
	sm_burkert_fit = burkert(x_fit, sm_r_s, sm_c200, sm_v200)


	sdm_yfit = np.stack(sdm_yfit)
	sdm_avg = np.mean(sdm_yfit, axis=0)
	params, cov = curve_fit(burkert, x_fit, sdm_avg)
	sdm_r_s, sdm_c200, sdm_v200 = params
	sdm_burkert_fit = burkert(x_fit, sdm_r_s, sdm_c200, sdm_v200)


	sd_yfit = np.stack(sd_yfit)
	sd_avg = np.mean(sd_yfit, axis=0)
	params, cov = curve_fit(burkert, x_fit, sd_avg)
	sd_r_s, sd_c200, sd_v200 = params
	sd_burkert_fit = burkert(x_fit, sd_r_s, sd_c200, sd_v200)

	BCD_yfit = np.stack(BCD_yfit)
	BCD_avg = np.mean(BCD_yfit, axis=0)
	params, cov = curve_fit(burkert, x_fit, BCD_avg)
	BCD_r_s, BCD_c200, BCD_v200 = params
	BCD_burkert_fit = burkert(x_fit, BCD_r_s, BCD_c200, BCD_v200)


	# Show the plot
	axs[0,0].plot(x_fit, sm_burkert_fit, linestyle = "--", label = "avg")
	axs[0,1].plot(x_fit, sdm_burkert_fit, linestyle = "--", label = "avg")
	axs[1,0].plot(x_fit, sd_burkert_fit, linestyle = "--", label = "avg")
	axs[1,1].plot(x_fit, BCD_burkert_fit, linestyle = "--", label = "avg")
	axs[0,0].set_xlabel("Radius/Rmax")
	axs[0,1].set_xlabel("Radius/Rmax")
	axs[1,0].set_xlabel("Radius/Rmax")
	axs[1,1].set_xlabel("Radius/Rmax")
	axs[0,0].set_ylabel("Velocity km/s")
	axs[0,1].set_ylabel("Velocity km/s")
	axs[1,0].set_ylabel("Velocity km/s")
	axs[1,1].set_ylabel("Velocity km/s")
	axs[0,0].set_title("sd type galaxy rotational velocity")
	axs[0,1].set_title("sdm type galaxy rotational velocity")
	axs[1,0].set_title("sm type galaxy rotational velocity")
	axs[1,1].set_title("BCD type galaxy rotational velocity")
	'''
	axs[0,0].legend()
	axs[0,1].legend()
	axs[1,0].legend()
	axs[1,1].legend()

	'''


	plt.tight_layout()

	plt.show()
runtime = time.time() - start_time
print("runtime: " + str(runtime))
