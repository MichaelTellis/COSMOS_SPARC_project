import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
import time 

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
arg 4: lower bound of the ML ratio you are checking (int)
arg 5: upper bound of the ML ratio you are checking (int)
arg 6: interval between ML ratio's you are checking (The lower this is, the longer the program will take)
'''
####Constants

start_time = time.time()
SUMMARY_FILE = sys.argv[1]
FILE_PATH = sys.argv[2]
START = int(sys.argv[3])
END = int(sys.argv[4])
INTERVAL = float(sys.argv[5])
COLUMN_NAMES = ["Galaxy Name","Hubble Type (1)", "Distance Mpc", "Mean error on D Mpc", 
"Distance Method (2)", "Inclination deg", "Mean error on Inc deg", "Total Luminosity at [3.6]_10+9solLum",
"Mean error on L[3.6]_10+9solLum", "Effective Radius at [3.6] kpc", "Effective Surface Brightness at [3.6]_solLum/pc2", 
"Disk Scale Length at [3.6]_kpc", "Disk Central Surface Brightness at [3.6]_solLum/pc2", "Total HI mass_10+9solMass", 
"HI radius at 1 Msun/pc2_kpc", "Asymptotically Flat Rotation Velocity km/s", "Mean error on Vflat km/s", "Quality Flag (3)",
"References for HI and Ha data (4)"]
HUBBLE_TYPE = {"0":"S0","1":"Sa","2":"Sab","3":"Sb", "4":"Sbc", "5":"Sc",
 "6":"Scd", "7":"Sd", "8":"Sdm", "9":"Sm", "10":"Im", "11":"BCD"}
COLSPECS = [(0,11),(12,14),(14,20),(20,25),(25,27),(27,31),(31,35),(35,42),(42,49),(49,54),(54,62),(62,67),(67,75),(75,82),(82,87),(87,92),(92,97),(97,100),(100,114)]
summary_df = pd.read_fwf(SUMMARY_FILE, header=None,names=COLUMN_NAMES, colspecs=COLSPECS)
galaxy_names = summary_df[COLUMN_NAMES[0]].tolist()
hubble_types = summary_df[COLUMN_NAMES[1]].tolist()
print(hubble_types)
for i in range(len(hubble_types)):
	hubble_types[i] = HUBBLE_TYPE[str(hubble_types[i])]
print(summary_df)
print(galaxy_names)


'''
R, V_obs, error_V_obs, V_gas, V_disk, V_bulge, SB_disk, SB_bulge = np.loadtxt(galaxy_file, unpack=True)

print(R)

'''

def Chi_test(galaxy_name,file_path, start, end, interval):
	
	galaxy_file = file_path + galaxy_name + "_rotmod.dat"
	R, V_obs, error_V_obs, V_gas, V_disk, V_bulge, SB_disk, SB_bulge = np.loadtxt(galaxy_file, unpack=True)

	V_visible = np.sqrt(V_gas**2 + V_disk**2 + V_bulge**2 )
	#V_DM = np.sqrt(V_obs**2 - V_visible**2)
	GNewton=4.43e-6   # kpc (km/s)2 / solar mass
	V_shell=[(4./3)*np.pi*R[0]**3]
	M_shell=[R[0]*V_obs[0]**2/GNewton]
	M_vis_shell=[ R[0]*V_visible[0]**2/GNewton ]
	R_midpt = [0.5*R[0]]
	for i in range(1,len(R)):
		V_shell.append( (4./3)*np.pi*(R[i]**3 - R[i-1]**3) )
		M_shell.append( R[i]*V_obs[i]**2/GNewton - M_shell[i-1])
		M_vis_shell.append(R[i]*V_visible[i]**2/GNewton - M_vis_shell[i-1] )
		R_midpt.append(0.5*(R[i]+R[i-1]))
	density = np.array(M_shell)/np.array(V_shell)
	vis_density = np.array(M_vis_shell)/np.array(V_shell)
	
	
	ML_ratio = []
	current = start
	while current <= end:
		ML_ratio.append(current)
		current += interval
	minChi = sys.float_info.max
	minML = 0

	for value in ML_ratio:
		V_stars_adjusted = np.sqrt(value)*V_disk
		V_bulge_adjusted = np.sqrt(value)*V_bulge
		V_visible_adjusted = np.sqrt(V_gas**2 + V_stars_adjusted**2 + V_bulge_adjusted**2 )
		#V_DM = np.sqrt(V_obs**2 - V_visible**2)
		chi2=0.
		for i in range(len(R)):
			chi2+=(V_obs[i]-V_visible_adjusted[i])**2/error_V_obs[i]**2
		deg_of_freedom=len(R)-1
		chi2_per_dof=chi2/deg_of_freedom
		if chi2_per_dof < minChi:
			minChi = chi2_per_dof
			minML = value
	print("\n" + galaxy_name)
	print("Chi^2/df: " + str(minChi))
	print("ML ratio: " + str(minML))
	print("data points: " + str(len(R)))
	return minChi, minML, R



chi_values = {"galaxy_name": [], "min_chi^2/df": [], "min_ML_ratio": [], "data_points": [], "galaxy_type": []}
for i in range(len(galaxy_names)):
	name = galaxy_names[i]
	galaxy_type = hubble_types[i]
	minChi, minML, R = Chi_test(name,FILE_PATH, START, END, INTERVAL)
	chi_values["galaxy_name"].append(name)
	chi_values["min_chi^2/df"].append(minChi)
	chi_values["min_ML_ratio"].append(minML)
	chi_values["data_points"].append(len(R))
	chi_values["galaxy_type"].append(galaxy_type)
	
chi_value = pd.DataFrame(chi_values)
chi_value.to_csv("Chi^2_SPARC.tsv", sep="\t", index=False)
print("successfully written to tsv: Chi^2_SPARC.tsv")
runtime = time.time() - start_time
print("runtime: " + str(runtime))
