import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd



SUMMARY_FILE = sys.argv[1]
PATH = sys.argv[2]


GNewton=4.43e-6   # kpc (km/s)2 / solar mass

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

fig, axs = plt.subplots(7, 25, figsize=(100, 28))
count = 0
for i in range(7):
	for j in range(25):
		galaxy_name = galaxy_names[count]
		R, V_obs, error_V_obs, V_gas, V_disk, V_bulge, SB_disk, SB_bulge = np.loadtxt(PATH + galaxy_name+"_rotmod.dat", unpack=True)
		V_visible = np.sqrt(V_gas**2 + V_disk**2 + V_bulge**2 )
		V_DM = np.sqrt(V_obs**2 - V_visible**2)
		V_shell=[(4./3)*np.pi*R[0]**3]
		M_shell=[R[0]*V_obs[0]**2/GNewton]
		M_vis_shell=[ R[0]*V_visible[0]**2/GNewton ]
		R_midpt = [0.5*R[0]]
		for k in range(1,len(R)):
			V_shell.append( (4./3)*np.pi*(R[k]**3 - R[k-1]**3) )
			M_shell.append( R[k]*V_obs[k]**2/GNewton - M_shell[k-1])
			M_vis_shell.append(R[k]*V_visible[k]**2/GNewton - M_vis_shell[k-1] )
			R_midpt.append(0.5*(R[k]+R[k-1]))
		density = np.array(M_shell)/np.array(V_shell)
		vis_density = np.array(M_vis_shell)/np.array(V_shell)
		ML_ratio= 0.5
		V_stars_adjusted = np.sqrt(ML_ratio)*V_disk
		V_bulge_adjusted = np.sqrt(ML_ratio)*V_bulge
		V_visible_adjusted = np.sqrt(V_gas**2 + V_stars_adjusted**2 + V_bulge_adjusted**2 )
		V_DM = np.sqrt(V_obs**2 - V_visible**2)


		
		axs[i,j].errorbar(R, V_obs , yerr=error_V_obs, marker='o', markersize = 1, color='chartreuse', label='Observed')
		axs[i,j].plot(R, V_gas, '--', color='green', linewidth=0.75, label='Gas')
		axs[i,j].plot(R, V_disk, '--', color='turquoise', linewidth=0.75, label='Disk (M/L=1)')
		axs[i,j].plot(R, V_stars_adjusted, '-', color='turquoise', linewidth=0.75, label='Disk (M/L='+repr(ML_ratio)+')')

		if sum(V_bulge)>0.:
			axs[i,j].plot(R, V_bulge, '--', color='purple', linewidth=0.75, label='Bulge (M/L=1)')
			axs[i,j].plot(R, V_bulge, '-', color='purple', linewidth=0.75, label='Bulge (M/L='+repr(ML_ratio)+')')

		axs[i,j].plot(R, V_visible_adjusted, '-', color='blue', label='Total visible matter (M/L='+repr(ML_ratio)+')')
		axs[i,j].plot(R, V_DM, '-', color='black', label='Inferred dark matter')
		axs[i,j].set_title(galaxy_name, fontsize = 6)
		
		print(galaxy_name)
		count+=1
plt.tight_layout()
plt.show()
	
