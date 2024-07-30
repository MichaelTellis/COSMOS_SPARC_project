import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import interp1d

BURKERT_PARAMETER = sys.argv[1]
SUMMARY_FILE = sys.argv[3]
PATH = sys.argv[2]
#burkert function
def burkert(r, rs, rho_s):
	DM_density = rho_s/ ((1 + (r / rs)) * (1 + (r / rs)**2))
	return DM_density

#read burkert doc 
with open(BURKERT_PARAMETER, 'rb') as file:
	content = file.read()[2221:]
	str(content)
	content = content.split()
	galaxy_names_list = []
	halo_radius_list = []
	halo_volume_list = []
	for i in range(0, len(content), 20):
		galaxy_names_list.append(content[i].decode('utf-8')) 
	for i in range(14, len(content), 20):
		halo_radius_list.append(content[i].decode('utf-8')) 
	for i in range(15, len(content), 20):
		halo_volume_list.append(content[i].decode('utf-8'))    
#read summary doc(make sure you delete all the text above actual data points)
with open (SUMMARY_FILE, 'rb') as file:
	content = file.read()
	str(content)
	content = content.split()
	galaxy_radius_list = []
	hubble_type_list = []
	for i in range(10, len(content), 19):
		galaxy_radius_list.append(content[i].decode('utf-8'))
	for i in range(1, len(content), 19):
		hubble_type_list.append(content[i].decode('utf-8'))


#make list of DM dominatn galaxies
galaxy_names_dominant_dm =[
	"DDO064", "DDO154", "DDO170",
	"IC2574", "NGC1705", "NGC2915", "NGC3741",
	"NGC4214", "UGC00731", "UGC00891", "UGC02259",
	"UGC04278", "UGC04325", "UGC05005", "UGC05716", "UGC05721",
	"UGC05750", "UGC05764", "UGC05829", "UGC05918", "UGC05986",
	"UGC05999", "UGC06399", "UGC06446", "UGC06667", "UGC06917",
	"UGC06930", "UGC06983", "UGC07151", "UGC07261",
	"UGC07399", "UGC07524", "UGC07603", "UGC07608",
	"UGC07866", "UGC08286", "UGC08490", "UGC08550", "UGC10310",
	"UGC11820", "UGC12632", "UGC12732", "UGCA442",
	"UGCA444"
]

#set other variabels
density_sd = []
density_sdm = []
density_sm = []
density_bcd = []
density_scd = []
density_im = []

max_sd = 0
max_sdm = 0
max_sm = 0
max_bcd = 0
max_scd = 0
max_im = 0



#get plots ready
fig, axs = plt.subplots(3,2, figsize=(50,50))

#go through all galaxies
for i in range(len(galaxy_names_list)):
	galaxy_name = galaxy_names_list[i]
	galaxy_file = PATH + galaxy_name + "_rotmod.dat"
	R, V_obs, error_V_obs, V_gas, V_disk, V_bulge, SB_disk, SB_bulge = np.loadtxt(galaxy_file, unpack=True)
	
	#turn into array
	halo_radius = np.array(halo_radius_list, dtype=float)
	halo_volume_density = np.array(halo_volume_list, dtype=float)
	galaxy_radius = np.array(galaxy_radius_list, dtype=float)
	hubble_type = np.array(hubble_type_list, dtype=float)
	
	#unlogify the density
	halo_volume_density = np.exp(halo_volume_density)

	#generate variables
	x_fit = np.linspace(0,galaxy_radius[i], 100) 
	y_fit = burkert(x_fit, list(halo_radius)[i], list(halo_volume_density)[i])
	x_fit = np.divide(x_fit, galaxy_radius[i])
	#graphy graphy and pruney pruney
	if (hubble_type[i] == 7) and (galaxy_name in galaxy_names_dominant_dm):
		axs[0,0].plot(x_fit, y_fit, label="Dark matter densisty for" + galaxy_name)
		axs[0,0].set_title("Sd Galaxies")
		density_sd.append(y_fit)
		
		
		if list(galaxy_radius)[i] > max_sd:
			max_sd = list(galaxy_radius)[i]
	elif (hubble_type[i] == 8) and (galaxy_name in galaxy_names_dominant_dm):

		axs[0,1].plot(x_fit, y_fit, label="Dark matter densisty for" + galaxy_name)
		axs[0,1].set_title("Sdm Galaxies")
		density_sdm.append(y_fit)
		
		if list(galaxy_radius)[i] > max_sdm:
			max_sdm = list(galaxy_radius)[i]
	elif (hubble_type[i] == 9) and (galaxy_name in galaxy_names_dominant_dm):
		axs[1,0].plot(x_fit, y_fit, label="Dark matter densisty for" + galaxy_name)
		axs[1,0].set_title("Sm Galaxies")
		density_sm.append(y_fit)
		
		if list(galaxy_radius)[i] > max_sm:
			max_sm = list(galaxy_radius)[i]
	elif (hubble_type[i] == 11) and (galaxy_name in galaxy_names_dominant_dm):
		axs[1,1].plot(x_fit, y_fit, label="Dark matter densisty for" + galaxy_name)  
		axs[1,1].set_title("BCD Galaxies")
		density_bcd.append(y_fit)
		
		if list(galaxy_radius)[i] > max_bcd:
			max_bcd = list(galaxy_radius)[i]
	elif (hubble_type[i] == 10) and (galaxy_name in galaxy_names_dominant_dm):
		axs[2,0].plot(x_fit, y_fit,)
		axs[2,0].set_title("Im Galaxies")
		density_im.append(y_fit)
		
		if list(galaxy_radius)[i] > max_im:
			max_im = list(galaxy_radius)[i]
	elif (hubble_type[i] == 6) and (galaxy_name in galaxy_names_dominant_dm):
		axs[2,1].plot(x_fit, y_fit)
		axs[2,1].set_title("Scd Galaxies")
		density_scd.append(y_fit)
		
		if list(galaxy_radius)[i] > max_scd:
			max_scd = list(galaxy_radius)[i]


threshold = 3

def get_avg(all_density, threshold):
	all_density= np.vstack(all_density)
	median = np.nanmedian(all_density, axis=0)
	abs_deviation = np.abs(all_density - median)

	# Calculate the median absolute deviation (MAD)
	mad = np.nanmedian(abs_deviation, axis=0)

	# Define a threshold for outliers
	# Identify outliers
	outliers = abs_deviation > (threshold * mad)
	# Replace outliers with NaN
	all_density[outliers] = np.nan
	average_density = np.nanmean(all_density, axis=0)
	return average_density

print(max_sm)
#logify axis
density_sd= get_avg(density_sd, threshold)
density_sdm= get_avg(density_sdm, threshold)
density_sm= get_avg(density_sm, threshold)
density_bcd= get_avg(density_bcd, threshold)
density_im= get_avg(density_im, threshold)
density_scd= get_avg(density_scd, threshold)




for ax in axs.flat:
   ax.set_yscale('log')
   ax.set_xscale('log')


#average functions

axs[0, 0].plot(x_fit, density_sd, color='black', linestyle='--', label='Average')


axs[0, 1].plot(x_fit, density_sdm, color='black', linestyle='--', label='Average')


axs[1, 0].plot(x_fit, density_sm, color='black', linestyle='--', label='Average')


axs[1, 1].plot(x_fit, density_bcd, color='black', linestyle='--', label='Average')


axs[2, 0].plot(x_fit, density_im, color='black', linestyle='--', label='Average')


axs[2, 1].plot(x_fit, density_scd, color='black', linestyle='--', label='Average')




'''
x_fit = np.linspace(0,max_sdm, 100)
avg_sdm = np.mean(average_density_sdm, axis=0)
axs[0, 1].plot(x_fit, avg_sdm, color='black', linestyle='--', label='Average')

x_fit = np.linspace(0,max_sm, 100)
avg_sm = np.mean(average_density_sm, axis=0)
axs[1, 0].plot(x_fit, avg_sm, color='black', linestyle='--', label='Average')

x_fit = np.linspace(0,max_bcd, 100)
avg_bcd = np.mean(average_density_bcd, axis=0)
axs[1, 1].plot(x_fit, avg_bcd, color='black', linestyle='--', label='Average')

x_fit = np.linspace(0,max_im, 100)
avg_im = np.mean(average_density_im, axis=0)
axs[2, 0].plot(x_fit, avg_im, color='black', linestyle='--', label='Average')

x_fit = np.linspace(0,max_scd, 100)
avg_scd = np.mean(average_density_scd, axis=0)
axs[2, 1].plot(x_fit, avg_scd, color='black', linestyle='--', label='Average')
'''
#name graphs
axs[2, 0].set_xlabel('Radius kpc')
axs[2, 1].set_xlabel('Radius kpc')
axs[0, 0].set_ylabel(r'Density M$_\odot$/pc$^3$')
axs[1, 0].set_ylabel(r'Density M$_\odot$/pc$^3$')
axs[0, 1].set_ylabel(r'Density M$_\odot$/pc$^3$')
axs[1, 1].set_ylabel(r'Density M$_\odot$/pc$^3$')
axs[2, 0].set_ylabel(r'Density M$_\odot$/pc$^3$')
axs[2, 1].set_ylabel(r'Density M$_\odot$/pc$^3$')

#give room on plots
plt.subplots_adjust(hspace=0.4)
plt.show()
