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
	 "NGC1705", "NGC2915", "NGC3741",
	"NGC4214", "UGC00731", "UGC00891", "UGC02259",
	"UGC04325", "UGC05005", "UGC05716", "UGC05721",
	"UGC05764", "UGC05829", "UGC05918", "UGC05986",
	"UGC05999", "UGC06399", "UGC06446", "UGC06667", "UGC06917",
	"UGC06930", "UGC06983", "UGC07151", "UGC07261",
	"UGC07399", "UGC07524", "UGC07603", "UGC07608",
	"UGC07866", "UGC08286", "UGC08490", "UGC08550", "UGC10310",
	 "UGC12632", "UGC12732", "UGCA442",
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
	halo_volume_density = np.power(10,halo_volume_density)

	#generate variables
	x_fit = np.linspace(0,galaxy_radius[i], 100) 
	y_fit = burkert(x_fit, list(halo_radius)[i], list(halo_volume_density)[i])
 
	#graphy graphy and pruney pruney
	if (hubble_type[i] == 7) and (galaxy_name in galaxy_names_dominant_dm):
		axs[0,0].plot(x_fit, y_fit, label="Dark matter densisty for" + galaxy_name)
		axs[0,0].set_title("Sd Galaxies", fontsize = 20)
		
		
		if list(galaxy_radius)[i] > max_sd:
			max_sd = list(galaxy_radius)[i]
	elif (hubble_type[i] == 8) and (galaxy_name in galaxy_names_dominant_dm):

		axs[0,1].plot(x_fit, y_fit, label="Dark matter densisty for" + galaxy_name)
		axs[0,1].set_title("Sdm Galaxies", fontsize = 20)
		
		if list(galaxy_radius)[i] > max_sdm:
			max_sdm = list(galaxy_radius)[i]
	elif (hubble_type[i] == 9) and (galaxy_name in galaxy_names_dominant_dm):
		axs[1,0].plot(x_fit, y_fit, label="Dark matter densisty for" + galaxy_name)
		axs[1,0].set_title("Sm Galaxies", fontsize = 20)
		
		if list(galaxy_radius)[i] > max_sm:
			max_sm = list(galaxy_radius)[i]
	elif (hubble_type[i] == 11) and (galaxy_name in galaxy_names_dominant_dm):
		axs[1,1].plot(x_fit, y_fit, label="Dark matter densisty for" + galaxy_name)  
		axs[1,1].set_title("BCD Galaxies", fontsize = 20)
		
		if list(galaxy_radius)[i] > max_bcd:
			max_bcd = list(galaxy_radius)[i]
	elif (hubble_type[i] == 10) and (galaxy_name in galaxy_names_dominant_dm):
		axs[2,0].plot(x_fit, y_fit,)
		axs[2,0].set_title("Im Galaxies", fontsize = 20)
		
		if list(galaxy_radius)[i] > max_im:
			max_im = list(galaxy_radius)[i]
	elif (hubble_type[i] == 6) and (galaxy_name in galaxy_names_dominant_dm):
		axs[2,1].plot(x_fit, y_fit)
		axs[2,1].set_title("Scd Galaxies", fontsize = 20)
		
		if list(galaxy_radius)[i] > max_scd:
			max_scd = list(galaxy_radius)[i]
print(max_sm)
#logify axis

extended_x_sd = np.linspace(0, max_sd, 100)
extended_x_sdm = np.linspace(0, max_sdm, 100)
extended_x_sm = np.linspace(0, max_sm, 100)
extended_x_bcd = np.linspace(0, max_bcd, 100)
extended_x_scd = np.linspace(0, max_scd, 100)
extended_x_im = np.linspace(0, max_im, 100)


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
	func = interp1d(x_fit, y_fit, bounds_error=False, fill_value=np.nan)
	#graphy graphy and pruney pruney
	if (hubble_type[i] == 7) and (galaxy_name in galaxy_names_dominant_dm):
		y_extended = func(extended_x_sd)
		density_sd.append(y_extended)
	
	elif (hubble_type[i] == 8) and (galaxy_name in galaxy_names_dominant_dm):
		y_extended = func(extended_x_sdm)
		density_sdm.append(y_extended)

	elif (hubble_type[i] == 9) and (galaxy_name in galaxy_names_dominant_dm):
		y_extended = func(extended_x_sm)
		density_sm.append(y_extended)

	elif (hubble_type[i] == 11) and (galaxy_name in galaxy_names_dominant_dm):
		y_extended = func(extended_x_bcd)
		density_bcd.append(y_extended)

	elif (hubble_type[i] == 10) and (galaxy_name in galaxy_names_dominant_dm):
		y_extended = func(extended_x_im)
		density_im.append(y_extended)

	elif (hubble_type[i] == 6) and (galaxy_name in galaxy_names_dominant_dm):
		y_extended = func(extended_x_scd)
		density_scd.append(y_extended)


for ax in axs.flat:
   ax.set_yscale('log')
   ax.set_xscale('log')

threshold = 3

def get_avg(max_r, all_density, threshold):
	x_fit = np.linspace(0,max_r, 100)
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
	return average_density, x_fit

#average functions
'''
avg_sd, x_fit_sd = get_avg(max_sd, density_sd, threshold)
axs[0, 0].plot(x_fit_sd, avg_sd, color='black', linestyle='--', label='Average')

avg_sdm, x_fit_sdm = get_avg(max_sdm, density_sdm, threshold)
axs[0, 1].plot(x_fit_sdm, avg_sdm, color='black', linestyle='--', label='Average')

avg_sm, x_fit_sm = get_avg(max_sm, density_sm, threshold)
axs[1, 0].plot(x_fit_sm, avg_sm, color='black', linestyle='--', label='Average')

avg_bcd, x_fit_bcd = get_avg(max_bcd, density_bcd, threshold)
axs[1, 1].plot(x_fit_bcd, avg_bcd, color='black', linestyle='--', label='Average')

avg_im, x_fit_im = get_avg(max_im, density_im, threshold)
axs[2, 0].plot(x_fit_im, avg_im, color='black', linestyle='--', label='Average')

avg_scd, x_fit_scd = get_avg(max_scd, density_scd, threshold)
axs[2, 1].plot(x_fit_scd, avg_scd, color='black', linestyle='--', label='Average')

'''


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
axs[2, 0].set_xlabel('Radius kpc', fontsize = 16)
axs[2, 1].set_xlabel('Radius kpc', fontsize = 16)
axs[0, 0].set_ylabel(r'Density M$_\odot$/pc$^3$', fontsize = 16)
axs[1, 0].set_ylabel(r'Density M$_\odot$/pc$^3$', fontsize = 16)
axs[0, 1].set_ylabel(r'Density M$_\odot$/pc$^3$', fontsize = 16)
axs[1, 1].set_ylabel(r'Density M$_\odot$/pc$^3$', fontsize = 16)
axs[2, 0].set_ylabel(r'Density M$_\odot$/pc$^3$', fontsize = 16)
axs[2, 1].set_ylabel(r'Density M$_\odot$/pc$^3$', fontsize = 16)
#give room on plots
plt.subplots_adjust(hspace=0.4)
plt.show()
