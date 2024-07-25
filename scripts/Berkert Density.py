import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#burkert function
def burkert(r, rs, rho_s):
	DM_density = rho_s/ ((1 + (r / rs)) * (1 + (r / rs)**2))
	return DM_density

#read burkert doc 
with open('parameter_Burkert.mrt', 'rb') as file:
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
with open ('galaxy_sample_summary.txt', 'rb') as file:
	content = file.read()
	str(content)
	content = content.split()
	galaxy_radius_list = []
	hubble_type_list = []
	for i in range(10, len(content), 19):
		galaxy_radius_list.append(content[i].decode('utf-8'))
	for i in range(1, len(content), 19):
		hubble_type_list.append(content[i].decode('utf-8'))

#get plots ready
fig, axs = plt.subplots(2,2, figsize=(50,50))

#go through all galaxies
for i in range(len(galaxy_names_list)):
	galaxy_name = galaxy_names_list[i]
	galaxy_file = galaxy_name + "_rotmod.dat"
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
 
	#divide for ratio
	x_ratio = np.divide(x_fit, 1) #(galaxy_radius[i])
	y_ratio = np.divide(y_fit, 1)
 
	#graphy graphy and pruney pruney
	if (hubble_type[i] == 7) and (galaxy_name!= "IC4202"):
		axs[0,0].plot(x_ratio, y_ratio, label="Dark matter densisty for" + galaxy_name)
		#axs[0,0].legend()
		axs[0,0].set_title("Sd Galaxies")
	elif (hubble_type[i] == 8) and (galaxy_name != "IC4202"):
		axs[0,1].plot(x_ratio, y_ratio, label="Dark matter densisty for" + galaxy_name)
		#axs[0,1].legend()
		axs[0,1].set_title("Sdm Galaxies")
	elif (hubble_type[i] == 9) and (galaxy_name != "IC4202"):
		axs[1,0].plot(x_ratio, y_ratio, label="Dark matter densisty for" + galaxy_name)
		#axs[1,0].legend()
		axs[1,0].set_title("Sm Galaxies")
	elif (hubble_type[i] == 11) and (galaxy_name != "IC4202"):
		axs[1,1].plot(x_ratio, y_ratio, label="Dark matter densisty for" + galaxy_name)  
		#axs[1,1].legend()
		axs[1,1].set_title("BCD Galaxies")
		
  
#logify y-axis
for ax in axs.flat:
   ax.set_yscale('log')
 
#name graphs
axs[1, 0].set_xlabel('R / Rmax')
axs[1, 1].set_xlabel('R / Rmax')
axs[0, 0].set_ylabel('Density SM/pc^3')
axs[1, 0].set_ylabel('Density SM/pc^3')
axs[0, 1].set_ylabel('Density SM/pc^3')
axs[1, 1].set_ylabel('Density SM/pc^3')
plt.show()
