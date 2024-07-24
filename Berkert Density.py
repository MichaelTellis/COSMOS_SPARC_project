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
    for i in range (2, len(content), 19):
        hubble_type_list.append(content[i].decode('utf-8'))
count = 0

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

    #calculate DM density
    x_fit = np.linspace(0,galaxy_radius[i], 100) 
 
    #generate y
    print(list(halo_radius)[i])
    print(list(halo_volume_density)[i])
    y_fit = burkert(x_fit, list(halo_radius)[i], list(halo_volume_density)[i])

    #divide radii for ratio
    x_ratio = np.divide(x_fit, galaxy_radius[i])

    #graphy graphy
    plt.plot(x_ratio, y_fit, label="Dark matter densisty for" + galaxy_name)
    count +=1
    if count > 10: #CHANGE THIS ONCE IT CAN SORT TO GO THROUGH ALL
        break
plt.xlabel("R/Rmax")
plt.ylabel("Density SM/pc^3")
plt.legend()
plt.show()