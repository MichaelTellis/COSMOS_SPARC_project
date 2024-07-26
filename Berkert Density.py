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



#get plots ready
fig, axs = plt.subplots(3,2, figsize=(50,50))

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
 
    '''
    #divide for ratio
    x_ratio = np.divide(x_fit, 1)
    y_ratio = np.divide(y_fit, 1)
    '''
    #graphy graphy and pruney pruney
    if (hubble_type[i] == 7) and (galaxy_name in galaxy_names_dominant_dm):
        axs[0,0].plot(x_fit, y_fit, label="Dark matter densisty for" + galaxy_name)
        #axs[0,0].legend()
        axs[0,0].set_title("Sd Galaxies")
    elif (hubble_type[i] == 8) and (galaxy_name in galaxy_names_dominant_dm):
        axs[0,1].plot(x_fit, y_fit, label="Dark matter densisty for" + galaxy_name)
        #axs[0,1].legend()
        axs[0,1].set_title("Sdm Galaxies")
    elif (hubble_type[i] == 9) and (galaxy_name in galaxy_names_dominant_dm):
        axs[1,0].plot(x_fit, y_fit, label="Dark matter densisty for" + galaxy_name)
        #axs[1,0].legend()
        axs[1,0].set_title("Sm Galaxies")
    elif (hubble_type[i] == 11) and (galaxy_name in galaxy_names_dominant_dm):
        axs[1,1].plot(x_fit, y_fit, label="Dark matter densisty for" + galaxy_name)  
        #axs[1,1].legend()
        axs[1,1].set_title("BCD Galaxies")
    elif (hubble_type[i] == 10) and (galaxy_name in galaxy_names_dominant_dm):
     axs[2,0].plot(x_fit, y_fit,)
     axs[2,0].set_title("Im Galaxies")
    elif (hubble_type[i] == 6) and (galaxy_name in galaxy_names_dominant_dm):
     axs[2,1].plot(x_fit, y_fit)
     axs[2,1].set_title("Scd Galaxies")


#logify axis
for ax in axs.flat:
   ax.set_yscale('log')
   ax.set_xscale('log')
 

#name graphs
axs[2, 0].set_xlabel('Radius kpc')
axs[2, 1].set_xlabel('Radius kpc')
axs[0, 0].set_ylabel(r'Density M$_\odot$/pc$^3$')
axs[1, 0].set_ylabel(r'Density M$_\odot$/pc$^3$')
axs[0, 1].set_ylabel(r'Density M$_\odot$/pc$^3$')
axs[1, 1].set_ylabel(r'Density M$_\odot$/pc$^3$')
axs[2, 0].set_ylabel(r'Density M$_\odot$/pc$^3$')
axs[2, 1].set_ylabel(r'Density M$_\odot$/pc$^3$')
plt.subplots_adjust(hspace=0.4)  # Increase horizontal and vertical space
plt.show()
