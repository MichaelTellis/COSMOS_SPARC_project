import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
import time 
from scipy.optimize import curve_fit
import sympy as sp







BURKERT_FITS = sys.argv[1]
SUMMARY_FILE = sys.argv[3]
PATH = sys.argv[2]
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

DM_GALAXY = [
    "DDO064", "DDO154", "DDO170",
    "IC2574", "NGC1705", "NGC2915", "NGC3741",
    "NGC4214", "UGC00731", "UGC00891", "UGC02259",
    "UGC04278", "UGC04325", "UGC05005", "UGC05716", "UGC05721",
    "UGC05750", "UGC05764", "UGC05829", "UGC05918", "UGC05986",
     "UGC06399", "UGC06446", "UGC06667", "UGC06917",
    "UGC06930", "UGC06983", "UGC07151", "UGC07261",
    "UGC07399", "UGC07524", "UGC07603", "UGC07608",
    "UGC07866", "UGC08286", "UGC08490", "UGC08550", "UGC10310",
    "UGC11820", "UGC12632", "UGC12732", "UGCA442",
    "UGCA444"]
kpc2_velocity = [46.4075,34.3000,28.8541,19.6559,72.9000,66.1145,31.9649,76.9630,46.0460,31.6342,70.4645,41.8290,82.1583,28.2901,48.6585,82.3932,21.8092,51.9355,35.1244,34.9675,70.6854,50.8928,62.5362,49.1223,62.0855,61.1612,71.5038,63.9000,63.3743,85.7456,41.0131,60.3479,47.1735,30.6000,70.3527,74.9438,49.2000,50.2211,48.6236,41.9749,50.3580,39.1545,34.2179]
kpc2_velocity = np.array(kpc2_velocity)
outer_velocity = [46.9000,45.5000,62.2000,67.5000,71.5000,86.5000,51.6000,80.6000,73.9000,63.7500,90.0000,92.8000,91.5000,99.1000,74.7000,79.5000,78.9000,49.9000,68.6000,44.5000,107.0000,87.6000,80.1000,85.7000,111.0000,108.0000,109.0000,76.2000,76.1000,106.0000,79.0000,64.0000,69.3000,33.1000,84.3000,77.6000,57.5000,73.2000,84.4500,73.1000,98.0000,56.5000,38.3000]
outer_velocity = np.array(outer_velocity)
kpc2_value = 2



burkert_2kpc = []
burkert_outer = []

x, r_s, C_200, V_200= sp.symbols('x r_s C_200 V_200')
burkert = sp.sqrt(C_200/(x/r_s))*sp.sqrt((0.5*sp.log(1 + (x/r_s)**2) + sp.log(1 + (x/r_s)) - sp.atan(x/r_s))/(0.5*sp.log(1 + C_200**2) + sp.log(1 + C_200) - sp.atan(C_200)))*(V_200)


for galaxy_name in galaxy_names:
    if galaxy_name in DM_GALAXY:
        galaxy_file = PATH + galaxy_name + "_rotmod.dat"
        R, V_obs, error_V_obs, V_gas, V_disk, V_bulge, SB_disk, SB_bulge = np.loadtxt(galaxy_file, unpack=True)
        print(len(R))
        y_value = burkert.subs({x: kpc2_value, r_s: burkert_df.loc[galaxy_name,"rs"], C_200: burkert_df.loc[galaxy_name,"c200"], V_200: burkert_df.loc[galaxy_name,"v200"]})
        burkert_2kpc.append(y_value)

        y_value = burkert.subs({x: np.max(R), r_s: burkert_df.loc[galaxy_name,"rs"], C_200: burkert_df.loc[galaxy_name,"c200"], V_200: burkert_df.loc[galaxy_name,"v200"]})
        burkert_outer.append(y_value)
#find difference vectors
difference_x = burkert_outer - outer_velocity 
difference_y = burkert_2kpc - kpc2_velocity
difference_x = np.array(difference_x, dtype = float)
difference_y = np.array(difference_y, dtype = float)
# find magnitudes
magnitude_diff = np.sqrt(difference_x**2 + difference_y**2)
#find mean values
mean_x_diff = np.mean(difference_x)
mean_y_diff = np.mean(difference_y)
mean_mag_diff = np.mean(magnitude_diff)
print(difference_x)
print("mean x diff: " + str(mean_x_diff))
print("mean y diff: " + str(mean_y_diff))
print("mean magnitude difference: " + str(mean_mag_diff))
#find standard dev 
print("std x diff: " + str(np.std(difference_x)))
print("std y diff: " +  str(np.std(difference_y)))
print("std mag diff: " + str(np.std(magnitude_diff)))
#find median
print("median x diff: " + str(np.median(difference_x)))
print("median y diff: " + str(np.median(difference_y)))
print("median magnitude diff: " + str(np.median(magnitude_diff)))
plt.scatter(kpc2_velocity, outer_velocity, label = "Observed Velocities", marker = "o", s = 150)
plt.scatter(burkert_2kpc, burkert_outer, label = "Burkert Fit Velocities", marker = "*", s = 150)
plt.quiver(kpc2_velocity,outer_velocity, difference_y,difference_x, angles='xy', scale_units='xy', scale=1, headlength = 5, headwidth= 3, width = 0.004, color='r')
plt.title("Velocity at 2kpc and Flat, Burkert vs Empirical", fontsize = 30)
plt.xlabel("Velocity at 2 kpc [km/s]", fontsize = 27)
plt.ylabel("Velocity Flat [km/s]", fontsize = 27)
plt.tick_params(axis='both', labelsize=20) 
plt.legend(fontsize = 25, markerscale = 2)
plt.show()
