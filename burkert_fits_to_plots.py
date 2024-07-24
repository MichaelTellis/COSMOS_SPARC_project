import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
import time 

start_time = time.time()

# use full path for burkert table
# sparc data's path

BURKERT_FITS = sys.argv[1]
SPARC_DATA_PATH = sys.argv[2]
GALAXY_NAME = sys.argv[3]

burkert_df = pd.read_csv(BURKERT_FITS, sep="\t",index_col = "galaxy name")
print(burkert_df)

galaxy_file = SPARC_DATA_PATH + GALAXY_NAME + "_rotmod.dat"
R, V_obs, error_V_obs, V_gas, V_disk, V_bulge, SB_disk, SB_bulge = np.loadtxt(galaxy_file, unpack=True)




def burkert(r, r_s, C_200, V_200):
	x = r / r_s
	return np.sqrt(C_200/x)*np.sqrt((0.5*np.log(1 + x**2) + np.log(1 + x) - np.arctan(x))/(0.5*np.log(1 + C_200**2) + np.log(1 + C_200) - np.arctan(C_200)))*(V_200)


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
runtime = time.time() - start_time
print("runtime: " + str(runtime))