import pandas as pd 
import sys
#Burkert table file path
BURKERT_FILE = sys.argv[1]
#### Constants
COLUMN_NAMES = ["Name", "Ydisk - Stellar mass-to-light ratio for disks", "Error on Ydisk", 
"Ybul - Stellar mass-to-light ratio for bulges", "Error on Ybul", "Galaxy distance - mpc", 
"Error on galaxy distance - mpc", "Disk inclination - deg", "Error on disk inclination - deg",
 "v200 - Rotation velocity at r200 - km/s", "Error on V200 - km/s", "c200 - Halo concentration", "Error on c200", 
 "rs - Halo scale radius - kpc", "Error on rs - kpc", "log(rhos) - Volume density of dark matter halos - solMass/pc^3",
 "Error on log(rhos) - solMass/pc^3", "log(M200) - Halo mass - solMass", "e_log(M200) - Error on halo mass- solMass", 
 "Reduced Chi^2"]
#range of bytes for each column
COLSPECS = [(0, 14), (14,20),(20,26),(26,32),(32,38),(38,45),(45,51),(51,57),(57,63),(63,71),(71,78),(78,86),(86,94),(94,102),(102,110),(110,117),(117,124),(124,131),(131,137),(137,144)]
burkert_df = pd.read_fwf(BURKERT_FILE, header = None, names = COLUMN_NAMES, colspecs = COLSPECS)
burkert_df.to_csv("burkert_data.tsv", sep = "\t", index = False)
print("burkert data successfully transformed into tsv file")
