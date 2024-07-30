import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
rho_s = [
    57543993.73, 26915348.04, 30902954.33, 3019951.72, 660693448, 
    141253754.5, 22908676.53, 25118864.32, 40738027.78, 14791083.88, 
    199526231.5, 12589254.12, 346736850.5, 6165950.019, 7943282.347, 
    380189396.3, 7762471.166, 223872113.9, 17782794.1, 57543993.73, 
    81283051.62, 27542287.03, 45708818.96, 87096359, 60255958.61, 
    58884365.54, 52480746.02, 72443596.01, 123026877.1, 208929613.1, 
    213796209, 30199517.2, 134896288.3, 47863009.23, 58884365.54, 
    162181009.7, 107151930.5, 87096359, 79432823.47, 1288249.552, 
    39810717.06, 9549925.86, 33113112.15, 54954087.39
]

r_s = [1.85 ,2.65  ,  2.99  ,  15.69 ,  0.78 ,   2.02  ,  2.93   , 3.17  ,  2.8, 5.14  ,  1.62  ,  9.19  ,  1.3 ,10.64 ,  5.23  ,  1.12 ,   7.71   , 0.99   , 3.1 ,1.49  ,  3.5 ,7.14 ,   3.49 ,   2.23 ,   3.09  ,  3.56  ,  3.69 ,   3.33  ,  1.57 ,   1.39,    1.83 ,   3.64 ,   1.51  ,  3.04 ,   0.96  ,  1.76  ,  1.96  ,  1.58  ,  2.04  ,  16.84 ,  2.93  ,  5.41 ,   2.81  ,  1.24]

r_s = np.array(r_s)
r_s_log = np.log10(r_s)
rho_s = np.array(rho_s)
rho_s_log = np.log10(rho_s)
slope, intercept, r_value, p_value, std_err = linregress(r_s_log, rho_s_log)
best_fit_line = slope * r_s_log + intercept
print(slope)
print(intercept)
plt.scatter(r_s_log, rho_s_log, marker = "o", color = "blue", label = "Data Points")
plt.plot(r_s_log, best_fit_line, color = "red", label = "Best Fit Function" )
plt.xlabel("log\u2081\u2080(rₛ) [kpc]", fontsize = 30)
plt.ylabel("log\u2081\u2080(\u03C1s) [M\u2609/kpc\u00B3]", fontsize =30)
plt.title("Relationship Between \u03C1s and rₛ", fontsize = 30)
# Adjust tick label size    
plt.tick_params(axis='both', which='major', labelsize=18)  
plt.tick_params(axis='both', which='minor', labelsize=10)
plt.legend(fontsize=28, markerscale=2)
plt.show()
