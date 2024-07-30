#finding chi square

import numpy as np
from scipy.stats import chi2_contingency
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score


FUNCTION_TYPE = ["negative_exponential_func", "logarithmic_func", "exponential_func", "polynomial_func", "power_func", "radical_func"]


v_at_two_kpc = [17.2322,   33.8225,   21.6483,   32.4165,   46.4075,   34.3000,   32.9000,   41.4431,   28.8541,   53.0165, 70.3322,   54.1471,   69.9266,   29.7193,   48.2503,   16.7772,   78.4961,   32.0591,   31.8894,   62.9576, 36.8201,   69.1000,   64.7476,   32.0000,   53.5891,   12.6000,   93.1166,   35.2271,   42.7941,   19.6559, 98.9343,   32.2957,   99.1005,   39.2935,   49.4189,   56.2946,  187.5085,   55.9004,  206.4323,  194.3535, 58.9072,   91.0159,   72.9000,   40.2842,   97.2000,  190.5191,   66.1145,  220.0200,   84.6802,  147.9555, 36.5970,   78.0804,  207.4677,   31.9649,   91.5617,   92.4839,  148.7782,   55.5916,  131.2661,   78.4782, 57.8992,  119.8448,   37.2553,   88.3993,   94.9303,  105.8178,   69.1229,   76.9630,  118.9248,   57.0914, 73.0833,  233.2469,  210.5301,  171.9580,  215.6024,  51.3072,  108.0000,  226.3872,  106.7364,  131.0000, 84.3909,  243.0261,   19.0503,   41.4280,   69.5308,   46.0460,   31.6342,   48.5938,   38.4274,   33.7116, 70.4645,   32.6993,  305.0000,  212.2679,  244.0000,  178.0038,  220.9567,   77.9709,   41.8290,   37.2329, 82.1583,   47.7365,   28.2901,  244.0054,   45.1478,   48.6585,   82.3932,   21.8092,   51.9355,   35.1244, 34.9675,   70.6854,   50.8928,   62.5362,   36.7738,   49.1223,  188.0772,  220.7253,  30.0200,   62.0855, 54.1575,   61.1612,  171.5414,   71.5038,   38.9052,   35.0876,   63.9000,   63.3743,   46.5256,   85.7456, 41.0131,   30.4105,   15.8995,   60.3479,   47.1735,   60.2997,   30.6000,   70.3527,   74.9438,   49.2000, 202.0000,   26.8177,   62.2000,  280.0000,   32.4467,   50.2211,   89.0155,   37.8048,   48.6236,  291.0000, 138.3354,   41.9749,   50.3580,   39.1545,   34.2179]
v_at_two_kpc = np.array(v_at_two_kpc)

all_endpoints = [20.1000,   35.9000,   25.0000,  57.3000,   46.9000,   45.5000,   66.1000,   52.0000,   62.2000,  178.0000,	112.0000,   62.7000,  312.0000,   50.4000,  106.0000,   27.3000,  118.0000,   83.1000,   52.2000,  142.0000, 120.0000,  118.0000,  144.0000,   83.9000,   99.7000,   40.0000,  114.0000,   85.8000,   69.9000,   67.5000, 247.0000,   34.2000,  110.0000,   86.5000,   91.2000,  107.0000,  165.0000,  93.5000,  216.0000,  208.0000, 115.0000,  160.0000,   71.5000,   49.4000,  134.0000,  180.0000,   86.5000,  227.0000,   85.3000,  203.0000,	67.3000,  149.0000,  206.0000,   51.6000,  113.0000,  169.0000,  167.0000,  137.0000,  169.0000,  134.0000, 122.0000,  153.0000,   41.9000,  136.0000,  174.0000,  185.0000,  113.0000,  80.6000,  178.0000,  110.0000, 119.0000,  265.0000,  196.0000,  172.0000,  213.0000,   89.4000,  152.0000,  246.0000,  115.0000,  154.0000,	90.8000,  214.0000,   18.3000,  125.0000,   83.8500,   73.9000,   63.7500,  103.0000,   56.9000,   58.8000, 90.0000,   61.0000,  298.0000,  181.0000,  272.0000,  220.0000,  193.0000,  124.0000,   92.8000,   33.0000,	91.5000,   74.3000,   99.1000,  218.0000,   61.4000,   74.7000,   79.5000,   78.9000,   49.9000,   68.6000, 44.5000,  107.0000,   87.6000,   80.1000,   42.3000,   85.7000,  211.0000,  255.0000,   74.4000,  111.0000,	81.1000,  108.0000,  180.0000,  109.0000,   79.1000,   64.9000,   76.2000,   76.1000,   85.6000,  106.0000, 79.0000,   32.1000,   17.8000,   64.0000,   69.3000,   55.9000,   33.1000,   84.3000,   77.6000,   57.5000,	183.0000,   48.0000,  152.0000,  229.0000,   34.3000,   73.2000,  266.0000,   84.5000,   84.4500,  305.0000,	225.0000,   73.1000,   98.0000,   56.5000,   38.3000]
all_endpoints = np.array(all_endpoints)

fig=plt.figure(figsize=(7,5))
ax=fig.add_subplot(111)
#ax.set_yscale('log')
ax.plot(v_at_two_kpc,all_endpoints,'o', color='red', label='Object velocities')


ax.set_xlabel('velocity at 2kpc')
ax.set_ylabel('velocity at farthest radius')
ax.set_title('comparison of velocities')
ax.legend()


def negative_exponential_func(x, a, b, c):
	return -a * np.exp(-b * x) + c

def exponential_func(x, a, b, c):
	return a * np.exp(b*x) + c


def logarithmic_func(x, a, b):
    return a * np.log(x) + b

def polynomial_func(x, a, b, c):
    return a * x**2 + b * x + c

def power_func(x,a, b, c):
	return a * x ** b + c

def radical_func(x,a, b):
	return a + b*np.sqrt(x)

print("hi1")

def find_best_func(v_at_two_kpc, all_endpoints):
	max_r_squared = 0
	best_fit = ""
	finala = 0
	finalb = 0
	finalc = 0
	for function in FUNCTION_TYPE:
		r_squared= 0
		a = 0
		b= 0
		c = 0
		#please find a better way to do this 
		if function == "negative_exponential_func":
			try:
				params, covariance = curve_fit(negative_exponential_func, v_at_two_kpc, all_endpoints)
				a, b , c= params
				y_pred = negative_exponential_func(v_at_two_kpc, a, b, c)
				r_squared = r2_score(all_endpoints, y_pred)
			except RuntimeError as e:
				pass
		elif function == "logarithmic_func":
			try:
				params, covariance = curve_fit(logarithmic_func, v_at_two_kpc, all_endpoints)
				a, b = params
				y_pred = logarithmic_func(v_at_two_kpc, a, b)
				r_squared = r2_score(all_endpoints, y_pred)
			except RuntimeError as e:
				pass
		elif function == "exponential_func":
			try:
				params, covariance = curve_fit(exponential_func, v_at_two_kpc, all_endpoints)
				a, b, c = params
				y_pred = exponential_func(v_at_two_kpc, a, b,c)
				r_squared = r2_score(all_endpoints, y_pred)
			except RuntimeError as e:
				pass
		elif function == "polynomial_func":
			try:
				params, covariance = curve_fit(polynomial_func, v_at_two_kpc, all_endpoints)
				a, b, c = params
				y_pred = polynomial_func(v_at_two_kpc, a, b,c)
				r_squared = r2_score(all_endpoints, y_pred)
			except RuntimeError as e:
				pass
		elif function == "power_func":
			try:
				params, covariance = curve_fit(power_func, v_at_two_kpc, all_endpoints)
				a, b, c = params
				y_pred = power_func(v_at_two_kpc, a, b,c)
				r_squared = r2_score(all_endpoints, y_pred)
			except RuntimeError as e:
				pass
		elif function == "radical_func":
			try:
				params, covariance = curve_fit(radical_func, v_at_two_kpc, all_endpoints)
				a, b = params
				y_pred = radical_func(v_at_two_kpc, a, b)
				r_squared = r2_score(all_endpoints, y_pred)
			except RuntimeError as e:
				pass
		else:
			print("IDK how you got here man")
		if r_squared > max_r_squared:
			max_r_squared = r_squared
			best_fit = function
			finala = a
			finalb = b
			finalc = c
	return best_fit, max_r_squared, finala, finalb, finalc

print("hi")

best_fit, max_r_squared, finala, finalb, finalc = find_best_func(v_at_two_kpc, all_endpoints)


x_fit = np.linspace(min(v_at_two_kpc), max(v_at_two_kpc), 100)
y_fit = None

if best_fit == "negative_exponential_func":
	y_fit = negative_exponential_func(x_fit, finala,finalb,finalc)
elif best_fit == "logarithmic_func":
	y_fit = logarithmic_func(x_fit, finala, finalb)
elif best_fit == "exponential_func":
	y_fit = exponential_func(x_fit, finala, finalb, finalc)
elif best_fit == "polynomial_func":
	y_fit = polynomial_func(x_fit, finala, finalb, finalc)
elif best_fit == "power_func":
	y_fit = power_func(x_fit, finala, finalb, finalc)
elif best_fit == "radical_func":
	y_fit = radical_func(x_fit, finala, finalb)
else:
	print("I ALSO HAVE NO IDEA HOW YOU GOT HERE")


#plt.scatter(R, V_obs, label='Data points')

#plt.plot(x_fit, y_fit, label=best_fit, color='chartreuse')
#plt.legend()
#plt.show()
#break 

fit = power_func(v_at_two_kpc, finala, finalb, finalc)
arr = np.array([list(v_at_two_kpc), list(fit)])

ax.plot(v_at_two_kpc,all_endpoints,'o', color='red', label='Object velocities')

ax.plot(y_fit,x_fit,label=best_fit)

chi_2, p_value, dof, expected = chi2_contingency(arr)

chi_2 = chi_2/dof

print(chi_2)

print(p_value)

#print(finala)
#print(finalb)
#print(finalc)

plt.show()
