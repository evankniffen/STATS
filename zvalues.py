from scipy.stats import norm

# find area of -inf to z using z=1
print(norm.cdf(1))
# find area of -z to z using z=1
print(norm.cdf(1) - norm.cdf(-1))

# find z for area -inf to z using area=0.5
print(norm.ppf(0.5))
# find z for area -z to z using area=0.5
p = 0.5
alpha = 1 - p
print(norm.ppf(1 - alpha / 2))
