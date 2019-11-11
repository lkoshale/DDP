import numpy as np 
import scipy.stats as stats


N = int(input())

mu = 100
sigma  = 25
lower, upper = 60, 140

B = stats.truncnorm( (lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma )


x = B.rvs(N)
y = B.rvs(N)

file = open("Cord.txt","w")

for i in range(0,N-1):
    file.write("{} {} \n".format(x[i],y[i]))

file.write("{} {} \n".format(mu,mu))
file.close()