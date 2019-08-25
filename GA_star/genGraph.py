import random
import math

N = int(input())
p = float(input())
esc = float(input())

file = open("pygraph","w")
count = 0
for i in range(0,N):
    ern = random.randint(0,100000)
    if ern > esc*100000:
        continue
    for j in range(0,N):
        if i == j:
            continue
        rn = random.randint(0,100000)
        if rn > p*100000:
            continue
        w = random.randint(1,100)
        file.write("{} {} {}\n".format(i,j,w))
        count+=1 

print(count)
print(N)
