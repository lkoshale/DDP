
import random

N = int(input())
E = int(input())

file = open("pygraph","w")

for i in range(0,E):
    a = random.randint(0,N-1)
    b = random.randint(0,N-1)
    c = random.randint(0,100)
    file.write("{} {} {}\n".format(a,b,c))