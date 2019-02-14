from random import random
from sys import stdout


def gen_rand_point(n, centre, disp=1):
    res = centre[:]
    for i in range(n):
        res[i] += (random() - 0.5) * disp
    return res

def print_cluster(n, centre, number, file=stdout, disp=1):
    for i in range(number):
        print(' '.join(map(str, gen_rand_point(n, centre, disp=disp))), file=file)
        
f = open('data_19_dim_2.txt', 'w')
n = 100
nc = 9
print(300, 2, file=f)
"""
print_cluster(15, [0] * 15 + [1], 50, f)
print_cluster(15, [120] * 15 + [2], 50, f)
print_cluster(15, [240] * 15 + [3], 50, f)
print_cluster(15, [450] * 15 + [4], 50, f)
print_cluster(15, [-340] * 15 + [5], 50, f)
for x in [-4, 0, 4]:
    for y in [-4, 0, 4]:
        print_cluster(2, [x, y], n, f, random() * 1.4 + 0.05)
"""

print_cluster(2, [-0.5, -0.5], 150, f)
print_cluster(2, [0, 0], 150, f)

f.close()
