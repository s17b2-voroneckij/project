from random import random
from sys import stdout


def gen_rand_point(n, centre, disp=1):
    res = centre[:]
    for i in range(n):
        res[i] += (random() - 0.5) * disp
    return res


def print_cluster(n, centre, number, file=stdout, disp=1):
    for i in range(number):
        print(' '.join(map(str, gen_rand_point(n, centre, disp=disp))), (i != 0), number, file=file)


def print_rubbish(n, number, size, file=stdout):
    for i in range(number):
        print(' '.join(map(str, gen_rand_point(n, [0] * n, disp=size))), 'Rubbish', number, file=file)


f = open('data_10_dim_2.txt', 'w')
nc = 9
n = 100
num_rub = 100
print(n * nc + num_rub, 2, file=f)
for x in [-3, 0, 3]:
    for y in [-3, 0, 3]:
        print_cluster(2, [x, y], n, file=f)
print_rubbish(2, num_rub, 7, file=f)
f.close()
