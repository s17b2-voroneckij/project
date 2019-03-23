from os import system
"""
for K in [1, 2, 3, 4]:
    for H in [0, 0.007, 0.1, 0.5, 0, 7, 1, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 3.5, 4, 4.5]:
        f = open('params.txt', 'w')
        print(K, H, file=f, end='')
        f.close()
        print(K, H)
        system('/home/dima/anaconda3/bin/python main.py < data_1_dim_2.txt > /dev/null')
"""

while True:
    inp = input().split()
    f = open('params.txt', 'w')
    f.close()
    system('/home/dima/anaconda3/bin/python main.py %s %s < data_%s_dim_2.txt > /dev/null' % (inp[0], inp[1], inp[2]))

