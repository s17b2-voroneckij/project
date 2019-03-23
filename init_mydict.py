import numpy as np
from copy import deepcopy
   
    

V = np.diag(np.load('ru_10000_freq_v.npy'))
S = np.load('ru_10000_freq_s.npy')
vectors = S.dot(V)
raw = open('unique_ru_dict.txt').read().split('\n')
raw.pop()
indices = dict()
words = []
norms = [np.linalg.norm(a) for a in vectors]
for i in range(len(raw)):
    elem = raw[i].split()[0].replace('?', '')
    indices[elem] = i
    words.append(elem)

f = open('mydict5.txt', 'w')
for i in range(len(words)):
    print(words[i], ' '.join(list(map(str, vectors[i][:5]))), file=f)
    if i % 1000 == 0: print(i)
f.close()