import numpy as np
from math import factorial as fact, pi
from sys import setrecursionlimit, argv, stdin, stderr
from random import random, seed
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import


class LSADict:
    def _closest_words(self, i, num, comp=-1):
        if comp == -1:
            comp = len(self.vectors[0])
        min_dists = [[0, 0] for i in range(num)]
        for j in range(len(self.vectors)):
            if self.norms[j] != 0:
                if self.norms[i] == 0 or self.norms[j] == 0:
                    continue
                d = np.dot(self.vectors[i][:comp], self.vectors[j][:comp]) / np.sqrt(
                    np.linalg.norm(self.vectors[i][:comp])) / np.sqrt(np.linalg.norm(self.vectors[j][:comp]))
                if d > min(min_dists)[0]:
                    min_dists.remove(min(min_dists))
                    min_dists.append([d, j])
        return [(a[1], self.words[a[1]]) for a in min_dists]

    def closest_words(self, i, num):
        min_dists = [[0, 0] for i in range(num)]
        for j in range(len(self.vectors)):
            if self.norms[j] != 0:
                if self.norms[i] == 0 or self.norms[j] == 0:
                    continue
                d = np.dot(self.vectors[i], self.vectors[j]) / np.sqrt(self.norms[i]) / np.sqrt(self.norms[j])
                if d > min(min_dists)[0]:
                    min_dists.remove(min(min_dists))
                    min_dists.append([d, j])
        return [(a[1], self.words[a[1]]) for a in min_dists]

    def closest_words2(self, word, num):
        return self.closest_words(self.words.index(word), num)

    def __init__(self, file, dim=10 ** 9):
        self.words = []
        self.vectors = []
        self.mydict = dict()
        self.norms = []
        s = file.readline()
        while s:
            s = s.split()
            if len(s) == 2:
                s = file.readline()
                continue
            word = s[0].split('_')[0]
            vector = np.array(list(map(np.float32, s[1:])))[:dim]
            self.words.append(word)
            self.vectors.append(vector)
            self.mydict[word] = vector
            self.norms.append(np.linalg.norm(vector))
            s = file.readline()


seed(179)
setrecursionlimit(1000000)
COLORS = ['blue', 'red', 'black', 'yellow', 'purple', 'orange', 'green', 'pink', 'violet', 'brown'] + \
         [[[random(), random(), random()]] for i in range(1000)]
g_all = []
g = []
prob = []
points = []
points_copy = []
n = 0
dim = 0
K = 0
H = 0

def dfact(n):
    if n <= 1:
        return 1
    else:
        return n * dfact(n - 2)


def volume(r):
    if dim % 2 == 0:
        k = dim // 2
        return (pi ** k) / fact(k) * (r ** dim)
    else:
        k = dim // 2
        return (2 ** (k + 1)) * (pi ** k) / dfact(dim) * (r ** dim)


def d_k(vec):
    v_sorted = sorted([(g_all[vec[1]][i], i) for i in range(n)])
    for i in range(K + 1):
        g[vec[1]].append(v_sorted[i][1])
    prob[vec[1]] = K / n / volume(v_sorted[K][0])
    return v_sorted[K][0]


class Cluster:
    def __init__(self, i=None):
        self.points = []
        self.max_prob = 0
        self.min_prob = 0
        self.is_deleted = False
        self.is_formed = False
        self.is_meaningful = (self.max_prob - self.min_prob) >= H
        if i is not None:
            self.add_point(i)

    def add_point(self, i):
        if len(self.points) == 0:
            self.max_prob = prob[i]
            self.min_prob = prob[i]
        else:
            self.max_prob = max(self.max_prob, prob[i])
            self.min_prob = min(self.min_prob, prob[i])
        self.points.append(i)
        self.is_meaningful = (self.max_prob - self.min_prob) >= H


def union_clusters(a, b, cluster, clusters):
    for point in clusters[b].points:
        cluster[point] = a
        clusters[a].add_point(point)
    clusters[b].points = []
    clusters[b].is_deleted = True


def Wishart():
    global g_all, g, prob
    prob = [0] * n
    g_all = [[np.float16(np.linalg.norm(points[i][0] - points[j][0])) for i in range(n)] for j in range(n)]
    g = [[] for i in range(n)]
    points.sort(key=d_k)
    clusters = [Cluster()]
    cluster = [-1] * n
    for i in range(n):
        connected_clusters = set()
        v = points[i][1]
        for u in g[v]:
            if cluster[u] != -1:
                connected_clusters.add(cluster[u])
        if not connected_clusters:
            clusters.append(Cluster(v))
            cluster[v] = len(clusters) - 1
        else:
            connected_clusters = list(connected_clusters)
            if len(connected_clusters) == 1:
                if clusters[connected_clusters[0]].is_formed:
                    clusters[0].add_point(v)
                    cluster[v] = 0
                else:
                    clusters[connected_clusters[0]].add_point(v)
                    cluster[v] = connected_clusters[0]
            else: # 3.3
                if min([clust.is_formed for clust in clusters]): # 3.3.1
                    cluster[v] = 0
                    clusters[0].add_point(v)
                else:
                    z_h = 0
                    for clust_i in connected_clusters:
                        if clusters[clust_i].is_meaningful: z_h += 1
                    if z_h > 1 or 0 in connected_clusters:
                        cluster[v] = 0
                        clusters[0].add_point(v)
                        for clust in connected_clusters:
                            if clusters[clust].is_meaningful:
                                clusters[clust].is_formed = True
                            elif clust:
                                union_clusters(0, clust, cluster, clusters)
                    else:
                        connected_clusters.sort()
                        clusters[connected_clusters[0]].add_point(v)
                        cluster[v] = connected_clusters[0]
                        for clust in connected_clusters[1:]:
                            union_clusters(connected_clusters[0], clust, cluster, clusters)
    return clusters, cluster, g


def print_clusters(clusters):
    nc = 0
    for cl in clusters:
        if not cl.is_deleted:
            print("Cluster %s" % str(nc), "p_min:", cl.min_prob, "p_max", cl.max_prob)
            nc += 1
            for i in cl.points:
                print(points_copy[i])
            print()
            print()


def plot(clusters, a=0, b=1):
    nc = 0
    for cl in clusters:
        if not cl.is_deleted:
            for i in cl.points:
                curr = points_copy[i]
                plt.scatter(curr[0][a], curr[0][b], color=COLORS[nc])
            nc += 1
    plt.title("K = %s, H = %s" % (str(K), str(H)))
    plt.show()


def plot3D(clusters, dim):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for a in range(dim):
        for b in range(a + 1, dim):
            for c in range(b + 1, dim):
                nc = 0
                for cl in clusters:
                    if not cl.is_deleted:
                        for i in cl.points:
                            curr = points_copy[i]
                            ax.scatter(curr[0][a], curr[0][b], curr[0][c], color=COLORS[nc])
                        nc += 1
    plt.title("K = %s, H = %s" % (str(K), str(H)))
    plt.show()


def plot_dim(clusters, dim):
    for a in range(dim - 1):
        for b in range(a + 1, dim):
            plot(clusters, a, b)


def quality_honest(clusters, cluster):
    all_points = set([i for i in range(n)])
    for i in range(len(g)):
        g[i] = set(g[i])
    for i in range(len(clusters)):
        clusters[i].points = set(clusters[i].points)
    try:
        glob_intern_density_up = 0
        glob_intern_density_down = 0
        for cl in clusters:
            if not cl.is_deleted:
                for v in cl.points:
                    for u in cl.points:
                        glob_intern_density_up += (u in g[v])
                glob_intern_density_down += (len(cl.points)) ** 2
                if len(cl.points) == 0:
                    cl.is_deleted = True

        glob_extern_density_up = 0
        glob_extern_density_down = 0
        for cl_i in range(len(clusters)):
            if not clusters[cl_i].is_deleted:
                for cl_j in range(len(clusters)):
                    if cl_i != cl_j and not clusters[cl_j].is_deleted:
                        for v in clusters[cl_i].points:
                            for u in clusters[cl_j].points:
                                glob_extern_density_up += (u in g[v])
                glob_extern_density_down += len(clusters[cl_i].points) * (n - len(clusters[cl_i].points))
        res_glob_dens = 0.5 * (1 + glob_intern_density_up / glob_intern_density_down - glob_extern_density_up / glob_extern_density_down)
    except:
        res_glob_dens = -1
    try:
        loc_dens = 0
        for cl_i in range(len(clusters)):
            if not clusters[cl_i].is_deleted:
                intern = 0
                extern = 0
                for v in clusters[cl_i].points:
                    intern += len(g[v] & clusters[cl_i].points)
                    extern += len(g[v] & (all_points - clusters[cl_i].points))
                intern /= len(clusters[cl_i].points) ** 2
                extern /= len(clusters[cl_i].points) * (n - len(clusters[cl_i].points))
                loc_dens += len(clusters[cl_i].points) / (2 * n) * (intern - extern + 1)
        res_loc_dens = loc_dens
    except:
        res_loc_dens = -1

    adj = np.array([[int(u in g[v]) for u in range(n)] for v in range(n)])
    A_v = np.array([[int(cluster[v] == cluster[u]) for u in range(n)] for v in range(n)])
    adj -= A_v
    dist_based = 0
    for row in adj:
        for elem in row:
            dist_based += abs(elem)
    dist_based /= n * n
    res_dist_based = 1 - dist_based

    try:
        node_memb = 0
        for v in range(n):
            node_memb += 1
            node_memb += len(g[v] & clusters[cluster[v]].points) / len(clusters[cluster[v]].points)
            node_memb -= len(g[v] & (all_points - clusters[cluster[v]].points)) / (n - len(clusters[cluster[v]].points))
        node_memb /= 2 * n
        res_node_memb = node_memb
    except:
        res_node_memb = -1

    # print(str(round(len(clusters[0].points) / n * 100, 2)) + '% points were not classified', file=stderr)
    return res_loc_dens, res_node_memb, res_glob_dens, res_dist_based


def readpoints(file=stdin):
    global points_copy, n, dim, points
    n, dim = map(int, file.readline().split())
    points = []
    for i in range(n):
        inp = file.readline().split()
        points.append((np.array(list(map(np.float16, inp[:dim]))), i, inp[dim:]))
    points_copy = points[:]


def set_default_params():
    global K, H
    if len(argv) >= 3:
        K = int(argv[1])
        H = float(argv[2])
    else:
        K = 5
        H = 1


def set_params(k, h):
    global K, H
    K = int(k)
    H = float(h)


def quality_cheat(cluster):
    wrong = 0
    i = 1
    while i < n:
        if points_copy[i][2][0] == 'True' and cluster[i] != cluster[i - 1]:
            wrong += 1
        elif points_copy[i][2][0] == 'False' and (cluster[i] == cluster[i - 1] or cluster[i] == 0):
            wrong += int(points_copy[i][2][1])
            i += int(points_copy[i][2][1]) - 1
        elif points_copy[i][2][0] == 'Rubbish' and cluster[i] != 0:
            wrong += 1
        i += 1
    return (n - wrong) / n


if __name__ == '__main__' and False:
    set_default_params()
    readpoints()
    clusters, cluster, g = Wishart()
    if dim == 2 and True:
        plot(clusters)
    if dim == 3:
        plot3D(clusters, 3)
    if dim >= 4:
        plot3D(clusters, dim)
    print_clusters(clusters)
    print(quality_cheat(cluster), file=stderr)


def make_paths(LSADict, filenames, z=3, _dim=5, k=5, h=0.5):
    global points_copy, n, dim, points
    dim = _dim
    set_params(k, h)
    _points = []
    points_copy = []
    for filename in filenames:
        f = open(filename)
        words = f.readlines()
        for s in words:
            if s:
                try:
                    s1, s2 = s.split()[:2]
                except ValueError:
                    pass
                if s2 not in ['pr', 'conj', 'apro', 'spro']:
                    s1 = s1.replace('?', '')
                    try:
                        _points.append(list(LSADict.mydict[s1][:dim]))
                    except KeyError:
                        pass
    points = []
    for i in range(len(_points) - z + 1):
        curr = []
        for j in range(z):
            curr += _points[i + j]
        points.append((np.array(curr), len(points)))
    n = len(points)
    points_copy = points[:]
    return Wishart()


if __name__ == '__main__' and True:
    LSA_DICT = LSADict(open(input('Enter dict name: ')))
