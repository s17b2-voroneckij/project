import numpy as np
from math import factorial as fact, pi
from sys import setrecursionlimit, argv, stdin
from random import random, seed
import matplotlib.pyplot as plt


seed(179)
setrecursionlimit(1000000)
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
    global g_all, g
    g_all = [[np.linalg.norm(points[i][0] - points[j][0]) for i in range(n)] for j in range(n)]
    g = [[] for i in range(n)]
    points_copy = points[:]
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

def plot(clusters):
    COLORS = ['blue', 'red', 'black', 'yellow', 'purple', 'orange', 'green', 'pink', 'violet', 'brown'] + [[[random(), random(), random()]] for i in range(1000)]
    nc = 0
    for cl in clusters:
        if not cl.is_deleted:
            for i in cl.points:
                curr = points_copy[i]
                plt.scatter(curr[0][0], curr[0][1], color=COLORS[nc])
            nc += 1
    plt.title("K = %s, H = %s" % (str(K), str(H)))
    plt.show()

def quality(clusters, cluster):
    from sys import stderr

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
        print('Global density quality function:', file=stderr)
        print(0.5 * (1 + glob_intern_density_up / glob_intern_density_down - glob_extern_density_up / glob_extern_density_down), '\n', file=stderr)
    except:
        pass
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
        print('Local weighted density quality function:', file=stderr)
        print(loc_dens, '\n', file=stderr)
    except:
        pass

    adj = np.array([[int(u in g[v]) for u in range(n)] for v in range(n)])
    A_v = np.array([[int(cluster[v] == cluster[u]) for u in range(n)] for v in range(n)])
    adj -= A_v
    dist_based = 0
    for row in adj:
        for elem in row:
            dist_based += abs(elem)
    dist_based /= n * n
    print('Distance based quality function:', 1 - dist_based, '\n', sep='\n', file=stderr)

    try:
        node_memb = 0
        for v in range(n):
            node_memb += 1
            node_memb += len(g[v] & clusters[cluster[v]].points) / len(clusters[cluster[v]].points)
            node_memb -= len(g[v] & (all_points - clusters[cluster[v]].points)) / (n - len(clusters[cluster[v]].points))
        node_memb /= 2 * n
        print('Node membership quality function', node_memb, '\n', file=stderr, sep='\n')
    except:
        pass

    print(str(round(len(clusters[0].points) / n * 100, 2)) + '% points were not classified', file=stderr)


def readpoints(file=stdin):
    global prob, points_copy, n, dim, points
    n, dim = map(int, file.readline().split())
    points = []
    for i in range(n):
        inp = file.readline().split()
        points.append((np.array(list(map(np.float32, inp[:dim]))), i, inp[dim:]))
    prob = [0] * n
    points_copy = points[:]


try:
    def set_params(k=argv[1], h=argv[2]):
        global K, H
        K = int(k)
        H = float(h)

except:
    def set_params(k, h):
        global K, H
        K = int(k)
        H = float(h)

if __name__ == '__main__':
    set_params()

    readpoints()

    clusters, cluster, g = Wishart()
    if dim == 2:
        plot(clusters)
    print_clusters(clusters)
    quality(clusters, cluster)