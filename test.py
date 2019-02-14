from main import *




f = open('data_1_dim_2.txt')
readpoints(f)
set_params(3, 1.5)
clusters, cluster, g = Wishart()
quality(clusters, cluster)