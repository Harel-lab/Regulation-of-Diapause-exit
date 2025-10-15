#### Packages ####
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np
from pyclustering.cluster.kmeans import kmeans
from pyclustering.utils.metric import type_metric, distance_metric
from pyclustering.cluster.center_initializer import kmeans_plusplus_initializer
from scipy import stats

### Upload and modify data ####
df_path = r"singe&double_sigmoid_fitted_data.csv"
df = pd.read_csv(df_path)  
df.set_index("gene", inplace = True)
sample = df.to_numpy()

#### Kmean clustering ####
k = 8
# Pearson corelation
def pearson_dist(x, y):
    r = stats.pearsonr(x, y)[0]
    return (1 - r) / 2
# define a empty metric
metric = distance_metric(type_metric.USER_DEFINED, func=pearson_dist)
# carry out a km++ init
initial_centers = kmeans_plusplus_initializer(sample, k).initialize()
# kmeans
kmeans_instance = kmeans(sample, initial_centers, metric=metric)
kmeans_instance.process()
clusters = kmeans_instance.get_clusters()

#### Gene list per cluster ####
# Create df for each cluster
for cluster in range(k):
    clusterdf = pd.DataFrame()
# Check if gene belong to cluster
    for i in clusters[cluster]:  
        clusterdf = pd.concat([clusterdf,df.iloc[i]],axis=1)
# Output
    report = pd.DataFrame(clusterdf.columns)
    file_name = "cluster " + str(cluster+1) + "genes.csv"
    report.to_csv(file_name)

