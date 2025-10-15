#kmean
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np
from pyclustering.cluster.kmeans import kmeans
from pyclustering.utils.metric import type_metric, distance_metric
from pyclustering.cluster.center_initializer import kmeans_plusplus_initializer
from scipy import stats

def pearson_dist(x, y):
    r = stats.pearsonr(x, y)[0]
    return (1 - r) / 2

def geneNameTransformer(query):
    df = pd.read_csv(r"C:\Users\omer9\Documents\לימודים\עבודה במעבדה\Itamar's lab\RNA seq project\Tools\Gene names transformer\genes_names.csv") 
    if type(query.iloc[0,0]) == str:
        mode = "Human->Fish"
        lookfor = 0
        lookin = "Human"
        match = "NCBI"
        fill = 1
    elif type(query.iloc[0,1]) == str:
        mode = "Fish->Human"
        lookfor = 1
        lookin = "NCBI"
        match = "Human"
        fill = 0
    #Overcome key-sensative 
    lookinForTest =[]
    for i in list(df[lookin]):
        lookinForTest.append(str(i).upper())
  
    #Find the index of the genes
    indexes = []
    for i in query.index:
        if query.iloc[i,lookfor].upper() in lookinForTest:
            indexes.append(list(lookinForTest).index(query.iloc[i,lookfor].upper()))
        else:
            indexes.append(None)       
    #Match the names
    for i in range(len(indexes)):
        if indexes[i] != None:
            query.iloc[i,fill] = df.loc[indexes[i],match]
    return query   
  
  
  
df = pd.read_csv(r"singe&double_sigmoid_fitted_data.csv")  
df.set_index("gene", inplace = True)



#kmean
k = 8
# change your df to numpy arr
sample = df.to_numpy()
# define a custom metric
metric = distance_metric(type_metric.USER_DEFINED, func=pearson_dist)
# carry out a km++ init
initial_centers = kmeans_plusplus_initializer(sample, k).initialize()
# execute kmeans
kmeans_instance = kmeans(sample, initial_centers, metric=metric)
# run cluster analysis
kmeans_instance.process()
# get clusts
clusters = kmeans_instance.get_clusters()


all_clusters = []

#foldchange from start
for i in range(len(df)):
    df.iloc[i,:] = df.iloc[i,:]-df.iloc[i,0]
    
#creat df for each cluster
for cluster in range(k):
    clusterdf = pd.DataFrame()
    

#check if gene belong to cluster
    for i in clusters[cluster]:  
        clusterdf = pd.concat([clusterdf,df.iloc[i]],axis=1)



#get cluster genes
    report = pd.DataFrame(clusterdf.columns)
    file_name = "cluster " + str(cluster+1) + "genes.csv"
    report.to_csv(file_name)


#preper for GO
    #geneNameTransformer parser
    report.columns = ["Fish"]
    report["Hum"] = None
    report = report.iloc[:,[1,0]]
    report = geneNameTransformer(report)
    file_name = "Cluster" + str(cluster+1) + " genes (for GO).xlsx"
    report["Hum"].to_excel(file_name,index = False)


#calc the mean of the cluster
    clusterdf["mean"] = clusterdf.mean(axis=1)

#calc the std of the cluster
    clusterdf["SD"] = clusterdf.iloc[:,:-2].std(axis=1)
    
    all_clusters.append(clusterdf.copy()) 
    
#make cluster graph
    x = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
    #x = [1,2,3,4,5,6,7]
    y1 = clusterdf["mean"] + clusterdf["SD"]
    y2 = clusterdf["mean"] - clusterdf["SD"]
    fig, ax = plt.subplots()
    ax.fill_between(x, y1, y2, alpha=.5, linewidth=0)
    ax.plot(x, (y1 + y2)/2, linewidth=2)
    ax.set(xticks=x)
    plt.grid(visible = True, axis ="x")
    title = "Cluster " + str(cluster+1) + "; n= " + str(clusterdf.shape[1])
    ax.set_title(title)
    plt.savefig(title+".jpg", format = "jpg")



# # plot actual expretion
#     fig, ax = plt.subplots()        
#     ax.plot(clusterdf,color = "blue", alpha =0.1)
#     title = "Cluster " + str(cluster+1) + "; n= " + str(clusterdf.shape[1])


#plot all clusters on one fig
    x = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
    y = clusterdf["mean"]
    plt.plot(x,y,label = "Cluster "+str(cluster+1))
    plt.grid(visible = True, axis ="x")
    plt.xticks(x)
    #plt.ylim(-1,1)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
