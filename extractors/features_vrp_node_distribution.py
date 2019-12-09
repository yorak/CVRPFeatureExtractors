# -*- coding: utf-8 -*-
"""
Created on Mon Mar 07 14:47:59 2016

@author: juherask

[Hu14] Hutter, F., Xu, L., Hoos, H.H., Leyton-Brown, K.: Algorithm runtime
prediction: Methods & evaluation. Artificial Intelligence206(2014) 79–111

[SM10] Smith-Miles, K., van Hemert, J., & Lim, X. Y. (2010). Understanding
TSP difficulty by learning from evolved instances. In Proceedings of the
4th Workshop on Learning and Intelligent Optimization (LION’10),
volume 6073 of LNCS, (pp. 266–280). Springer-Verlag.
"""

from helpers.rkm_stats import calculate_distribution_statistics
from helpers.rkm_util import d

import numpy as np
from sklearn.cluster import DBSCAN # Used to check for clusters
from sklearn import metrics

from collections import Counter
from math import sqrt
    
X = 0; Y = 1 # ~DEFINE for readablity
PT_NORMALIZATION_TO_RANGE = (0.0,400.0)

def cluster_nodes(pts, D, rangex = None, rangey = None):
    n = len(pts)
    if not rangex or not rangey:
        xl, yl = zip(*pts)
        rangex = max(xl)-min(xl)
        rangey = max(yl)-min(xl)
        
    # Use the eps estimation method suggested by Steinhous 2015
    #  assume cities placed uniformly on a b x b -grid.
    nn_4th_estimated_eps = sqrt(rangex*rangey)/(sqrt(n)-1)
    
    # TODO: Alternative way to determine eps for DBSCAN (from Steinhous 2015)
    # Run DBSCAN multiple times using a range of values between [median,85 th % ] for Eps,
    #  and choose an Eps value that corresponds to the most frequent number of clusters found
    #  across all runs of DBSCAN.
    
#    print D
#    print D.shape
    
    db = DBSCAN(eps=nn_4th_estimated_eps, metric='precomputed', min_samples=4).fit( D )
    return db.labels_, db.core_sample_indices_
   
def node_distribution_features(npts, nD, clusters=None, cluster_core_point_idxs=None):
    features = []
    
    ## 2-4. & 51. Shape of the Cost Matrix Distribution ## 
    
    flattened_D = [item for i, row in enumerate(nD[:-1]) for item in row[i+1:]]
    # in addition to the st.dv. in [SM10], include also other statistical
    #  moments: mean, sd, skew, kurthosis
    ## Desbribe Cost Matrix Distribution ##
    features.append( ("ND1", "SotD of Cost Matrix [SM10, Hu14]",
                      calculate_distribution_statistics(flattened_D, level=3)) )
    
    ## 52-55. Fraction of Distinct Distances
    dd_frac = []
    np_fD = np.array(flattened_D)
    tot_nbr_D = float(len(np_fD))
    for precision in range(1,5):
        dd_frac.append( len(set(np_fD.round(precision)))/tot_nbr_D ) 
    features.append( ("ND2", "Fraction of Rounded Distinct Distances [SM10]", tuple(dd_frac) ) )

    ## 56-58. Centroid
    xl, yl = zip(*npts)
    centroid = (sum(xl)/len(npts), sum(yl)/len(npts))
    features.append( ("ND3", "Centroid of (0,400) Normalized Vertex Positions (x,y) [SM10, Hu14]", centroid) )
    d_to_centroid = [d(centroid,npt) for npt in npts]
    features.append( ("ND4", "SotD of Distance to Centroid [SM10]",
                      calculate_distribution_statistics(d_to_centroid, level=1)) )
    
    # distance to depot(s) (VRP specific feature)
    features.append( ("DC2", "Normalized Depot Location (x,y) [Sh15]", npts[0]) )
                      
    # distance to depot(s) (VRP specific feature)
    features.append( ("DC3", "Distance between Depot and Centroid [RKM16]",
                      d(centroid, npts[0])) )
                      
    d_to_depot = [d(npts[0],npt) for npt in npts[1:]]
    features.append( ("DC4", "SotD of Distance to Depot [RKM16]",
                      calculate_distribution_statistics(d_to_depot)) )
    
    # 59. Area of the Enclosing Rectangle [SM10]
    max_x_bb = max(xl)
    max_y_bb = max(yl)
    rel_area = (max_x_bb*max_y_bb) / ((max(PT_NORMALIZATION_TO_RANGE)-min(PT_NORMALIZATION_TO_RANGE))**2)
    features.append( ("G1","Area of the Enclosing Rectangle [SM10]", rel_area ) )
    
    ## 56-58. Clusters
    # TODO: GDBSCAN? Check difference between it and DBSCAN
    
    if clusters is None or cluster_core_point_idxs is None:
        labels, cluster_core_point_idxs = cluster_nodes(npts, nD, rangex = max_x_bb,rangey = max_y_bb)
    else:
        labels = clusters
        
    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_outliers_ = sum(labels == -1)
    n_core_ = len(cluster_core_point_idxs)
    
    features.append( ("ND5-1","Number of clusters (absolute, DBSCAN) [SM10]", n_clusters_) )
    features.append( ("ND5-2","Number of clusters (relative to N, DBSCAN) [SM10]", float(n_clusters_)/len(npts)) )
    features.append( ("ND6-1","Number of Nodes are Outliers (relative, DBSCAN) [SM10]",
                      float(n_outliers_) / len(npts) ) )
    features.append( ("ND6-2","Number of Nodes are Edge Points (relative, DBSCAN) [SM10]",
                      float(len(npts)-n_outliers_-n_core_) / len(npts) ) )
    features.append( ("ND6-3","Number of Nodes are Core Points (relative, DBSCAN) [SM10]",
                      float(n_core_) / len(npts) ) )
    
    if not 1 < n_clusters_ < len(npts):
        #only outliers
        ss = 0.0 #only points "on the cluster border"
    else:
        # silhouette_score does not like outlier labels (-1), so
        # r eplace all outliers with valid cluster IDs (clusters of 1 element)
        ss_labels = np.array(labels)
        ss_outlier_idx = max(labels)+1
        for sslidx in range(len(ss_labels)):
            if ss_labels[sslidx]==-1:
                ss_labels[sslidx] = ss_outlier_idx
                ss_outlier_idx+=1
        # the silhouette_score will throw a "Mean of an empty slice"
        #  RuntimeWarning for 1 member clusters,
        #  
        old_settings = np.seterr(all='ignore')  #seterr to known value
        ss = metrics.silhouette_score(nD, ss_labels, metric='precomputed')
        np.seterr(**old_settings)  # reset to default
        #print "stilhouette score : ", ss
        
    features.append( ("ND9", "Silhouette Coefficient (DBSCAN) [RKM16]", ss ) )
    
    nnodes_in_clusters = Counter(labels)
    # ignore outliers
    del nnodes_in_clusters[-1]
    nnode_count_list = nnodes_in_clusters.values()
    if len(nnode_count_list)==0:
        nnsm = (1.0, 0.0, 0.0, 0.0, 0.0) #only outliers
    else:
        nnsm = calculate_distribution_statistics(nnode_count_list)
    features.append( ("ND8", "SotD of Nodes in a Cluster [SM10]", nnsm) )
    
    npts = np.array(npts)
    ds_to_cluster_centroid = []
    for cluster_id, cluster_size in nnodes_in_clusters.items():
        cluster_pts = npts[labels==cluster_id]
        cluster_sum_x, cluster_sum_y = np.sum(cluster_pts , axis=0)
        cluster_centroid = (cluster_sum_x/cluster_size, cluster_sum_y/cluster_size)
        for cpt in cluster_pts:
            ds_to_cluster_centroid.append(d(cluster_centroid, cpt))
    if len(ds_to_cluster_centroid)==0:
        cdsm = (0.0, 0.0, 0.0, 0.0, 0.0) #only outliers
    else:
        cdsm = calculate_distribution_statistics(nnode_count_list)
    features.append( ("ND7","SotD of Distance to Cluster Centroid [SM10]", cdsm) )
    
        
    return features

def smoke_test():
    #dummy_pts = [(14.3, 4.4, 3), (-123.3, 12.4, 2), (0.0, -22.2, 7)]    
    #pprint( _normalize_to_rect(dummy_pts, (0,400), True ) )
    #from helpers.rkm_util import calculate_D
    #pprint( calculate_D(dummy_pts) )
    
    from helpers.cvrp_rkm16_io import generate_CVRP
    (N, pts,demands, D, C) = generate_CVRP(50, 100.0, 5.0, 2.0)
    from pprint import pprint
    pprint(node_distribution_features(pts, D))
    
if __name__=="__main__":
    smoke_test()