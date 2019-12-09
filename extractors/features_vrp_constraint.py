# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 10:28:51 2016

@author: juherask

Based on:
[Sh15] Steinhaus, M. "The Application of the Self Organizing Map to the Vehicle
Routing Problem.". Open Access Dissertations. Paper 383. (2015). University of
Rhode Island.
"""


from helpers.rkm_stats import calculate_distribution_statistics, STAT_LEVEL_EXTRA
from tsp_solver import solve_tsp_acotsp

import numpy as np

from math import ceil
from collections import Counter

def constraint_features(demands_with_depot, C, clusters, D, L=None, k=None):
    """
    pts is list of tuples (x,y,d)
    C is the truck capacity
    k is the min. nbr. of trucks (if known)
    """
    features = []
    # ignore the depot (it has demand of 0)
    
    if demands_with_depot:
        demands = demands_with_depot[1:]
        #normalize with vehicle capacity
        nds = np.array([d/float(C) for d in demands])
        sum_nds = sum(nds)
    else:
        C = 1.0
        nds = [0.0]*(len(D)-1)
        sum_nds = 0.0
    
    # if k is not set, estimate the number of trucks needed as in [Sh15]
    if not k:
        k = ceil(sum_nds)
    
    normalized_demand_stats = calculate_distribution_statistics(nds)
    features.append( ("DC5", "SotD of Normalized Node Demand [Sh15]", normalized_demand_stats) )
    features.append( ("DC6", "Ratio of Total Demand to Total Capacity [Sh15]", sum_nds/float(k)) )
    features.append( ("DC9", "Ratio of Single Biggest Demand to Capacity [Sh15]", max(nds) ) )
    how_many_average_sized_request_per_car = 1.0/(sum_nds/float(len(nds))) if sum_nds else float("inf")
    features.append( ("DC10", "Number of Average Sized Customers fits on a Vehicle [Sh15]", how_many_average_sized_request_per_car ) )
    features.append( ("DC11", "Minimum Number of Vehicles due to Capacity [Sh15]", ceil(sum_nds) ) )
    
    # clusters
    if clusters is not type(np.array):
        clusters = np.array(clusters)
    # ignore the depot
    node_clusters = clusters[1:] 
    
    cluster_ids = set(node_clusters)
    cluster_demands = []
    outlier_demand = 0.0
    for cluster_id in cluster_ids:
        cluster_demand = sum(nds[ node_clusters==cluster_id ]) if sum_nds else 0.0
        
        
        # outliers        
        if cluster_id==-1:
            outlier_demand = cluster_demand
        else:
            cluster_demands.append(cluster_demand)
    ratio_cluster_demands_to_C = (0.0, 0.0, 0.0, 0.0, float('inf'), 0.0, 0.0, 0.0)
    if cluster_demands:
        ratio_cluster_demands_to_C = calculate_distribution_statistics(cluster_demands, level=STAT_LEVEL_EXTRA)
            
    features.append( ("DC7", "SotD of Ratio of Cluster Demands to Vehicle Capacity [Sh15]",
                      ratio_cluster_demands_to_C ) )
    features.append( ("DC8", "Ratio of Outlier Demand to Total Demand [Sh15]",
                      outlier_demand/float(sum_nds) if sum_nds else 0.0 ) )
    
    # avoid NaN
    ratio_TSP_to_L = 1.0
    diff_CKlb_to_LKlb = 0.0
    ratio_CKlb_to_LKlb = 1.0
    ratio_clusterTSP_to_L = (1.0, 0.0, 0.0, 0.0, float('inf'))
    
    if L:
        # route all customers using acotsp NN init + local search
        _, tsp_f = solve_tsp_acotsp(D, range(len(D)))
        
        ratio_TSP_to_L = tsp_f/float(L) # ~ lower bound on how many vehicles
        diff_CKlb_to_LKlb = ceil(sum_nds)-ceil(ratio_TSP_to_L)
        ratio_CKlb_to_LKlb  = sum_nds/ratio_TSP_to_L
        
        
        cluster_tsp_fs = []
        for cluster_id in cluster_ids:
            cluster_node_ids = [0]+np.flatnonzero( node_clusters==cluster_id).tolist()
            _, cluster_tsp_f = solve_tsp_acotsp(D, cluster_node_ids)
            cluster_tsp_fs.append(cluster_tsp_f/float(L))
        ratio_clusterTSP_to_L = calculate_distribution_statistics(cluster_tsp_fs)
        
    
    features.append( ("DC12", "Lower Bound for Vehicles due to Route Duration [RMK19]", ratio_TSP_to_L) )
    features.append( ("DC13", "Difference between Capacity and Route Duration Constraint Estimated Vehicle Counts [RMK19]",diff_CKlb_to_LKlb) )
    features.append( ("DC14", "Ratio in Vehicle Count Lower Bound between Capacity and Route Duration Constraints [RMK19]",ratio_CKlb_to_LKlb) )
    features.append( ("DC15", "SotD of Cluster Routing to Route Duration Constraints [RMK19]",                
                      ratio_clusterTSP_to_L) )

    return features
    
def test():
    from cvrp_rkm16_io import read_TSPLIB_CVRP
    from pprint import pprint
    from math import sqrt
    from sklearn.cluster import DBSCAN
    
    (N, pts, dd, demands, D, C) = read_TSPLIB_CVRP("../instances/A-n32-k5.vrp")
    k = 5
    
    # Cluster
    n = len(pts)
    xs, ys = zip(*pts)
    nn_4th_estimated_eps = sqrt((max(xs)-min(xs))*(max(ys)-min(ys)))/(sqrt(n)-1)
    node_clusters = DBSCAN(eps=nn_4th_estimated_eps, metric='precomputed', min_samples=4).fit( D ).labels_

    pprint( demand_features(demands, C, node_clusters, k) )# 0 is the depot in VRP
    
if __name__=="__main__":
    test()