# -*- coding: utf-8 -*-
"""
Created on Wed Mar 09 11:22:51 2016

@author: juherask

Based on:
[Hu14] Hutter, F., Xu, L., Hoos, H.H., Leyton-Brown, K.: Algorithm runtime
prediction: Methods & evaluation. Artificial Intelligence206(2014) 79–111
"""

from prim_mst import minimum_spanning_tree
from dijkstra import dijkstra
from helpers.rkm_stats import calculate_distribution_statistics
from collections import Counter, defaultdict

def mst_features(npts, nD, root_note_idx=None, return_MST_graph=False):
    features = []
    
    # Use the Prin's MST algorithm as implemented by Andreas Mueller
    mst_edges = minimum_spanning_tree(nD)    
    
    mds_ds = [nD[i][j] for i,j in mst_edges]
    node_degrees = Counter((ij for edge in mst_edges for ij in edge)).values()
    
    features.append( ("MST1", "SotD of Minimum Spanning Tree Edge Costs [Hu14,Me13]",
                      calculate_distribution_statistics(mds_ds, level=1)) )
    features.append( ("MST2", "SotD of Minimum Spanning Tree Node Degrees [Hu14]",
                      calculate_distribution_statistics(node_degrees)) )

    ## Depth from the root node
    if root_note_idx==None:
        root_note_idx = mst_edges[0][0]
        
    # convert to graph
    G = defaultdict(dict)
    for i,j in mst_edges:
        G[i][j]=1
        G[j][i]=1
    # calculate node depth from root to every other node
    depths = dijkstra(G,root_note_idx)
    features.append( ("MST3", "SotD of Minimum Spanning Tree Node Depth from the Root [Me13]", 
                      calculate_distribution_statistics(depths.values(), level=1)) )
    
    features.append( ("MST4", "Normalized Minimum Spanning Tree Cost Sum [Me13]",
                      sum(mds_ds)/sum(nD.flatten()) ) )
    
    if return_MST_graph:
        return features, G
    else:
        return features
    
def smoke_test():
    from cvrp_rkm16_io import generate_CVRP
    (N, pts,demands, D, C) = generate_CVRP(50, 100.0, 5.0, 2.0)
    
    #print pts[0]
    print mst_features(pts, D, 0) # 0 is the depot in VRP
    
if __name__=="__main__":
    smoke_test()