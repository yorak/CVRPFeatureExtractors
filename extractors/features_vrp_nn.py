# -*- coding: utf-8 -*-
"""
Created on Wed Mar 09 11:22:51 2016

@author: juherask

Based on:
[Hu14] Hutter, F., Xu, L., Hoos, H.H., Leyton-Brown, K.: Algorithm runtime
prediction: Methods & evaluation. Artificial Intelligence206(2014) 79–111
"""

from helpers.rkm_stats import calculate_distribution_statistics
from helpers.rkm_util import produce_nn_list

from collections import Counter
from math import atan2, pi
from scipy.spatial import distance
import sys

try:    
    import networkx as nx
except:
    print >> sys.stderr, "WARNING: no networkx, only some graph features are calculated", 
    nx = None

X = 0; Y = 1

def _angle_in_rads(p1, p2):
    dx = p2[X]-p1[X]
    dy = p2[Y]-p1[Y]
    rads = atan2(dy,dx)
    rads %= (2*pi)
    return rads
    
def nn_features(npts, nD, n_NN_D=None):
    """ n_NN_D is the normalized (if you want), nearest neighbour distances,
    where each row contains a list of NN's for that point. The NN's are 
    listed as pairs with (node_id, distance). If this is not given, it
    is calcluated using util.produce_nn_list(nD)
    """
    n = len(nD)
    features = []
    
    if n_NN_D is None:
        n_NN_D = produce_nn_list(nD)
        
     # 60. nNNd Distribition of distances to the nearest neighbor distances
    nns = []
    angles_2nn = []
    cosine_similarity_2nn = []
    for i in range(n):
        # take distance to nearest neighbour (1=first, 1=value)
        nns.append( n_NN_D[i][1][1] )
    
        from_pt = npts[i]
        to_pt_nn_idx = 1
        to_pt = None
        # Do not accept overlapping nodes for calculating angle (not defined)
        while (to_pt_nn_idx<n) and ((to_pt is None) or (tuple(to_pt)==tuple(from_pt))):
            to_idx = n_NN_D[i][to_pt_nn_idx][0]
            to_pt = npts[to_idx]
            to_pt_nn_idx+=1
        v =  (to_pt[X]-from_pt[X], to_pt[Y]-from_pt[Y])
        
        to_pt = None
        # Do not accept overlapping nodes for calculating angle (not defined)
        while (to_pt_nn_idx<n) and ((to_pt is None) or (tuple(to_pt)==tuple(from_pt))):
            to_idx = n_NN_D[i][to_pt_nn_idx][0]
            to_pt = npts[to_idx]
            to_pt_nn_idx+=1
        w =  (to_pt[X]-from_pt[X], to_pt[Y]-from_pt[Y])
        
        # n_NN_D[i][j][0] = index of the node that is the j_th nearest neighbour of i
        angles_2nn.append(  _angle_in_rads(v, w) )
        cosine_similarity_2nn.append(distance.cosine(v,w))
        
    features.append( ("NN1","SotD of 1st Nearest Neighbour Distance [SM10]", 
                      calculate_distribution_statistics(nns)) )

    #Angle of edges connecting 2 nearest neighbours   
    features.append( ("NN2","SotD of Angle (in radians) Between Edges to Two Nearest Neighbours [Me13]",
                      calculate_distribution_statistics(angles_2nn, level=1) ) )
    #Cosinus of the angles connecting 2-nn    
    features.append( ("NN3","SotD of Cosine Similarity Between Edges to Two Nearest Neighbours [PM14]",
                      calculate_distribution_statistics(cosine_similarity_2nn, level=1) ) )
        
    # PM14 kNN features for different k
    ks = sorted(list(set([3, 5, 7])))
                #int(round(n**(1.0/3))), int(round(2*n**(1.0/3))),
                #int(round(0.5*n**(0.5))), int(round(n**(0.5)))])))
    #print ks
    nn_feature_idx = 4
    for k in ks:
        if nx:
            diG = nx.DiGraph()
            diG.add_nodes_from(range(n))
        
        #Input degree in Directed kNN graph (DkNNG) [PM14]
        incoming = []
        for i in range(n):
            knn, _ = zip(*n_NN_D[i][1:k+1])
            incoming+=knn
            
            # build networkx digraph            
            if diG is not None:
                diG.add_edges_from( zip([i]*k, knn) )
            
        
        in_knn_node_deg = Counter(incoming).values()
		# level 2 stats include min, max, median, and quantiles
        features.append( ("NN%d"%nn_feature_idx,
                          "SotD of Input Degrees in Directed %dNN Graph [PM14]"%k,
                          calculate_distribution_statistics(in_knn_node_deg, level=3) ) )
        nn_feature_idx+=1
        
        if nx:
            #udG = diG.to_undirected()
            
            sccs = nx.strongly_connected_components(diG)
            wccs = nx.weakly_connected_components(diG)
            
            scc_sizes = [len(scc) for scc in sccs ] 
            nscc_sizes = [s/float(n) for s in scc_sizes]
            wcc_sizes = [len(wcc) for wcc in wccs ]
            nwcc_sizes = [s/float(n) for s in wcc_sizes]
            
            #print k, scc_sizes, wcc_sizes
            
            #Strongly connected component count in DkNNG            
            features.append( ("NN%d"%nn_feature_idx,
                              "Number of Strongly Connected %dNN Digraph Components [PM14]"%k,
                              len(scc_sizes)) )
            nn_feature_idx+=1
            features.append( ("NN%d"%nn_feature_idx,
                              "Number of Strongly Connected %dNN Digraph Components (normalized) [PM14]"%k,
                              len(scc_sizes)/float(n)) )
            nn_feature_idx+=1
            #Strongly connected component size in DkNNG
            features.append( ("NN%d"%nn_feature_idx,
                              "SotD of Strongly Connected %dNN Digraph Component Size [PM14]"%k,
                              calculate_distribution_statistics(scc_sizes, level=1)) )
            nn_feature_idx+=1
            features.append( ("NN%d"%nn_feature_idx,
                              "SotD of Strongly Connected %dNN Digraph Component Size (normalized) [PM14]"%k,
                              calculate_distribution_statistics(nscc_sizes, level=1)) )
            nn_feature_idx+=1
                              
            #Weakly connected component count in kNNG
            features.append( ("NN%d"%nn_feature_idx,
                             "Number of Weakly Connected %dNN Graph Components [PM14]"%k,
                             len(wcc_sizes)) )
            nn_feature_idx+=1
            features.append( ("NN%d"%nn_feature_idx,
                              "Number of Weakly Connected %dNN Graph Components (normalized) [PM14]"%k,
                              len(wcc_sizes)/float(n)) )
            nn_feature_idx+=1
            #Weakly connected component size in kNNG
            features.append( ("NN%d"%nn_feature_idx,
                              "SotD of Weakly Connected %dNN Graph Component Size [PM14]"%k,
                              calculate_distribution_statistics(wcc_sizes, level=1)) )
            nn_feature_idx+=1
            features.append( ("NN%d"%nn_feature_idx,
                              "SotD of Weakly Connected %dNN Graph Component Size (normalized) [PM14]"%k,
                              calculate_distribution_statistics(nwcc_sizes, level=1)) )
            nn_feature_idx+=1
            
            #Ratio of number of strongly and weakly connected cmp.
            features.append( ("NN%d"%nn_feature_idx,
                              "Ratio of Strongly and Weakly Connected %dNN (Di)Graph Components [PM14]"%k,
                              len(scc_sizes)/float(len(wcc_sizes))) )
            nn_feature_idx+=1
            
    
    return features
    
def smoke_test():
    from cvrp_rkm16_io import generate_CVRP
    (N, pts,demands, D, C) = generate_CVRP(50, 100.0, 5.0, 2.0)
    #print pts[0]
    print nn_features(pts, D)
    
if __name__=="__main__":
    smoke_test()