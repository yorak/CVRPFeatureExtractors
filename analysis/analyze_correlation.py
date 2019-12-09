# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 21:29:48 2016

@author: juherask
"""

from analyze_process import *
from analyze_clusters import plot_clusters

import numpy as np
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import MinMaxScaler

from math import isnan
from scipy import stats
import sys


def _robust_correlation(a,b):
    s_a = np.std(a)
    s_b = np.std(b)
    np_a = np.array(a)
    np_b = np.array(b)
        
    # other series is constant, correlation is undefined, but to
    #  us it servers well that we say there is NO CORRELATION
    if s_a==0.0 or s_b==0.0:
        return (0.0, 1.0)
    # if there is nan values in the list, only use rows that have valid
    #  values
    if isnan(s_a) or isnan(s_b):
        non_nan_a_selector = np.array([not isnan(v_a) for v_a in np_a], dtype=bool)
        non_nan_b_selector = np.array([not isnan(v_b) for v_b in np_b], dtype=bool)
        #print np_a
        #print np_b
        #print non_nan_a_selector & non_nan_b_selector 
        
        np_a = np_a[non_nan_a_selector & non_nan_b_selector]
        np_b = np_b[non_nan_a_selector & non_nan_b_selector]
        
        # If there is no values to use -> no correlation
        if len(np_a)<=1:
            return (0.0, 1.0)
    return stats.pearsonr( np_a, np_b )

def correlation_distance_matrix(float_cols):
    N = len(float_cols)
    D = np.zeros( (N, N) )
    for i in range(N):
        for j in range(i, N):
            d = 1.0-abs(_robust_correlation(float_cols[i],float_cols[j])[0])
            if isnan(d):
                d = 1.0
            D[i][j]=d
            D[j][i]=d
            #print d
    return D
        
    
    
def print_promising_correlations(float_headers, float_cols, how_many_best=10):
    print "Correlation targets!"
    print "*",float_headers[12] # neg and pos separately?
    print "*",float_headers[85]
    print "*",float_headers[90]
    print
    
    
    best_corr = []
    all_corr = []
    
    TOP_K = 10
    selpos_bc_values = float_cols[12]>=0
    selneg_bc_values = float_cols[12]<0
    correlations = []
    print "top %d correlates with B&C pos quality (feasible solutions)" % TOP_K
    for fid, fcol in enumerate(float_cols):
        if fid == 12:
            continue
        correlations.append( ( 12+1, fid+1, _robust_correlation( float_cols[12][selpos_bc_values], fcol[selpos_bc_values] )))
    correlations.sort(key=lambda corr: abs(corr[2][0]), reverse=True) # (id, (corr, p-value))
    all_corr+=correlations
    for i in range(TOP_K):
        print correlations[i]
        best_corr.append(correlations[i])
    print 
    
    correlations = []    
    print "top %d correlates with B&C neg quality (infeasible solutions -> lower bound)" % TOP_K
    for fid, fcol in enumerate(float_cols):
        if fid == 12 or fid == 11:
            continue
        correlations.append( ( 12+1, fid+1, _robust_correlation( float_cols[12][selneg_bc_values], fcol[selneg_bc_values] )))
    correlations.sort(key=lambda corr: abs(corr[2][0]), reverse=True) # (id, (corr, p-value))
    all_corr+=correlations
    for i in range(TOP_K):
        print correlations[i]
        best_corr.append(correlations[i])
    print 
    
    correlations = []    
    print "top %d correlates with C&W quality" % TOP_K
    for fid, fcol in enumerate(float_cols):
        if fid == 85:
            continue
        correlations.append( ( 85+1, fid+1, _robust_correlation( float_cols[85], fcol )))
    correlations.sort(key=lambda corr: abs(corr[2][0]), reverse=True) # (id, (corr, p-value))
    all_corr+=correlations
    for i in range(TOP_K):
        print correlations[i]
        best_corr.append(correlations[i])
    print 
    
    correlations = []
    print "top %d correlates with LS quality" % TOP_K
    for fid, fcol in enumerate(float_cols):
        if fid == 90:
            continue
        correlations.append( ( 90+1, fid+1, _robust_correlation( float_cols[90], fcol )))
    correlations.sort(key=lambda corr: abs(corr[2][0]), reverse=True) # (id, (corr, p-value))
    all_corr+=correlations
    for i in range(TOP_K):
        print correlations[i]
        best_corr.append(correlations[i])
    print 

    candidates = all_corr
    
    # 
    from pprint import pprint
    print
    for c in candidates:
        if isnan(c[2][0]):
            print c
            print float_cols[c[1]-1]
    print
    
    from collections import defaultdict
    c_score = defaultdict(float)
    for c in candidates:
        c_score[c[1]]+=abs(c[2][0])
    score_list = c_score.items()
    score_list.sort(key=lambda sc: sc[1], reverse=True)
    
    print 
    print "The best"
    pprint(score_list[:how_many_best])
    print
    
    return score_list[:how_many_best]

                    
#features_file = "vrp_features_Augerat-A-14.csv"
#features_file = "vrp_features_Augerat-B-14.csv"
#features_file = "vrp_features_Augerat-AandB-14.csv"
features_file = 'vrp_all_features_FINAL.csv'

header, cols, float_headers, float_cols = read_features_from_file_and_validate(features_file)

print_best_features = False
if print_best_features:
    best_features = print_promising_correlations(float_headers, float_cols, how_many_best=50)
    
    sys.stdout.write("instance")
    for col_idx, (bc_id, bc_corr) in enumerate(best_features):
        sys.stdout.write(",f%02d"%bc_id)
    print
    
    for i in range(14):
        sys.stdout.write(cols[1][i]+".vrp") #instance name
        for col_idx, (bc_id, bc_corr) in enumerate(best_features):
            sys.stdout.write(","+str(float_cols[bc_id-1][i]))
        print

plot_correlated_features = True
if plot_correlated_features:
    var_col = np.array([np.std(col)!=0.0 for col in float_cols],dtype=bool)
    float_cols = np.array(float_cols)
    
    import os.path
    if os.path.isfile("MDS_feature_pts.txt"):
        pts = np.loadtxt("MDS_feature_pts.txt")
        D = np.loadtxt("D.txt")
    else:
        D = correlation_distance_matrix(float_cols[var_col])
        np.savetxt("D.txt",D)
        
        from sklearn import manifold
        mds = manifold.MDS(n_components=2, dissimilarity='precomputed')
        results = mds.fit(D)
        pts = results.embedding_
        
        np.savetxt("MDS_feature_pts.txt", pts)
    
    #print pts
    
    
    
    kth_nns = []
    for i in range(len(D)):
        fkth = sorted(D[i,:])[3]
        kth_nns.append(fkth)
    kth_eps = np.mean(kth_nns)
    
    #for seekeps in np.linspace(0.07, 0.11, 40):
    from sklearn.cluster import DBSCAN
    
    # Got 0.085 experimentally
    db = DBSCAN(eps=0.085, min_samples=4, metric='precomputed').fit( D  )
    labels = db.labels_
    core_samples_mask = np.zeros_like(labels, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True        

    from collections import Counter
    counts = Counter(labels)
    #print
    #print seekeps, len(counts)
    #exit()
    
    onlyclustered = np.zeros( len(labels), dtype=bool )
    for cluster_id in sorted(counts.keys()):
        if cluster_id == -1:
            continue

        print
        print "cluster", cluster_id
        for i in range(len(labels)):
            if labels[i]==cluster_id:
                print float_headers[i]
                onlyclustered[i]=True
     
    for cluster_id in sorted(counts.keys()):      
        print "cluster", cluster_id, counts[cluster_id] 
    
    

    import matplotlib.pyplot as plt
    
    #ptsD = squareform( pdist(pts) )    
    #expDeltaD = 1/(0.1+np.abs(D-ptsD))/10
    #z = np.sum(expDeltaD, axis=0)
    
    mdsClustered = True
    if mdsClustered:
        
        ocD = (D[onlyclustered,:])[:,onlyclustered]
        
        print ocD.shape
        
        from sklearn import manifold
        mds = manifold.MDS(n_components=2, dissimilarity='precomputed')
        results = mds.fit(ocD)
        ocPts = results.embedding_
        
        sqsumD = np.sum(np.square(ocD[np.triu_indices(len(ocD),1)]))
        print mds.stress_/sqsumD
    
        pts = ocPts
        x,y = zip(*ocPts )
        x = np.array(x)
        y = np.array(y)
        labels = labels[onlyclustered]    
        
        core_samples_mask = core_samples_mask[onlyclustered]
    
    plot_clusters(None, labels, core_samples_mask=core_samples_mask, pts=pts)    
        
    #    x,y = zip(*pts )
    #    x = np.array(x)[onlyclustered]
    #    y = np.array(y)[onlyclustered]
    #plt.scatter(x,y, c=labels)
    #plt.colorbar()
    #colors = np.random.rand(N)
        #area = np.pi * (15 * np.random.rand(N))**2  # 0 to 15 point radiuses
        #plt.(x, y, s=area, c=colors, alpha=0.5)
    #plt.show()
        
    