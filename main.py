# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 13:01:02 2016

@author: juherask
"""

from extractors.features_vrp_autocorrelation import autocorrelation_features
from extractors.features_vrp_bcp import bcp_features
from extractors.features_vrp_bottleneck import bottleneck_features_MST
from extractors.features_vrp_constraint import constraint_features
from extractors.features_vrp_geometric import geometric_features
from extractors.features_vrp_lsp import lsp_features
from extractors.features_vrp_mst import mst_features
from extractors.features_vrp_node_distribution import node_distribution_features, cluster_nodes
from extractors.features_vrp_nn import nn_features

from helpers.cvrp_rkm16_io import read_TSPLIB_CVRP, calculate_D, read_TSBLIB_additional_constraints
from helpers.cvrp_rkm16_ops import D2D_c
from helpers.rkm_stats import get_distribution_statistics_labels
from helpers.rkm_util import normalize_to_rect, produce_nn_list

import numpy as np

from glob import glob
from os import path
from math import ceil
import time
import re

#import warnings
#warnings.simplefilter("error")

regexp_k_from_filename = re.compile("[_-]k([0-9]+)")

CALC_AUTOCORRELATION = 1
CALC_BRANCH_AND_CUT_PROBING = 2
CALC_BOTTLENECK_COST = 4
CALC_CONSTAINT =  8
CALC_GEOMETRIC = 16
CALC_LOCAL_SEACH_PROBING = 32
CALC_MINIMUM_SPANNING_TREE = 64
CALC_NEAREST_NEIGHBOR = 128
CALC_NODE_DISTRIBUTION = 256


def calculate_features(vrp_file, task_falgs, flatten_feature_list=False):
    #name+".vrp"
    N, pts, dd, demands, D, C, _ = read_TSPLIB_CVRP(vrp_file)
    K, L, service_time = read_TSBLIB_additional_constraints(vrp_file)
    
    # Embed service times to D (if needed)
    D = D if service_time is None else D2D_c(D, service_time)
    
    D_and_pts_match = False

    # pts and D come from the file, everything is A-OK
    if pts is not None and D is not None:
        D_and_pts_match = True
    # We do not have point coodrinates, but we have D!
    # -> Use MDS to get an approximation for the point locations  
    if ((pts is None and dd is None) and D is not None):
        from sklearn import manifold
        mds = manifold.MDS(n_components=2, dissimilarity='precomputed')
        results = mds.fit(D)
        pts = results.embedding_
        D_and_pts_match = True
    # The pts are not given, but display points are, use them as pts
    #  however, we must make the D match later on...
    if pts is None and dd is not None:
        D_and_pts_match = False
        pts = dd
        
    # normalize as in [SM10], but keep the scale
    npts = normalize_to_rect(pts, (0,400), keep_aspect_ratio=True)
    
    if (D is None):
        D_and_pts_match = True
        nD = calculate_D(npts, None, 'EXACT_2D') 
        
    # We have normalized the points, but now the predefined distance matrix
    #  does not match. It may have weird rounding rules etc. so we just scale
    #  it by same amount we scaled the points in euclidian space.
    else:
        refnD = calculate_D(npts, None, 'EXACT_2D')
        mean_distance_value_refnD = np.mean( refnD[np.triu_indices(len(refnD),1)] )
        if D_and_pts_match:
            origD = calculate_D(pts, None, 'EXACT_2D')                       
            mean_distance_value_orginalD = np.mean( origD[np.triu_indices(len(origD),1)] )
            scale_factor = mean_distance_value_refnD/mean_distance_value_orginalD
        else:
            mean_distance_value_orginalD = np.mean( D[np.triu_indices(len(D),1)] )
            scale_factor = mean_distance_value_refnD/mean_distance_value_orginalD
            #print scale_factor, mean_distance_value_refnD, mean_distance_value_orginalD
        
        nD = D.copy().astype(float)
        if type(nD) is not np.ndarray:
            nD = np.array(nD)
        # Scale the distances in the distance matrix
        nD*=scale_factor 
    
    ## Precalculate some data structures ##    
    # TODO: time these?
    if task_falgs & CALC_CONSTAINT or\
       task_falgs & CALC_NODE_DISTRIBUTION:   
        clusters, core_pt_idxs = cluster_nodes(npts, nD)
    if task_falgs & CALC_BOTTLENECK_COST or\
       task_falgs & CALC_NEAREST_NEIGHBOR:  
        n_NN_D = produce_nn_list(nD)
    k=None
    soft_k = False
    if task_falgs & CALC_BRANCH_AND_CUT_PROBING or\
       task_falgs & CALC_CONSTAINT:
        m = regexp_k_from_filename.search( path.basename(vrp_file) )
        k = None
        if m:
            k = int(m.group(1)) # from filename
        elif K:
            k = K # from file VEHICLES / NUMBER_OF_TRUCKS field
        else:
            k = ceil(sum(demands)*1.05/C) # estimate k by giving 5 % slack capacity
            soft_k = True
    
    all_features = []
    timing_features = []
    upper_bound  = None
    
    all_features.append( ("DC1", "Size of the Problem [SM10]", N) )
    
    if task_falgs & CALC_AUTOCORRELATION:
        start = time.clock()
        all_features += autocorrelation_features(npts, nD, demands, C)
        elapsed = time.clock()-start
        timing_features.append( ("T8", "Timing for Autocorrelation Features [Hu14]", elapsed ) )
    
    if task_falgs & CALC_LOCAL_SEACH_PROBING:
        start = time.clock()
        if task_falgs & CALC_BRANCH_AND_CUT_PROBING:    
            calculated_lsp_features, upper_bound, upper_bound_k = \
                lsp_features(vrp_file, pts, nD, demands, C, repeats=20,
                             return_upper_bound=True,
                             prioritize_upper_bound_k=None if soft_k else k)
            if soft_k and upper_bound_k is not None:
                k = upper_bound_k
            all_features += calculated_lsp_features
        else:
            all_features += lsp_features(vrp_file, pts, nD, demands, C, repeats=20)    
        elapsed = time.clock()-start
        timing_features.append( ("T3", "Timing for Local Search Probing Features [Hu14]", elapsed ) )

    if task_falgs & CALC_BRANCH_AND_CUT_PROBING:
        start = time.clock()
        if task_falgs & CALC_LOCAL_SEACH_PROBING:
            all_features += bcp_features(vrp_file, N, k, time_limit=3.0, upper_bound=upper_bound)
        else:
            all_features += bcp_features(vrp_file, N, k, time_limit=3.0)
        elapsed = time.clock()-start
        timing_features.append( ("T4", "Timing for Branch-and-Cut Probing Features [Hu14]", elapsed ) )

    if task_falgs & CALC_CONSTAINT:    
        start = time.clock()
        # note! it is important to use the non-scaled, non-normalized D for accuracy!
        all_features += constraint_features(demands, C, clusters, D, L, k=k)    
        elapsed = time.clock()-start
        timing_features.append( ("T7", "Timing for Constraint Features [Hu14]", elapsed ) )
        
    if task_falgs & CALC_NODE_DISTRIBUTION:
        start = time.clock()
        all_features += node_distribution_features(npts, nD, clusters, core_pt_idxs)    
        elapsed = time.clock()-start
        timing_features.append( ("T1", "Timing for Node Distribution Features [Hu14]", elapsed ) )
                
    if task_falgs & CALC_GEOMETRIC:
        start = time.clock()    
        all_features += geometric_features(npts)
        elapsed = time.clock()-start
        timing_features.append( ("T5", "Timing for Geometric Features [Hu14]", elapsed ) )
                 
    MST_graph = None
    if task_falgs & CALC_MINIMUM_SPANNING_TREE:
        start = time.clock()        
        calculated_mst_features, MST_graph = mst_features(npts, nD, root_note_idx=0, return_MST_graph=True)
        all_features += calculated_mst_features
        elapsed = time.clock()-start
        timing_features.append( ("T2", "Timing for Minimum Spanning Tree Features [Hu14]", elapsed ) )
    
    if task_falgs & CALC_BOTTLENECK_COST:    
        start = time.clock()
        # reuse the MST graph from MST feature computation if available
        all_features += bottleneck_features_MST(nD, MST_graph)
        elapsed = time.clock()-start
        timing_features.append( ("T9", "Timing for Bottleneck Features [Hu14]", elapsed ) )

    if task_falgs & CALC_NEAREST_NEIGHBOR:
        start = time.clock()        
        all_features += nn_features(npts, nD, n_NN_D)        
        elapsed = time.clock()-start
        timing_features.append( ("T6", "Timing for Nearest Neighbour Features [Hu14]", elapsed ) )
    
    

    #TODO: Compare to [Hu14] timing features (and regroup?)
    all_features += timing_features
        
    if not flatten_feature_list:
        return all_features
    
    FEATURE_CODE = 0    
    FEATURE_NAME = 1
    FEATURE_VALUE = 2
    flattened_features = []
    statistics_labels = get_distribution_statistics_labels()
    frac_labels = ["1 decimal","2 decimals","3 decimals","4 decimals",]
    for feature in all_features:
        try:
            num_f = len(feature[FEATURE_VALUE])
        except:
            num_f = 1
        if num_f>1:
            labels = statistics_labels
            if "(x,y)" in feature[FEATURE_NAME]:
                labels = ["x","y"]
            elif "Fraction of Rounded Distinct Distances" in feature[FEATURE_NAME]:
                labels = frac_labels
            elif not "SotD" in feature[FEATURE_NAME]:
                print  feature[FEATURE_NAME]
                raise ValueError("all features with many values should be Distribution statistics")
            for i in range(num_f):
                flattened_features.append(
                    (feature[FEATURE_CODE]+"-%d"%(i+1),
                    feature[FEATURE_NAME]+" "+labels[i], #name
                    feature[FEATURE_VALUE][i]) #value
                )
        else:
            flattened_features.append(feature)
            
    return flattened_features

def main():
    # All of the .vrp files in these folders will be processed
    folders = [("Examples", "./instances/")]
    
    # If any filenames (not paths!) are listed here, only those are processed
    #  use it to process, e.g., just one file.
    whitelist = [] #["E-n13-k4.vrp"]
    
    blacklist = [#"att-n48-k4.vrp"
        #"A-n37-k5.vrp", #Degrees of freedom <= 0 for slice -> only one cluster, gets values still
        #"B-n31-k5.vrp", #Mean of empty slice. -> Calculation in Silhoulette, produces a value nontheless
        #"E-n22-k4.vrp", # Degrees of freedom <= 0 for slice
        #"P-n19-k2.vrp", # RuntimeWarning: invalid value encountered in double_scalars
    ]

    # Select which features to calculate for each of the problem instances    
    #features_to_calculate = CALC_LOCAL_SEACH_PROBING | CALC_BRANCH_AND_CUT_PROBING    
    features_to_calculate = \
        CALC_BOTTLENECK_COST |\
        CALC_CONSTAINT |\
        CALC_GEOMETRIC |\
        CALC_MINIMUM_SPANNING_TREE |\
        CALC_NODE_DISTRIBUTION |\
        CALC_AUTOCORRELATION |\
        CALC_LOCAL_SEACH_PROBING |\
        CALC_BRANCH_AND_CUT_PROBING |\
        CALC_NEAREST_NEIGHBOR
                
    with open("results/features.csv", "w") as outfile:
        first = True
        for pset, pfolder in folders:
            search_path = path.join(pfolder,r"*.vrp")
            #print search_path
            for vrp_file in glob(search_path):
                
                file_basename = vrp_file.split('\\')[-1]
                if whitelist:
                    if file_basename not in whitelist:
                        continue
                elif file_basename in blacklist:
                    continue
                
                print "processing", vrp_file
                #try:
                fl = calculate_features(vrp_file, features_to_calculate, flatten_feature_list=True)
                print "processing done!"
                #except IOError as e:
                #    print str(e)
                #    continue
                    
                flc, fln, flv = zip(*fl)
                if first:
                    outfile.write("Problem Set; Problem Name;")
                    outfile.write( ";".join(fln) )
                    outfile.write("\n")
                    outfile.write("Problem Set; Problem Name;")
                    outfile.write( ";".join(flc) )
                    outfile.write("\n")                
                    first=False
                outfile.write(pset+";"+file_basename.replace(".vrp", "")+";")
                outfile.write( ";".join((str(v) for v in flv)) )
                outfile.write("\n")
                
                print ""
        
    print "DONE!"
    
    
if __name__=="__main__":
    np.seterr(all='raise')
    main()
    
        
            
