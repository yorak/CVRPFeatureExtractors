# -*- coding: utf-8 -*-
"""
Created on Wed Mar 09 11:22:51 2016

@author: juherask

Based on:
[Hu14] Hutter, F., Xu, L., Hoos, H.H., Leyton-Brown, K.: Algorithm runtime
prediction: Methods & evaluation. Artificial Intelligence206(2014) 79–111
"""

from helpers.rkm_stats import calculate_distribution_statistics
from scipy.spatial import ConvexHull
from math import sqrt
import numpy as np

X=0;Y=1; #~Defines for readablilty

# Shoelace formula
def _poly_area(corners):
    n = len(corners) # of corners
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area
    
def _dist_sqd(v, w):
    return (v[X] - w[X])**2 + (v[Y] - w[Y])**2 

# a basic python implementation of distance to squared. 
#  vectorization and changing the modus operandi to point to sengment_list
#  should make this faster.
#  using a module like shapely could be more computationally efficient
#  http://stackoverflow.com/questions/24415806/coordinate-of-the-closest-point-on-a-line
def _dist_to_segment_squared(p, v, w):   
    l2 = _dist_sqd(v, w)
    if l2 == 0:
        return _dist_sqd(p, v)
    t = ((p[X] - v[X]) * (w[X] - v[X]) + (p[Y] - v[Y]) * (w[Y] - v[Y])) / l2
    t = max(0, min(1, t))
    p2 = (v[X] + t * (w[X] - v[X]), v[Y] + t * (w[Y] - v[Y]))
    return _dist_sqd(p, p2)

def geometric_features(npts):
    features = []

    # Convex hull features    
    
    xy = np.array(npts)
    hull = ConvexHull(xy)    
    hull_corner_points = xy[hull.vertices]
    
    features.append( ("G2", "Area for Convex Hull containing All Nodes [Me13]", _poly_area(hull_corner_points) ) )
    
    # How many nodes on the hull line and what are the distances of inner nodes
    nodes_on_hull = 0
    distance_to_hull_for_internal_nodes = []    
    for p in xy:        
        end_segment = hull_corner_points[-1]
        dists = []
        for start_segment in hull_corner_points:
            dists.append( _dist_to_segment_squared(p, start_segment, end_segment) )
            end_segment = start_segment
        dist = sqrt(min(dists))
        if dist < 0.00001:
            nodes_on_hull+=1
        else:
            distance_to_hull_for_internal_nodes.append(dist)
    features.append( ("G3", "Ratio of Nodes on the Convex Hull [PM14]",
                      float(nodes_on_hull)/len(xy) ) )
    features.append( ("G4", "SotD of Distance of Inner Nodes to the Convex Hull [PM14]", 
                      calculate_distribution_statistics(distance_to_hull_for_internal_nodes, level=1) ) )
    
    # edge lenghts of the convex hull
    end_segment = hull_corner_points[-1]
    convex_hull_elenghts = []
    for start_segment in hull_corner_points:
        convex_hull_elenghts.append( sqrt(_dist_sqd(start_segment, end_segment) ) )
        end_segment = start_segment
    features.append( ("G5", "SotD of Edge Lenghts of the Convex Hull [PM14]",
                      calculate_distribution_statistics(convex_hull_elenghts, level=1) ) )
          
    return features
  
    
    
def test():
    from cvrp_rkm16_io import read_TSPLIB_CVRP
    (N, pts, dd, demands, D, C) = read_TSPLIB_CVRP("../instances/A-n32-k5.vrp")

    #print pts[0]
    print geometric_features(pts)

    #Synthetic
    spts = [[2.0,2.0],
            [2.0,5.0],#~on segment
            [1.0,4.0],
            [4.0,2.0],#on segment
            [4.0,7.0],
            [4.0,4.0],#dead center d=2?
            [7.0,4.0],
            [5.0,5.0],#quite close to a segm.
            [5.0,2.0]]
    print geometric_features(spts)    
    
    #Random
    ra = np.random.rand(2000,2)
    print geometric_features(ra)
    
    
if __name__=="__main__":
    test()