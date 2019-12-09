# -*- coding: utf-8 -*-
"""
Created on Wed Mar 09 11:22:51 2016

@author: juherask

Based on:
[Hu14] Hutter, F., Xu, L., Hoos, H.H., Leyton-Brown, K.: Algorithm runtime
prediction: Methods & evaluation. Artificial Intelligence206(2014) 79–111

The feature introduced in
Neto, D.: Efficient Cluster Compensation for Lin-Kernighan Heuristics. PhD thesis,
Computer Science, University of Toronto (1999)

"""


from helpers.rkm_stats import calculate_distribution_statistics
from helpers.rkm_util import produce_nn_list
from dijkstra import dijkstra_path

import numpy as np

from heapq import heappop, heappush
from collections import defaultdict

def calc_bottlenecks_using_nn_lists(to_node, D, NNl):
    """
    This algorithm borrows some ideas from Dijkstra's algorithm and calculates 
    the minimum bottleneck cost between target node and every other node.
    The bottleneck cost of a path is defined as the largest cost along the path.
    We get the mimimal cost by considering any path between them.
    
    This implementation is O(n^3)
    
    TODO: For large instances (n>100) this takes quite long. Consider 
     some optmizations such as reusing / caching caclucations.
     
    >>> D = [[ 0.,  4.,  4.,  6.],
    ...      [ 4.,  0.,  7.,  5.],
    ...      [ 4.,  7.,  0.,  4.],
    ...      [ 6.,  5.,  4.,  0.]]
    >>> NNl = produce_nn_list(D)
    >>> for i in range(len(D)):
    ...   print calc_bottlenecks_using_nn_lists(i, D, NNl)
    ...
    [0.0, 4.0, 4.0, 4.0]
    [4.0, 0.0, 4.0, 4.0]
    [4.0, 4.0, 0.0, 4.0]
    [4.0, 4.0, 4.0, 0.0]

    >>> D = [[  0.,   6.,   4.,   5.],
    ...      [  6.,   0.,   4.,  10.],
    ...      [  4.,   4.,   0.,   6.],
    ...      [  5.,  10.,   6.,   0.]]
    >>> NNl = produce_nn_list(D)
    >>> for i in range(len(D)):
    ...   print calc_bottlenecks_using_nn_lists(i, D, NNl)
    ...
    [0.0, 4.0, 4.0, 5.0]
    [4.0, 0.0, 4.0, 5.0]
    [4.0, 4.0, 0.0, 5.0]
    [5.0, 5.0, 5.0, 0.0]
    
    
    """
    n = len(NNl)
    bottleneck_table = [D[to_node][i] for i in range(n)]
    # prev_table makes sure we leave from the tgt at least once, and that 
    #  if there are 0.0 length edges, we do not loop!
    prev_table = [None]*n
    prev_table[to_node]=to_node
    bottlenecks_sorted = sorted(bottleneck_table)
    max_bottleneck = bottlenecks_sorted.pop()

    Q = [(0, to_node)]
    while Q:
        (bottleneck_cost, v1) = heappop(Q)
        
        print "max bn", max_bottleneck 
        print "edge", (bottleneck_cost, v1) 
        print "table", bottleneck_table
        print "... "
        
        for v2, edge_cost in NNl[v1][1:]:
            if edge_cost>max_bottleneck:
                break #everything from now on are bigger than this
            
            if edge_cost>bottleneck_cost:
                current_bn  = edge_cost
            else:
                current_bn = bottleneck_cost
            v2_bn = bottleneck_table[v2]
            
            # makes sure we visit at least once +
            #  avoids looping with edges of 0.0 length
            if prev_table[v2]==None or \
                (current_bn<v2_bn and prev_table[v2]!=v1): 
                    
                bottleneck_table[v2] = current_bn
                prev_table[v2] = v1
                
                # replaced the node with largest bottleneck value
                #  lower the limit of max edge weight worth looking into
                if v2_bn==max_bottleneck:
                    max_bottleneck = bottlenecks_sorted.pop()
                
                if (edge_cost>0.0):
                    heappush(Q, (current_bn, v2))
    
    print "... EOL ... "
    return bottleneck_table

def calc_bottlenecks_using_mst(to_node, D, MBST_graph):
    """
    We use the property that MST is MBST
    the CVRP is symmetric, so it makes no difference if we start from the end
    
    >>> D = np.array([[ 0.,  4.,  4.,  6.],
    ...      [ 4.,  0.,  7.,  5.],
    ...      [ 4.,  7.,  0.,  4.],
    ...      [ 6.,  5.,  4.,  0.]])
    >>> G = _create_MST_graph(D)
    >>> for i in range(len(D)):
    ...   print calc_bottlenecks_using_mst(i, D, G)
    ...
    [0.0, 4.0, 4.0, 4.0]
    [4.0, 0.0, 4.0, 4.0]
    [4.0, 4.0, 0.0, 4.0]
    [4.0, 4.0, 4.0, 0.0]
    """
    N = len(D)
    bottleneck_table = [0.0]*N
    
    prevs = dijkstra_path(MBST_graph,to_node)
    for from_node in xrange(N):
        if from_node!=to_node:
            prev_node = prevs[from_node]
            bottleneck = D[prev_node,from_node]
            while prev_node!=to_node:
                next_node =  prevs[prev_node]
                edgewt = D[prev_node,next_node]
                if edgewt>bottleneck:
                    bottleneck = edgewt
                prev_node = next_node
            bottleneck_table[from_node]=bottleneck
    return bottleneck_table
    

def bottleneck_features_NN(nD, n_NN_D=None):
    """ n_NN_D is the normalized (if you want), nearest neighbour distances,
    where each row contains a list of NN's for that point. The NN's are 
    listed as pairs with (node_id, distance). If this is not given, it
    is calcluated using util.produce_nn_list(nD)
    
    Alternatively a minimum bottleneck spanning tree can be used. This is 
     handy especially if MST is calculated as each MST is also MBST.
    """
    n = len(nD)
    
    if n_NN_D is None:
        n_NN_D = produce_nn_list(nD)
            
    min_bottleneck_costs=[]
    # Assume symmetric D. If not, run to n
    for to_node in range(n-1):
        #print "to node i", i
        costs = calc_bottlenecks_using_nn_lists(to_node,nD,n_NN_D)
        # upper corner of the MBNC -matrix
        min_bottleneck_costs += costs[to_node+1:] # is symmetric, we need only half
    
    #print min_bottleneck_costs
    return [ ("ND10", "SotD of Minimum Bottleneck Cost [Hu14]", calculate_distribution_statistics(min_bottleneck_costs)) ]

def _create_MST_graph(D):
    from prim_mst import minimum_spanning_tree
    MST_edges = minimum_spanning_tree(D)  
    # convert to graph
    MST_graph = defaultdict(dict)
    for i,j in MST_edges:
        MST_graph[i][j]=D[i][j]
        MST_graph[j][i]=D[i][j]
    return MST_graph
    
def bottleneck_features_MST(nD, MST_graph=None):
    """ n_NN_D is the normalized (if you want), nearest neighbour distances,
    where each row contains a list of NN's for that point. The NN's are 
    listed as pairs with (node_id, distance). If this is not given, it
    is calcluated using util.produce_nn_list(nD)
    
    Alternatively a minimum bottleneck spanning tree can be used. This is 
     handy especially if MST is calculated as each MST is also MBST.
    """
    n = len(nD)
    
    if MST_graph is None:
        MST_graph = _create_MST_graph(nD)
            
    min_bottleneck_costs=[]
    # Assume symmetric D. If not, run to n
    for i in range(n-1):
        #print "to node i", i
        costs = calc_bottlenecks_using_mst(i,nD,MST_graph )
        # upper corner of the MBNC -matrix
        min_bottleneck_costs += costs[i+1:] # is symmetric, we need only half
    
    #print min_bottleneck_costs
    return [ ("ND10", "SotD of Minimum Bottleneck Cost Between All Node Pairs [Hu14]", calculate_distribution_statistics(min_bottleneck_costs)) ]
                       

def test_timeit_random(r, test_MST=True, test_NN=False):
    for i in range(100):
        from cvrp_rkm16_io import generate_CVRP
        (N, pts,demands, D, C) = generate_CVRP(r, 100.0, 5.0, 2.0, R=5)
        if test_MST:
            bottleneck_features_MST(pts, D)
        if test_NN:
            bottleneck_features_MST(pts, D)

               
def smoke_test():
    from cvrp_rkm16_io import read_TSPLIB_CVRP, generate_CVRP
    (N, pts, dd, demands, D, C) = read_TSPLIB_CVRP("../instances/A-n32-k5.vrp")
    print "bottleneck_features for A-n32-k5"
    print bottleneck_features_NN(D)
    print
    print "bottleneck_features of the first 3 points of A-n32-k5 with fixed D"
    print bottleneck_features_MST(#A B C D
                       np.array([[0.0,3.0,4.0,1.0],
                        [3.0,0.0,1.0,5.0],
                        [4.0,1.0,0.0,2.0],
                        [1.0,5.0,2.0,0.0]]))

    TEST_SAMENESS = True
    if TEST_SAMENESS:
        print
        print "MST and NN are equvalen test:"
        for i in range(100):
            
            (N, pts,demands, D, C) = generate_CVRP(3, 100.0, 5.0, 2.0, R=5)
            MST_stat =  bottleneck_features_NN(D)[0][2]
            NN_stat =  bottleneck_features_MST(D)[0][2]
            if MST_stat!=NN_stat:
                from pprint import pprint
                print D
                pprint(pts)
                print MST_stat
                print NN_stat
                print "DIFFERENT!"      
                exit()
            else:
                print "SAME"
                
def timeit_test(): 
    TEST_TIMEIT = False
    if TEST_TIMEIT:
        print 
        print "timeit test:"
        import timeit
        for r in [2**p for p in range(3,8)]:
            #test_timeit_random(r)
            timeit_report = timeit.timeit("""test_timeit_random(%d)"""%r, "from __main__ import test_timeit_random, bottleneck_features_MST, bottleneck_features_NN", number=1)
            print r, timeit_report
    
    
    
if __name__=="__main__":
    import doctest
    doctest.testmod()
    
    #smoke_test()
    
    