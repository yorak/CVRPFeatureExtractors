# -*- coding: utf-8 -*-
"""
Created on Wed Mar 09 11:22:51 2016

@author: juherask

Local search probing with VRPH
"""

from helpers.rkm_stats import calculate_distribution_statistics
from scipy import spatial
from scipy.stats.mstats import mquantiles
from subprocess import check_output, CalledProcessError
from random import uniform 
import numpy as np
import sys
import os # to remove tmp file


SOLUTION_FILE = 'tmp_solution.txt'
SOLVER = os.path.join(os.path.dirname(__file__), '../solvers/vrp_init')
X=0;Y=1; # ~define

DEBUG_LSP_TO_FILE = False

# VRPH can fail with some random lambdas
#  this is very rare, so just retry if this happens
RETRY_ON_FAILURE_TIMES = 10

#NOTE: the vrph_init should be compiled with the CLEAN_DEBUG that prints
# extra output on local search step
vprh_solver = SOLVER
if sys.platform == 'win32':
    vprh_solver+='.exe'
    
    
def _get_tour_edges(tour, D):
    """ The D is already normalized, just return the edge lenghts.
    """
    return [D[tour[i]][tour[i+1]] for i in range(len(tour)-1)]
    
def _normalize_tour_edges(tour, D):
    """ Normalize the tour edge lengths to total tour length. The tour can be 
    a route or a giant tour encoded solution.
    """
    tour_edges = [D[tour[i]][tour[i+1]] for i in range(len(tour)-1)]
    tour_len = float(sum(tour_edges))
    # special case that is eg. in F-n45-k4 where there is a tour with just one 
    #  node overlapping the depot
    if tour_len==0.0:
        return tour_edges
    else:
        return [e_len/tour_len for e_len in tour_edges]  
    
def _normalize_solution_values(vals, bks, floorsol=0):
    return [(val-floorsol)/float(bks) for val in vals]

# Taken from http://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
#  might have a problem with colinarity
def _ccw(A,B,C):
    return (C[Y]-A[Y]) * (B[X]-A[X]) > (B[Y]-A[Y]) * (C[X]-A[X])
def _intersect(A,B,C,D):
    return _ccw(A,C,D) != _ccw(B,C,D) and _ccw(A,B,C) != _ccw(A,B,D)
    
def _analyze_tour_for_intersections( tour, pts ):
    """ Find if any two edges in the tour intersect on the plane. The tour node
    coodrinates are given in <pts> indexed by the tour indices. The tour can be 
    a route or a giant tour encoded solution.
    """
    intersections = 0
    for i in range(1, len(tour)):
        for j in range(1, len(tour)):
            edge1 = (pts[tour[i-1]],pts[tour[i]])
            edge2 = (pts[tour[j-1]],pts[tour[j]])
            
            shared_endpoint = \
              (edge1[0][X]==edge2[0][X] and edge1[0][Y]==edge2[0][Y]) or\
              (edge1[1][X]==edge2[0][X] and edge1[1][Y]==edge2[0][Y]) or\
              (edge1[0][X]==edge2[1][X] and edge1[0][Y]==edge2[1][Y]) or\
              (edge1[1][X]==edge2[1][X] and edge1[1][Y]==edge2[1][Y])
            
            if not shared_endpoint and\
                _intersect(edge1[0],edge1[1], edge2[0],edge2[1]):
               intersections+=1
    return intersections
        

def _analyze_tour_edges( 
    #in 
    ntour_edge_lens,
    #out
    quantile_edges,
    tour_segment_edge_counts,
    tour_segment_lengths,
    tour_segment_edge_lengths ):
    """ Processes a tour (a route)  and
    distributes edges into quantiles based on their edge lenghts/weights. 
    Also splits the tour into segments that are left behind when the long edges
    are cut away from the tour accorting to following rules:
    
      1) Top 5% of the longest edges of the tour are removed
      2) Everything over 1.5 x the shortest edge in 4th quartile are removed
    
    The quantile edges and tour segment edges, counts and edge lenghts are 
    added into the output parameters.
    """
    
    # Gather data for:                
    # - LSP Edge lenghts in quartiles
    # - LSP Tour segments
    #print tour
    
    e_cnt = len(ntour_edge_lens)
    ntour_edge_lens_sorted = sorted( zip(ntour_edge_lens, range(e_cnt)) )
    
    q_ends = [ntour_edge_lens_sorted[0][0],]+\
             list(mquantiles(ntour_edge_lens))+\
             [ntour_edge_lens_sorted[-1][0],]

    edge_idx = 0
    quantile_idx = 0
    first_eq_dist_quantile_idx = 0
    shortest_q4_edge_wt = None
    while edge_idx<e_cnt:
        quantile_start = q_ends[quantile_idx]
        quantile_end = q_ends[quantile_idx+1]
        edge_wt = ntour_edge_lens_sorted[edge_idx][0]
        # record shortest of the 75% quantile for segment computation below
        if quantile_idx==3 and shortest_q4_edge_wt is None:
            shortest_q4_edge_wt = ntour_edge_lens_sorted[edge_idx][0]
            cut_rule_shortest_q4_edge_wt = shortest_q4_edge_wt*1.5
        # record the edge to a active quantile
        if edge_wt<=quantile_end:
            quantile_edges[quantile_idx].append( edge_wt )
            edge_idx+=1
            if edge_idx==e_cnt:
                break
            next_edge_wt = ntour_edge_lens_sorted[edge_idx][0]           
            # if it could belong to any quantile, distribute evenly
            if edge_wt==next_edge_wt:
                if quantile_start==next_edge_wt or quantile_end==next_edge_wt:
                    quantile_idx+=1
                    if quantile_idx==4 or next_edge_wt>q_ends[quantile_idx+1]:
                        quantile_idx = first_eq_dist_quantile_idx
                
            else:
                first_eq_dist_quantile_idx = quantile_idx
        else:
            quantile_idx+=1
            first_eq_dist_quantile_idx = quantile_idx
            
    
    # These are the edges to cut away from tour to get segements
    #  use two rules:
    #  1) Top 5% of the longest edges are removed
    #  2) Everything over 1.5 x the shortest edge in 4th quartile is removed
    cut_5p_rule_from_idx = int(round(e_cnt*0.95))
    cut_away_edge_indices = set( i for wt,i in ntour_edge_lens_sorted[cut_5p_rule_from_idx:] )
    if shortest_q4_edge_wt:
        cut_away_edge_indices.update( set(i for wt, i in ntour_edge_lens_sorted
                                          if wt>=cut_rule_shortest_q4_edge_wt) )
    
    segment_length = 0
    segment_edge_count = 0
    for ntour_edge_idx, ntour_edge_len in enumerate(ntour_edge_lens):
        if ntour_edge_idx in cut_away_edge_indices:
            if segment_edge_count>0:
                tour_segment_edge_counts.append(segment_edge_count)
                tour_segment_lengths.append(segment_length)
            segment_edge_count = 0
            segment_length = 0
        else:
            tour_segment_edge_lengths.append(ntour_edge_len)
            segment_edge_count += 1
            segment_length += ntour_edge_len
    if segment_edge_count>0:
        tour_segment_edge_counts.append(segment_edge_count)
        tour_segment_lengths.append(segment_length)
    
def lsp_features(problem_file_name, pts, D, qs, Q, repeats=20, normalize=True,
                 return_upper_bound=False, prioritize_upper_bound_k=None):
    """ qs are the demands (with 0 being the depot with demand 0f 0)
        Q is the capacity"""
    
    features = []
        
    initial_solution_values = []
    ls_LO_solution_values = [] # LO = local optima
    ls_improvements = []
    step_compositions = []
    ls_steps_to_LO = []
    best_known_solution_value = None
    problem_size = None
    N = None
    local_optima = None
    upper_bound = None
    upper_bound_k = None    
    
    LO_solutions = []

    quantile_edges = [[],[],[],[]]
    tour_segment_edge_counts = []
    tour_segment_lengths = []
    tour_segment_edge_lengths = []
    normalized_route_demands = []
    normalized_nodes_per_route = []

    intersections = []    
    
    for r in range(repeats):
        
        retries_left = RETRY_ON_FAILURE_TIMES
        while retries_left>0:
            random_lambda = str(uniform(0.5,2.0))
            sover_cmd_w_args = [vprh_solver,
                  '-f', problem_file_name,
                  '-m', '0', # Clarke-Wright construction heuristic
                  '-l', random_lambda, #randomize scaling factor %lambda for CW
                  '-c', #use local search to do intra-route LS with ONE_POINT_MOVE+TWO_POINT_MOVE+TWO_OPT
                  '-t', SOLUTION_FILE, # write routes to a file
                  '-h', "ONE_POINT_MOVE",
                  '-h', "TWO_POINT_MOVE",
                  '-h', "TWO_OPT",
                  '-h', "THREE_OPT"
                  ]
            
            try:
                vrph_output = check_output(sover_cmd_w_args)
            except CalledProcessError:
                retries_left-=1;
            else:
                break # 
        
        #print " ".join(sover_cmd_w_args)
        if DEBUG_LSP_TO_FILE:
            print " ".join(sover_cmd_w_args)
            problem_id = os.path.basename(problem_file_name).replace(".vrp", "")
            debug_output_file = file(os.path.join(os.path.dirname(__file__), '../solvers/debug_vrph_out_%s_(%s).txt'%(problem_id,random_lambda)  ), 'w')
            debug_output_file.writelines(vrph_output)
        
              
        # Process VRPH output to gather statistics
        #  remember that CLEAN_DEBUG have to be enabled!
        
        current_ls_value = None 
        nbr_improving_steps = 0
        route_count = 0
        
        for line in vrph_output.splitlines():
            if r"CLEAN::" in line:
                if r"end_rlen" in line:
                    ls_value = float(line.split("=")[1])
                    if ls_value<=current_ls_value:
                        ls_improvements.append(current_ls_value-ls_value)
                    if ls_value<current_ls_value:
                        current_ls_value=ls_value
                        nbr_improving_steps+=1
                elif r"start_val" in line:
                    current_ls_value = float(line.split("=")[1]) 
                elif r"tried_num" in line:
                    parts = line.split()
                    tried_num = int(parts[1].split("=")[1])
                    better_num = int(parts[2].split("=")[1])
                    best_num = int(parts[3].split("=")[1])
                    moves_num = int(parts[4].split("=")[1])
                    step_compositions.append( (tried_num, better_num, best_num, moves_num) )
            elif r"Total route length before clean up:" in line:
                initial_solution_values.append( float(line.split(":")[1]) )
            elif r"Total route length:" in line:
                local_optima = float(line.split(":")[1])
                ls_LO_solution_values.append( local_optima )
            elif best_known_solution_value==None and "Best known solution: " in line:
                best_known_solution_value = float(line.split(":")[1])
            elif r"Number of nodes visited" in line:
                problem_size = int(line.split(":")[1])
            elif r"Route " in line[:7]:
                route_count +=1
                
        ls_steps_to_LO.append(nbr_improving_steps)
        
        # only allow those local_optimas to be a UB that satisfy the number
        #  of routes condition
        if prioritize_upper_bound_k==None:
            if upper_bound is None or local_optima<upper_bound:
                upper_bound = local_optima
                upper_bound_k = route_count
        else:
            # prioritize_upper_bound_k is set. Thus, set a new upper bound if ...
            # 1) no previous upper bound is set OR
            # 2) first or better with the preferred upper number of routes is found
            # 3) not 
            no_upper_bound_set = upper_bound_k is None
            best_preferred_upper_bound = (prioritize_upper_bound_k==route_count and
                                          (upper_bound_k!=prioritize_upper_bound_k
                                           or local_optima<upper_bound) )
            best_non_preferred_upper_bound = (prioritize_upper_bound_k!=upper_bound_k
                                              and local_optima<upper_bound)
            if no_upper_bound_set or\
               best_preferred_upper_bound or\
               best_non_preferred_upper_bound:
               
                upper_bound = local_optima
                upper_bound_k = route_count
                
    
        # read the solution file
        N = problem_size+1
        solution = np.zeros((N, N), dtype=int)
        with open(SOLUTION_FILE, 'r') as solfile:
            # first value is the size of the problem, ignore it, but leave 
            #  0 at the end, which is the "end solution" 
            solfile_routes = [int(p) for p in solfile.readline().split()[1:]]
            
            from_pt_idx = 0
            route = None
            gt_solution = []
            
            route_demand = 0.0
            nodes_per_route = 0
            for pt_idx in solfile_routes:
                # negative index is new route, 0 is the end of the solution
                if pt_idx<=0:
                    gt_solution.append(0)
                    
                    # close and analyze previous route
                    if route:
                        if Q:
                            normalized_route_demands.append(route_demand/float(Q))
                        else:
                            normalized_route_demands.append(0)
                            
                        normalized_nodes_per_route.append(nodes_per_route/float(N))
                        route.append(0)
                        
                        #ntour_edge_lens = _normalize_tour_edges(tour, D)
                        # we assume that the D is already noremalized, no need to normalize again
                        tour_edge_lens = _get_tour_edges(route, D)
                        
                        _analyze_tour_edges( tour_edge_lens, #in
                           quantile_edges, tour_segment_edge_counts, #out
                           tour_segment_lengths, tour_segment_edge_lengths)
    
                    # start a new tour
                    if pt_idx<0:
                        route_demand = 0.0
                        nodes_per_route = 0
                        solution[from_pt_idx, 0] = 1
                        solution[0, from_pt_idx] = 1
                        from_pt_idx = 0
                        route = [0]
                        
                if pt_idx!=0:                
                    a_pt_idx = abs(pt_idx)
                    if qs:
                        route_demand+=qs[a_pt_idx]
                    nodes_per_route+=1
                    gt_solution.append(a_pt_idx)
                    route.append(a_pt_idx)
                    solution[from_pt_idx, a_pt_idx] = 1
                    solution[a_pt_idx, from_pt_idx] = 1
                    from_pt_idx = a_pt_idx
                    
            intersections.append( _analyze_tour_for_intersections(gt_solution, pts) )
            LO_solutions.append(solution.flatten())
        
    os.remove(SOLUTION_FILE)
    
    # Postprocess feature values
    if normalize:
        gap_initial = _normalize_solution_values(initial_solution_values, best_known_solution_value, best_known_solution_value)
        gap_after_ls = _normalize_solution_values(ls_LO_solution_values, best_known_solution_value, best_known_solution_value)
        gap_improvements = _normalize_solution_values(ls_improvements, best_known_solution_value)
    else:
        gap_initial = initial_solution_values
        gap_after_ls = ls_LO_solution_values
        gap_improvements = ls_improvements
    
    # No improvements made. Mark the feature as 0
    if len(gap_improvements)==0:
        gap_improvements=[0.0]
    
    features.append( ("LSP1", "SotD of Solution Values after Heuristic Construction [RKM]", calculate_distribution_statistics(gap_initial)) )
    features.append( ("LSP2", "SotD of Local Optimum after Intra-route Heuristic Improvement [Hu14]", calculate_distribution_statistics(gap_after_ls)) ) 
    features.append( ("LSP3", "SotD of Improvement per Local Search Step [Hu14]", calculate_distribution_statistics(gap_improvements)) )
    features.append( ("LSP4", "SotD of Number of Local Search Steps to Local Optimum [Hu14]", calculate_distribution_statistics(ls_steps_to_LO)) )
    
    
    sol_to_sol_normalized_distance = []
    probability_of_edge_in_solution = np.zeros(N*N)
    for i in range(repeats):
        probability_of_edge_in_solution+=LO_solutions[i]
        for j in range(i+1,repeats):
            manhattan_distance = spatial.distance.cityblock(LO_solutions[i], LO_solutions[j])
            # calculate similarity (the diagonal is always equal, do not count those)
            sol_to_sol_normalized_distance.append( manhattan_distance/(float(N**2-N ) ) )
    probability_of_edge_in_solution/=float(repeats)
    features.append( ("LSP5", "SotD of Normalized Manhattan Distance Between Local Optima [Hu14]", calculate_distribution_statistics(sol_to_sol_normalized_distance)) )
    features.append( ("LSP6", "SotD of Probability of Edges in Local Optima [Hu14]", calculate_distribution_statistics(probability_of_edge_in_solution)) )
    
    #TODO: LSP Number of improving steps / LSP Number of best improving steps? -> would require best_accept, first_accept
    
    for quartile_idx, quartile_edge_lens in enumerate(quantile_edges): 
        #print "###", quartile_idx, len(quartile_edge_lens), sorted(quartile_edge_lens), calculate_distribution_statistics(quartile_edge_lens)
        features.append( ("LSP7-%d"%(quartile_idx+1), "SotD of Intra-route Normalized Edge Lenghts in %d. Quartile [PM14]" % (quartile_idx+1),
                          calculate_distribution_statistics(quartile_edge_lens)) )
    
        
    features.append( ("LSP8", "SotD of Intra-route Normalized Segment Lenghts [PM14]",
                          calculate_distribution_statistics(tour_segment_lengths)) )
    features.append( ("LSP9", "SotD of Intra-route Segment Edge Count[PM14]",
                          calculate_distribution_statistics(tour_segment_edge_counts)) )
    features.append( ("LSP10", "SotD of Intra-route Segment Normalized Edge Lenghts [PM14]",
                          calculate_distribution_statistics(tour_segment_edge_lengths)) )                      

    features.append( ("LSP11", "SotD of Edge Intersections [PM14]",
                          calculate_distribution_statistics(intersections)) )                      

    features.append( ("LSP13", "SotD of Local Optimum Route Fill Ratio [RKM17]", calculate_distribution_statistics(normalized_route_demands)) )
    features.append( ("LSP14", "SotD of Normalized Number of Nodes on Local Optimum Routes [RKM17]", calculate_distribution_statistics(normalized_nodes_per_route)) )
    
    if step_compositions:
        tried_nums, better_nums, best_nums, moves_nums = zip(*step_compositions)
    else:
        tried_nums, better_nums, best_nums, moves_nums = [],[],[],[]
    features.append( ("LSP15", "SotD of Number of Tried Moves per Step [PM14]",
                      calculate_distribution_statistics(tried_nums)) )
    features.append( ("LSP16", "SotD of Number of Improving Moves per Step [PM14]",
                      calculate_distribution_statistics(better_nums)) )
    features.append( ("LSP17", "SotD of Number of Best Moves per Step [PM14]",
                      calculate_distribution_statistics(best_nums)) )
    features.append( ("LSP18", "SotD of Number of Made Moves per Step [RKM17]",
                      calculate_distribution_statistics(moves_nums)) )
                      
    if return_upper_bound:
        return features, upper_bound, upper_bound_k
    else:
        return features
    
def test():
    vrpfile = "../instances/A-n32-k5.vrp"
    from cvrp_rkm16_io import read_TSPLIB_CVRP
    (N, pts, _, demands, D, C) = read_TSPLIB_CVRP(vrpfile)
    #import re
    #k = int(re.match("k(\d+)", vrpfile).group(1))
    from pprint import pprint
    pprint(lsp_features(vrpfile, pts, D, demands, C, repeats=5))    
    
if __name__=="__main__":
    np.seterr(all='raise')
    test()