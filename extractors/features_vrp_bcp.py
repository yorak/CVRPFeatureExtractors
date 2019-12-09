# -*- coding: utf-8 -*-
"""
Created on Wed Mar 09 11:22:51 2016

@author: juherask

Local search probing with VRPH
"""

from helpers.rkm_stats import calculate_distribution_statistics, STAT_LEVEL_QUANTILES
from subprocess import Popen, PIPE
from os import path
import numpy as np
import sys


DEBUG_BCP_TO_FILE = False

SOLVER_CWD =  path.abspath(path.join(path.dirname(__file__), r'..', 'solvers'))
branch_and_cut_solver = path.join(SOLVER_CWD,r'vrp')
if sys.platform == 'win32':
    branch_and_cut_solver+='.exe'
#./vrp.exe 

class InfeasibleException(Exception):
    pass

def _normalize_solution_values(vals, bks, floorsol=0):
    if not hasattr(vals, '__iter__'):
        return _normalize_solution_values([vals], bks, floorsol)[0]
    return [(val-floorsol)/float(bks) for val in vals]

def invoke_symphony(problem_file, k, time_limit, verbosity=5,
                    upper_bound=None, solve_as_TSP=False):
    #vrp.exe -F sample.vrp -N 5 -t 10 -v 3 > solve_sample_v3.log
    cmd_w_args = [branch_and_cut_solver,
          '-F', problem_file,
          '-p', '1',
          '-N', str(k), 
          '-t', str(time_limit),
          '-v', str(verbosity)
          #'-o', 'bcp_probe.sol' 
          #TODO: Find out what could be read from the b/c - tree (sol) file 
    ]
    
    ub_at_end = None
    if upper_bound:
        ub_at_end = upper_bound
        cmd_w_args.append('-u')
        cmd_w_args.append(str(upper_bound+0.1))
    if solve_as_TSP:
        cmd_w_args.append('-R')
        
    #print "REMOVEME", " ".join(cmd_w_args)

    lb_at_end = None
    lb_at_root = None
        
    lower_bound_improvements_per_cut = []
    improvements_between_feasible_solutions = []
    
    prev_rlp_val = None
    cuts_added = None 
    feasible_solution_value = None
    optimal_found  = False
    result_infeasible = False
    non_integer_lp_variable_values = []
    relaxed_lp_val = None
    total_wallclock_time = 0.0
    
    p = Popen(cmd_w_args, stdout=PIPE, bufsize=1, cwd=SOLVER_CWD)
    
    #print " ".join(cmd_w_args)
    if DEBUG_BCP_TO_FILE :
        print " ".join(cmd_w_args)
        debug_output = file(path.join(SOLVER_CWD,"debug_symphony_out.txt"), "w")
    #for line in vrph_output.splitlines():
    read_mode = None
    with p.stdout:
        for line in iter(p.stdout.readline, b''):
            #print line,
            if DEBUG_BCP_TO_FILE :
                debug_output.write(line)
                
            
            # check if we need to change the mode
            if line and line[0]=="*":
                if r"* Now processing NODE" in line:
                    prev_rlp_val = None
                    cuts_added = None
                    read_mode = 'node_processing'                
                if r"* Time Limit Reached " in line:
                    read_mode = 'stats'
                if r"* Problem Infeasible " in line:
                    #note: this also mean that the time ran out!
                    # check it later
                    result_infeasible = True
                    read_mode = 'stats'
                if r"* Optimal Solution Found " in line:
                    optimal_found = True
                    non_integer_lp_variable_values = []
                    read_mode = 'stats'
                    
            # process lines according to the read mode

            if read_mode == 'variable_values':
                if line.strip()=="":
                    read_mode='node_processing'
                elif not "+" in line:
                    str_idx, str_value = line.split()
                    idx, value = int(str_idx), float(str_value)
                    # only take nonintegers 
                    if value!=1.0:
                        non_integer_lp_variable_values.append( value )
            
            elif read_mode == 'node_processing':
                if r"The LP value is:" in line:
                    relaxed_lp_val = float(line.split()[4])
                    if cuts_added and prev_rlp_val:
                        lower_bound_improvements_per_cut.append( (relaxed_lp_val-prev_rlp_val)/cuts_added )
                    prev_rlp_val = relaxed_lp_val
                elif r"cuts added altogether" in line:
                    #Continue with this node. 29 cuts added altogether in iteration 11
                    cuts_added = int(line.split()[4])
                elif r"Solution Cost:" in line:
                    new_feasible_solution_value = float(line.split(":")[1])
                    if feasible_solution_value:
                        d_fsv = feasible_solution_value-new_feasible_solution_value
                        improvements_between_feasible_solutions.append(d_fsv)
                    feasible_solution_value = new_feasible_solution_value
                    non_integer_lp_variable_values = []
                elif r" User indices and values of nonzeros in the solution" in line:
                    read_mode='variable_values'
                    non_integer_lp_variable_values = []
                    
            elif read_mode == 'stats':
                if r"Lower Bound in Root: " in line:
                    lb_at_root = float(line.split(":")[1])
                elif r"Current Lower Bound: " in line:
                    lb_at_end = float(line.split(":")[1])
                elif r"Current Upper Bound: " in line:
                    ub_at_end = float(line.split(":")[1])
                elif r"Total Wallclock Time" in line:
                    total_wallclock_time = float(line.split()[-1])
        
    # make sure subprocess finishes
    p.wait()
    
    if DEBUG_BCP_TO_FILE :
        debug_output.close()
    
    # if for some reason this is not set (eg. incorrect
    #  infeasibility) use the latest relaxed lp value.
    if not lb_at_end and relaxed_lp_val:
        lb_at_end = relaxed_lp_val   
    if not lb_at_end and ub_at_end:
        ub_at_end  = ub_at_end
    if optimal_found:
        lb_at_end = feasible_solution_value
    elif result_infeasible:
        # it says it is infeasible, but it used all of the time, do not 
        #  believe it!
        if abs(total_wallclock_time-time_limit)<0.1:
            result_infeasible = False # -> TIMEOUT
            ub_at_end = upper_bound
        # it seems the optimization was terminated early, return nothing
        else:
            lb_at_end, ub_at_end, lb_at_root, feasible_solution_value = None, None, None, None
    
    #print "REMOVEME:", p.returncode, upper_bound, lb_at_end, feasible_solution_value
    
    #print "lb_at_end, ub_at_end, lb_at_root, feasible_solution_value, result_infeasible, optimal_found"
    #print lb_at_end, ub_at_end, lb_at_root, feasible_solution_value, result_infeasible, optimal_found
    return lb_at_end, ub_at_end, lb_at_root, feasible_solution_value, \
      result_infeasible, optimal_found,\
      lower_bound_improvements_per_cut, improvements_between_feasible_solutions,\
      non_integer_lp_variable_values
   
def bcp_features(problem_file_name, N, k, time_limit=3.0, upper_bound=None):
    
    # Get BKS from the .vrp file    
    best_known_solution_value = None            
    is_exact_2d = False
    vrpfile = open(problem_file_name, "r")
    has_demands = False
    for fline in vrpfile.readlines():
        if "BEST_KNOWN"  in fline:
            best_known_solution_value = float(fline.split(":")[1])
        if "EXACT_2D" in fline:
            is_exact_2d = True
        if "DEMAND_SECTION" in fline:
            has_demands = True
            
    vrpfile.close()
    
    if upper_bound!=None and best_known_solution_value!=None and upper_bound<best_known_solution_value:
        print "WARNING: Upper bound is smaller than BKS! This is impossible. Ignoring UB."
        print upper_bound, "vs", best_known_solution_value
        upper_bound = None        
        #raise UserWarning("Upper bound is smaller than BKS! Ignoring UB.")
    
    # QUICKHACK: The support for EXACT_2D coordinates in SYPHONY is a hack, 
    #  adjust the upper bound and BKS to the same 5 decimal accuracy.
    if is_exact_2d:
        if upper_bound:
            upper_bound *= 10000;
            upper_bound = int(upper_bound)
        if best_known_solution_value:
            best_known_solution_value *= 10000
            best_known_solution_value = int(best_known_solution_value)

    features = []
    retrys_total = 5
    retry_count = retrys_total
    
    #invoke_symphony(problem_file_name, k, time_limit, upper_bound=upper_bound)
    
    while retry_count>0:
        try:
            lb_at_end, ub_at_end, lb_at_root, feasible_solution_value, \
                result_infeasible, is_optimal,\
                lower_bound_improvements_per_cut,\
                improvements_between_feasible_solutions,\
                non_integer_lp_variable_values = invoke_symphony(
                  problem_file_name, k, time_limit,
                  upper_bound=upper_bound,
                  solve_as_TSP=not has_demands)
            
            #print "NILPV", non_integer_lp_variable_values
            
            # HACK: SYMPHONY sometimes (probably) incorrectly generates cuts that
            #  make the subproblem infeasible resulting to incorrect infeasibility
            #  wiggling the time limit a bit sometimes fixes this.
            if result_infeasible:
                raise InfeasibleException("INFEASIBLE! retrying (%d of %d retries)" % (retrys_total-retry_count+1, retrys_total))
                                  
            if len(lower_bound_improvements_per_cut)==0:
                lower_bound_improvements_per_cut.append(0.0)
                #raise ValueError("Something is wrong, did not get improvements per cuts")
            if not lb_at_end:
                lb_at_end = lb_at_root
                
            # the UB may be None
            ub_lb_rate = float('Inf')
            if lb_at_end and ub_at_end:
                ub_lb_rate = ub_at_end/lb_at_end
                
            
            n_lb_at_end = _normalize_solution_values(lb_at_end, best_known_solution_value)
            n_lb_improvements = _normalize_solution_values(lower_bound_improvements_per_cut, best_known_solution_value)
            #n_feasible_improvements = _normalize_solution_values(improvements_between_feasible_solutions, best_known_solution_value)
                
            best_solution = lb_at_end
            if feasible_solution_value:
                best_solution = feasible_solution_value
            
            # negative -> infeasble
            n_best_solution = _normalize_solution_values(best_solution, best_known_solution_value, best_known_solution_value)
            
            features.append( ("BCP1", "SotD Normalized Branch-and-Cut Improvement per Cut [Hu14]", calculate_distribution_statistics(n_lb_improvements)) )   
            features.append( ("BCP2", "Ratio of Branch-and-Cut Upper and Lower Bound at the End [Hu14]", ub_lb_rate) ) 
            features.append( ("BCP3", "Best Found Branch-and-Cut Solution (negative is relaxed solution) [Hu14]", n_best_solution) ) 
            features.append( ("BCP4", "Normalized Branch-and-Cut Lower Bound at the End [Hu14]", n_lb_at_end) )                         
            
            N_w_depot = N+1
            decision_variable_count = (N_w_depot*N_w_depot-N_w_depot)/2
            nonint_to_int_dcvs = len(non_integer_lp_variable_values)/float(decision_variable_count)
            features.append( ("BCP5", "Ratio of LP relaxation non-integer to integer decision variable values [Hu14]", nonint_to_int_dcvs) )                         
            
            # avoid nans
            if len(non_integer_lp_variable_values)==0:
                non_integer_lp_variable_values = [0.0]
            dcsv_stats = calculate_distribution_statistics(non_integer_lp_variable_values, STAT_LEVEL_QUANTILES)
            features.append( ("BCP6", "SotD of LP relaxation non-integer decision variable values [Hu14]", dcsv_stats ) ) 
            break
        except Exception as e:
            #raise e
            if type(e) is InfeasibleException:
                print e.message
            else:
                print "SOLVER CRASH! retrying (%d of %d retries)" % (retrys_total-retry_count+1, retrys_total)
            # "wiggle" it a little: try to give more time and more slack to the upper bound
            time_limit+=0.01
            if upper_bound:
                upper_bound=upper_bound*1.01
            retry_count-=1
            # for the last retry, ignore the UB altoghter and increase K
            if retry_count==1:
                k += 1
                upper_bound = None
        
    if retry_count==0:
        features.append( ("BCP1", "SotD Normalized Branch-and-Cut Improvement per Cut [Hu14]", [np.nan]*5 ) )   
        features.append( ("BCP2", "Ratio of Branch-and-Cut Upper and Lower Bound at the End [Hu14]", np.nan) ) 
        features.append( ("BCP3", "Best Found Branch-and-Cut Solution (negative is relaxed solution) [Hu14]", np.nan) ) 
        features.append( ("BCP4", "Normalized Branch-and-Cut Lower Bound at the End [Hu14]", np.nan) )         
        features.append( ("BCP5", "Ratio of LP relaxation non-integer to integer decision variable values [Hu14]", np.nan) )                         
        features.append( ("BCP6", "SotD of LP relaxation non-integer decision variable values [Hu14]", [np.nan]*11 ) ) 
            
    return features
    
def test():
    global DEBUG_BCP_TO_FILE 
    DEBUG_BCP_TO_FILE = True
    #print bcp_features("../instances/test_instance.vrp", 3.0)
    print bcp_features(r"../instances/A-n33-k6.vrp",
                       32, 6, time_limit=3.0, upper_bound=774.0)
    print 
    
    print bcp_features(r"../instances/X-n1001-k43.vrp",
                       1000, 43, time_limit=3.0, upper_bound=76766.0)
    print 
    
    print bcp_features(r"../instances/F-n135-k7.vrp",
                       134, 7, time_limit=3.0, upper_bound=1191.0)
    print 
    
    
    print bcp_features(r"../instances/CMT06.vrp",
                       51, 6, time_limit=3.0, upper_bound=560.09)
    print 
    
    
    
if __name__=="__main__":
    test()