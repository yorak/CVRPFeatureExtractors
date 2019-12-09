# -*- coding: utf-8 -*-
"""
Created on Wed Mar 09 11:22:51 2016

@author: juherask

Based on:
[Hu14] Hutter, F., Xu, L., Hoos, H.H., Leyton-Brown, K.: Algorithm runtime
prediction: Methods & evaluation. Artificial Intelligence206(2014) 79–111
"""

from helpers.rkm_stats import calculate_distribution_statistics
from cvrp_solver import LS_RANDOM_POINT_MOVE, solve
from math import log, sqrt


def _autocorrelation(z, i):
    """Calculate autocorrelation of a trace produced by a search algirithm.
    It is assumed that the fitness landscape is statistically isomorphoc
    (same everywhere, regardless the staring point)
    """

    T = len(z) # length of the time series
    
    z_mean = sum(z) / float(T)
    z_var_est =  sum( ((z[t]-z_mean)**2 for t in xrange(0,T)) )  
    
    # Variance 0, but do not let it be
    if z_var_est < 1e-10:
        return 1.0

    # From Hordjik1996 
    r_i = sum( (((z[t]-z_mean)*(z[t+i]-z_mean)) for t in xrange(0,T-i)) ) / z_var_est
    return r_i
    
def _correlation_length(trace, method = "exp"):
    """ Calculates the correlation lenght, that is a estimate for 
    ruggedenss for the solution space based on the random walk.
    
    if method is "exp":
    The method from Stadler1996 via Richter2006 is used
    lower the lambda more rugged the landscape is. 
    
    if method is "sig":
    The method from Hordijk1996 (that follows the Granger and Newbold1986)
    required that T>>i
    """
    
    solutions, z, f = zip(*trace)
    if method=="exp":
        r_1 = _autocorrelation(z, 1)
        cl_lambda = -1.0 / log(abs(r_1))        
        return cl_lambda
    elif method=="sig":
        T = len(z)
        significance_lvl = 2/sqrt(T) 
        cl_i = 1
        while cl_i<T:
            r_i = _autocorrelation(z, cl_i)
            if abs(r_i)<significance_lvl:
                # The last time lag cl_i was bigger than significance_lvl
                return cl_i-1
            cl_i+=1

def autocorrelation_features(npts, nD, demands, C, walk_count=10, walk_length=None, cl_method = "exp"):
    N = len(npts)
    if not walk_length:
        walk_length = 2*N
    
    ls_ops = [LS_RANDOM_POINT_MOVE]
    
    cl_rw_probs = []
    for r in xrange(walk_count): 
        #solve(points, demands, C, D, solutions, steps, local_search = [LS_ONE_POINT_MOVE]):
        trace = solve(npts, demands, C, nD, None, walk_length, ls_ops)
        lo_sol, lo_q, lo_f = trace[-1]
        cl_rw_probs.append( _correlation_length(trace, cl_method) )
    
    #if print_trace_out:
    #    print "LS to (local) optimum"
    #    for sol in trace:
    #        print sol
    #    print "Optimum is", lo_q, "at", lo_sol
    #    print 
    #    
    #    print "Average correlation length is ", cl_rw_prob
    #    print
    
    
    features= []
    features.append( ("LSP12", "SotD of Autocorrelation Lenghts [Me13]", calculate_distribution_statistics(cl_rw_probs)) )
    
    return features
    
def smoke_test():
    from cvrp_rkm16_io import generate_CVRP
    (N, pts,demands, D, C) = generate_CVRP(50, 100.0, 5.0, 2.0)
    print autocorrelation_features(pts,D, demands, C)
    
    
    
if __name__=="__main__":
    smoke_test()