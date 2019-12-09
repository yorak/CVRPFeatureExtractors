# -*- coding: utf-8 -*-

import helpers.cvrp_rkm16_ops as cvrp_ops
from random import choice, randint, shuffle, randrange

LS_ONE_POINT_MOVE = 1
LS_RANDOM_POINT_MOVE = 2

def _random_point_move(N, points, demands, C,  D, from_solution):
    new_solution = from_solution
    
    while (from_solution==new_solution):
        
        lsol = len(from_solution)
        
        # Find the point to move (!=0)
        while 1:
            move_from_idx = randint(0, lsol-1)
            point_to_move = from_solution[move_from_idx]
            if point_to_move!=0:
                break
        
        # Prevent removed point to leave 0, 0 behind
        pt_inc = 1
        while from_solution[move_from_idx-1]==0 and\
              move_from_idx+pt_inc<=lsol and\
              (move_from_idx+pt_inc>=lsol or from_solution[move_from_idx+pt_inc]==0):
            pt_inc+=1
                
        new_solution = from_solution[:move_from_idx] + from_solution[move_from_idx+pt_inc:] 
    
        while 1:
            move_to_idx = randint(1, lsol-1)
            if move_to_idx != move_from_idx:
                break;
                
        new_solution.insert(move_to_idx, point_to_move)
    
        if new_solution[-1]!=0:
            new_solution.append(0)
            
    new_solution_q = cvrp_ops.calculate_objective(new_solution, D)
    new_solution_f = all(cvrp_ops.check_solution_feasibility(new_solution, D, demands, C))     
    #print move_from_idx, from_solution, move_to_idx, new_solution
        
    return new_solution, new_solution_q, new_solution_f
        
        
    
def _one_point_move(N, points, demands, C,  D, from_solution):
    
    lsol = len(from_solution)

    best_solution = from_solution
    best_solution_q = cvrp_ops.calculate_objective(N, from_solution, D)
    
    for point_to_move_idx in xrange(1, lsol-1):
        
        pt_to_move = from_solution[point_to_move_idx]
        if pt_to_move==0: continue #Skip moving route delimeters
        
        
        # Prevent removed point to leave 0, 0 behind
        pt_inc = 1
        while from_solution[point_to_move_idx-1]==0 and\
              point_to_move_idx+pt_inc<=lsol and\
              (point_to_move_idx+pt_inc>=lsol or from_solution[point_to_move_idx+pt_inc]==0):
            pt_inc+=1
            
        #print from_solution, point_to_move_idx, pt_inc
    
        # Create new candidates
        for move_after_idx in xrange(0, lsol):
            new_solution = from_solution[:point_to_move_idx] + from_solution[point_to_move_idx+pt_inc:] 
            
            new_solution.insert(move_after_idx, pt_to_move)
            
            if new_solution[-1]!=0:
                new_solution.append(0)
                
            if all(cvrp_ops.check_solution_feasibility(new_solution, D, demands, C)):
               q = cvrp_ops.calculate_objective(N, new_solution, D)
               if q < best_solution_q:
                   best_solution = new_solution
                   best_solution_q = q
    
    return best_solution, best_solution_q, True

def _choose_random_initial_solution(N, points, demands, C, D, solutions, require_feasibility):
    while 1:    
        initial_sol = choice(solutions)
        if not require_feasibility or \
           all(cvrp_ops.check_solution_feasibility(initial_sol, D, demands, C)):
               
           initial_q =  cvrp_ops.calculate_objective(N, initial_sol, D)
           initial_f =  all(cvrp_ops.check_solution_feasibility(initial_sol,D,demands,C)) 
           break
       
    return (initial_sol, initial_q, initial_f)
    
    
def _build_random_initial_solution(N, points, demands, C, D, require_feasibility):
    
    ids = range(1,N+1)
    shuffle(ids)
    
    solution = [0]+ids+[0] 
    new_route_candidates = range(1,N)
    
    num_routes = randint(1,N-1)
    t = 1
    
    # Split to random number of routes, but make sure the solution is feasible if needed
    while t<num_routes or (require_feasibility and \
          not all(cvrp_ops.check_solution_feasibility(solution, D, demands, C))):
              
        t+=1
    
        lcand = len(new_route_candidates)
        new_route_pointer = randrange(0,lcand)
        new_route_idx = new_route_candidates[new_route_pointer]
        
        new_route_candidates.pop(new_route_pointer)
        for i in xrange(new_route_pointer, lcand-1):
            new_route_candidates[i]+=1
    
        solution.insert(new_route_idx+1, 0)        
        
        #print "t, num_routes", t, num_routes 
        #print "feasible:", cvrp.check_solution_feasibility(solution,D,demands,C)
        #print "picked:", new_route_idx
        #print "candidates:", new_route_candidates
        #print "solution:", solution
    
    
    solution_q =  cvrp_ops.calculate_objective(solution, D)
    solution_f =  all(cvrp_ops.check_solution_feasibility(solution, D, demands, C))
       
    return (solution, solution_q, solution_f)


def solve(points, demands, C, D, solutions, steps, local_search = [LS_ONE_POINT_MOVE]):
    """
    Solve a CVRP with points representing the problem customer points with
    demands being the capacity, C the capacity limit, D the distance matrix
    between the points (including the depot with idx 0), steps is the
    number of steps to do local search (unless the local optimum is found)
    and local_search being the local search algorithm to use.
    local_search can be a list of operators for each step. Alternatives are
        LS_ONE_POINT_MOVE
        RANDOM_POINT_MOVE        
    """    
    
    # One is depot
    N = len(points)-1
    
    # solutions visited (with qualities)
    trace = []
    
    # Pick initial solution (must be feasible) 
    require_feasible_init = True
    if len(local_search)==1 and local_search[0] == LS_RANDOM_POINT_MOVE:
        require_feasible_init = False
    
    if solutions:
        current_sol = _choose_random_initial_solution(N, points, demands, C, D, solutions, require_feasible_init)
    else:
        current_sol = _build_random_initial_solution(N, points, demands, C, D, require_feasible_init)
        
    trace.append( current_sol )
    
    at_lo = False
    for step in xrange(steps):
        for ls_op in local_search:
            if ls_op == LS_ONE_POINT_MOVE:
                next_sol = _one_point_move(N, points, demands, C, D, current_sol[0] )
            elif ls_op == LS_RANDOM_POINT_MOVE:
                next_sol = _random_point_move(N, points, demands, C, D, current_sol[0] )
            
            #print next_sol, current_sol
            if (next_sol[1]==current_sol[1] and next_sol[0]==current_sol[0]):
                #print ("Local optimum reached")
                at_lo = True
            
            trace.append(next_sol)
            current_sol = next_sol
        if at_lo:
            break
    return trace
    
def smoke_test():
    from cvrp_io import generate_CVRP
    (N, pts,demands, D, C) = generate_CVRP(50, 100.0, 5.0, 2.0)
    print solve(pts, demands, C, D, None, 100)
    #from pprint import pprint
    #pprint( zip(["solution", "solution_q", "solution_f"], list(solve(pts, demands, C, D, None, 100))) )
    
if __name__=="__main__":
    smoke_test()
    