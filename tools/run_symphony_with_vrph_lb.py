# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 15:41:01 2016

@author: juherask
"""

from extractors.features_vrp_lsp import lsp_features
from extractors.features_vrp_bcp import invoke_symphony
from helpers.cvrp_rkm16_io import read_TSPLIB_CVRP
from math import ceil
from main import regexp_k_from_filename
from os.path import basename
from glob import glob

folders = [
("VRPLIB-A", r"../../Benchmarks/VRPLIB/A"),
("VRPLIB-B", r"../../Benchmarks/VRPLIB/B"),
("VRPLIB-E", r"../../Benchmarks/VRPLIB/E"),
("VRPLIB-F", r"../../Benchmarks/VRPLIB/F"),
("VRPLIB-G", r"../../Benchmarks/VRPLIB/G"),
("VRPLIB-M", r"../../Benchmarks/VRPLIB/M"),
("VRPLIB-P", r"../../Benchmarks/VRPLIB/P"),
("VRPLIB-V", r"../../Benchmarks/VRPLIB/V"),

("CMT",      r"../../Benchmarks/Christofides"),
("Golden",   r"../../Benchmarks/Golden"),
("Li",       r"../../Benchmarks/Li"),
("Taillard", r"../../Benchmarks/Taillard"),
]

whitelist =  None # ["01-eil7.vrp"] #None#"A-n45-k7.vrp"]
blacklist = []

outfile = open("symphony_vrp_results.csv", "w")

first = True
for pset, pfolder in folders:
    for vrp_file in glob(pfolder+r"\*.vrp"):
        file_basename = vrp_file.split('\\')[-1]
        if whitelist is not None:
            if file_basename not in whitelist:
                continue
        elif file_basename in blacklist:
            continue
        
        
        print "processing", vrp_file
        
        # Get the necessary background information
        (N, pts, dd, demands, D, C) = read_TSPLIB_CVRP(vrp_file)
        m = regexp_k_from_filename.search( basename(vrp_file) )
        k = None
        if m:
            k = int(m.group(1))            
        else:
            print "WARNING: using k estimate for", vrp_file
            k = ceil(sum(demands)*1.05/C) # 10 % slack
         
        ub_with_heurs = None
        feasible_solution_value = None
        elapsed_t = 0.0
        ub_t = 0.0
        try:
            import time
                     
            start_t = time.time()
            lspf = dict(lsp_features(vrp_file, pts, D, demands, C, repeats=1, normalize=False))
            ub_t = time.time()-start_t
            #from pprint import pprint
            #pprint(lspf)
            ub_with_heurs = lspf['SotD of Local Optimum after Intra-route Improvement Phase (LS Probe) [Hu14]'][0]
                        
            lb_at_end, ub_at_end, lb_at_root, feasible_solution_value, is_optimal,\
                lower_bound_improvements_per_cut,\
                improvements_between_feasible_solutions =\
                    invoke_symphony(vrp_file, k, time_limit=60, verbosity=0, upper_bound=ub_with_heurs)
            elapsed_t = time.time() - start_t
            
        except IOError as e:
            print str(e)
            continue
            
            #TODO: RUNTIME! of SYMPHONY! Q tightness etc.
        if first:   
            outfile.write("Problem Set; Problem Name;Q;tightness;UB;OPT;ub_time;tot_time")
            outfile.write("\n")
            first=False
        outfile.write(pset+";"+file_basename.replace(".vrp", "")+";")

        outfile.write( str(C)+";")
        outfile.write( str(sum(demands)/float(k*C))+";")
        
        outfile.write( str(ub_with_heurs)+";")
        outfile.write( str(feasible_solution_value)+";")
        outfile.write( str(ub_t)+";")
        outfile.write( str(elapsed_t))#+";")
        
        outfile.write("\n")
    
outfile.close()
print "DONE!"
