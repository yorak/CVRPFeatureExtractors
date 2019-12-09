# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 13:19:58 2016

@author: juherask

This script finds the problem instance folders for similar problem instances /
possible duplicates.
"""

from helpers.cvrp_rkm16_io import read_TSPLIB_CVRP
from os.path import basename
from glob import glob
import numpy as np


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
match_from_folders = folders
match_from_folders = [("Examples", r"../instances/")]

first = True

for vrp_folder_t in match_from_folders:
    print 
    print vrp_folder_t[0]
    for vrp_file_t in glob(vrp_folder_t[1]+r"\*.vrp"):
        
        (tN, tpoints, _, tdemands, tD, tC) = read_TSPLIB_CVRP(vrp_file_t)
        
        tD = np.array(tD)
        n_tD = tD / np.max(tD)
        tdemands.sort()
        
        
        match = None
        match_lvl = 0
        for pset, pfolder in folders:        
            if pfolder==vrp_folder_t[1]:
                continue
        
            for vrp_file_c in glob(pfolder+r"\*.vrp"):
                #print "processing", vrp_file
                (cN, cpoints, _, cdemands, cD, cC) = read_TSPLIB_CVRP(vrp_file_c)
                
                cD = np.array(cD)
                n_cD = cD / np.max(cD)
                cdemands.sort()
                
                if cN==tN:
                    #print basename(vrp_file_t), "and", basename(vrp_file_c), "same N"
                    if not match:
                        match = vrp_file_c
                        match_lvl = 1
                    
                    distance_of_DMs = np.sum(n_tD-n_cD)/(tN*cN)
                    #print distance_of_DMs
                    if distance_of_DMs<0.01:
                        if match_lvl<=1:
                            match = vrp_file_c
                            match_lvl = 2
                        elif match_lvl == 2:
                            print "Warning: already lvl 2 match with", basename(match)
                            match = vrp_file_c 
                            
                        if tdemands==cdemands:
                            #print basename(vrp_file_t), "and", basename(vrp_file_c), "same demands"
                            
                            if match_lvl<=2:
                                match = vrp_file_c
                                match_lvl = 3
                            elif match_lvl == 3:
                                print "Warning: already lvl 3 match with", basename(match)
                                match = vrp_file_c

                            if cC==tC:
                                if match_lvl<=3:
                                    match = vrp_file_c
                                    match_lvl = 4
                                elif match_lvl == 4:
                                    print "Warning: already lvl 4 match with", basename(match)
                                    match = vrp_file_c
                            
                            
                                
        if match:
            print "MATCHED level", match_lvl, "for", basename(vrp_file_t), "with", basename(match)
        else:
            print "NO MATCH for", basename(vrp_file_t)