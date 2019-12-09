# -*- coding: utf-8 -*-
"""
Created on Sat Apr 09 00:58:00 2016

@author: juherask
"""

from helpers.cvrp_rkm16_io import read_TSPLIB_CVRP
from sklearn import manifold
from glob import glob
from os import path

folders = [
    ("VRPLIB-A", r"../../Benchmarks/VRPLIB/A"),
    ("VRPLIB-B", r"../../Benchmarks/VRPLIB/B"),
    ("VRPLIB-E", r"../../Benchmarks/VRPLIB/E"),
    ("VRPLIB-F", r"../../Benchmarks/VRPLIB/F"),
    ("VRPLIB-G", r"../../Benchmarks/VRPLIB/G"),
    ("VRPLIB-M", r"../../Benchmarks/VRPLIB/M"),
    ("VRPLIB-P", r"../../Benchmarks/VRPLIB/P"),
    ("VRPLIB-V", r"../../Benchmarks/VRPLIB/V"),

    ("CMT", r"../../Benchmarks/Christofides"),
    ("Golden", r"../../Benchmarks/Golden"),
    ("Li", r"../../Benchmarks/Li"),
    ("Taillard", r"../../Benchmarks/Taillard"),
    ]

#folders = [folders[1]]

for pset, pfolder in folders:
    print "generating plots for folder", pfolder
    for vrp_file in glob(pfolder+r"/*.vrp"):
        basename = path.basename(vrp_file)
        #print basename
        
        (N, points, dd_points, demands, D, C) = read_TSPLIB_CVRP(vrp_file)
        if points == None and dd_points==None:
            print
            print "WARNING:", basename, "has no coordinates"
            
            mds = manifold.MDS(n_components=2, dissimilarity='precomputed')
            results = mds.fit(D)
            pts = results.embedding_
            print "DISPLAY_DATA_SECTION"
            for i, pts in enumerate(pts):
                print "%i %f %f" % (i+1, pts[0], pts[1])
            print
            print