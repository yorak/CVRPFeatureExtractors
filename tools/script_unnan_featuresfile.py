# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 13:57:07 2016

@author: juherask
"""

write_output = False

# open featurefile
iff = open(r"../results/vrp_all_features_FINAL.csv", 'r')
if write_output:
    off = open(r"./results/vrp_all_features_FINAL_nonnan.csv", 'w')
header = None
for i, l in enumerate(iff.readlines()):
    if not header:
        header = l.split(";")
    if "nan" in l:
        cols = l.split(";")
        for j,c in enumerate(cols):
            if "nan" in c:
                print i, j, header[j]
                # no upper bound found
                if j==13:
                    cols[j]="2.0"
                # no improving LS steps found
                if 98<=j<=102:
                    cols[j]="0.0"
        l = ";".join(cols)
    
    if write_output:                
        off.write(l)
        