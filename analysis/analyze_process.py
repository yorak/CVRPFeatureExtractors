# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 21:29:48 2016

@author: juherask
"""
import numpy as np
from scipy import stats
from sklearn.preprocessing import MinMaxScaler

from math import isnan
import csv

def read_features_from_file_and_validate(features_file):
    with open(features_file, 'rb') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=';')
        #print spamreader.fieldnames
        header = None
        data = []
        for row in spamreader:
            if not header:
                header = row
            else:
                data.append(tuple(row))
                
        cols = zip(*data)
        
        float_headers = []
        float_cols = []
        for col_id, col in enumerate(cols):
            if col_id<2:
                continue
            
            float_headers.append(header[col_id])
            col_values = np.array([float(sv) for sv in col])
            float_cols.append(col_values)
            #print header[col_id], stats.describe(col_values)
            for row_id, value in enumerate(col_values):            
                if isnan(value):
                    print "WARNING: nan in", (row_id, col_id)
        return header, cols, float_headers, float_cols

def preprocess_feature_data(float_cols):
    M_f = np.array(float_cols)
    M_f = M_f.transpose()

    # Fix missing values
    for i in range(len(M_f)):
        for j in range(len(M_f[i])):
            if isnan(M_f[i][j]):
                #print i, j
                if j==11: # BCP lb/ub features
                    M_f[i][j] = 2.0
                    #print "replace", j
                elif j in range(96,101): # LSP features
                    M_f[i][j] = 0.0
                    
    np.savetxt("Mf.txt", M_f)
    # scale to 0.0-1.0
    mms = MinMaxScaler()
    return mms.fit_transform(M_f)
    
#features_file_A = "vrp_features_Augerat-A-14.csv"
#features_file_B = "vrp_features_Augerat-B-14.csv"
#features_file = "vrp_features_Augerat-AandB-14.csv"

def test():
    features_file = 'vrp_all_features_nn_updated.csv'    
    header, cols, float_headers, float_cols = read_features_from_file_and_validate(features_file)
    
    for fcol_idx, fcol_values in enumerate(float_cols):
        print float_headers[fcol_idx], stats.describe(fcol_values)
    
if __name__=='__main__':
    test()