# -*- coding: utf-8 -*-

from analyze_process import read_features_from_file_and_validate

from scipy import stats
import matplotlib.pyplot as plt
import numpy as np

from math import sqrt

def plot_hist(data, titles):
    
    #data=np.random.random((4,10))
    #xaxes = ['x1','x2','x3','x4']
    #yaxes = ['y1','y2','y3','y4']
    #titles = ['t1','t2','t3','t4'] 
    
    Nf = len(titles)
    dim = int(sqrt(Nf)+0.5)
    
    f,a = plt.subplots(dim,dim)
    a = a.ravel()
    for idx,ax in enumerate(a):
        print idx, titles[idx]
        ax.hist(data[idx])
        ax.set_title(titles[idx])
        #ax.set_xlabel(xaxes[idx])
        #ax.set_ylabel(yaxes[idx])
    plt.tight_layout()
    plt.savefig('full_feature_histog.png')

def test():
    features_file = '../vrplib_features_A_070317.csv'    
    header, cols, float_headers, float_cols = read_features_from_file_and_validate(features_file)
    
    for fcol_idx, fcol_values in enumerate(float_cols):
        print float_headers[fcol_idx], stats.describe(fcol_values)
        if np.nan in fcol_values:
            print "ALERT: Nan"
        
    for
    plot_hist(of float_cols, float_headers)
    
if __name__=='__main__':
    test()