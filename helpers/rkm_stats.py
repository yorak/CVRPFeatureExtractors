# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 10:49:11 2017

@author: juherask
"""

from scipy.stats import describe
from scipy.stats.mstats import mquantiles
import numpy as np
from collections import Counter

STAT_LEVEL_BASIC = 0
STAT_LEVEL_EXTRA = 1
STAT_LEVEL_QUANTILES = 2
STAT_LEVEL_MODE = 3

def get_distribution_statistics_labels():
    return ["mean", "CV", "stdev", "skew", "kurt",
    "min", "max", "median",
    "25p_qtr", "50p_qtr", "75p_qtr",
    "#modes", "#mode_values", "mode(avg)"]
    
def calculate_distribution_statistics(l, level=STAT_LEVEL_BASIC):
    """
    level = 0, STAT_LEVEL_BASIC, "mean", "CV", "stdev", "skew", "kurthosis",
    level = 1, STAT_LEVEL_EXTRA, also "min", "max", "median",
    level = 2, STAT_LEVEL_QUANTILES, also "25/50/72p_quantile"
    level = 3, STAT_LEVEL_MODE, also "#modes", "#mode_values", "mode(avg)",
    """
    
    counts = Counter(l)
    
    # handle special cases here to avoid RuntimeWarnings from the stats and numpy
    if len(l)==0:
        statcnt = 5
        if level>=STAT_LEVEL_EXTRA:
            # min max median
            statcnt+=3        
        if level>=STAT_LEVEL_QUANTILES:
            # quantiles
            statcnt+=3
        if level>=STAT_LEVEL_MODE:
            # modes
            statcnt+=3
        return [np.nan]*statcnt
    # there is only one measurement, or the measurements are constant
    elif len(l)==1 or len(counts)==1:
        oneval_stats = [l[0], 0.0, 0.0, 0.0, float('Inf')]
        
        if level>=STAT_LEVEL_EXTRA:
            # min max median
            oneval_stats +=[l[0],l[0],l[0]]        
        if level>=STAT_LEVEL_QUANTILES:
            # quantiles
            oneval_stats +=[l[0],l[0],l[0]]        
        if level>=STAT_LEVEL_MODE:
            # modes
            oneval_stats +=[1,1,l[0]]        
            
        return oneval_stats
    # Ok, there seems to be some measurements to describe, good
    else:  
        statmom = list(describe(l))    
        
        stdev = np.std(l)
        mean = statmom[2]
        if mean==0 and stdev==0:
            coefficient_of_variation = 0
        else:
            coefficient_of_variation = stdev/mean
                
        #leave out 0=N, 1=(min,max)
        statistics = [statmom[2]]+[coefficient_of_variation, stdev]+statmom[4:]
        
        if level>=STAT_LEVEL_EXTRA:    
            # min, max
            statistics+=list(statmom[1])
            # median
            statistics+=[ np.median(l) ]
        if level>=STAT_LEVEL_QUANTILES:    
            statistics+=list(mquantiles(l))
        if level>=STAT_LEVEL_MODE:
            # mode    
            mode = counts.most_common(1)[0]
            modes = [m for m in counts.most_common() if m[1]==mode[1]]
            modevals,_ = zip(*modes)
            # qmodes, fmodes, mode_avg
            statistics+=[ len(modes), mode[1], np.mean(modevals) ] 
            
        return statistics