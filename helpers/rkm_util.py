# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 10:35:54 2011

@author: juherask
"""

from math import radians, cos, sin, asin, sqrt
from operator import itemgetter
from itertools import groupby
X = 0; Y = 1 # ~DEFINE for readablity
    
def normalize_to_rect(pts, to_range, keep_aspect_ratio=False):
    """ pts elements are (x,y,d), where d is demand
    to_range is a 2-tuple for the desired output range min and max"""
    
    xl, yl  = zip(*pts)
    minx = min(xl)
    maxx = max(xl)
    rangex = maxx-minx
    miny = min(yl)
    maxy = max(yl)
    rangey = maxy-miny
    minr = min(to_range)
    maxr = max(to_range)
    ranger = maxr-minr
    
    if keep_aspect_ratio:
        rangex = max(rangex, rangey)
        rangey = rangex
    
    new_xl = [(x-minx)/float(rangex)*(ranger)+minr for x in xl]
    new_yl = [(y-miny)/float(rangey)*(ranger)+minr for y in yl]
    
    return zip(new_xl, new_yl)    

def d(pt1,pt2):
    return sqrt( (pt1[X]-pt2[X])**2 + (pt1[Y]-pt2[Y])**2 )


def produce_nn_list(D):
    n = len(D)
    # preprocess D to sorted_per_line_D
    NN_D = [None]*n
    for i in range(n):
        # sort each row
        NN_D[i] = sorted(enumerate(D[i][:]), key=itemgetter(1))
    return NN_D
    

def haversine(pt1, pt2):
    """
    # from http://stackoverflow.com/questions/4913349/
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    The distance should be within ~0.3% of the correct value.
    """    
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [pt1[0], pt1[1], pt2[0], pt2[1]])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    km = 6367 * c
    return km

def list_trim(l, e):
    """ works like string trimming but for lists """
    trimmed = list(l)
    while len(trimmed)>0 and trimmed[0]==e:
        del trimmed[0]
    while len(trimmed)>0 and trimmed[-1]==e:
        del trimmed[-1]
    return trimmed
    
def list_split(l, e):
    """ works like string splitting but for lists """
    return [list(group) for k, group in groupby(l, lambda x:x==e) if not k]    