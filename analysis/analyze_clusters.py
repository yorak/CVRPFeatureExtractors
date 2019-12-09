# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 21:29:48 2016

@author: juherask
"""

from analyze_process import *
from helpers.rkm_util import d

from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import MinMaxScaler
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN, k_means
import matplotlib.pyplot as plt

import numpy as np

from os import path
from random import random

        
name_translator_table = {
"Golden_01":"G01",  "Golden_05":"G05",  "Golden_09":"G09",  "Golden_13":"G13",  "Golden_17":"G17",
"Golden_02":"G02",  "Golden_06":"G06",  "Golden_10":"G10",  "Golden_14":"G14",  "Golden_18":"G18",
"Golden_03":"G03",  "Golden_07":"G07",  "Golden_11":"G11",  "Golden_15":"G15",  "Golden_19":"G19",
"Golden_04":"G04",  "Golden_08":"G08",  "Golden_12":"G12",  "Golden_16":"G16",  "Golden_20":"G20",

"Li_21":"L21",  "Li_23":"L23", "Li_25":"L25",  "Li_27":"L27",  "Li_29":"L29",  "Li_31":"L31",
"Li_22":"L22",  "Li_24":"L24", "Li_26":"L26",  "Li_28":"L28",  "Li_30":"L30",  "Li_32":"L32",

"Christofides_01":"C01", "Christofides_02":"C02", "Christofides_03":"C03", "Christofides_04":"C04", "Christofides_05":"C05",
"Christofides_06":"C06", "Christofides_07":"C07", "Christofides_08":"C08", "Christofides_09":"C09", "Christofides_10":"C10",
"Christofides_11":"C11", "Christofides_12":"C12", "Christofides_13":"C13", "Christofides_14":"C14", 

"Taillard_100A":"T01",  "Taillard_150A":"T05",  "Taillard_385":"T09",  "Taillard_75D":"T13",
"Taillard_100B":"T02",  "Taillard_150B":"T06",  "Taillard_75A":"T10",
"Taillard_100C":"T03",  "Taillard_150C":"T07",  "Taillard_75B":"T11",
"Taillard_100D":"T04",  "Taillard_150D":"T08",  "Taillard_75C":"T12",

"A-n32-k5":"A01",  "A-n37-k6":"A07",  "A-n45-k7":"A13",  "A-n60-k9":"A19",   "A-n65-k9":"A25",
"A-n33-k5":"A02",  "A-n38-k5":"A08",  "A-n46-k7":"A14",  "A-n61-k9":"A20",   "A-n69-k9":"A26",
"A-n33-k6":"A03",  "A-n39-k5":"A09",  "A-n48-k7":"A15",  "A-n62-k8":"A21",   "A-n80-k10":"A27",
"A-n34-k5":"A04",  "A-n39-k6":"A10",  "A-n53-k7":"A16",  "A-n63-k10":"A22",
"A-n36-k5":"A05",  "A-n44-k6":"A11",  "A-n54-k7":"A17",  "A-n63-k9":"A23",
"A-n37-k5":"A06",  "A-n45-k6":"A12",  "A-n55-k9":"A18",  "A-n64-k9":"A24",

"B-n31-k5":"B01",  "B-n41-k6":"B06",  "B-n50-k7":"B11",  "B-n57-k7":"B16",   "B-n67-k10":"B21",
"B-n34-k5":"B02",  "B-n43-k6":"B07",  "B-n50-k8":"B12",  "B-n57-k9":"B17",   "B-n68-k9":"B22",
"B-n35-k5":"B03",  "B-n44-k7":"B08",  "B-n51-k7":"B13",  "B-n63-k10":"B18",  "B-n78-k10":"B23",
"B-n38-k6":"B04",  "B-n45-k5":"B09",  "B-n52-k7":"B14",  "B-n64-k9":"B19",
"B-n39-k5":"B05",  "B-n45-k6":"B10",  "B-n56-k7":"B15",  "B-n66-k9":"B20",

"E-n101-k14":"E01",  "E-n22-k4":"E04",  "E-n31-k7":"E07",  "E-n76-k10":"E10",  "E-n76-k8":"E13",
"E-n101-k8":"E02",   "E-n23-k3":"E05",  "E-n33-k4":"E08",  "E-n76-k14":"E11",
"E-n13-k4":"E03",    "E-n30-k3":"E06",  "E-n51-k5":"E09",  "E-n76-k7":"E12",

"F-n135-k7":"F1",  "F-n45-k4":"F2", "F-n72-k4":"F3", "G-n262-k25":"G1",

"M-n101-k10":"M1",  "M-n121-k7":"M2",  "M-n151-k12":"M3",  "M-n200-k16":"M4",  "M-n200-k17":"M5",

"P-n101-k4":"P01",  "P-n22-k2":"P06",  "P-n50-k10":"P11",  "P-n55-k15":"P16",  "P-n65-k10":"P21",
"P-n16-k8":"P02",   "P-n22-k8":"P07",  "P-n50-k7":"P12",   "P-n55-k7":"P17",   "P-n70-k10":"P22",
"P-n19-k2":"P03",   "P-n23-k8":"P08",  "P-n50-k8":"P13",   "P-n55-k8":"P18",   "P-n76-k4":"P23",
"P-n20-k2":"P04",   "P-n40-k5":"P09",  "P-n51-k10":"P14",  "P-n60-k10":"P19",  "P-n76-k5":"P24",
"P-n21-k2":"P05",   "P-n45-k5":"P10",  "P-n55-k10":"P15",  "P-n60-k15":"P20",

"att-n48-k4":"V01",      "fri-n26-k3":"V05",  "gr-n48-k3":"V09",       "ulysses-n22-k4":"V13",
"bayg-n29-k4":"V02",     "gr-n17-k3":"V06",   "hk-n48-k4":"V10",
"bays-n29-k5":"V03",     "gr-n21-k3":"V07",   "swiss-n42-k5":"V11",
"dantzig-n42-k4":"V04",  "gr-n24-k4":"V08",   "ulysses-n16-k3":"V12"}


def _closest_node(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = np.sum((nodes - node)**2, axis=1)
    return np.argmin(dist_2)
    
def plot_clusters_with_images(pts, labels, imagefiles, img_size = (800,600), img_file='out.tiff', bg_colors=False):
    from PIL import Image, ImageDraw, ImageOps, ImageFont
    background = Image.new('RGBA', img_size, (255, 255, 255, 255))
    bg_w, bg_h = img_size
    
    x,y = zip(*pts)
    minx = min(x)
    miny = min(y)
    rangex = max(x)-minx
    rangey = max(y)-miny
    
    imagespace_pts = []
    for i, pt in enumerate(pts):
        posx = (pt[0]-minx)/rangex * bg_w
        posy = (1-(pt[1]-miny)/rangey) * bg_h
        imagespace_pts.append( (posx, posy) )        
    
    if bg_colors:
        pixels = background.load()
        imagespace_pts = np.array(imagespace_pts)
        unique_labels = set(labels)
        colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
        for x in range(img_size[0]):
            print x, "of", img_size[0]
            for y in range(img_size[1]):
                cpi = _closest_node( (x,y), imagespace_pts)
                #print cpi
                pixel_label = labels[cpi]
                
                
                pixel_color = (200,200,200, 0)
                if pixel_label!=-1:
                    pixel_color = (
                        int(colors[pixel_label][0]*255),
                        int(colors[pixel_label][1]*255),
                        int(colors[pixel_label][2]*255),
                        0)
                
                #print pixel_color
                pixels[x,y] = pixel_color
        
        #background.save('outbg.tiff')
    
    font = ImageFont.truetype("ANTQUAB.TTF", 14)
    #font = ImageFont.load_default().font
    for i, pt in enumerate(pts):
        img = Image.open(imagefiles[i], 'r')
        img_w, img_h = img.size
        
        mask=Image.new('L', img.size, color=255)
        draw=ImageDraw.Draw(mask) 
        #draw.ellipse((2*img_w/10,img_h/10,img_w*9.5/10, img_h*9/10), fill=0)
        #imask = ImageOps.invert(mask)
        #img.putalpha(mask)
    
        img_w, img_h = img.size
        posx, posy = imagespace_pts[i]
        offset = (int(posx - (img_w)/2), int(posy - (img_h)/2))
        background.paste(img, offset)#, imask)
        
        fbname = path.basename(imagefiles[i])
        draw=ImageDraw.Draw(background) 
        draw.text((offset[0]+7,offset[1]+7),fbname, 255, font=font)
               
    background.save(img_file)


def plot_clusters(M, labels, outfiles=None, core_samples_mask=None, names=None, pts=None, xy_labels=None):    
    
    if pts==None:
        pca2d = PCA(n_components=2)
        PCA_2f = pca2d.fit_transform(M)
        xy_labels = [
            "PC 1 (%.2f%%)" % (pca2d.explained_variance_ratio_[0]*100), 
            "PC 2 (%.2f%%)" % (pca2d.explained_variance_ratio_[1]*100)
        ]
    else:
        PCA_2f = pts
        
    #Plot result
    
    
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    
    if core_samples_mask is None:
        core_samples_mask = np.zeros_like(labels, dtype=bool)
    
    # Black removed and is used for noise instead.
    unique_labels = set(labels)
    colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
    for k, col in zip(unique_labels, colors):
        if k == -1:
            # Black used for noise.
            col = 'k'
    
        class_member_mask = (labels == k)
    
        xy = PCA_2f[class_member_mask & core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
                 markeredgecolor='k', markersize=10)
    
        xy = PCA_2f[class_member_mask & ~core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
                 markeredgecolor='k', markersize=5)
    
    frame1 = plt.gca()
    ax = frame1.axes
    if names:
        for i, (x,y) in enumerate(PCA_2f):
            #r = 0.2+random()
            #angle
            dx = random()/5.0
            dy = random()/5.0
            annotate = True
            for j, (ox,oy) in enumerate(PCA_2f):
                
                if i!=j and d((x,y), (ox,oy))<0.2:
                    if random()<0.4:
                       annotate = False
            if annotate:
                ax.annotate(names[i], xy=(x,y), xytext=(x+dx,y+dy))
            
    
    #frame1.axes.get_xaxis().set_visible(False)
    #frame1.axes.get_yaxis().set_visible(False)
    
    frame1.axes.get_xaxis().set_ticklabels([])
    frame1.axes.get_yaxis().set_ticklabels([])
    frame1.axes.set_xlabel( xy_labels[0] )
    frame1.axes.set_ylabel( xy_labels[1] )
    
    plt.title('Number of clusters: %d' % n_clusters_)
    if outfiles:
        for outfile in outfiles:
            plt.savefig(outfile)
    plt.show()

#def plot3d:
def main():                    
    features_file_A = "vrp_features_Augerat-A-14.csv"
    features_file_B = "vrp_features_Augerat-B-14.csv"
    #features_file = "vrp_features_Augerat-AandB-14.csv"
    features_file = 'vrp_all_features_FINAL.csv'
    
    header, cols, float_headers, float_cols = read_features_from_file_and_validate(features_file)
    
    pca_and_cluster = True
    if pca_and_cluster:
        
        nM_f = preprocess_feature_data(float_cols)
        this is wrong! Should normalize first!
        pca7d = PCA(n_components=7)
        PCA_7f = pca7d.fit_transform(nM_f)
        
        print pca7d.explained_variance_ratio_ 
        print sum(pca7d.explained_variance_ratio_)
        exit()
        #print pca7d.components_.shape
        #pca7d.components_.tofile("7dcomponents_ndarray.csv")
        
        this is wrong! Should apply column by column!
        mms = MinMaxScaler()
        normalized_7df = mms.fit_transform(PCA_7f)
        
        DO_DBSCAN = True
        min_samples_aka_k = 3
        USE_ADAPTIVE_EPS = False
        ADAPTIVE_EPS_TARGET_CLUSTERS = 3
        PLOT_HUGE_IMAGE = False
        if USE_ADAPTIVE_EPS:
            
            #print len(squareform())
            nndist = np.sort( squareform(pdist(normalized_7df)), axis=1)
            #print nndist
            kdist = sorted(nndist[:,min_samples_aka_k])
            start_eps = np.median( kdist )
            end_eps = kdist[int(0.85*len(kdist))]
            eps_steps = 10
            
            start_eps = 0.50
            end_eps = 0.65
    
            for eps in np.linspace(start_eps, end_eps, eps_steps):
                db = DBSCAN(eps=eps, min_samples=min_samples_aka_k).fit( normalized_7df  )
                labels = db.labels_
                n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
                print eps, n_clusters_, np.sum( labels==-1 )
    
        elif DO_DBSCAN:
            eps_for_all = 0.20
            eps_for_Aset = 0.85
            eps_for_ABset = 0.60
            db = DBSCAN(eps=eps_for_all, min_samples=3).fit( normalized_7df  )
            labels = db.labels_
            
            core_samples_mask = np.zeros_like(labels, dtype=bool)
            core_samples_mask[db.core_sample_indices_] = True        
    
            pca2d = PCA(n_components=2)
            pts = pca2d.fit_transform(nM_f)
            
            if PLOT_HUGE_IMAGE:
                imagefiles = [r"C:\MyTemp\TEMP\Benchmarks_with_images\Benchmarks\all\%s.png" % iname for iname in cols[1]]
                plot_clusters_with_images(pts, labels, imagefiles, img_size = (16384,12288), img_file='out.tiff', bg_colors=True)
                #plot_clusters_with_images(pts, labels, imagefiles, img_file='out.tiff')
            else:
                shortnames = [name_translator_table[n] for n in cols[1]]
                plot_clusters(nM_f, labels, ["clusters.png", "clusters.ps", "clusters.pdf"], core_samples_mask, shortnames)
            
            
    
        else:
            centroids, labels, _ = k_means(normalized_7df, ADAPTIVE_EPS_TARGET_CLUSTERS)
            plot_clusters(nM_f, labels)
            
        
        
        
        print "Clustering result for A"
        instance_names = np.array(cols[1])
        for label in set(labels):
            print label, instance_names[labels==label]
    
        exit()
        
        header_B, cols_B, float_headers_B, float_cols_B = read_features_from_file_and_validate(features_file_B)
        
        nM_B_f = preprocess_feature_data(float_cols_B)
        PCA_B_7fN = pca7d.transform(nM_B_f)
        normalized_7df_B = mms.fit_transform(PCA_B_7fN)
        
        if DO_DBSCAN:
            db = DBSCAN(eps=eps_for_Aset, min_samples=3).fit( normalized_7df_B  )
            labelsB = db.labels_
        else:
            centroidsB, labelsB, _ = k_means(normalized_7df, ADAPTIVE_EPS_TARGET_CLUSTERS)
            
            
        plot_clusters(nM_B_f, labelsB)
        
        print "Clustering result for B using hte PCA transform of A"
        instance_names = np.array(cols_B[1])
        for label in set(labelsB):
            print label, instance_names[labels==label]
        
        AB = np.concatenate((normalized_7df,normalized_7df_B))
        D_AB = squareform(pdist(AB))
        A_cluster_labels_for_B = []
        for b_idx in range(14,28):
            closest_a_idx = np.argmin(D_AB[b_idx][0:14])
            #print np.min(D_AB[b_idx][0:14])
            A_cluster_labels_for_B.append(labels[closest_a_idx])
            
        print "Clustering B using the clusters for A"
        instance_names = np.array(cols_B[1])
        for label in set(A_cluster_labels_for_B):
            print label, instance_names[A_cluster_labels_for_B==label]
            
if __name__=="__main__":
    main()    
            
    
    
        
        
    