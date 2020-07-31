# THIS is from https://github.com/yorak/CVRPFeatureExtractors/blob/78af90e8b44e25000cf6d0fe8bbe3e6574d6ee9b/helpers/rkm_util.py#L13
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
    
from random import uniform
import matplotlib.pyplot as plt
import matplotlib.patches as patches

pts = []
for i in range(100):
    pts.append( ( uniform(-400,400), uniform(-200,200) ) )
pts_rect = list(normalize_to_rect(pts, (0,400), False))
pts_kar = list(normalize_to_rect(pts, (0,400), True))

plt.figure(figsize=(9, 3))
colors = ['b','r', 'g']
for i, p in enumerate([pts, pts_rect, pts_kar]):
  ax = plt.subplot(130+(i+1))
  plt.ylim((-500,500)); plt.xlim((-500,500));
  x, y = zip(*p) ; plt.scatter( x,y, c=colors[i] );
  if i>0:
    rect = patches.Rectangle((0,0),400,400, linewidth=1, 
      edgecolor='gray', facecolor='none')
    ax.add_patch(rect)

plt.suptitle('Original / Rect scaling / Shape preserving scaling')
plt.show()

