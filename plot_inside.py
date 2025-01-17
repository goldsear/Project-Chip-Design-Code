import matplotlib.pyplot as plt
from shapely.geometry import box
import itertools
import json
import sys
import read
import numpy as np
# inst_001.json

def plot(instance, file = True):
 if file:
  data = read.load_instance(instance)
 else:
  data = instance

 ax = plt.axes()
 # Fix aspect ratio to 1:1
 ax.set_aspect('equal')

 rectangles = []

 for rect in data["rectangleData"]:
   min_x, min_y, max_x, max_y = rect['coords']
   rectangles.append(box(min_x, min_y, max_x, max_y))

 for rect in rectangles:
   ax.fill(*rect.exterior.xy, alpha = 0.5)
   #ax.plot(*rect.exterior.xy)

 # Calculate and visualize overlap between rectangles
 total_inter = 0
 for a, b in itertools.combinations(rectangles, 2):
     if a.intersects(b):
       intersection = a.intersection(b)
       total_inter += intersection.area
       #ax.fill(*intersection.exterior.xy, color='black')

 xgrid = data['bb'][0] + np.arange(0, data['bb'][2], data['bb'][2])
 ax.set_xticks(xgrid)
 ygrid = data['bb'][1] + np.arange(0, data['bb'][3], data['bb'][3])
 ax.set_yticks(ygrid)
 plt.grid()

 print("Total area of intersection is: ", total_inter, "square units")

 ax.set_xlim(data['bb'][0], data['bb'][2])
 ax.set_ylim(data['bb'][1], data['bb'][3])

 plt.show()


