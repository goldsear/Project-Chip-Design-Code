import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import colorsys
import json
import sys

cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
def load_instance(filename):
  with open(filename) as input_file:
    return json.load(input_file)

def adjust_lightness(color, amount):
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

def plot_data(data):
  plt.figure()
  ax = plt.axes()
  ax.clear()
  bb = data['bb']   # bounding box

  ax.set_aspect('equal')
  ax.get_xaxis().set_visible(False)
  ax.get_yaxis().set_visible(False)
  ax.set_xlim(bb[0], bb[2])
  ax.set_ylim(bb[1], bb[3])

  c_num = 0

  for rect in data['rectangleData']:
    r = rect['coords']    # rectangle
    m = rect['margins']   # margin

    col = cycle[c_num]

    # Plot "outer" darker margin
    #ax.add_patch(mpatches.Rectangle((r[0]-m[0], r[1]-m[1]),
    #                       width=((r[2]+m[2])-(r[0]-m[0])),
    #                       height=((r[3]+m[3])-(r[1]-m[1])), alpha=0.8, facecolor=col))
    
    # Plot "inner" lighter rectangle (just plots over the margin, covering it up)
    ax.add_patch(mpatches.Rectangle((r[0], r[1]), width=r[2]-r[0], height=r[3]-r[1], facecolor=adjust_lightness(col, amount=1.3)))

    c_num += 1
    c_num %= len(cycle)
  plt.show()

def plot_final_data(init_data, final_data):
  ax = plt.axes()


  bb = final_data['bb']   # bounding box

  ax.set_aspect('equal')
  ax.get_xaxis().set_visible(False)
  ax.get_yaxis().set_visible(False)
  ax.set_xlim(bb[0], bb[2])
  ax.set_ylim(bb[1], bb[3])

  c_num = 0

  for rect in final_data['rectangleData']:
    r = rect['coords']    # rectangle
    m = rect['margins']   # margin
    col = cycle[c_num]

    # Plot "outer" darker margin
    ax.add_patch(mpatches.Rectangle((r[0]-m[0], r[1]-m[1]), 
                           width=((r[2]+m[2])-(r[0]-m[0])), 
                           height=((r[3]+m[3])-(r[1]-m[1])), alpha=0.8, facecolor=col))
    
    # Plot "inner" lighter rectangle (just plots over the margin, covering it up)
    ax.add_patch(mpatches.Rectangle((r[0], r[1]), width=r[2]-r[0], height=r[3]-r[1], facecolor=adjust_lightness(col, amount=1.5)))

    c_num += 1
    c_num %= len(cycle)

  for i in range(len(final_data['rectangleData'])):
    # Plot lines for macro movement
    start_box = init_data['rectangleData'][i]['coords']
    start_x, start_y = ((start_box[2] + start_box[0])/2, (start_box[3] + start_box[1])/2)
    end_box = final_data['rectangleData'][i]['coords']
    end_x, end_y = ((end_box[2] + end_box[0]) / 2, (end_box[3] + end_box[1]) / 2)
    ax.plot((start_x, end_x), (start_y, end_y), c="black")
  plt.show()


plt.rcParams["figure.figsize"] = (10,10)
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']   # color_cycle  

if __name__ == "__main__":
    if len(sys.argv) < 3:
        data = load_instance(sys.argv[1])
        plot_data(data)
    else:
        init_data = load_instance(sys.argv[1])
        final_data = load_instance(sys.argv[2])
        plot_final_data(init_data, final_data)
