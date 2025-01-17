import read
import os
from sys import argv

instance = read.load_instance(argv[1])

folder_name = argv[1]
folder_name = folder_name[0:len(folder_name)-5] #.json has length 5
os.makedirs(folder_name)
os.chdir(folder_name)

#create DEF file
with open(folder_name+'.def', 'w') as f:
 x = f.write(f'VERSION 5.7 ;\nDESIGN {folder_name}\nUNITS DISTANCE MICRONS 1 ;\n\n') #discard result x
 bb = instance['bb']
 x = f.write(f'DIEAREA ( {bb[0]} {bb[1]} ) ( {bb[2]} {bb[3]} ) ;\n\n')
 c = len(instance['rectangleData'])
 x = f.write(f'COMPONENTS {c} ;\n')
 i = 0
 for macro in instance['rectangleData']:
  x = f.write(f"   - o442432 block_{i}\n      + PLACED ( macro['coords'][0] macro['coords'][1] ) N ;\n")
  i += 1
 x = f.write(f'END COMPONENTS\n\n\n\nEND DESIGN\n\n')
 f.close()

#create LEF file
with open(folder_name+'.lef', 'w') as f:
 i = 0
 for macro in instance['rectangleData']:
  coords = macro['coords']
  x = f.write(f"MACRO block_{i}\n   SIZE {coords[2]-coords[0]} BY {coords[3]-coords[1]} ;\nEND block_{i}\n\n")
  i += 1
 x = f.write('END LIBRARY\n')
f.close()
 

#create TXT constraint file
#calculate constants, arbitrarily really
powerplan, spacing = 0.0, 0.0
for macro in instance['rectangleData']:
 margins = macro['margins']
 powerplan += margins[0]+margins[2]
 spacing += margins[1]+margins[3]
len = len(margins)
powerplan /= len
spacing /= len
with open(folder_name+'.txt', 'w') as f:
 x = f.write(f'powerplan_width_constraint {powerplan}\nminimum_channel_spacing_between_macros_constraint {spacing}\nbuffer_area_reservation_extended_distance_constraint 80\nweight_alpha 1\nweight_beta 0')

#powerplan_width_constraint 50 #universal x margin
#minimum_channel_spacing_between_macros_constraint 4 #universal y margin
#buffer_area_reservation_extended_distance_constraint 80 #not used
#weight_alpha 1 #cost multiplier
#weight_beta 4 #some other cost multiplier, for this use set to 0


