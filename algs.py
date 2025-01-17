from mip import *
from random import random
from random import choice

import read
from verify import Rectangle
from verify import verify
from verify import EPSILON
from verify import verify_macro
from verify import score
from verify import directions_of_illegality

from copy import deepcopy
from itertools import combinations

from fractions import Fraction as Frac
from math import exp
from math import sqrt

#import networkx as nx

from time import sleep
from time import time

def violates_bounding_box(macro, bounding_box):
 if macro['coords'][0]-macro['margins'][0] < bounding_box[0]-EPSILON:
  return 'left'
 if macro['coords'][2]+macro['margins'][2] > bounding_box[2]+EPSILON:
  return 'right'
 if macro['coords'][1]-macro['margins'][1] < bounding_box[1]-EPSILON:
  return 'bottom'
 if macro['coords'][3]+macro['margins'][3] > bounding_box[3]+EPSILON:
  return 'top'

 return None

#shuffles the order of the rectangleData of an instance
#order seems to influence the performance of the mip solver, very strange
#rectData is assumed to be a list of macro info, like in the json files
def shuffle(rectData):
 sh_rectData = {} #need to remember original order
 l = list(range(len(rectData))) #remaining indices
 while l != []:
  r = choice(l)
  sh_rectData[r] = rectData[r]
  for item in range(len(l)):
   if l[item] == r:
    del l[item]
    break

 return sh_rectData

#solves an initial instance to a legal instance using an exact mip_solver
#only for the two linear metrics
#instance is assumed to be a dict in the Python version of the json format
#metric is 'd1' or 'd1w', the latter for the weighted version
#MAX_DURATION in seconds, default 1800 (30 min)
def mip_solver(instance, metric, MAX_DURATION = 1800):
 m = Model()
 original = deepcopy(instance)
 instance = deepcopy(instance) #need the original later on for the comparison
 instance = read.to_int(instance) #everything to int for grid alignment
 instance = read.Problem.load(instance) #convert to read.Problem object, for easy properties

 epsx = EPSILON/instance.grid[2]/10.0
 epsy = EPSILON/instance.grid[3]/10.0

 bleft = str(instance.bb[0])
 bleft = Frac(bleft)*(1+EPSILON)
 bbottom = str(instance.bb[1])
 bbottom = Frac(bbottom)*(1+EPSILON)
 bright = str(instance.bb[2])
 bright = Frac(bright)*(1+EPSILON)
 btop = str(instance.bb[3])
 btop = Frac(btop)*(1+EPSILON)
 M = (bright-bleft)+(btop-bbottom) #large number in the implementation

 #add the variables
 x = {}
 y = {}
 absx = {}
 absy = {}
 for macro in instance.macros:
  mleft = str(macro.left)
  mleft = Frac(mleft)
  #mright = str(macro.right) #not used apparently
  #mright = Frac(bright)
  #mtop = str(macro.top) #not used apparently
  #mtop = Frac(btop) #not used apparently
  mbottom = str(macro.bottom)
  mbottom = Frac(mbottom)
  mperimeter = str(macro.perimeter)
  mperimeter = Frac(mperimeter)
  b_x = Frac(str(instance.xspacing_b))
  c_x = Frac(str(instance.xspacing_c))
  b_y = Frac(str(instance.yspacing_b))
  c_y = Frac(str(instance.yspacing_c))
  #macro coordinates: int because of grid alignment
  if not macro.blockage:
   x[macro] = m.add_var(var_type = INTEGER) #left coo of macro
   y[macro] = m.add_var(var_type = INTEGER) #bottom coo of macro

   #absolute values in the objective
   absx[macro] = m.add_var()
   absy[macro] = m.add_var()

   #absolute value constraints
   if metric == 'd1':
    m += absx[macro] >= x[macro]-mleft
    m += absx[macro] >= mleft-x[macro]
    m += absy[macro] >= y[macro]-mbottom
    m += absy[macro] >= mbottom-y[macro]
   elif metric == 'd1w':
    m += absx[macro] >= mperimeter*(x[macro]-mleft)
    m += absx[macro] >= mperimeter*(mleft-x[macro])
    m += absy[macro] >= mperimeter*(y[macro]-mbottom)
    m += absy[macro] >= mperimeter*(mbottom-y[macro])

   #bounding box constraints
   m += x[macro] >= bleft+macro.margins.left-epsx
   m += y[macro] >= bbottom+macro.margins.bottom-epsy
   m += x[macro]+macro.width <= bright-macro.margins.right+epsx
   m += y[macro]+macro.height <= btop-macro.margins.top+epsy


 d = {}

 #combinations only gives the 'upper triangle' of combinations!
 #x = (0,1,2,3)
 #combinations(x,2) = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]
 for macro1, macro2 in combinations(instance.macros, 2):
  if macro1 not in d.keys():
   d[macro1] = {}

  #assume macro1, macro2 are not blockages
  if not macro1.blockage and not macro2.blockage:
   #spacing rule vars
   d[macro1][macro2] = {}
   d[macro1][macro2]['rightb'] = m.add_var(var_type = BINARY)
   d[macro1][macro2]['rightc'] = m.add_var(var_type = BINARY)
   d[macro1][macro2]['leftb'] = m.add_var(var_type = BINARY)
   d[macro1][macro2]['leftc'] = m.add_var(var_type = BINARY)
   d[macro1][macro2]['belowb'] = m.add_var(var_type = BINARY)
   d[macro1][macro2]['belowc'] = m.add_var(var_type = BINARY)
   d[macro1][macro2]['aboveb'] = m.add_var(var_type = BINARY)
   d[macro1][macro2]['abovec'] = m.add_var(var_type = BINARY)

   #spacing constraints
   #macro1 right of macro2:
   m += x[macro1]-x[macro2]-macro2.width+M*(1-d[macro1][macro2]['rightb']) >= b_x-epsx
   m += x[macro1]-x[macro2]-macro2.width-M*(1-d[macro1][macro2]['rightb']) <= b_x+epsx
   m += x[macro1]-x[macro2]-macro2.width+M*(1-d[macro1][macro2]['rightc']) >= c_x-epsx

   #macro1 left of macro2:
   m += -x[macro1]-macro1.width+x[macro2]+M*(1-d[macro1][macro2]['leftb']) >= b_x-epsx
   m += -x[macro1]-macro1.width+x[macro2]-M*(1-d[macro1][macro2]['leftb']) <= b_x+epsx
   m += -x[macro1]-macro1.width+x[macro2]+M*(1-d[macro1][macro2]['leftc']) >= c_x-epsx

   #macro1 below macro2:
   m += -y[macro1]-macro1.height+y[macro2]+M*(1-d[macro1][macro2]['belowb']) >= b_y-epsy
   m += -y[macro1]-macro1.height+y[macro2]-M*(1-d[macro1][macro2]['belowb']) <= b_y+epsy
   m += -y[macro1]-macro1.height+y[macro2]+M*(1-d[macro1][macro2]['belowc']) >= c_y-epsy

   #macro1 above macro2:
   m += y[macro1]-y[macro2]-macro2.height+M*(1-d[macro1][macro2]['aboveb']) >= b_y-epsy
   m += y[macro1]-y[macro2]-macro2.height-M*(1-d[macro1][macro2]['aboveb']) <= b_y+epsy
   m += y[macro1]-y[macro2]-macro2.height+M*(1-d[macro1][macro2]['abovec']) >= c_y-epsy

   #at least one of the d[m1][m2][x/y][b/c] has to be 1
   m += xsum(d[macro1][macro2][i] for i in ['rightb','leftb','rightc','leftc','belowb','aboveb','belowc','abovec']) >= 1

   #margin vars
   d[macro1][macro2]['margin'] = [m.add_var(var_type = BINARY) for i in range(4)] #right, left, bottom, top

   #macro1 to the right, left, below and above macro2 respectively
   m += x[macro1]-macro1.margins.left-macro2.margins.right-macro2.width-x[macro2]+M*(1-d[macro1][macro2]['margin'][0]) >= 0
   m += -x[macro1]-macro1.width-macro1.margins.right-macro2.margins.left+x[macro2]+M*(1-d[macro1][macro2]['margin'][1]) >= 0
   m += -y[macro1]-macro1.height-macro1.margins.top-macro2.margins.bottom+y[macro2]+M*(1-d[macro1][macro2]['margin'][2]) >= 0
   m += y[macro1]-macro1.margins.bottom-macro2.margins.top-macro2.width-y[macro2]+M*(1-d[macro1][macro2]['margin'][3]) >= 0

   #at least one of the constraints true
   m += xsum(d[macro1][macro2]['margin'][i] for i in range(4)) >= 1

  #blockages; macro1 is the blockage
  elif macro1.blockage and not macro2.blockage:
   d[macro1][macro2] = [m.add_var(var_type = BINARY) for i in range(4)] #right, left, bottom, top

   #note: x[macro1],y[macro1] don't exist! (not variables, because of blockage)
   m += macro1.left-macro2.margins.right-macro2.width-x[macro2]+M*(1-d[macro1][macro2][0]) >= 0 #blockage right of macro
   m += -macro1.right-macro2.margins.left+x[macro2]+M*(1-d[macro1][macro2][1]) >= 0 #left
   m += -macro1.top-macro2.margins.bottom+y[macro2]+M*(1-d[macro1][macro2][2]) >= 0 #below
   m += macro1.bottom-macro2.margins.top-macro2.height-y[macro2]+M*(1-d[macro1][macro2][3]) >= 0 #above

   #at least one of the constraints true
   m += xsum(d[macro1][macro2][i] for i in range(4)) >= 1

  #blockages; macro2 is the blockage
  elif not macro1.blockage and macro2.blockage:
   #blockages; macro2 is the blockage
   d[macro1][macro2] = [m.add_var(var_type = BINARY) for i in range(4)] #right, left, bottom, top

   #note: x[macro2],y[macro2] don't exist! (not variables, because of blockage)
   m += macro2.left-macro1.margins.right-macro1.width-x[macro1]+M*(1-d[macro1][macro2][0]) >= 0 #blockage right of macro
   m += -macro2.right-macro1.margins.left+x[macro1]+M*(1-d[macro1][macro2][1]) >= 0 #left
   m += -macro2.top-macro1.margins.bottom+y[macro1]+M*(1-d[macro1][macro2][2]) >= 0 #below
   m += macro2.bottom-macro1.margins.top-macro1.height-y[macro1]+M*(1-d[macro1][macro2][3]) >= 0 #above

   #at least one of the constraints true
   m += xsum(d[macro1][macro2][i] for i in range(4)) >= 1

  #macro1, macro2 both blockages: do nothing, or continue
 
 #metric included in def of absx, absy
 #xsum is mip lingo for sum of variables
 m.objective = original['grid'][2]*xsum(absx[macro] for macro in instance.macros if not macro.blockage)+original['grid'][3]*xsum(absy[macro] for macro in instance.macros if not macro.blockage)

 #make sure the algorithm terminates after 30 minutes
 status = m.optimize(max_seconds = MAX_DURATION)
 
 print(status)
 if status == OptimizationStatus.OPTIMAL:
  print('optimal solution found')
 elif status == OptimizationStatus.FEASIBLE:
  print('feasible solution found')
 elif status == OptimizationStatus.ERROR:
  print('an error of unknown type occurred')
 elif status == OptimizationStatus.INFEASIBLE:
  print('problem not solvable')
 elif status == OptimizationStatus.UNBOUNDED:
  print('problem unbounded: forgot a constraint?')
 elif status == OptimizationStatus.INT_INFEASIBLE:
  print('problem only solvable in the lp case')
 elif status == OptimizationStatus.NO_SOLUTION_FOUND:
  print('no integer solution found')

 #load bb from original
 #load spacing_rules from original
 #load margins from original
 #load grid from original

 rectangleData = []

 #x = [A, B, C, D]
 #y = [0, 1]
 #zip(x, y) = [(A, 0), (B, 1)]

 #can only fetch if a feasible solution is found
 if status in [OptimizationStatus.OPTIMAL, OptimizationStatus.FEASIBLE]:
  for macro, oldmacro in zip(instance.macros, original['rectangleData']):
   #rectangleData[oldmacro] = {'coords': [x[macro].x, y[macro].x, x[macro.x]+macro.width, y[macro].x+macro.height], 'blockage': original['rectangleData'][oldmacro]['blockage'], 'margins': original['rectangleData'][oldmacro]['margins']}
   rectangle = {}

   #convert back to original coordinates
   if macro.blockage:
    rectangle['coords'] = oldmacro['coords']
   else:
    rectangle['coords'] = [x[macro].x*original['grid'][2], y[macro].x*original['grid'][3], (x[macro].x+macro.width)*original['grid'][2], (y[macro].x+macro.height)*original['grid'][3]]

   rectangle['blockage'] = oldmacro['blockage']
   rectangle['margins'] = oldmacro['margins']
   rectangleData.append(rectangle)

 new_instance = {}
 new_instance['spacing_rule_widths'] = original['spacing_rule_widths']
 new_instance['spacing_rule_heights'] = original['spacing_rule_heights']
 new_instance['rectangleData'] = rectangleData #this thing has changed (empty if no feasible solution)
 new_instance['bb'] = original['bb']
 new_instance['grid'] = original['grid']

 return new_instance

#it appears the mip solver's result depends on the order in which the rectangles are loaded
#this function shuffles an instance randomly, in order to try different orderings
def shuffle_instance(instance):
 instance = deepcopy(instance) #make sure instance isn't modified
 new_instance = {}
 new_instance['spacing_rule_widths'] = instance['spacing_rule_widths']
 new_instance['spacing_rule_heights'] = instance['spacing_rule_heights']

 rectangleData = []
 lis = list(range(len(instance['rectangleData'])))
 order = []
 while lis != []:
  r = choice(lis)
  order.append(r)
  rectangleData.append(instance['rectangleData'][r])
  for item in range(len(lis)):
   if lis[item] == r:
    del lis[item]
    break

 new_instance['rectangleData'] = rectangleData #this thing has changed (empty if no feasible solution)

 new_instance['bb'] = instance['bb']
 new_instance['grid'] = instance['grid']

 return new_instance, order

#dual function to above
#brings the macros back in their original order
def restore_order(instance, order):
 new_instance = {}
 new_instance['spacing_rule_widths'] = instance['spacing_rule_widths']
 new_instance['spacing_rule_heights'] = instance['spacing_rule_heights']

 rectangleData = []
 for index in range(len(order)):
  for index2 in range(len(order)):
   if order[index2] == index:
    rectangleData.append(instance['rectangleData'][index2])

 new_instance['rectangleData'] = rectangleData

 new_instance['bb'] = instance['bb']
 new_instance['grid'] = instance['grid']

 return new_instance

#version of the mip solver that shuffles the instance
#it appears the mip solver's result depends on the order in which the rectangles are loaded
#hopefully gives a better result
def mip_shuffle(instance, metric, maxiter, MAX_DURATION = 1800):
 instance = deepcopy(instance)
 best_score = float('inf')
 best_solution = []
 for index in range(maxiter):
  sh_instance, order = shuffle_instance(instance)
  new_instance = mip_solver(sh_instance, metric, MAX_DURATION/maxiter)
  if verify(new_instance, sh_instance):
   s = score(new_instance, sh_instance)
   sc = float('inf')
   if metric == 'd1':
    sc = s[0]
   elif metric == 'd1w':
    sc = s[1]
   if sc < best_score:
    best_solution = new_instance
    best_solution = restore_order(best_solution, order)
    best_score = sc

 return best_solution


#implements the greedy algorithm from 4.3.2
#metric is 'd1' or 'd1w'
#maxiter determines the number of tries upon no improvement
#large values of maxiter are more likely to produce good results, but take longer to compute
#again, MAX_DURATION in seconds, default 1800
def greedy(instance, metric, maxiter, MAX_DURATION = 1800):
 instance = deepcopy(instance) #need the original later on for the comparison
 _original = read.Problem.load(instance) #again, need another copy of the original in read.Problem format
 instance = read.to_int(instance) #everything to int for grid alignment
 original = read.Problem.load(instance) #copy, to revert after an iteration of the algorithm
 instance = read.Problem.load(instance) #convert to read.Problem object, for easy properties
 M = (instance.bb[2]-instance.bb[0])+(instance.bb[3]-instance.bb[1]) #large number in the implementation

 S = {}
 #need ordering, to be updated after each round
 for i in range(len(instance.macros)):
  S[i] = i #initially, set up as the original ordering

 n = 0
 rectangleData = []
 currentBest = float('inf') #current best solution value, initialized as infinity
 while n < maxiter:
  solution_value = {} #need individual values for reordering

  #fixed macros in current round of the algorithm
  SF = []

  #restore original macro coords for repeating of the algorithm
  #additionally, add the blockages to SF, as every macro has to deal with them
  for macro, oldmacro in zip(instance.macros, original.macros):
   for i in range(4):
    macro.coords[i] = oldmacro.coords[i] #per value, to prevent pointer-entanglement
   if macro.blockage:
    SF.append(macro)

  rectDataCurrent = {} #rectangle data for current best solution
  for j in range(len(instance.macros)):
   macro1 = instance.macros[S[j]]
   oldmacro = original.macros[S[j]]
   _oldmacro = _original.macros[S[j]]

   if not macro1.blockage:
    m = Model()
    x = m.add_var(var_type = INTEGER)
    y = m.add_var(var_type = INTEGER)
    absx = m.add_var()
    absy = m.add_var()

    #absolute value constraints
    if metric == 'd1':
     m += absx >= x-macro1.left
     m += absx >= macro1.left-x
     m += absy >= y-macro1.bottom
     m += absy >= macro1.bottom-y
    elif metric == 'd1w':
     m += absx >= macro1.perimeter*(x-macro1.left)
     m += absx >= macro1.perimeter*(macro1.left-x)
     m += absy >= macro1.perimeter*(y-macro1.bottom)
     m += absy >= macro1.perimeter*(macro1.bottom-y)

    #bounding box constraints
    m += x >= instance.bb[0]+macro1.margins.left
    m += y >= instance.bb[1]+macro1.margins.bottom
    m += x+macro1.width <= instance.bb[2]-macro1.margins.right
    m += y+macro1.height <= instance.bb[3]-macro1.margins.top

    d = {}
    for macro2 in SF:
     if macro2 not in d.keys():
      d[macro2] = {}

     if not macro2.blockage:
      #spacing rule vars
      d[macro2]['rightb'] = m.add_var(var_type = BINARY)
      d[macro2]['rightc'] = m.add_var(var_type = BINARY)
      d[macro2]['leftb'] = m.add_var(var_type = BINARY)
      d[macro2]['leftc'] = m.add_var(var_type = BINARY)
      d[macro2]['belowb'] = m.add_var(var_type = BINARY)
      d[macro2]['belowc'] = m.add_var(var_type = BINARY)
      d[macro2]['aboveb'] = m.add_var(var_type = BINARY)
      d[macro2]['abovec'] = m.add_var(var_type = BINARY)

      #constraints
      #macro1 right of macro2:
      m += x-macro2.right+M*(1-d[macro2]['rightb']) >= instance.xspacing_b #b_x
      m += x-macro2.right-M*(1-d[macro2]['rightb']) <= instance.xspacing_b #b_x
      m += x-macro2.right+M*(1-d[macro2]['rightc']) >= instance.xspacing_c #c_x

      #macro1 left of macro2:
      m += -x-macro1.width+macro2.left+M*(1-d[macro2]['leftb']) >= instance.xspacing_b #b_x
      m += -x-macro1.width+macro2.left-M*(1-d[macro2]['leftb']) <= instance.xspacing_b #b_x
      m += -x-macro1.width+macro2.left+M*(1-d[macro2]['leftc']) >= instance.xspacing_c #c_x

      #macro1 below macro2:
      m += -y-macro1.height+macro2.bottom+M*(1-d[macro2]['belowb']) >= instance.yspacing_b #b_y
      m += -y-macro1.height+macro2.bottom-M*(1-d[macro2]['belowb']) <= instance.yspacing_b #b_y
      m += -y-macro1.height+macro2.bottom+M*(1-d[macro2]['belowc']) >= instance.yspacing_c #c_y

      #macro1 above macro2:
      m += y-macro2.top+M*(1-d[macro2]['aboveb']) >= instance.yspacing_b #b_y
      m += y-macro2.top-M*(1-d[macro2]['aboveb']) <= instance.yspacing_b #b_y
      m += y-macro2.top+M*(1-d[macro2]['abovec']) >= instance.yspacing_c #c_y

      #at least one of the d[m2][x/y][b/c] has to be 1
      m += xsum(d[macro2][i] for i in ['rightb','leftb','rightc','leftc','belowb','aboveb','belowc','abovec']) >= 1

      #margin vars
      d[macro2]['margin'] = [m.add_var(var_type = BINARY) for i in range(4)] #right, left, bottom, top

      #macro1 to the right, left, below and above macro2 respectively
      m += x-macro1.margins.left-macro2.margins.right-macro2.right+M*(1-d[macro2]['margin'][0]) >= 0
      m += -x-macro1.width-macro1.margins.right-macro2.margins.left+macro2.left+M*(1-d[macro2]['margin'][1]) >= 0
      m += -y-macro1.height-macro1.margins.top-macro2.margins.bottom+macro2.bottom+M*(1-d[macro2]['margin'][2]) >= 0
      m += y-macro1.margins.bottom-macro2.margins.top-macro2.top+M*(1-d[macro2]['margin'][3]) >= 0

      #at least one of the constraints true
      m += xsum(d[macro2]['margin'][i] for i in range(4)) >= 1

     #macro2 is a blockage
     else:
      d[macro2] = [m.add_var(var_type = BINARY) for i in range(4)] #right, left, bottom, top

      m += macro2.left-macro1.margins.right-macro1.width-x+M*(1-d[macro2][0]) >= 0 #blockage right of macro
      m += -macro2.right-macro1.margins.left+x+M*(1-d[macro2][1]) >= 0 #left
      m += -macro2.top-macro1.margins.bottom+y+M*(1-d[macro2][2]) >= 0 #below
      m += macro2.bottom-macro1.margins.top-macro1.height-y+M*(1-d[macro2][3]) >= 0 #above

      #at least one of the constraints true
      m += xsum(d[macro2][i] for i in range(4)) >= 1

    #macro1, macro2 both blockages: do nothing, or continue
 
    #metric included in def of absx, absy
    m.objective = original.grid[2]*absx+original.grid[3]*absy

    #make sure the algorithm terminates after 30 minutes
    status = m.optimize(max_seconds = MAX_DURATION/maxiter) #needs to be updated
    print(status)
    if status in [OptimizationStatus.OPTIMAL, OptimizationStatus.FEASIBLE]:
     solution_value[S[j]] = m.objective_value
     rectData = {'coords': [x.x*original.grid[2], y.x*original.grid[3], (x.x+macro1.width)*original.grid[2], (y.x+macro1.height)*original.grid[3]], 'blockage': _oldmacro.blockage, 'margins': _oldmacro.margins}
     rectDataCurrent[S[j]] = rectData
     print(rectData)
     #edit macro1 to be appended to SF with the new coordinates according to rectData
     #note: to_inted data needed!
     macro1.coords = Rectangle([x.x, y.x, x.x+macro1.width, y.x+macro1.height])

    else:
     solution_value[S[j]] = float('inf') #no solution, so infinity
    print('solution_value = ', solution_value)

   #macro1 is blockage: still need rectangle data and solution value zero (as immovable)
   else:
    rectDataCurrent[S[j]] = {'coords': [oldmacro.coords[0]*original.grid[2], oldmacro.coords[1]*original.grid[3], oldmacro.coords[2]*original.grid[2], oldmacro.coords[3]*original.grid[3]], 'blockage': _oldmacro.blockage, 'margins': _oldmacro.margins}
    solution_value[S[j]] = 0

   SF.append(macro1)

  if sum(solution_value.values()) < currentBest:
   currentBest = sum(solution_value.values())
   print('current best = ', currentBest)
   n = 0
   rectangleData = []
   for i in range(len(instance.macros)):
    rectangleData.append(rectDataCurrent[i]) #keep order for comparison
  else:
   n += 1

  #reorder S
  #sorted by value, keeping keys
  #sorted low to high
  solution_value_sorted = {k: v for k, v in sorted(solution_value.items(), key = lambda item: item[1])}

  #need high to low, hence reverse
  #item picks the key rather than value
  i = len(instance.macros)-1
  for item in solution_value_sorted:
   S[i] = item
   i -= 1

 solution = {}
 solution['spacing_rule_widths'] = instance['spacing_rule_widths']
 solution['spacing_rule_heights'] = instance['spacing_rule_heights']
 solution['rectangleData'] = rectangleData
 solution['bb'] = instance['bb']
 solution['grid'] = instance['grid']
 return solution

#relaxation of the mip solver, using continuous variables and a slight extension of the margins
#only for the two linear metrics
#instance is assumed to be a dict in the Python version of the json format
#metric is 'd1' or 'd1w', the latter for the weighted version
#MAX_DURATION in seconds, default 1800 (30 min)
#relax_percent applies to outer spacing rules, for now
def mip_relaxed(instance, metric, relax_percent, MAX_DURATION = 1800):
 m = Model()
 original = deepcopy(instance)
 instance = deepcopy(instance) #need the original later on for the comparison
 instance = read.Problem.load(instance) #convert to read.Problem object, for easy properties

 epsx = EPSILON
 epsy = EPSILON

 bleft = str(instance.bb[0])
 bleft = Frac(bleft)*(1+EPSILON)
 bbottom = str(instance.bb[1])
 bbottom = Frac(bbottom)*(1+EPSILON)
 bright = str(instance.bb[2])
 bright = Frac(bright)*(1+EPSILON)
 btop = str(instance.bb[3])
 btop = Frac(btop)*(1+EPSILON)
 M = (bright-bleft)+(btop-bbottom) #large number in the implementation

 #add the variables
 x = {}
 y = {}
 absx = {}
 absy = {}
 for macro in instance.macros:
  mleft = str(macro.left)
  mleft = Frac(mleft)
  #mright = str(macro.right) #not used apparently
  #mright = Frac(bright)
  #mtop = str(macro.top) #not used apparently
  #mtop = Frac(btop) #not used apparently
  mbottom = str(macro.bottom)
  mbottom = Frac(mbottom)
  mperimeter = str(macro.perimeter)
  mperimeter = Frac(mperimeter)
  b_x = Frac(str(instance.xspacing_b))
  c_x = Frac(str(instance.xspacing_c))*Frac(1+relax_percent)
  b_y = Frac(str(instance.yspacing_b))
  c_y = Frac(str(instance.yspacing_c))*Frac(1+relax_percent)
  #macro coordinates: continuous in the relaxation
  if not macro.blockage:
   x[macro] = m.add_var() #left coo of macro
   y[macro] = m.add_var() #bottom coo of macro

   #absolute values in the objective
   absx[macro] = m.add_var()
   absy[macro] = m.add_var()

   #absolute value constraints
   if metric == 'd1':
    m += absx[macro] >= x[macro]-mleft
    m += absx[macro] >= mleft-x[macro]
    m += absy[macro] >= y[macro]-mbottom
    m += absy[macro] >= mbottom-y[macro]
   elif metric == 'd1w':
    m += absx[macro] >= mperimeter*(x[macro]-mleft)
    m += absx[macro] >= mperimeter*(mleft-x[macro])
    m += absy[macro] >= mperimeter*(y[macro]-mbottom)
    m += absy[macro] >= mperimeter*(mbottom-y[macro])

   #bounding box constraints
   m += x[macro] >= bleft+macro.margins.left-epsx
   m += y[macro] >= bbottom+macro.margins.bottom-epsy
   m += x[macro]+macro.width <= bright-macro.margins.right+epsx
   m += y[macro]+macro.height <= btop-macro.margins.top+epsy


 d = {}

 #combinations only gives the 'upper triangle' of combinations!
 #x = (0,1,2,3)
 #combinations(x,2) = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]
 for macro1, macro2 in combinations(instance.macros, 2):
  if macro1 not in d.keys():
   d[macro1] = {}

  #assume macro1, macro2 are not blockages
  if not macro1.blockage and not macro2.blockage:
   #spacing rule vars
   d[macro1][macro2] = {}
   d[macro1][macro2]['rightb'] = m.add_var(var_type = BINARY)
   d[macro1][macro2]['rightc'] = m.add_var(var_type = BINARY)
   d[macro1][macro2]['leftb'] = m.add_var(var_type = BINARY)
   d[macro1][macro2]['leftc'] = m.add_var(var_type = BINARY)
   d[macro1][macro2]['belowb'] = m.add_var(var_type = BINARY)
   d[macro1][macro2]['belowc'] = m.add_var(var_type = BINARY)
   d[macro1][macro2]['aboveb'] = m.add_var(var_type = BINARY)
   d[macro1][macro2]['abovec'] = m.add_var(var_type = BINARY)

   #spacing constraints
   #macro1 right of macro2:
   m += x[macro1]-x[macro2]-macro2.width+M*(1-d[macro1][macro2]['rightb']) >= b_x-epsx
   m += x[macro1]-x[macro2]-macro2.width-M*(1-d[macro1][macro2]['rightb']) <= b_x+epsx
   m += x[macro1]-x[macro2]-macro2.width+M*(1-d[macro1][macro2]['rightc']) >= c_x-epsx

   #macro1 left of macro2:
   m += -x[macro1]-macro1.width+x[macro2]+M*(1-d[macro1][macro2]['leftb']) >= b_x-epsx
   m += -x[macro1]-macro1.width+x[macro2]-M*(1-d[macro1][macro2]['leftb']) <= b_x+epsx
   m += -x[macro1]-macro1.width+x[macro2]+M*(1-d[macro1][macro2]['leftc']) >= c_x-epsx

   #macro1 below macro2:
   m += -y[macro1]-macro1.height+y[macro2]+M*(1-d[macro1][macro2]['belowb']) >= b_y-epsy
   m += -y[macro1]-macro1.height+y[macro2]-M*(1-d[macro1][macro2]['belowb']) <= b_y+epsy
   m += -y[macro1]-macro1.height+y[macro2]+M*(1-d[macro1][macro2]['belowc']) >= c_y-epsy

   #macro1 above macro2:
   m += y[macro1]-y[macro2]-macro2.height+M*(1-d[macro1][macro2]['aboveb']) >= b_y-epsy
   m += y[macro1]-y[macro2]-macro2.height-M*(1-d[macro1][macro2]['aboveb']) <= b_y+epsy
   m += y[macro1]-y[macro2]-macro2.height+M*(1-d[macro1][macro2]['abovec']) >= c_y-epsy

   #at least one of the d[m1][m2][x/y][b/c] has to be 1
   m += xsum(d[macro1][macro2][i] for i in ['rightb','leftb','rightc','leftc','belowb','aboveb','belowc','abovec']) >= 1

   #margin vars
   d[macro1][macro2]['margin'] = [m.add_var(var_type = BINARY) for i in range(4)] #right, left, bottom, top

   #macro1 to the right, left, below and above macro2 respectively
   m += x[macro1]-macro1.margins.left-macro2.margins.right-macro2.width-x[macro2]+M*(1-d[macro1][macro2]['margin'][0]) >= 0
   m += -x[macro1]-macro1.width-macro1.margins.right-macro2.margins.left+x[macro2]+M*(1-d[macro1][macro2]['margin'][1]) >= 0
   m += -y[macro1]-macro1.height-macro1.margins.top-macro2.margins.bottom+y[macro2]+M*(1-d[macro1][macro2]['margin'][2]) >= 0
   m += y[macro1]-macro1.margins.bottom-macro2.margins.top-macro2.width-y[macro2]+M*(1-d[macro1][macro2]['margin'][3]) >= 0

   #at least one of the constraints true
   m += xsum(d[macro1][macro2]['margin'][i] for i in range(4)) >= 1

  #blockages; macro1 is the blockage
  elif macro1.blockage and not macro2.blockage:
   d[macro1][macro2] = [m.add_var(var_type = BINARY) for i in range(4)] #right, left, bottom, top

   #note: x[macro1],y[macro1] don't exist! (not variables, because of blockage)
   m += macro1.left-macro2.margins.right-macro2.width-x[macro2]+M*(1-d[macro1][macro2][0]) >= 0 #blockage right of macro
   m += -macro1.right-macro2.margins.left+x[macro2]+M*(1-d[macro1][macro2][1]) >= 0 #left
   m += -macro1.top-macro2.margins.bottom+y[macro2]+M*(1-d[macro1][macro2][2]) >= 0 #below
   m += macro1.bottom-macro2.margins.top-macro2.height-y[macro2]+M*(1-d[macro1][macro2][3]) >= 0 #above

   #at least one of the constraints true
   m += xsum(d[macro1][macro2][i] for i in range(4)) >= 1

  #blockages; macro2 is the blockage
  elif not macro1.blockage and macro2.blockage:
   #blockages; macro2 is the blockage
   d[macro1][macro2] = [m.add_var(var_type = BINARY) for i in range(4)] #right, left, bottom, top

   #note: x[macro2],y[macro2] don't exist! (not variables, because of blockage)
   m += macro2.left-macro1.margins.right-macro1.width-x[macro1]+M*(1-d[macro1][macro2][0]) >= 0 #blockage right of macro
   m += -macro2.right-macro1.margins.left+x[macro1]+M*(1-d[macro1][macro2][1]) >= 0 #left
   m += -macro2.top-macro1.margins.bottom+y[macro1]+M*(1-d[macro1][macro2][2]) >= 0 #below
   m += macro2.bottom-macro1.margins.top-macro1.height-y[macro1]+M*(1-d[macro1][macro2][3]) >= 0 #above

   #at least one of the constraints true
   m += xsum(d[macro1][macro2][i] for i in range(4)) >= 1

  #macro1, macro2 both blockages: do nothing, or continue
 
 #metric included in def of absx, absy
 #xsum is mip lingo for sum of variables
 m.objective = xsum(absx[macro] for macro in instance.macros if not macro.blockage)+xsum(absy[macro] for macro in instance.macros if not macro.blockage)

 #make sure the algorithm terminates after 30 minutes
 status = m.optimize(max_seconds = MAX_DURATION)
 
 print(status)

 #load bb from original
 #load spacing_rules from original
 #load margins from original
 #load grid from original

 rectangleData = []

 #x = [A, B, C, D]
 #y = [0, 1]
 #zip(x, y) = [(A, 0), (B, 1)]

 #can only fetch if a feasible solution is found
 if status in [OptimizationStatus.OPTIMAL, OptimizationStatus.FEASIBLE]:
  for macro, oldmacro in zip(instance.macros, original['rectangleData']):
   #rectangleData[oldmacro] = {'coords': [x[macro].x, y[macro].x, x[macro.x]+macro.width, y[macro].x+macro.height], 'blockage': original['rectangleData'][oldmacro]['blockage'], 'margins': original['rectangleData'][oldmacro]['margins']}
   rectangle = {}

   #convert back to original coordinates
   if macro.blockage:
    rectangle['coords'] = oldmacro['coords']
   else:
    #round to grid
    left_coo = round(x[macro].x/original['grid'][2])*original['grid'][2]
    bottom_coo = round(y[macro].x/original['grid'][3])*original['grid'][3]
    rectangle['coords'] = [left_coo, bottom_coo, left_coo+macro.width, bottom_coo+macro.height]

   rectangle['blockage'] = oldmacro['blockage']
   rectangle['margins'] = oldmacro['margins']
   rectangleData.append(rectangle)

 new_instance = {}
 new_instance['spacing_rule_widths'] = original['spacing_rule_widths']
 new_instance['spacing_rule_heights'] = original['spacing_rule_heights']
 new_instance['rectangleData'] = rectangleData #this thing has changed (empty if no feasible solution)
 new_instance['bb'] = original['bb']
 new_instance['grid'] = original['grid']

 return new_instance
 


#relaxation of the mip solver, using continuous variables and a slight extension of the margins
#relax_percent applies to outer (c) spacing rules, for now
#keep original orientation for faster processing
#only for the two linear metrics
#instance is assumed to be a dict in the Python version of the json format
#metric is 'd1' or 'd1w', the latter for the weighted version
#MAX_DURATION in seconds, default 1800 (30 min)
def flexible_greedy(instance, metric, relax_percent, MAX_DURATION = 1800):
 m = Model()
 original = deepcopy(instance)
 instance = deepcopy(instance) #need the original later on for the comparison
 rectangleData = instance['rectangleData']

 epsx = EPSILON
 epsy = EPSILON

 bleft = Frac(str(original['bb'][0]))
 bbottom = Frac(str(original['bb'][1]))
 bright = Frac(str(original['bb'][2]))*Frac(1+EPSILON)
 btop = Frac(str(original['bb'][3]))*Frac(1+EPSILON)

 b_x = Frac(str(original['spacing_rule_widths'][0]))
 c_x = Frac(str(original['spacing_rule_widths'][1]))*Frac(1+relax_percent)
 b_y = Frac(str(original['spacing_rule_heights'][0]))
 c_y = Frac(str(original['spacing_rule_heights'][1]))*Frac(1+relax_percent) 

 #add the variables
 x = {}
 y = {}
 absx = {}
 absy = {}
 for macro in range(len(rectangleData)):
  if not rectangleData[macro]['blockage']:
   #coordinates as Fractions
   mleft = Frac(str(rectangleData[macro]['coords'][0]))
   mright = Frac(str(rectangleData[macro]['coords'][2]))
   mtop = Frac(str(rectangleData[macro]['coords'][3]))
   mbottom = Frac(str(rectangleData[macro]['coords'][1]))
   mwidth = mright-mleft
   mheight = mtop-mbottom
   mperimeter = mwidth+mheight

   #margins as Fractions
   mmarginsleft = Frac(str(rectangleData[macro]['margins'][0]))
   mmarginsbottom = Frac(str(rectangleData[macro]['margins'][1]))
   mmarginsright = Frac(str(rectangleData[macro]['margins'][2]))
   mmarginstop = Frac(str(rectangleData[macro]['margins'][3]))

   x[macro] = m.add_var() #left coo of macro
   y[macro] = m.add_var() #bottom coo of macro

   #absolute values in the objective
   absx[macro] = m.add_var()
   absy[macro] = m.add_var()

   #absolute value constraints
   if metric == 'd1':
    m += absx[macro] >= x[macro]-mleft
    m += absx[macro] >= mleft-x[macro]
    m += absy[macro] >= y[macro]-mbottom
    m += absy[macro] >= mbottom-y[macro]
   elif metric == 'd1w':
    m += absx[macro] >= mperimeter*(x[macro]-mleft)
    m += absx[macro] >= mperimeter*(mleft-x[macro])
    m += absy[macro] >= mperimeter*(y[macro]-mbottom)
    m += absy[macro] >= mperimeter*(mbottom-y[macro])

   #bounding box constraints
   m += x[macro] >= bleft+mmarginsleft-epsx
   m += y[macro] >= bbottom+mmarginsbottom-epsy
   m += x[macro]+mwidth <= bright-mmarginsright+epsx
   m += y[macro]+mheight <= btop-mmarginstop+epsy

 #combinations only gives the 'upper triangle' of combinations!
 #x = (0,1,2,3)
 #combinations(x,2) = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]
 for macro1, macro2 in combinations(range(len(rectangleData)), 2):
  #coordinates as Fractions
  m1left = Frac(str(rectangleData[macro1]['coords'][0]))
  m1right = Frac(str(rectangleData[macro1]['coords'][2]))
  m1top = Frac(str(rectangleData[macro1]['coords'][3]))
  m1bottom = Frac(str(rectangleData[macro1]['coords'][1]))
  m2left = Frac(str(rectangleData[macro2]['coords'][0]))
  m2right = Frac(str(rectangleData[macro2]['coords'][2]))
  m2top = Frac(str(rectangleData[macro2]['coords'][3]))
  m2bottom = Frac(str(rectangleData[macro2]['coords'][1]))

  m1width = m1right-m1left
  m1height = m1top-m1bottom
  m2width = m2right-m2left
  m2height = m2top-m2bottom

  #margins as Fractions
  m1marginsleft = Frac(str(rectangleData[macro1]['margins'][0]))
  m1marginsbottom = Frac(str(rectangleData[macro1]['margins'][1]))
  m1marginsright = Frac(str(rectangleData[macro1]['margins'][2]))
  m1marginstop = Frac(str(rectangleData[macro1]['margins'][3]))
  m2marginsleft = Frac(str(rectangleData[macro2]['margins'][0]))
  m2marginsbottom = Frac(str(rectangleData[macro2]['margins'][1]))
  m2marginsright = Frac(str(rectangleData[macro2]['margins'][2]))
  m2marginstop = Frac(str(rectangleData[macro2]['margins'][3]))

  #macro1 and macro2 are not blockages
  if not rectangleData[macro1]['blockage'] and not rectangleData[macro2]['blockage']:

   #skip the binaries

   #spacing constraints
   #macro1 right of macro2:
   if m1left+m1right > m2left+m2right:
    #margins need to fit b to be allowed to be placed at a distance b
    if m1top > m2bottom:
     if m2top > m1bottom:
      if m1left <= m2right+(b_x+c_x)/2 and m1marginsleft+m2marginsright <= b_x:
       m += x[macro1]-x[macro2]-m2width >= b_x-epsx
       m += x[macro1]-x[macro2]-m2width <= b_x+epsx
      else:
       m += x[macro1]-x[macro2]-m2width >= c_x-epsx
      #margin constraint
      m += x[macro1]-m1marginsleft-m2marginsright-m2width-x[macro2] >= 0

   #macro1 left of macro2:
   elif m1left+m1right <= m2left+m2right:
    #margins need to fit b to be allowed to be placed at a distance b
    if m1top > m2bottom:
     if m2top > m1bottom:
      if m2left <= m1right+(b_x+c_x)/2 and m2marginsleft+m1marginsright <= b_x:
       m += -x[macro1]-m1width+x[macro2] >= b_x-epsx
       m += -x[macro1]-m1width+x[macro2] <= b_x+epsx
      else:
       m += -x[macro1]-m1width+x[macro2] >= c_x-epsx
      #margin constraint
      m += -x[macro1]-m1width-m1marginsright-m2marginsleft+x[macro2] >= 0

   #macro1 below macro2:
   if m1top+m1bottom < m2top+m2bottom:
    #margins need to fit b to be allowed to be placed at a distance b
    if m1right > m2left:
     if m2right > m1left:
      if m1top <= m2bottom+(b_x+c_x)/2 and m1marginstop+m2marginsbottom <= b_y:
       m += -y[macro1]-m1height+y[macro2] >= b_y-epsy
       m += -y[macro1]-m1height+y[macro2] <= b_y+epsy
      else:
       m += -y[macro1]-m1height+y[macro2] >= c_y-epsy
      #margin constraint
      m += -y[macro1]-m1height-m1marginstop-m2marginsbottom+y[macro2] >= 0

   #macro1 above macro2:
   elif m1top+m1bottom >= m2top+m2bottom:
    #margins need to fit b to be allowed to be placed at a distance b
    if m1right > m2left:
     if m2right > m1left:
      if m2top <= m1bottom+(b_x+c_x)/2 and m2marginstop+m1marginsbottom <= b_y:
       m += y[macro1]-y[macro2]-m2height >= b_y-epsy
       m += y[macro1]-y[macro2]-m2height <= b_y+epsy
      else:
       m += y[macro1]-y[macro2]-m2height >= c_y-epsy
      #margin constraint
      m += y[macro1]-m1marginsbottom-m2marginstop-m2width-y[macro2] >= 0

      
  #blockages; macro1 is the blockage
  elif rectangleData[macro1]['blockage'] and not rectangleData[macro2]['blockage']:

   #note: x[macro1],y[macro1] don't exist! (not variables, because of blockage)
   #macro1 to the right of macro2
   if m1left+m1right > m2left+m2right:
    if m1top > m2bottom:
     if m2top > m1bottom:
      m += m1left-m2marginsright-m2width-x[macro2] >= 0

   #macro1 to the left of macro2
   else:
    if m1top > m2bottom:
     if m2top > m1bottom:
      m += -m1right-m2marginsleft+x[macro2] >= 0

   #macro1 below macro2
   if m1top+m1bottom < m2top+m2bottom:
    if m1right > m2left:
     if m2right > m1left:
      m += -m1top-m2marginsbottom+y[macro2] >= 0

   #macro1 above macro2
   else:
    if m1right > m2left:
     if m2right > m1left:
      m += m1bottom-m2marginstop-m2height-y[macro2] >= 0

  #blockages; macro2 is the blockage
  elif not rectangleData[macro1]['blockage'] and rectangleData[macro2]['blockage']:
   #blockages; macro2 is the blockage

   #note: x[macro2],y[macro2] don't exist! (not variables, because of blockage)
   #blockage right of macro
   if m1left+m1right > m2left+m2right:
    if m1top > m2bottom:
     if m2top > m1bottom:
      m += m2left-m1marginsright-m1width-x[macro1] >= 0

   #blockage left of macro
   else:
    if m1top > m2bottom:
     if m2top > m1bottom:
      m += -m2right-m1marginsleft+x[macro1] >= 0

   #blockage below macro
   if m1top+m1bottom < m2top+m2bottom:
    if m1right > m2left:
     if m2right > m1left:
      m += -m2top-m1marginsbottom+y[macro1] >= 0

   #blockage above macro
   else:
    if m1right > m2left:
     if m2right > m1left:
      m += m2bottom-m1marginstop-m1height-y[macro1] >= 0

  #macro1, macro2 both blockages: do nothing, or continue
 
 #metric included in def of absx, absy
 #xsum is mip lingo for sum of variables
 m.objective = xsum(absx[macro] for macro in range(len(rectangleData)) if not rectangleData[macro]['blockage'])+xsum(absy[macro] for macro in range(len(rectangleData)) if not rectangleData[macro]['blockage'])

 #make sure the algorithm terminates after 30 minutes
 status = m.optimize(max_seconds = MAX_DURATION)
 
 print(status)

 rectangleDataNew = []

 #x = [A, B, C, D]
 #y = [0, 1]
 #zip(x, y) = [(A, 0), (B, 1)]

 #can only fetch if a feasible solution is found
 if status in [OptimizationStatus.OPTIMAL, OptimizationStatus.FEASIBLE]:
  for macro, oldmacro in zip(range(len(rectangleData)), original['rectangleData']):
   #rectangleData[oldmacro] = {'coords': [x[macro].x, y[macro].x, x[macro.x]+macro.width, y[macro].x+macro.height], 'blockage': original['rectangleData'][oldmacro]['blockage'], 'margins': original['rectangleData'][oldmacro]['margins']}
   rectangle = {}

   mleft = Frac(str(rectangleData[macro]['coords'][0]))
   mright = Frac(str(rectangleData[macro]['coords'][2]))
   mtop = Frac(str(rectangleData[macro]['coords'][3]))
   mbottom = Frac(str(rectangleData[macro]['coords'][1]))
   mwidth = mright-mleft
   mheight = mtop-mbottom

   #if blockage, load original coordinates
   if rectangleData[macro]['blockage']:
    rectangle['coords'] = oldmacro['coords']
   else:
    #round to grid
    left_coo = round(x[macro].x/original['grid'][2])*original['grid'][2]
    bottom_coo = round(y[macro].x/original['grid'][3])*original['grid'][3]
    rectangle['coords'] = [left_coo, bottom_coo, left_coo+mwidth, bottom_coo+mheight]

   rectangle['blockage'] = oldmacro['blockage']
   rectangle['margins'] = oldmacro['margins']
   rectangleDataNew.append(rectangle)

 new_instance = {}
 new_instance['spacing_rule_widths'] = original['spacing_rule_widths']
 new_instance['spacing_rule_heights'] = original['spacing_rule_heights']
 new_instance['rectangleData'] = rectangleDataNew #this thing has changed (empty if no feasible solution)
 new_instance['bb'] = original['bb']
 new_instance['grid'] = original['grid']

 return new_instance



#energy function: equals metric
#energy according to rectDataReference
def energy(rectDataNew, rectDataReference, metric):
 sum = 0
 if metric == 'd1':
  for macro, oldmacro in zip(rectDataNew,rectDataReference):
   sum += abs(macro['coords'][0]-oldmacro['coords'][0])+abs(macro['coords'][1]-oldmacro['coords'][1])
 elif metric == 'd1w':
  for macro, oldmacro in zip(rectDataNew,rectDataReference):
   w = oldmacro['coords'][2]-oldmacro['coords'][0]+oldmacro['coords'][3]-oldmacro['coords'][1] #half the perimeter, but as everything is scaled no need to multiply by 2
   sum += w*(abs(macro['coords'][0]-oldmacro['coords'][0])+abs(macro['coords'][1]-oldmacro['coords'][1]))
 elif metric == 'd1w2':
  for macro, oldmacro in zip(rectDataNew,rectDataReference):
   w = oldmacro['coords'][2]-oldmacro['coords'][0]+oldmacro['coords'][3]-oldmacro['coords'][1] #half the perimeter, but as everything is scaled no need to multiply by 2
   sum += w*(abs(macro['coords'][0]-oldmacro['coords'][0])+abs(macro['coords'][1]-oldmacro['coords'][1]))*(abs(macro['coords'][0]-oldmacro['coords'][0])+abs(macro['coords'][1]-oldmacro['coords'][1]))
 elif metric == 'd2w':
  for macro, oldmacro in zip(rectDataNew,rectDataReference):
   w = oldmacro['coords'][2]-oldmacro['coords'][0]+oldmacro['coords'][3]-oldmacro['coords'][1] #half the perimeter, but as everything is scaled no need to multiply by 2
   sum += w*(macro['coords'][0]-oldmacro['coords'][0])*(macro['coords'][0]-oldmacro['coords'][0])+w*(macro['coords'][1]-oldmacro['coords'][1])*(macro['coords'][1]-oldmacro['coords'][1])

 return sum



#brownian_motion algorithm
#works best for one of the weighted metrics
#metric is 'd1', 'd1w', 'd1w2' or 'd2w'
#T is the 'temperature' of the simulated annealing
#under construction: T omitted atm
#in theory, low T gives a better chance of a good solution at the expense of more running time
def simulated_annealing(instance, metric, MAX_DURATION = 1800):
 start_time = time()
 solution = deepcopy(instance)
 original = deepcopy(instance) #needed for the energy function below
 rectangleData = solution['rectangleData']

 #return best found legal solution
 #fill up the rest later, just to be sure
 best_solution = {}

 b_x = Frac(str(original['spacing_rule_widths'][0]))
 c_x = Frac(str(original['spacing_rule_widths'][1]))
 b_y = Frac(str(original['spacing_rule_heights'][0]))
 c_y = Frac(str(original['spacing_rule_heights'][1]))
 grid_x = Frac(str(original['grid'][2]))
 grid_y = Frac(str(original['grid'][3]))
 bb_x = Frac(str(original['bb'][2]))-Frac(str(original['bb'][0]))
 bb_y = Frac(str(original['bb'][3]))-Frac(str(original['bb'][1]))

 scale_x = grid_x/bb_x
 scale_y = grid_y/bb_y

 #quickly verify overlap of blockage in the algorithm, not having to search for them every time
 #quickly access the macros that are not blockages, not having to check that every time
 blockages = [] #indices of macros
 macros = [] #indices of macros
 for macro in range(len(rectangleData)):
  if rectangleData[macro]['blockage']:
   blockages.append(macro)
  else:
   macros.append(macro)

 #first, come up with a legal solution
 #later, improve the solution
 while True: #gets broken at the empty illegal_macros check
  #determine illegally placed macros
  illegal_macros = [] #indices of macros
  for macro in macros:
    if not verify_macro(rectangleData[macro], solution):
     illegal_macros.append(macro)

  #instead of halting, continue to the next phase of the algorithm
  if illegal_macros == []:
   break

  #pick a macro
  W = 0
  w = {}
  for macro in illegal_macros:
   w[macro] = rectangleData[macro]['coords'][3]-rectangleData[macro]['coords'][1]+rectangleData[macro]['coords'][2]-rectangleData[macro]['coords'][0] #again, half the perimeter
   W += 1/w[macro]
  r = random()
  sum = 0
  for macro in illegal_macros:
   sum += 1/(W*w[macro])
   if sum >= r:
    chosen_macro = macro
    break

  #determine directions of illegality
  directions = directions_of_illegality(rectangleData[chosen_macro], solution)
  #print('dir of ill =',directions)

  #define probabilities
  prob = {}
  remaining = 0.8 #prob of doing nothing is 0.2
  all = ['left','down','right','up']
  for direction in directions:
   prob[direction] = 0.05
   remaining -= 0.05
   all.remove(direction)
  if all != []:
   for direction in all:
    prob[direction] = remaining/len(all)
  
  #choose a direction according to prob
  r = random()
  all = ['left','down','right','up']
  c = 0
  chosen_direction = None
  for direction in all:
   c += prob[direction]
   if r < c:
    chosen_direction = direction
    break

  #determine magnitude of movement
  r = random()
  if r < 0.15:
   magnitude_of_movement = c_x-b_x
  elif r < 0.3:
   magnitude_of_movement = c_y-b_y
  elif r < 0.4:
   magnitude_of_movement = grid_x
  elif r < 0.5:
   magnitude_of_movement = grid_y
  else:
   magnitude_of_movement = r*(bb_x+bb_y)*0.01

  #move the cluster (or leave it)
  #print(chosen_direction)
  if chosen_direction == 'left':
   mwidth = original['rectangleData'][chosen_macro]['coords'][2]-original['rectangleData'][chosen_macro]['coords'][0]
   for i in [0,2]:
    rectangleData[chosen_macro]['coords'][i] -= magnitude_of_movement
   rectangleData[chosen_macro]['coords'][0] = round(rectangleData[chosen_macro]['coords'][0]/grid_x)*grid_x #round to grid
   rectangleData[chosen_macro]['coords'][2] = rectangleData[chosen_macro]['coords'][0]+mwidth
  if chosen_direction == 'down':
   mheight = original['rectangleData'][chosen_macro]['coords'][3]-original['rectangleData'][chosen_macro]['coords'][1]
   for i in [1,3]:
    rectangleData[chosen_macro]['coords'][i] -= magnitude_of_movement
   rectangleData[chosen_macro]['coords'][1] = round(rectangleData[chosen_macro]['coords'][1]/grid_y)*grid_y #round to grid
   rectangleData[chosen_macro]['coords'][3] = rectangleData[chosen_macro]['coords'][1]+mheight
  if chosen_direction == 'right':
   mwidth = original['rectangleData'][chosen_macro]['coords'][2]-original['rectangleData'][chosen_macro]['coords'][0]
   for i in [0,2]:
    rectangleData[chosen_macro]['coords'][i] += magnitude_of_movement
   rectangleData[chosen_macro]['coords'][0] = round(rectangleData[chosen_macro]['coords'][0]/grid_x)*grid_x #round to grid
   rectangleData[chosen_macro]['coords'][2] = rectangleData[chosen_macro]['coords'][0]+mwidth
  if chosen_direction == 'up':
   mheight = original['rectangleData'][chosen_macro]['coords'][3]-original['rectangleData'][chosen_macro]['coords'][1]
   for i in [1,3]:
    rectangleData[chosen_macro]['coords'][i] += magnitude_of_movement
   rectangleData[chosen_macro]['coords'][1] = round(rectangleData[chosen_macro]['coords'][1]/grid_y)*grid_y #round to grid
   rectangleData[chosen_macro]['coords'][3] = rectangleData[chosen_macro]['coords'][1]+mheight
  #print('new coo of chosen_macro', chosen_macro, 'is', float(rectangleData[chosen_macro]['coords'][0]), float(rectangleData[chosen_macro]['coords'][1]))
  #if chosen_direction is None, do nothing

  #in case of the new position being outside of the chip or intersecting a blockage
  #instead of reverting, do nothing, trust the probabilities
  #solution['rectangleData'] = rectangleData
  #sleep(0.5)

 print('exited legalizing loop: assume solution is legal, try to improve in the remaining time')
 if not verify(solution, original):
  print('something is wrong: investigate')
  raise AssertionError

 current_time = time()
 time_diff = current_time-start_time
 print('time remaining:', MAX_DURATION-time_diff,'s')
 if time_diff < MAX_DURATION:
  best_solution = sa_legal(original, solution, metric, MAX_DURATION-time_diff)
  return best_solution
 else:
  return solution


#improve solution using simulated annealing
#assumes instance is legal
#original is the original (non-legal) solution, towards which the optimization takes place
#time_left is the remaining time given in seconds, passed on from the legalizer
def sa_legal(original, instance, metric, time_left):
 solution = deepcopy(instance)

 #return best found legal solution
 #fill up the rest later, just to be sure
 best_solution = deepcopy(instance)
 best_energy = energy(instance['rectangleData'], original['rectangleData'], metric)

 b_x = Frac(str(original['spacing_rule_widths'][0]))
 c_x = Frac(str(original['spacing_rule_widths'][1]))
 b_y = Frac(str(original['spacing_rule_heights'][0]))
 c_y = Frac(str(original['spacing_rule_heights'][1]))
 grid_x = Frac(str(original['grid'][2]))
 grid_y = Frac(str(original['grid'][3]))
 bb_x = Frac(str(original['bb'][2]-original['bb'][0]))
 bb_y = Frac(str(original['bb'][3]-original['bb'][1]))

 scale_x = grid_x/bb_x
 scale_y = grid_y/bb_y

 #determine initial temp
 #exp(-typical_change/temp) should be ~1-scale in the beginning
 #so, typical_change/temp should be ~scale
 scale = scale_x+scale_y
 typical_change = 0
 if metric == 'd1':
  typical_change = grid_x+grid_y
 elif metric == 'd1w':
  typical_change = (grid_x+grid_y)*(bb_x+bb_y)/len(original['rectangleData'])
 elif metric == 'd1w2':
  typical_change = (grid_x+grid_y)*(bb_x+bb_y)/len(original['rectangleData'])/10 #unknown distance, so assume low
 elif metric == 'd2w':
  typical_change = (grid_x+grid_y)*(bb_x+bb_y)/len(original['rectangleData'])/10 #again, unknown distance so assume low
 temp = typical_change/scale

 current_temp = temp
 while True:
  start_time = time()

  #perturb solution
  macro = choice(range(len(solution['rectangleData']))) #index of random macro
  rectDataOriginal = original['rectangleData'][macro] #for the energy function, which in this case only depends on the moved macro
  rectDataOld = deepcopy(solution['rectangleData'][macro]) #be able to revert if necessary
  bc_gap = {'x':c_x-b_x, 'y': c_y-b_y}
  r = random()
  if r < 1/5:
   #move left
   solution['rectangleData'][macro]['coords'][0] -= grid_x
   solution['rectangleData'][macro]['coords'][2] -= grid_x
  elif r < 2/5:
   #move down
   solution['rectangleData'][macro]['coords'][1] -= grid_y
   solution['rectangleData'][macro]['coords'][3] -= grid_y
  elif r < 3/5:
   #move right
   solution['rectangleData'][macro]['coords'][0] += grid_x
   solution['rectangleData'][macro]['coords'][2] += grid_x
  elif r < 4/5:
   #move up
   solution['rectangleData'][macro]['coords'][1] += grid_y
   solution['rectangleData'][macro]['coords'][3] += grid_y

  #sometimes jump between b and c alignment instead of grid
  else:
   for othermacro in solution['rectangleData']:

    #macro c_x-spaced right of othermacro: move to b_x
    if c_x-EPSILON <= solution['rectangleData'][macro]['coords'][0]-othermacro['coords'][2] <= c_x+EPSILON:
     solution['rectangleData'][macro]['coords'][0] -= bc_gap['x']
     solution['rectangleData'][macro]['coords'][2] -= bc_gap['x']
     break
    #macro b_x-spaced right of othermacro: move to c_x
    elif b_x-EPSILON <= solution['rectangleData'][macro]['coords'][0]-othermacro['coords'][2] <= b_x+EPSILON:
     solution['rectangleData'][macro]['coords'][0] += bc_gap['x']
     solution['rectangleData'][macro]['coords'][2] += bc_gap['x']
     break
    #macro b_x-spaced left of othermacro: move to c_x
    elif b_x-EPSILON <= -solution['rectangleData'][macro]['coords'][0]+othermacro['coords'][2] <= b_x+EPSILON:
     solution['rectangleData'][macro]['coords'][0] -= bc_gap['x']
     solution['rectangleData'][macro]['coords'][2] -= bc_gap['x']
     break
    #macro c_x-spaced left of othermacro: move to b_x
    elif c_x-EPSILON <= -solution['rectangleData'][macro]['coords'][0]+othermacro['coords'][2] <= c_x+EPSILON:
     solution['rectangleData'][macro]['coords'][0] += bc_gap['x']
     solution['rectangleData'][macro]['coords'][2] += bc_gap['x']
     break
    #macro c_y-spaced above othermacro: move to b_y
    if c_y-EPSILON <= solution['rectangleData'][macro]['coords'][1]-othermacro['coords'][3] <= c_y+EPSILON:
     solution['rectangleData'][macro]['coords'][1] -= bc_gap['y']
     solution['rectangleData'][macro]['coords'][3] -= bc_gap['y']
     break
    #macro b_y-spaced above othermacro: move to c_y
    elif b_y-EPSILON <= solution['rectangleData'][macro]['coords'][1]-othermacro['coords'][3] <= b_y+EPSILON:
     solution['rectangleData'][macro]['coords'][1] += bc_gap['y']
     solution['rectangleData'][macro]['coords'][3] += bc_gap['y']
     break
    #macro b_y-spaced below othermacro: move to c_y
    elif b_y-EPSILON <= -solution['rectangleData'][macro]['coords'][1]+othermacro['coords'][3] <= b_y+EPSILON:
     solution['rectangleData'][macro]['coords'][1] -= bc_gap['y']
     solution['rectangleData'][macro]['coords'][3] -= bc_gap['y']
     break
    #macro c_y-spaced below othermacro: move to b_y
    elif c_y-EPSILON <= -solution['rectangleData'][macro]['coords'][1]+othermacro['coords'][3] <= c_y+EPSILON:
     solution['rectangleData'][macro]['coords'][1] += bc_gap['y']
     solution['rectangleData'][macro]['coords'][3] += bc_gap['y']
     break

  #check energy and legality
  #if not legal, revert
  if not verify_macro(solution['rectangleData'][macro], original):
   solution['rectangleData'][macro] = rectDataOld
  else:
   #if energy is worse, revert with prob 1-exp(-energy/temp)
   energyOld = energy([rectDataOld], [rectDataOriginal], metric)
   energyNew = energy([solution['rectangleData'][macro]], [rectDataOriginal], metric)
   if energyNew-energyOld > 0:
    r = random()
    if r >= exp((energyOld-energyNew)/current_temp):
     solution['rectangleData'][macro] = rectDataOld

  #if found new best solution: update
  energyTotal = energy(solution['rectangleData'], original['rectangleData'], metric)
  if energyTotal < best_energy:
   best_solution = deepcopy(solution)
   best_energy = energyTotal

  #decrease temp
  #multiply temp by alpha, with alpha depending on time
  current_time = time()
  time_factor = (current_time-start_time)/time_left
  #alpha^number of steps = TERMINATE_CONSTANT
  #number of steps ~ 1/time_factor
  TERMINATE_CONSTANT = 1/1000
  alpha = TERMINATE_CONSTANT**time_factor
  current_temp *= alpha

  #restart condition
  #if energy is worse than a number of times the typical energy change
  FACTOR = 20
  if energyTotal > best_energy+typical_change*FACTOR:
   solution = deepcopy(best_solution) #revert to best solution, continue from there

  #terminate condition
  if current_temp < temp*TERMINATE_CONSTANT:
   return best_solution

'''
#brownian_motion algorithm
#works best for one of the weighted metrics
#metric is 'd1', 'd1w', 'd1w2' or 'd2w'
def graph_prepare(instance):
 solution = deepcopy(instance)
 original = deepcopy(instance) #needed for the energy function below
 rectangleData = solution['rectangleData']

 #cut ties to solution['rectangleData'] as want to sync only afterwards
 #additionally, convert to dict for easy labeling for the clusters
 rectangleData = shuffle(rectangleData)

 b_x = original['spacing_rule_widths'][0]
 c_x = original['spacing_rule_widths'][1]
 b_y = original['spacing_rule_heights'][0]
 c_y = original['spacing_rule_heights'][1]
 grid_x = Frac(str(original['grid'][2]))
 grid_y = Frac(str(original['grid'][3]))

 #quickly verify overlap of blockage in the algorithm, not having to search for them every time
 #quickly access the macros that are not blockages, not having to check that every time
 blockages = []
 macros = []
 for macro in rectangleData:
  if rectangleData[macro]['blockage']:
   blockages.append(macro)
  else:
   macros.append(macro)
   G.add_node(macro)

 #first, come up with a legal solution
 #later, improve the solution
 while True: #gets broken at the empty illegal_clusters check
  #make Gh, Gv
  Gh = nx.DiGraph()
  Gv = nx.DiGraph()
  for macro in macros:
   Gh.add_edge('s',macro)
   Gh.add_edge(macro,'t')
   Gv.add_edge('s',macro)
   Gv.add_edge(macro,'t')
  for macro1, macro2 in combinations(macros, 2):
   if rectangleData[macro1][0] <= rectangleData[macro2][0]:
    freespace_horizontal = rectangleData[macro2][0]-rectangleData[macro1][2]
   else:
    freespace_vertical = rectangleData[macro2][0]-rectangleData[macro1][2]
   if freespace_horizontal > freespace_vertical:
    Gh.add_edge(macro1, macro2)
   else:
    Gv.add_edge(macro1, macro2)

  #add weights
  for macro in macros:
   Gh['s'][macro]['weight'] = 0
   Gv['s'][macro]['weight'] = 0
  for macro1, macro2 in combinations(macros, 2):
   if (macro1, macro2) in Gh.edges:
    #width of source of edge+closest_spacing (i.e. spacing b_x or c_x, the one requiring minimal movement)
   if (macro1, macro2) in Gv.edges:
    #height of source of edge+closest_spacing (note: b_y and c_y for this one)
  for macro in macros:
   Gh[macro]['t']['weight'] = #width of macro
   Gv[macro]['t']['weight'] = #height of macro

  #determine longest path in both graphs
  #longest path, so weights negative
  #sadly, NetworkX doesn't allow to return both the path and it's length in one go
  dijkstra_h = nx.dijkstra_path(Gh, 's','t','weight')
  dijkstra_h_length = nx.shortest_path_length(Gh,'s','t','weight')
  dijkstra_v = nx.dijkstra_path(Gv, 's','t','weight')
  dijkstra_v_length = nx.shortest_path_length(Gv,'s','t','weight')

  #if longest path is too long, adjust, mostly by swapping between the graphs

 #TODO: time remaining
 if not verify(solution, original):
  print('something is wrong: not legal yet')

 return solution

#move each macro to the closest current legal position, according to metric
#macro is moved to legal position, w.r.t. the macros and blockages it intersects with
#if there is no legal position, randomly move and restart
def greedy2(instance, metric):
 original = deepcopy(instance)
 instance = deepcopy(instance) #don't change the input instance

 legalCount = 0
 legalCountMax = len(original['rectangleData'])*(1+100/sqrt(len(original['rectangleData']))
 while True:
  #pick a macro
  macro = choice(instance['rectangleData']) #here, macro is the actual macro, not it's index

  if not macro['blockage']:
   #check if it's legally placed
   legal = verify_macro(macro, original)

   #keep track of number of consecutive encounters of legally placed macros
   #if over a certain number, which depends on the number of macros, check legality of the entire thing
   if legal:
    legalCount += 1
    if legalCount > legalCountMax:
     if verify(instance):
      return instance
     else:
      legalCount = 0 #keep checking

   #if not, determine nearby legal positions
   else:
    vbb = violates_bounding_box(macro)
    if vbb is not None:
     #move back in, directly, accounting for margins
     if vbb == 'left':
      mwidth = Frac(str(macro['coords'][2]))-Frac(str(macro['coords'][0]))
      macro['coords'][0] = Frac(str(original['bb'][0]))+Frac(str(macro['margins'][0]))
      macro['coords'][2] = Frac(str(original['bb'][0]))+Frac(str(macro['margins'][0]))+mwidth
     elif vbb == 'right':
      mwidth = Frac(str(macro['coords'][2]))-Frac(str(macro['coords'][0]))
      macro['coords'][0] = Frac(str(original['bb'][2]))-Frac(str(macro['margins'][2]))
      macro['coords'][2] = Frac(str(original['bb'][2]))-Frac(str(macro['margins'][2]))+mwidth
     elif vbb == 'bottom':
      mheight = Frac(str(macro['coords'][3]))-Frac(str(macro['coords'][1]))
      macro['coords'][1] = Frac(str(original['bb'][1]))+Frac(str(macro['margins'][1]))
      macro['coords'][3] = Frac(str(original['bb'][1]))+Frac(str(macro['margins'][1]))+mheight
     elif vbb == 'top':
      mheight = Frac(str(macro['coords'][3]))-Frac(str(macro['coords'][1]))
      macro['coords'][1] = Frac(str(original['bb'][3]))-Frac(str(macro['margins'][3]))
      macro['coords'][3] = Frac(str(original['bb'][3]))-Frac(str(macro['margins'][3]))+mheight


    othermacro = choice(instance['rectangleData']

   #determine which ones are the 'cheapest'

   #optionally, build in some randomness according to SA

   #move it

'''
