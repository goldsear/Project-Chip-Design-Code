from re import A
import sys
import json
import itertools

EPSILON = 10**-6 # Error margin (distance units) to account for floating point errors

class Rectangle(list):
    @property
    def left(self):
        return self[0]
    @property
    def bottom(self):
        return self[1]
    @property
    def right(self):
        return self[2]
    @property
    def top(self):
        return self[3]

    @property
    def width(self):
      assert(self.right - self.left >= 0)
      return self.right - self.left

    @property
    def height(self):
      assert(self.top - self.bottom >= 0)
      return self.top - self.bottom

    @property
    def area(self):
      return self.width * self.height

    @property
    def center_x(self):
      return (self.right + self.left / 2)

    @property
    def center_y(self):
      return (self.bottom + self.top / 2)

    def distance_to(self, other: "Rectangle"):
        if self.is_before_x(other):
            rx = other.left - self.right
        else:
            rx = self.left - other.right
        if self.is_before_y(other):
            ry = other.bottom - self.top
        else:
            ry = self.bottom - other.top
        return rx, ry
        #return (abs(max(other.left, self.left) - min(other.right, self.right)), abs(max(other.bottom, self.bottom) - min(other.top, self.top)))
    def is_before_x(self, other: "Rectangle"):
        return other.center_x - self.center_x >= 0
    def is_before_y(self, other: "Rectangle"):
        return other.center_y - self.center_y >= 0
    def intersects_with(self, other: "Rectangle"):
        if self.right <= other.left:
            return False
        elif self.left >= other.right:
            return False
        elif self.bottom >= other.top:
            return False
        elif self.top <= other.bottom:
            return False

        return True

    def move_to_x(self, x):
        width = self.width
        self[0] = x
        self[2] = x + width

    def move_to_y(self, y):
        height = self.height
        assert(self.bottom != self.top)
        self[1] = y
        self[3] = y + height
        assert(self.bottom != self.top)

    def intersection(self, other):
        if not self.intersects_with(other):
            return Rectangle([0, 0, 0, 0])
        return Rectangle([max(self[0], other[0]), max(self[1], other[1]), min(self[2], other[2]), min(self[3], other[3])])

    def add(self, x, y):
        return Rectangle([self[0] + x, self[1] + y, self[2] + x, self[3] + y])

    def expand(self, x, y):
        return Rectangle([self[0] - x, self[1] - y, self[2] + x, self[3] + y])



class Macro(object):
    def __init__(self,macro_dict):
        super().__init__()
        self.target = Rectangle(macro_dict['coords'])
        if "target" in macro_dict:
            self.target = Rectangle(macro_dict['target'])
        self.coords = Rectangle(macro_dict['coords'])

        self.blockage = macro_dict['blockage']
        self.margins = Rectangle(macro_dict['margins'])

        assert(self.coords.left <= self.coords.right)
        assert(self.coords.bottom <= self.coords.top)
        assert(not self.blockage or all(margin == 0 for margin in self.margins))

    @property
    def left(self):
        return self.coords.left

    @property
    def bottom(self):
        return self.coords.bottom

    @property
    def right(self):
        return self.coords.right

    @property
    def top(self):
        return self.coords.top

    @property
    def width(self):
      assert(self.right - self.left >= 0)
      return self.right - self.left

    @property
    def height(self):
      assert(self.top - self.bottom >= 0)
      return self.top - self.bottom

    @property
    def perimeter(self):
        return 2*(self.coords.top-self.coords.bottom) + 2*(self.coords.right-self.coords.left)

    def to_json(self):
        return {"coords": self.coords, "blockage": self.blockage, "margins": self.margins, "target": self.target}
    def __str__(self):
        return f"Macro(coords = {self.coords}, blockage = {self.blockage}, margins = {self.margins})"


def display_instance_information(instance):
    print("Instance information:")
    print(f"Spacing rule widths: {instance['spacing_rule_widths']}")
    print(f"Spacing rule heights: {instance['spacing_rule_heights']}")
    print(f"Number of macros: {len(instance['rectangleData'])}")
    print(f"Bounding Box: {instance['bb']}")
    print(f"Grid: {instance['grid']}")
    print('')

def load_instance(filename):
    with open(filename) as input_file:
        return json.load(input_file)

def verify(solution, input_instance = None):
    macros = [Macro(macro) for macro in solution['rectangleData']]

    if input_instance is None:
        print("Warning: cannot compare solution to input instance, because the input instance was not given.\n")
    else:
        if solution['spacing_rule_widths'] != input_instance['spacing_rule_widths']:
            print("Spacing rule widths of solution and input are unequal\n")
            return False
        if solution['spacing_rule_heights'] != input_instance['spacing_rule_heights']:
            print("Spacing rule heights of solution and input are unequal\n")
            return False
        if solution['bb'] != input_instance['bb']:
            print("Bounding box of solution and input are unequal\n")
            return False
        if solution['grid'] != input_instance['grid']:
            print("Grid of solution and input are unequal\n")
            return False
        
        input_macros = [Macro(macro) for macro in input_instance['rectangleData']]
        if len(input_macros) != len(macros):
            print(f"Solutions has a different amount of macros ({len(macros)}) from the input ({len(input_macros)})\n")
            return False
        for macro,input_macro in zip(macros,input_macros):
            if macro.margins != input_macro.margins:
                print(f"Macro margins of the following two macros are unequal:\n{macro} (solution)\n{input_macro} (input)\n")
                return False
            if macro.blockage != input_macro.blockage:
                print(f"Macro blockage status of the following two macros are unequal:\n{macro} (solution)\n{input_macro} (input)\n")
                return False
            if abs(macro.coords.top - macro.coords.bottom-(input_macro.coords.top - input_macro.coords.bottom)) > EPSILON or abs(macro.coords.right - macro.coords.left-(input_macro.coords.right - input_macro.coords.left)) > EPSILON:
                print(f"The following macros are not the same size:\n{macro} (solution)\n{input_macro} (input)\n")
                return False
            if macro.blockage and (abs(macro.bottom - input_macro.bottom) > EPSILON or abs(macro.left - input_macro.left) > EPSILON):
                print(f"The following blockage macro has moved:\n{macro} (solution)\n{input_macro} (input)\n")
                return False
    
    b_x,c_x = solution['spacing_rule_widths']
    b_y,c_y = solution['spacing_rule_heights']
    grid_offset_x,grid_offset_y,grid_spacing_x,grid_spacing_y = solution['grid']

    bounding_box = Rectangle(solution['bb'])

    def violates_bounding_box(macro):
        violation = False
        if macro.coords.left - macro.margins.left < bounding_box.left - EPSILON:
            violation = True
        elif macro.coords.right + macro.margins.right > bounding_box.right + EPSILON:
            violation = True
        elif macro.coords.bottom - macro.margins.bottom < bounding_box.bottom - EPSILON:
            violation = True
        elif macro.coords.top + macro.margins.top > bounding_box.top + EPSILON:
            violation = True

        return violation
    
    def satisfies_grid(macro):
        x_units = (macro.coords.left - grid_offset_x)/grid_spacing_x #int version of the x coo according to grid
        if abs(x_units - round(x_units)) * grid_spacing_x > EPSILON: #x coo needs to be EPSILON close to grid point in original coordinates
            return False
        
        y_units = (macro.coords.bottom - grid_offset_y)/grid_spacing_y #same for y
        if abs(y_units - round(y_units)) * grid_spacing_y > EPSILON:
            return False
        
        return True

    def valid_spacing(macro_1,macro_2):
        # EPSILON close to b or at least c (up to EPSILON) far apart
        if macro_1.coords.left - macro_2.right >= c_x - EPSILON or abs(macro_1.coords.left - macro_2.right - b_x) < EPSILON:
            return True
        if macro_2.coords.left - macro_1.right >= c_x - EPSILON or abs(macro_2.coords.left - macro_1.right - b_x) < EPSILON:
            return True
        if macro_1.bottom - macro_2.coords.top >= c_y - EPSILON or abs(macro_1.bottom - macro_2.coords.top - b_y) < EPSILON:
            return True
        if macro_2.bottom - macro_1.coords.top >= c_y - EPSILON or abs(macro_2.bottom - macro_1.coords.top - b_y) < EPSILON:
            return True
        
        return False
    
    def valid_margins(macro_1,macro_2):
        #no overlap of margins, so need to sum the margins
        if macro_1.coords.left   - macro_2.right >= macro_1.margins.left   + macro_2.margins.right - EPSILON:
            return True
        if macro_2.coords.left   - macro_1.right >= macro_2.margins.left   + macro_1.margins.right - EPSILON:
            return True
        if macro_1.bottom - macro_2.coords.top   >= macro_1.margins.bottom + macro_2.margins.top   - EPSILON:
            return True
        if macro_2.bottom - macro_1.coords.top   >= macro_2.margins.bottom + macro_1.margins.top   - EPSILON:
            return True

        return False

    for macro in macros:
        if violates_bounding_box(macro):
            print(f"The following macro violates the bounding box:\n{macro}\n")
            return False
        
        if not satisfies_grid(macro) and not macro.blockage:
            print(f"The following macro does not satisfy the grid:\n{macro}\n")
            return False

    for macro_1,macro_2 in itertools.combinations(macros,2):
        if not valid_spacing(macro_1,macro_2):
            print(f"The following macros violate the spacing rule or intersect:\n{macro_1}\n{macro_2}\n")
            return False
        if not valid_margins(macro_1,macro_2):
            print(f"The following macros violate a keep-out zone:\n{macro_1}\n{macro_2}\n")
            return False
    
    print("Solution is correct!")
    return True

def directions_of_illegality(macro, solution):
 macro = Macro(macro) #assume macro in dict format
 directions = set()
 for othermacro in solution['rectangleData']:
  othermacro = Macro(othermacro)

  #overlap with other macro or blockage, including margins
  #othermacro left of macro
  if othermacro.left < macro.left:
   if othermacro.right+othermacro.margins.right > macro.left-macro.margins.left:
    if othermacro.bottom-othermacro.margins.bottom < macro.top+macro.margins.top:
     if othermacro.top+othermacro.margins.top > macro.bottom-macro.margins.bottom:
      directions.add('left')

  #othermacro below macro
  if othermacro.bottom < macro.bottom:
   if othermacro.top+othermacro.margins.top > macro.bottom-macro.margins.bottom:
    if othermacro.left-othermacro.margins.left < macro.right+macro.margins.right:
     if othermacro.right+othermacro.margins.right > macro.left-macro.margins.left:
      directions.add('down')

  #othermacro right of macro
  if othermacro.left > macro.left:
   if othermacro.left-othermacro.margins.left < macro.right+macro.margins.right:
    if othermacro.bottom-othermacro.margins.bottom < macro.top+macro.margins.top:
     if othermacro.top+othermacro.margins.top > macro.bottom-macro.margins.bottom:
      directions.add('right')

  #othermacro above macro
  if othermacro.bottom > macro.bottom:
   if othermacro.bottom-othermacro.margins.bottom < macro.top+macro.margins.top:
    if othermacro.left-othermacro.margins.left < macro.right+macro.margins.right:
     if othermacro.right+othermacro.margins.right > macro.left-macro.margins.left:
      directions.add('up')

  #outside of bounding box, including margins
  if macro.left-macro.margins.left < solution['bb'][0]:
   directions.add('left')
  if macro.bottom-macro.margins.bottom < solution['bb'][1]:
   directions.add('down')
  if macro.right+macro.margins.right > solution['bb'][2]:
   directions.add('right')
  if macro.top+macro.margins.top > solution['bb'][3]:
   directions.add('up')
 return directions


#verifies if a single macro verifies the placement rules w.r.t. the bounding box and the other macros and blockages
def verify_macro(macro, solution):
    macro = Macro(macro) #assume macro in dict format
    b_x,c_x = solution['spacing_rule_widths']
    b_y,c_y = solution['spacing_rule_heights']
    grid_offset_x,grid_offset_y,grid_spacing_x,grid_spacing_y = solution['grid']

    bounding_box = Rectangle(solution['bb'])

    def violates_bounding_box(macro):
        violation = False
        if macro.coords.left - macro.margins.left < bounding_box.left - EPSILON:
            violation = True
        elif macro.coords.right + macro.margins.right > bounding_box.right + EPSILON:
            violation = True
        elif macro.coords.bottom - macro.margins.bottom < bounding_box.bottom - EPSILON:
            violation = True
        elif macro.coords.top + macro.margins.top > bounding_box.top + EPSILON:
            violation = True

        return violation
    
    def satisfies_grid(macro):
        x_units = (macro.coords.left - grid_offset_x)/grid_spacing_x #int version of the x coo according to grid
        if abs(x_units - round(x_units)) * grid_spacing_x > EPSILON: #x coo needs to be EPSILON close to grid point in original coordinates
            return False
        
        y_units = (macro.coords.bottom - grid_offset_y)/grid_spacing_y #same for y
        if abs(y_units - round(y_units)) * grid_spacing_y > EPSILON:
            return False
        
        return True

    def valid_spacing(macro_1,macro_2):
        # EPSILON close to b or at least c (up to EPSILON) far apart
        if macro_1.coords.left - macro_2.right >= c_x - EPSILON or abs(macro_1.coords.left - macro_2.right - b_x) < EPSILON:
            return True
        if macro_2.coords.left - macro_1.right >= c_x - EPSILON or abs(macro_2.coords.left - macro_1.right - b_x) < EPSILON:
            return True
        if macro_1.bottom - macro_2.coords.top >= c_y - EPSILON or abs(macro_1.bottom - macro_2.coords.top - b_y) < EPSILON:
            return True
        if macro_2.bottom - macro_1.coords.top >= c_y - EPSILON or abs(macro_2.bottom - macro_1.coords.top - b_y) < EPSILON:
            return True
        
        return False
    
    def valid_margins(macro_1,macro_2):
        #no overlap of margins, so need to sum the margins
        if macro_1.coords.left   - macro_2.right >= macro_1.margins.left   + macro_2.margins.right - EPSILON:
            return True
        if macro_2.coords.left   - macro_1.right >= macro_2.margins.left   + macro_1.margins.right - EPSILON:
            return True
        if macro_1.bottom - macro_2.coords.top   >= macro_1.margins.bottom + macro_2.margins.top   - EPSILON:
            return True
        if macro_2.bottom - macro_1.coords.top   >= macro_2.margins.bottom + macro_1.margins.top   - EPSILON:
            return True

        return False

    if violates_bounding_box(macro):
        return False

    if not satisfies_grid(macro):
        return False

    for othermacro in solution['rectangleData']:
        othermacro = Macro(othermacro)
        if othermacro.coords != macro.coords:
            if not valid_spacing(macro, othermacro):
                return False
            if not valid_margins(macro, othermacro):
                return False

    return True


def score_vectors(solution,input_instance):
    sum_of_absolute_values = []
    sum_of_absolute_values_weighted = []
    sum_of_absolute_values_squared_weighted = []
    sum_of_squares_weighted = []

    solution_macros = [Macro(macro) for macro in solution['rectangleData']]
    input_macros = [Macro(macro) for macro in input_instance['rectangleData']]

    for solution_macro,input_macro in zip(solution_macros,input_macros):
        weight = solution_macro.perimeter
        dx = abs(solution_macro.coords.left - input_macro.coords.left)
        dy = abs(solution_macro.coords.bottom - input_macro.coords.bottom)
        
        sum_of_absolute_values.append(dx + dy)
        sum_of_absolute_values_weighted.append(weight*(dx + dy))
        sum_of_absolute_values_squared_weighted.append(weight*((dx + dy)**2))
        sum_of_squares_weighted.append(weight*(dx**2 + dy**2))
    
    score_vectors = (sum_of_absolute_values,sum_of_absolute_values_weighted,sum_of_absolute_values_squared_weighted,sum_of_squares_weighted)

    return score_vectors

def score(solution,input_instance,verbose = False):
    sum_of_absolute_values = 0
    sum_of_absolute_values_weighted = 0
    sum_of_absolute_values_squared_weighted = 0
    sum_of_squares_weighted = 0

    solution_macros = [Macro(macro) for macro in solution['rectangleData']]
    input_macros = [Macro(macro) for macro in input_instance['rectangleData']]

    for solution_macro,input_macro in zip(solution_macros,input_macros):
        weight = solution_macro.perimeter
        dx = abs(solution_macro.coords.left - input_macro.coords.left)
        dy = abs(solution_macro.coords.bottom - input_macro.coords.bottom)
        
        sum_of_absolute_values += (dx + dy)
        sum_of_absolute_values_weighted += weight*(dx + dy)
        sum_of_absolute_values_squared_weighted += weight*((dx + dy)**2)
        sum_of_squares_weighted += weight*(dx**2 + dy**2)
    
    if verbose:
        print(f"Sum of absolute values: {sum_of_absolute_values}\nWeighted sum of absolute values: {sum_of_absolute_values_weighted}\nWeighted sum of (absolute values) squared: {sum_of_absolute_values_squared_weighted}\nWeighted sum of squares: {sum_of_squares_weighted}")

    return (sum_of_absolute_values,sum_of_absolute_values_weighted,sum_of_absolute_values_squared_weighted,sum_of_squares_weighted)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print('Please give the solution to verify as a command line argument.\n Optionally, give the reference instance as a second argument.\n')
        exit()

    solution = load_instance(sys.argv[1])
    input_instance = None
    
    if len(sys.argv) >= 3:
        input_instance = load_instance(sys.argv[2])
    
    verify(solution,input_instance)

    if input_instance is None:
        print("Cannot score solution, because the input instance was not given as the second command line argument.\n")
    else:
        score(solution,input_instance,verbose=True)
