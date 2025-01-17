import sys
from typing import List
try:
    from shapely.geometry import box
    from shapely.ops import nearest_points
except Exception:
    print("Warning: could not import shapely\n")
    pass
import json
import copy

from verify import Macro

from fractions import Fraction as Frac


class Problem():
    def __init__(self, macros: List[Macro], bb, spacing_rule_widths, spacing_rule_heights, grid):
        self.macros = macros
        self.bb = bb
        self.spacing_rule_widths = spacing_rule_widths
        self.spacing_rule_heights = spacing_rule_heights
        self.grid = grid

    @property
    def width(self):
        return self.bb[2]

    @property
    def height(self):
        return self.bb[3]

    @property
    def xspacing_b(self):
        return self.spacing_rule_widths[0]

    @property
    def xspacing_c(self):
        return self.spacing_rule_widths[1]

    @property
    def yspacing_b(self):
        return self.spacing_rule_heights[0]

    @property
    def yspacing_c(self):
        return self.spacing_rule_heights[1]

    @classmethod
    def load(cls, instance: dict):
        return Problem([Macro(macro) for macro in instance['rectangleData']], instance["bb"], instance['spacing_rule_widths'], instance['spacing_rule_heights'], instance["grid"])

    def to_json(self):
        return {"rectangleData": list(map(Macro.to_json, self.macros)), "bb": self.bb, "spacing_rule_widths": self.spacing_rule_widths, "spacing_rule_heights": self.spacing_rule_heights, "grid": self.grid}


def display_instance_information(instance):
    print("Instance information:")
    print(f"Spacing rule widths: {instance['spacing_rule_widths']}")
    print(f"Spacing rule heights: {instance['spacing_rule_heights']}")
    print(f"Number of rectangles: {len(instance['rectangleData'])}")
    print(f"Bounding Box: {instance['bb']}")
    print(f"Grid: {instance['grid']}")
    print('')


def load_instance(filename):
    with open(filename) as input_file:
        return json.load(input_file)

def write_instance(instance: Problem, filename):
    json.dump(instance.to_json(), open(filename, 'w'))

def write_instance_json(instance, filename):
    json.dump(instance, open(filename, 'w'))


def create_instance(spacing_rule_widths, spacing_rule_heights, rectangleData, bb, grid, do_sanity_checks=True):
    if do_sanity_checks:
        # spacing_rule_widths
        assert(len(spacing_rule_widths) == 2)
        assert(0 <= spacing_rule_widths[0] <= spacing_rule_widths[1])
        # spacing_rule_heights
        assert(len(spacing_rule_heights) == 2)
        assert(0 <= spacing_rule_heights[0] <= spacing_rule_heights[1])
        # rectangleData
        for rectangle in rectangleData:
            assert("coords" in rectangle)
            assert(len(rectangle["coords"]) == 4)
            assert("blockage" in rectangle)
            assert(type(rectangle["blockage"]) == bool)
            assert("margins" in rectangle)
            assert(len(rectangle["margins"]) == 4)
        # bb
        assert(len(bb) == 4)
        assert(bb[2] > bb[0])
        assert(bb[3] > bb[1])
        # grid
        assert(len(grid) == 4)

    instance = {}
    instance['spacing_rule_widths'] = spacing_rule_widths
    instance['spacing_rule_heights'] = spacing_rule_heights
    instance['rectangle_Data'] = rectangleData
    instance['bb'] = bb
    instance['grid'] = grid

    return instance


def preprocess(instance):
    for item in instance["rectangleData"]:
        min_x, min_y, max_x, max_y = item['coords']
        item['rect'] = box(min_x, min_y, max_x, max_y)


#converts an instance to integer grid coordinates, to be able to use the mip solver
def to_int(instance):
    instance = copy.deepcopy(instance)

    def convertx(v):
        return v / instance['grid'][2]

    def converty(v):
        return v / instance['grid'][3]
    for item in instance["rectangleData"]:
        item['coords'] = [convertx(c) if i % 2 == 0 else converty(
            c) for i, c in enumerate(item['coords'])]
        item['margins'] = [convertx(c) if i % 2 == 0 else converty(
            c) for i, c in enumerate(item['margins'])]
    instance['bb'] = [convertx(c) if i % 2 == 0 else converty(c)
                      for i, c in enumerate(instance['bb'])]
    instance['spacing_rule_widths'] = [
        convertx(c) for c in instance['spacing_rule_widths']]
    instance['spacing_rule_heights'] = [
        converty(c) for c in instance['spacing_rule_heights']]
    return instance

#inverse of above function
def to_int_inv(instance):
    instance = copy.deepcopy(instance)

    def convertx(v):
        return v * instance['grid'][2]

    def converty(v):
        return v * instance['grid'][3]
    for item in instance["rectangleData"]:
        item['coords'] = [convertx(c) if i % 2 == 0 else converty(
            c) for i, c in enumerate(item['coords'])]
        item['margins'] = [convertx(c) if i % 2 == 0 else converty(
            c) for i, c in enumerate(item['margins'])]
    instance['bb'] = [convertx(c) if i % 2 == 0 else converty(c)
                      for i, c in enumerate(instance['bb'])]
    instance['spacing_rule_widths'] = [
        convertx(c) for c in instance['spacing_rule_widths']]
    instance['spacing_rule_heights'] = [
        converty(c) for c in instance['spacing_rule_heights']]
    return instance



def distance_vector(a_obj, b_obj):
    """Get the vector of the horizonatal and vertical distance between rectangles a_obj and b_obj"""
    #point_a, point_b = nearest_points(a, b)
    # return [abs(point_a.x - point_b.x), abs(point_a.y - point_b.y)]
    a = a_obj['coords']
    b = b_obj['coords']
    return [abs(max(a[0], b[0]) - min(a[2], b[2])), abs(max(a[1], b[1]) - min(a[3], b[3]))]


if __name__ == "__main__":
    instance = load_instance(sys.argv[1])
    display_instance_information(instance)
