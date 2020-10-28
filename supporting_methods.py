import numpy as np
from copy import deepcopy

def rectangular_to_polar(complex_number):
    r = np.sqrt(complex_number.real*complex_number.real+complex_number.imag*complex_number.imag).real
    if complex_number.real < 0 and complex_number.imag > 0: # second quadrant
        angle = np.pi - np.arctan(abs(complex_number.imag / complex_number.real))
    elif complex_number.real < 0 and complex_number.imag < 0: # third quadrant
        angle = -(np.pi-np.arctan(abs(complex_number.imag / complex_number.real)))
    elif complex_number.real > 0 and complex_number.imag < 0: # fourth quadrant
        angle = -np.arctan(abs(complex_number.imag / complex_number.real))
    else: # first quadrant
        angle = np.arctan(abs(complex_number.imag / complex_number.real))
    return r, angle

def polar_to_rectangular(abs, angle):
    a = abs*np.cos(angle)
    b = abs*np.sin(angle)
    return complex(a, b)

def get_from_and_to_bus(input_string):
    """
    input string should be on the format "x-y"
    where x is the from bus and y is the to bus
    from and to bus are outputed as integers
    """
    from_bus, to_bus = input_string.split("-")
    return int(from_bus), int(to_bus)

def print_title1(input_string=None):
    """
    Prints a main title generally for a new task
    """
    symbol_count = 120
    symbol = "#"
    if input_string:
        print("\n", int(symbol_count / 2) * symbol, input_string, int(symbol_count / 2) * symbol)
    else:
        print("\n", symbol_count * symbol)

def print_title2(input_string=None):
    """
    Prints a undertitle generally for a new problem
    """
    symbol_count = 60
    symbol = "*"
    if input_string:
        print("\n", int(symbol_count / 2) * symbol, input_string, int(symbol_count / 2) * symbol)
    else:
        print("\n", symbol_count * symbol)

def print_title3(input_string=None):
    """
    Prints a sub undertitle generally for a new iteration
    """
    symbol_count = 50
    symbol = "-"
    if input_string:
        print("\n", int(symbol_count / 2) * symbol, input_string, int(symbol_count / 2) * symbol)
    else:
        print("\n", symbol_count * symbol)

def create_simplified_y_bus(B, lines, slack_bus_number):
    """
    Method for creating the simplified Y-bus used in fast decoupled power flow and for distribution factors
    """
    for line in lines:
        B[line.from_bus.bus_number - 1, line.to_bus.bus_number - 1] = -1 / line.reactance
        B[line.to_bus.bus_number - 1, line.from_bus.bus_number - 1] = B[
            line.from_bus.bus_number - 1, line.to_bus.bus_number - 1]
    # Get the sum of the rows
    diagonal_elements = np.sum(B, axis=1)  # axis 1 meaning the we sum each colomn along the rows
    for i, Y_ii in enumerate(diagonal_elements):
        B[i, i] = -Y_ii  # subracting because the off diagonal elements are negative (--=+)
    # Remove row and column of slack bus
    B = np.delete(B, slack_bus_number - 1, axis=0)
    B = np.delete(B, slack_bus_number - 1, axis=1)
    return B

