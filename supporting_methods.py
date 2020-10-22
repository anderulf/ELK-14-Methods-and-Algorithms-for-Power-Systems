import numpy as np

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

def complex_angle(complex_number):
    r, angle = rectangular_to_polar(complex_number)
    return angle

def complex_radius(complex_number):
    r, angle = rectangular_to_polar(complex_number)
    return r

def print_title1(input_string=None):
    """
    Prints a main title generally for a new task
    """
    symbol_count = 100
    symbol = "#"
    if input_string:
        print("\n", int(symbol_count / 2) * symbol, input_string, int(symbol_count / 2) * symbol)
    else:
        print("\n", symbol_count * symbol)

def print_title2(input_string=None):
    """
    Prints a undertitle generally for a new problem
    """
    symbol_count = 80
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