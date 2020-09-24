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