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

def run_newton_raphson(N_R, printing=True, total_iterations=None, convergence=False):
    """
    Runs the general newton raphson method
    N_R: Load_Flow, Contination or Fast_Decoupled object
    set printing=False to suppress printing
    Set input total_iterations if external iterations should be counted
    Set convergence=True to track convergence
    """
    if convergence:
        converging = True
    while N_R.power_error() > 0.0001:
        N_R.iteration += 1
        if total_iterations:
            total_iterations += 1
        N_R.reset_values()
        N_R.calc_new_power_injections()
        N_R.error_specified_vs_calculated()
        N_R.jacobian.create()
        N_R.find_x_diff()
        N_R.update_values()
        if printing:
            print_title3("Iteration: {}".format(N_R.iteration))
            N_R.print_matrices()
        if N_R.diverging():
            if printing:
                print_title3("No convergence")
                converging = False
            break
    if total_iterations and convergence:
        return total_iterations, converging
    elif convergence and not total_iterations:
        return converging
    elif total_iterations and not convergence:
        return  total_iterations

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

def calculate_distribution_factors(B_p, P_array, buses, lines, slack_bus_number, avoid_line=None, printing=True):
    """
    Runs a distribution factor calculation
    """
    right_hand_side_dict = {}
    a_dict = {}
    right_hand_side = np.zeros([len(P_array), 1])
    for line in lines:
        for bus in buses.values():
            if bus.bus_number == slack_bus_number:
                pass
            elif bus.bus_number == line.from_bus.bus_number:
                right_hand_side[int(bus.bus_number) - 1] = 1 / line.reactance
            elif bus.bus_number == line.to_bus.bus_number:
                right_hand_side[int(bus.bus_number) - 1] = -1 / line.reactance
            else:
                right_hand_side[int(bus.bus_number) - 1] = 0
        right_hand_side_dict[line.name] = deepcopy(right_hand_side)
        a = np.linalg.solve(B_p, right_hand_side)
        a_dict[line.name] = a
        line.p_power_flow = np.matmul(np.transpose(a), P_array)
        if line.name != avoid_line and printing:
            print_title3("{}:".format(line.name))
            print("\nRight hand side: \n", np.round(right_hand_side_dict[line.name], 4))
            print("\nDistribution factors: \n", np.round(a_dict[line.name], 4))
            print("\nPower flow on line:\n", "{}pu".format(round(line.p_power_flow[0][0], 4)))
    return right_hand_side_dict, a_dict