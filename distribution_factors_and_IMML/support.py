import numpy as np
from copy import deepcopy
from supporting_methods import print_title3

# Numpy printing options
np.set_printoptions(suppress=True)  # suppress scientific notations

def IMML_algorithm(P_array, buses, lines, slack_bus_number, from_bus, to_bus, h_modification=1, printing=True):
    """
    Outage should be a string in the format "x-y" where x is the from bus, and y is the to bus
    The buses are changed by the algorithm so the bus data can be used after the algorithm is run

    from_bus and to_bus are the buses on which outage line is connected to

    Printing is set to True by default

    h_modification is a float or integer used to modify a line in the network. Defaulted to 1 if everything is normal
    Set h_modification to 0.5 if a line is removed
    or h_modification to 2 if a line is added
    """
    H = np.zeros([len(buses), len(buses)], dtype=float)
    for line in lines:
        H[line.from_bus.bus_number - 1, line.to_bus.bus_number - 1] = -1 / line.reactance
        H[line.to_bus.bus_number - 1, line.from_bus.bus_number - 1] = H[
            line.from_bus.bus_number - 1, line.to_bus.bus_number - 1]
    # Get the sum of the rows
    diagonal_elements = np.sum(H, axis=1)  # axis 1 meaning the we sum each column along the rows
    for i, Y_ii in enumerate(diagonal_elements):
        H[i, i] = -Y_ii  # subtracting because the off diagonal elements are negative (--=+)
    # Remove slack bus from H-matrix
    H = np.delete(H, slack_bus_number - 1, axis=0)
    H = np.delete(H, slack_bus_number - 1, axis=1)
    delta_h = H[from_bus-1, to_bus-1]*h_modification
    M = np.zeros([len(buses)-1, 1])
    M[from_bus-1] = 1
    M[to_bus-1] = -1

    # X-vector = H_invers * M
    # x = H_invers[:,from_bus-1] - H_invers[:,to_bus-1]
    x = np.linalg.solve(H, M)

    # z = M_transpose * x = M_transpose *H_invers*M
    z = x[from_bus-1] - x[to_bus -1]

    #P_array = np.zeros([len(buses)-1, 1])
    #for bus in specified_active_powers:
        #if int(bus) == slack_bus_number:
            #pass
        #else:
            #P_array[int(bus)-1] = specified_active_powers[bus]
    delta_0 = np.linalg.solve(H, P_array)

    # M_transpose*H_invers*P
    angle_diff = (delta_0[from_bus-1] - delta_0[to_bus-1])[0]

    # temp_correction_2 = np.matmul(np.transpose(M), delta_0)[0][0] # Is a scalar value

    c_inverse = (1/delta_h + z)[0] # Is generally a scalar value
    # Utilizing scalar values is a lot more efficiency compared to using matrix multiplication for larger systems
    # Special Case occurs when running several modifaction in the system simultaneously
    c = 1/c_inverse
    #delta_correction_temp_1 = c * angle_diff
    delta_correction = -x * c * angle_diff
    delta = delta_0 + delta_correction
    index = 0
    for bus in buses.values():
        if bus.bus_number == slack_bus_number :
            pass
        else:
            bus.delta = delta[index][0]
            index += 1
    if printing:
        print("\nH-matrix: \n", H)
        print("\ndelta_h: ", delta_h)
        print("\nM:\n", M)
        print("\nP_array: \n", P_array)
        print("\nDelta_0:\n ", np.round(delta_0, 4))
        print("\nX-vector: \n", np.round(x, 4))
        print("\nc: \n", np.round(c, 4))
        print("\ndelta_corrections: \n", np.round(delta_correction, 4))
        print("\nNew angles: \n", np.round(delta, 4), "\n")

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