import numpy as np

def IMML_algorithm(specified_active_powers, buses, lines, slack_bus_number, outage, h_modification=1, printing=True):
    """
    Outage should be a string in the format "x-y" where x is the from bus, and y is the to bus
    The buses are changed by the algorithm so the bus data can be used after the algorithm is run

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
    # Remove slack bus
    H = np.delete(H, slack_bus_number - 1, axis=0)
    H = np.delete(H, slack_bus_number - 1, axis=1)
    from_bus, to_bus = outage.split("-")
    from_bus = int(from_bus)
    to_bus = int(to_bus)
    delta_h = H[from_bus-1, to_bus-1]*h_modification
    M = np.zeros([len(buses)-1, 1])
    M[from_bus-1] = 1
    M[to_bus-1] = -1

    # X-vector = H_invers * M
    #x = H_invers[:,from_bus-1] - H_invers[:,to_bus-1]
    x = np.linalg.solve(H, M)

    #z = M_transpose * x = M_transpose *H_invers*M
    z = x[from_bus-1] - x[to_bus -1]
    print("z = ", z)

    P_array = np.zeros([len(buses)-1, 1])
    for bus in specified_active_powers:
        if int(bus) == slack_bus_number:
            pass
        else:
            P_array[int(bus)-1] = specified_active_powers[bus]
    delta_0 = np.linalg.solve(H, P_array)

    #M_transpose*H_invers*P
    angle_diff = (delta_0[from_bus-1] - delta_0[to_bus-1])[0]

    #temp_correction_2 = np.matmul(np.transpose(M), delta_0)[0][0] # Is a scalar value

    c_inverse = (1/delta_h + z)[0] # Is generally a scalar value
    # Utilizing scalar values is a lot more efficiency compared to using matrix multiplication for larger systems
    # Special Case occurs when running several modifaction in the systems simustaneously
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
        print("\nTemp correction 2: \n", np.round(angle_diff, 4))
        print("\nc_inverse: \n", np.round(c_inverse, 4))
        print("\nc: \n", np.round(c, 4))
        print("\ndelta_corrections: \n", np.round(delta_correction, 4))
        print("\nNew angles: \n", np.round(delta, 4), "\n")
    return from_bus, to_bus