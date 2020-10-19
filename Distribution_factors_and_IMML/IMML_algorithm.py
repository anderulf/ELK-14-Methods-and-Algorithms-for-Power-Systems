import numpy as np

def IMML_algorithm(specified_active_powers, buses, lines, slack_bus_number, outage, printing=True):
    """
    Outage should be a string in the format "x-y" where x is the from bus, and y is the to bus
    The buses are changed by the algorithm so the bus data can be used after the algorithm is run

    Printing is set to True by default
    """
    H = np.zeros([len(buses), len(buses)], dtype=float)
    for line in lines:
        H[line.from_bus.bus_number - 1, line.to_bus.bus_number - 1] = -1 / line.reactance
        H[line.to_bus.bus_number - 1, line.from_bus.bus_number - 1] = H[
            line.from_bus.bus_number - 1, line.to_bus.bus_number - 1]
    # Get the sum of the rows
    diagonal_elements = np.sum(H, axis=1)  # axis 1 meaning the we sum each colomn along the rows
    for i, Y_ii in enumerate(diagonal_elements):
        H[i, i] = -Y_ii  # subracting because the off diagonal elements are negative (--=+)
    # Remove slack bus
    H = np.delete(H, slack_bus_number - 1, axis=0)
    H = np.delete(H, slack_bus_number - 1, axis=1)
    from_bus, to_bus = outage.split("-")
    from_bus = int(from_bus)
    to_bus = int(to_bus)
    delta_h = H[from_bus-1, to_bus-1]/2
    M = np.zeros([len(buses)-1, 1])
    M[from_bus-1] = 1
    M[to_bus-1] = -1
    # X-vector
    x = np.linalg.solve(H, M)
    P_array = np.zeros([len(buses)-1, 1])
    for bus in specified_active_powers:
        if int(bus) == slack_bus_number:
            pass
        else:
            P_array[int(bus)-1] = specified_active_powers[bus]
    delta_0 = np.linalg.solve(H, P_array)
    temp_correction_2 = np.matmul(np.transpose(M), delta_0)
    c_inverse = 1/delta_h + np.linalg.inv(H)[from_bus-1, from_bus-1] - np.linalg.inv(H)[from_bus-1, to_bus-1]
    c = 1/c_inverse
    delta_correction_temp_1 = c * temp_correction_2
    #delta_correction = -x * delta_correction_temp_1
    delta_correction = np.matmul(-x, delta_correction_temp_1)
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
        print("\nDelta_0:\n ", delta_0)
        print("\nX-vector: \n", x)
        print("\nTemp correction 2: \n", temp_correction_2)
        print("\nc_inverse: \n", c_inverse)
        print("\nc: \n", c)
        print("\ndelta_corrections: \n", delta_correction)

