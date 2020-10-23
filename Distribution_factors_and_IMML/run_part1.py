import numpy as np
from classes import Line, Bus
from supporting_methods import print_title1, print_title2, print_title3, create_simplified_y_bus, calculate_distribution_factors
import copy
"""
Settings
"""
avoid_line = "Line 2-3x"
"""
Input values
"""
slack_bus_number = 4
V = {"1": 1, "2": 1, "3": 1, "4": 1 }
delta = {"1": 0, "2": 0, "3": 0, "4": 0}
# Q values from assignment
Q = {"1": 0, "2": 0, "3": 0, "4": None}
# P values from project
P = {"1": -1.25, "2": -0.4, "3": -0.6, "4": None}
# line data
r = {"1-2": 0.0, "1-3": 0.0, "2-3": 0.0, "3-4": 0}
x = {"1-2": 0.2, "1-3": 0.1, "2-3": 0.25, "3-4": 0.25}

# Create buses
buses = {}
for bus_number in V:
    buses[int(bus_number)] = Bus(int(bus_number), P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])
# Add lines
line_12 = Line(buses[1], buses[2], r["1-2"], x["1-2"])
line_13 = Line(buses[1], buses[3], r["1-3"], x["1-3"])
line_23 = Line(buses[2], buses[3], r["2-3"], x["2-3"])
line_34 = Line(buses[3], buses[4], r["3-4"], x["3-4"])

lines = [line_12, line_13, line_23, line_34]

#Assignment 4 part 1

"""
X_12 = 0.2
X_13 = 0.1
X_23 = 0.25
X_34 = 0.25

P_1 = -1.25
P_2 = -0.4
P_3 = -0.6
"""

print_title1("Task 1")
print("In order to build the given equation system it is assumed that R = 0 and that the voltage magnitudes are constant\n"
      "and equal to one, which means the flow of reactive power is zero. The simplified Y-bus (the B-matrix) is built\n"
      "as normal, but with R = 0 and exluding the row and column corresponding to the slack bus. Here is the B-matrix:")
# B_p is created too big, but is skimmed down in the create_simplified_y_bus method
B_p = np.zeros([len(P), len(P)])
B_p = create_simplified_y_bus(B_p, lines, slack_bus_number)
print("\nB'-matrix:\n", B_p)


print_title1("Task 2")

P_array = np.zeros([len(P)-1, 1])
voltage_angles_labels = []
for index, p_spec in enumerate(P.values()):
    if p_spec:
        P_array[index] = p_spec
        voltage_angles_labels.insert(index, "\u03B4{}".format(index+1))
voltage_angles = np.linalg.solve(B_p, P_array)
print("\nSolving the given equation system for the voltage angles")
print("\nVoltage angles:\n", np.c_[voltage_angles_labels, np.round(voltage_angles, 4)])

print_title1("Task 3")

print("\nBy utilizing distribution factors there is no need to calculate the voltage angles directly like in task 2.\n"
      "In addition, fewer constraints are needed in time demanding parts of the OPF. The line violations in the OPF are\n"
      "checked iteratively and only the violated lines are considered in the next iterations. Uneccesarry calculations\n"
      "is avoided. Therefore the use of distribution factors is very effective in OPF calculations. \n\n"
      "Distribution factors functions as sensitivities with respect to lines which in congestion cases can be used as \n"
      "an indication to understand at which buses the power injections should be changed in order to uncongest the \n"
      "corresponding line. ")

right_hand_side_dict, a_dict = calculate_distribution_factors(B_p, P_array, buses, lines, slack_bus_number, avoid_line)
lines_from_task3 = copy.deepcopy(lines)

print_title1("Task 4")
bus_to_be_changed_task4 = "1"
p_diff_task4 = 0.5
P[bus_to_be_changed_task4] -= p_diff_task4
print_title2("Load on bus {} increased with {}pu".format(bus_to_be_changed_task4, p_diff_task4))
P_array[int(bus_to_be_changed_task4)-1] = P[bus_to_be_changed_task4]

for index, line in enumerate(lines):
    print_title3(line.name)
    # Calculate new power flow based on changed specified power
    line.p_power_flow = np.matmul(np.transpose(a_dict[line.name]), P_array)
    print("\nPower flow on line:".format(line.name), "{}pu".format(round(line.p_power_flow[0][0], 4)) )
    print("Change from Basecase: {}pu".format(round(line.p_power_flow[0][0] - lines_from_task3[index].p_power_flow[0][0], 4)))

print("\nWe can observe that the flow from the slack bus (bus 4) to bus 3 is equal to the increase in load on bus {}.\n"
      "This means the slack bus has to cover the change in load.".format(bus_to_be_changed_task4))
lines_from_task4 = lines

print_title1("Task 5")
bus_to_be_changed_task5 = "2"
p_diff_task5 = 0.3
P[bus_to_be_changed_task5] += p_diff_task5
print_title2("Load on bus {} increased with {}pu and load on bus {} decreased with {}pu".format(bus_to_be_changed_task4,p_diff_task4, bus_to_be_changed_task5, p_diff_task5))
P_array[int(bus_to_be_changed_task5)-1] = P[bus_to_be_changed_task5]

for index, line in enumerate(lines):
    print_title3(line.name)
    # Calculate new power flow based on changed specified power
    line.p_power_flow = np.matmul(np.transpose(a_dict[line.name]), P_array)
    print("\nPower flow on line:".format(line.name), "{}pu".format(round(line.p_power_flow[0][0], 4)) )
    print("Change from Basecase: {}pu".format(round(line.p_power_flow[0][0] - lines_from_task3[index].p_power_flow[0][0], 4)))

print("\nLike in task 4 the slack bus covers the net increase in the load in the system which is {}pu. Additionally..."
      " ".format(p_diff_task4-p_diff_task5))