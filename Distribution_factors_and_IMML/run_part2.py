from classes import Bus, Line
from Distribution_factors_and_IMML.IMML_algorithm import IMML_algorithm
import numpy as np
from supporting_methods import print_title1, print_title3

"""
Settings
"""
outage_task_2 = "1-2"
outage_task_3 = "1-3"

line_12_flow_from_part1 = -0.0455
line_13_flow_from_part1 = -1.2045
line_23_flow_from_part1 = -0.4455
line_34_flow_from_part1 = -2.25
lines_from_part1 = [line_12_flow_from_part1, line_13_flow_from_part1, line_23_flow_from_part1, line_34_flow_from_part1]
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
    buses[int(bus_number)] = Bus(int(bus_number), P[bus_number], Q[bus_number], V[bus_number], delta[bus_number] )
# Add lines
line_12 = Line(buses[1], buses[2], r["1-2"], x["1-2"])
line_13 = Line(buses[1], buses[3], r["1-3"], x["1-3"])
line_23 = Line(buses[2], buses[3], r["2-3"], x["2-3"])
line_34 =Line(buses[3], buses[4], r["3-4"], x["3-4"])

lines = [line_12, line_13, line_23, line_34]

"""
Program
"""
#Part 2
print_title1("Task 1")
print("\nIMML is a fast and general technique for simulating impacts of modifying a network topolgy, such as outages\n"
      "and compensations. The method is based on creating three sub-matrices, H, delta_h and M, which will be constant\n"
      "during the entire process. H describes the original network topology, delta_h reflects what to me manipulated\n"
      "and M reflects which lines/buses to be affected by the manipulation. The results after applying an IMML are\n"
      "the state of the new voltage angles as well as the new power flow for the given modification\n")

print_title1("Task 2")
#Find voltage angles and the power flow when the line 1-2 is disconnected by using the IMML
from_bus, to_bus = IMML_algorithm(P, buses, lines, slack_bus_number, outage_task_2)
# from_bus and to_bus is the line which is not considered
for index, line in enumerate(lines):
    if from_bus == line.from_bus.bus_number and to_bus == line.to_bus.bus_number:
        pass
    else:
        line.p_power_flow =  (line.from_bus.delta - line.to_bus.delta)/line.reactance
        print_title3(line.name)
        print("\nActive power flow on line:", round(line.p_power_flow, 3))
        print("Change from basecase: {}pu".format(round(line.p_power_flow - lines_from_part1[index], 3)))

print("\nWhen line {} is disconnected the flow which originally was flowing between bus 1 and 2 is flowing on line 1-3\n"
      "instead. Hence the change on line 1-3 and line 2-3 is opposite. Less power flows from bus 3 to bus 2 because it\n"
      "cannot flow over to bus 1.".format(outage_task_2))

print_title1("Task 3")

from_bus, to_bus  = IMML_algorithm(P, buses, lines, slack_bus_number, outage_task_3, h_modification=0.5)
# Removal of one line means the equivalent impedance on the remaining line is doubled
line_13.reactance *= 2
for index, line in enumerate(lines):
    line.p_power_flow = (line.from_bus.delta - line.to_bus.delta)/line.reactance
    print_title3(line.name)
    print("\nActive power flow on line:", round(line.p_power_flow, 3))
    print("Change from basecase: {}pu".format(round(line.p_power_flow - lines_from_part1[index], 3)))

print("\nBecause one of the lines are disconnected, the net reactance doubles to {} on line 1-3. The flow from bus 3\n"
      "to bus 1 has decreased because of this increase in reactance. Because of this, more power flows from bus 3 to \n"
      "bus 1 via bus 2. The flow on line 3-4 has not changed because the total load in the system is unchanged."
      "".format(line_13.reactance))
P_array_new = np.zeros([len(buses) - 1, 1])
for line in lines:
    if line.from_bus.bus_number == slack_bus_number:
        pass
    else:
        P_array_new[line.from_bus.bus_number-1][0] += (line.from_bus.delta - line.to_bus.delta)/line.reactance
    if line.to_bus.bus_number != slack_bus_number:
        P_array_new[line.to_bus.bus_number - 1][0] += -(line.from_bus.delta - line.to_bus.delta) / line.reactance
    else:
        pass

print("\nP_array calculated from flow on lines: \n", np.round(P_array_new, 4))
print("\nNote that the above P_array is equal to the specified in the input. This verifies the line flows and that the\n"
      "IMML is an alternative method to calculate the same result as by using the distribution factors. Task 2 and 3 \n"
      "shows that topology changes are easily calculated using the IMML method which is beneficial for contingency \n"
      "analysis. This is because there is no need to change the topology matrix H (B'). Only the M-matrix and \u0394h\n"
      "is changed which is a very fast modification.")

print_title1()

