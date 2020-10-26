from classes import Bus, Line
from supporting_methods import print_title1, print_title3, create_simplified_y_bus, calculate_distribution_factors
import numpy as np
from IMML_algorithm import IMML_algorithm

#import pyomo.environ as pyo
#from pyomo.opt import SolverFactory
import copy
"""
Settings
"""
outage_task_1 = "2-3"
error = 1e-5
basecase_dispatch = {"1": 0.8, "2": 0, "3": 1.3, "4": 1}
"""
Input values
"""
slack_bus_number = 4
V = {"1": 1, "2": 1, "3": 1, "4": 1 }
delta = {"1": 0, "2": 0, "3": 0, "4": 0}
# Q values from assignment
Q = {"1": 0, "2": 0, "3": 0, "4": None}
# P values from project
P = {"1": -1.6, "2": -0.9, "3": -0.6, "4": 0}
# From LTsolve
#P = {"1": -0.8, "2": -0.9, "3": 0.7, "4": 0}
# line data
r = {"1-2": 0.0, "1-3": 0.0, "2-3": 0.0, "3-4": 0}
x = {"1-2": 0.2, "1-3": 0.1, "2-3": 0.2, "3-4": 0.25}
gen_cost = {"Bus 1": 4, "Bus 2": 5, "Bus 3": 3, "Bus 4": 2}
trans_cap = {"Line 1-2": 1, "Line 1-3": 1, "Line 2-3": 1.5, "Line 3-4": 1}


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

for line in lines:
    line.transfer_capacity = trans_cap[line.name]

P_array = np.zeros([len(P)-1, 1])
voltage_angles_labels = []
for index, p_spec in enumerate(P.values()):
    if p_spec:
        P_array[index] = p_spec
        voltage_angles_labels.insert(index, "\u03B4{}".format(index+1))


print_title1("Task 1")

from_bus, to_bus  = IMML_algorithm(P, buses, lines, slack_bus_number, outage_task_1, h_modification=0.5, printing=False)
# Removal of one line means the equivalent impedance on the remaining line is doubled
line_23.reactance *= 2
line_23.transfer_capacity /=2

# Congestion check
congested = False
congested_lines = []
print_title3("Power flow on lines")
for line in lines:
    line.p_power_flow = (line.from_bus.delta - line.to_bus.delta) / line.reactance
    print("\n{} power flow: {}pu".format(line.name, round(line.p_power_flow, 4)))
    print("Transfer capacity {}pu".format(line.transfer_capacity))
    if line.transfer_capacity <= (abs(line.p_power_flow) - error):
        congested = True
        congested_lines.append(line)

if congested:
    print("\nCongested lines:")
    for line in congested_lines:
        print("", line.name)
    print("\nSince one or more lines are congested, a subproblem is formulated and solved.")
else:
    print("\nNo lines are congested.")

print("\nNote that only the constraints corresponding to the congested lines are included in the OPF using LPsolve.")



print_title1("Task 2")
#Behold restriksjonene som var fra før, line 1-3 and line 3-4.

# For convenience the distribution factors are calculated
B_p = np.zeros([len(P), len(P)])
B_p = create_simplified_y_bus(B_p, lines, slack_bus_number)
_, a_dict = calculate_distribution_factors(B_p, P_array, buses, lines, slack_bus_number, printing=False)
"""
# Calculate the flow on lines using the distribution factors
for line in lines:
    line.p_power_flow = 0
    for bus in buses.values():
        if bus.bus_number != slack_bus_number:
            line.p_power_flow += a_dict[line.name][bus.bus_number-1][0] * (basecase_dispatch[str(bus.bus_number)] + bus.p_spec)
    print("\n{} power flow: {}".format(line.name, round(line.p_power_flow, 4)), "pu")
"""
# The subproblem is formulated


# Create the objective function for the subproblem

print("\nObjective function for subproblem: ")
print("min: Ks = (PG1-PG10) + (PG2-PG20) + (PG3-PG30) + (PG4-PG40)- (PG1-PG10) - (PG2-PG20) - (PG3-PG30) - (PG4-PG40)")
# Create constraints based on IMML congestion analysis to be considered in the subproblem
print("\nConstraint(s) to be added to the subproblem formulation in addition to the load balance: \n")
for line in congested_lines:
    line.p_power_flow = line.p_power_flow[0][0]
    bus_numbers = "{}{}".format(line.from_bus.bus_number, line.to_bus.bus_number)
    s_up = "P{}: {} (PG1-PG10) + {} (PG2-PG20) + {} (PG3-PG30) >= {}".format(bus_numbers, round(a_dict[line.name][0][0], 2),
                                                round(a_dict[line.name][1][0], 2), round(a_dict[line.name][2][0],2),
                                                round(-line.transfer_capacity - line.p_power_flow, 2))
    s_down = "P{}: {} (PG1-PG10) + {} (PG2-PG20) + {} (PG3-PG30) <= {}".format(bus_numbers, round(-a_dict[line.name][0][0], 2),
                                                round(-a_dict[line.name][1][0], 2), round(-a_dict[line.name][2][0], 2),
                                                  round(line.transfer_capacity - line.p_power_flow,2))
    print(s_up)
    print(s_down)

print_title1("Task 3")

# Input this objective function and these constraints to LPsolve and get new optimal solution


