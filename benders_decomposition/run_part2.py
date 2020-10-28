from supporting_classes import Bus, Line
from supporting_methods import print_title1, print_title3, create_simplified_y_bus, get_from_and_to_bus
import numpy as np
from distribution_factors_and_IMML.support import IMML_algorithm, calculate_distribution_factors

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
gen_cost = {"1": 4, "2": 5, "3": 3, "4": 2}
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

B_p = np.zeros([len(P), len(P)])
B_p = create_simplified_y_bus(B_p, lines, slack_bus_number)
_, a_dict_old = calculate_distribution_factors(B_p, P_array, buses, lines, slack_bus_number, printing=False)

print_title1("Task 1")
from_bus, to_bus = get_from_and_to_bus(outage_task_1)
IMML_algorithm(P, buses, lines, slack_bus_number, from_bus, to_bus, h_modification=0.5, printing=False)
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
    print("\nSince one or more lines are congested, a subproblem is formulated and solved. In the subproblem, devitations\n"
          "from the basecase is treated as penalties in order to keep the solution close to the basecase. The output \n"
          "from the subproblem solution is a new constraint to be added to the basecase such that the outage does not\n"
          "give any line violations.")
else:
    print("\nNo lines are congested.")

print("\nNote that only the constraints corresponding to the congested lines are included in the OPF using LPsolve.")



print_title1("Task 2")
#Behold restriksjonene som var fra før, line 1-3 and line 3-4.

# For convenience the distribution factors are calculated
B_p = np.zeros([len(P), len(P)])
B_p = create_simplified_y_bus(B_p, lines, slack_bus_number)
_, a_dict = calculate_distribution_factors(B_p, P_array, buses, lines, slack_bus_number, printing=False)

print("\nDist factors before outage: \n", a_dict_old)
print("\nDist factors after outage \n", a_dict)
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

# Get the output from LPsolve
k = 0
new_dispatch = {"1": 0, "2": 0, "3": 0, "4": 0}
new_marginal_cost = {"1": 0, "2": 0, "3": 0, "4": 0}
for bus in buses.values():
    bus.marginal_cost = new_marginal_cost["{}".format(bus_number)]
# Print the results

print("\nOptimal objective function value, from LPsolve:", k)
for bus_number in new_dispatch.keys():
    print("Dispatch for bus {} : {} pu".format(int(bus_number)+1, new_dispatch[bus_number]))

print_title1("Task 4")

# Formulate the constraint to be added to master problem in order to make a feasible solution with the contingency

master_constraint = "Ks + C1 * MC1 * \u0394PG1 + C2 * MC2 * \u0394PG2 + C3 * MC3 * \u0394PG3 + C4 * MC4 * \u0394PG4 <= 0"
master_constraint_values = "{} ".format(k)
for bus in buses.values():
    master_constraint_values += "+ {}*(PG1-{})".format(round(bus.marginal_cost,3), basecase_dispatch["{}".format(bus.bus_number)])

master_constraint_values += " <= 0"
print_title3("Final constraint to be added to master problem: ")
print(master_constraint)
print(master_constraint_values)

print_title3("Final master formulation:")

print("min: Ks = {} PG1 + {} PG2 + {} PG3 + {} PG4".format(gen_cost["1"], gen_cost["2"], gen_cost["3"], gen_cost["4"]))
# Create constraints based on IMML congestion analysis to be considered in the subproblem

# Either get the constraints from task 2 or make new with a_dict before contingency, if so just calculate before
#  removing the line in the beginning of the code and create with code below
print("\nConstraints: \n")
for line in congested_lines:
    line.p_power_flow = line.p_power_flow
    bus_numbers = "{}{}".format(line.from_bus.bus_number, line.to_bus.bus_number)
    s_up = "P{}: {} (PG1-PG10) + {} (PG2-PG20) + {} (PG3-PG30) >= {}".format(bus_numbers, round(a_dict_old[line.name][0][0], 2),
                                                round(a_dict_old[line.name][1][0], 2), round(a_dict_old[line.name][2][0],2),
                                                round(-line.transfer_capacity - line.p_power_flow, 2))
    s_down = "P{}: {} (PG1-PG10) + {} (PG2-PG20) + {} (PG3-PG30) <= {}".format(bus_numbers, round(-a_dict_old[line.name][0][0], 2),
                                                round(-a_dict_old[line.name][1][0], 2), round(-a_dict_old[line.name][2][0], 2),
                                             round(line.transfer_capacity - line.p_power_flow,2))
    print(s_up)
    print(s_down)
    print("\n")

print(master_constraint_values)
print("\nLoad balance: PG1 + PG2 + PG3 + PG4 = {}".format(sum(P_array)[0]))