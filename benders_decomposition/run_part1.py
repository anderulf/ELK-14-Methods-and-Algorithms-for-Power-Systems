from supporting_classes import Bus, Line
from supporting_methods import print_title1, print_title3, create_simplified_y_bus
from distribution_factors_and_IMML.support import calculate_distribution_factors
import numpy as np
#import pyomo.environ as pyo
#from pyomo.opt import SolverFactory
import copy
"""
Settings
"""
error = 1e-5
# From LTsolve
#P = {"1": -0.8, "2": -0.9, "3": 0.7, "4": 0}
"""
Input values
"""
slack_bus_number = 4
V = {"1": 1, "2": 1, "3": 1, "4": 1 }
delta = {"1": 0, "2": 0, "3": 0, "4": 0}
# Specified active load at bus
P = {"1": -1.6, "2": -0.9, "3": -0.6, "4": 0}
# Specified reactive load at bus
Q = {"1": 0, "2": 0, "3": 0, "4": None}
# line data
r = {"1-2": 0.0, "1-3": 0.0, "2-3": 0.0, "3-4": 0}
x = {"1-2": 0.2, "1-3": 0.1, "2-3": 0.2, "3-4": 0.25}
gen_cost = {"1": 4, "2": 5, "3": 3, "4": 2}
trans_cap = {"Line 1-2": 1, "Line 1-3": 1, "Line 2-3": 1.5, "Line 3-4": 1}


# Create buses
bus1 = Bus(1, P["1"], Q["1"], V["1"], delta["1"])
bus2 = Bus(2, P["2"], Q["2"], V["2"], delta["2"])
bus3 = Bus(3, P["3"], Q["3"], V["3"], delta["3"])
bus4 = Bus(4, P["4"], Q["4"], V["4"], delta["4"])
buses = {1: bus1, 2: bus2, 3: bus3, 4: bus4}
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

B_p = np.zeros([len(P), len(P)])
B_p = create_simplified_y_bus(B_p, lines, slack_bus_number)

_, a_dict = calculate_distribution_factors(B_p, P_array, buses, lines, slack_bus_number, printing=False)

# Congestion check
congested = False
congested_lines = []
print_title3("Power flow on lines")
for line in lines:
    line.p_power_flow = np.matmul(np.transpose(a_dict[line.name]), P_array)[0][0]
    print("\n{} power flow: {}pu".format(line.name, round(line.p_power_flow, 4)))
    print("Transfer capacity {}pu".format(line.transfer_capacity))
    if line.transfer_capacity <= (abs(line.p_power_flow) - error):
        congested = True
        congested_lines.append(line)

if congested:
    print("\nCongested line(s):\n")
    for line in congested_lines:
        print(line.name)
    print("\nSince one or more lines are congested, linear programming is required.")
else:
    print("\nNo lines are congested.")

print("Only the constraints corresponding to the congested lines are included in the OPF using LPsolve.")

# Create LPsolve formulation
print_title3("Master formulation:")

print("min: {} PG1 + {} PG2 + {} PG3 + {} PG4;".format(gen_cost["1"], gen_cost["2"], gen_cost["3"], gen_cost["4"]))
# Create constraints based on IMML congestion analysis to be considered in the subproblem

# Either get the constraints from task 2 or make new with a_dict before contingency, if so just calculate before
#  removing the line in the beginning of the code and create with code below
print("\n/*Constraints: */\n")
for line in congested_lines:
    bus_numbers = "{}{}".format(line.from_bus.bus_number, line.to_bus.bus_number)
    line_constraint_basecase = "P{0}: {1} <= {2} PG1 + {3} + {4} PG2 + {5} + {6} PG3 + {7} <= {8};".format(bus_numbers, round(-line.transfer_capacity, 2),
            round(a_dict[line.name][0][0], 2), round(a_dict[line.name][0][0] * bus1.p_spec,3),round(a_dict[line.name][1][0], 2), round(a_dict[line.name][1][0] * bus2.p_spec,3), round(a_dict[line.name][2][0],2),
            round(a_dict[line.name][2][0] * bus3.p_spec,3), round(line.transfer_capacity, 2))

    print(line_constraint_basecase)

print("\nLoad_balance: PG1 + PG2 + PG3 + PG4 = {};\n".format(-sum(P_array)[0]))

for bus in buses.values():
    print("PG{1}_limit: PG{1} >= 0;".format(bus.bus_number, bus.bus_number))

print_title1("Task 5")

# Calculated values from LP-solve for part 1 task 1
dispatch = {"1": 0.8, "2": 0, "3": 1.3, "4": 1}
k = 9.1
dispatch_duals   = {"1": 0, "2": 0, "3": 0, "4": 0} # Dual values for dispatch limits set to prod > 0
for index, bus in enumerate(buses.values()):
    bus.dispatch = dispatch["{}".format(bus.bus_number)]
    bus.gen_cost = gen_cost["{}".format(bus.bus_number)]
    bus.p_gen = dispatch["{}".format(bus.bus_number)]
    bus.marginal_cost = dispatch_duals["{}".format(bus.bus_number)]



print("\nOptimal objective function value, from LPsolve:", k)
for bus_number in dispatch.keys():
    print("Dispatch for bus {} : {} pu".format(int(bus_number)+1, dispatch[bus_number]))

# New flow on line after optimization
print_title3("New flow on lines after optimization")
for line in lines:
    # Reset current power flow
    line.p_power_flow = 0
    for bus in buses.values():
        if bus.bus_number != slack_bus_number:
            line.p_power_flow += a_dict[line.name][bus.bus_number-1][0] * (bus.dispatch + bus.p_spec)
    print("{} power flow: {} pu".format(line.name, round(line.p_power_flow, 3)))

print_title1("Task 2")
#Task 2
#Check the marginal costs (reduced cost, dual variables) and check these against the operating cost at each bus.

duals = {"1": 0, "2": 1.5, "3": 0, "4": 0}

print_title3("Dual values")
for index, dual in enumerate(duals.values()):
    print("Bus {} dual value: {}".format(index+1, dual))


print("\nDual values gives indications on the marginal cost of production. All generators that produce will obtain an\n"
      "dual value equal to zero. This is because a upper bound for dispatch does not exist. On the other hand we \n"
      "observe that PG2 has a dual value equal to {}. This is because PG2 is limited by the lower bound, zero, which\n"
      "means the objective function would be decreased by 1.5 if PG2 were allowed to produce -1 pu. As this is not\n"
      "a realistic case all producers will have marginal costs equal to their production cost, C. Line constraints\n"
      "with dual values different from zero are limited, which means the objective function will decrease if the \n"
      "transfer capacity on these lines are increased.".format(duals["2"]))
print_title3("Marginal costs")
for bus in buses.values():
    print("Bus {} marginal cost: {}".format(bus.bus_number, bus.gen_cost))