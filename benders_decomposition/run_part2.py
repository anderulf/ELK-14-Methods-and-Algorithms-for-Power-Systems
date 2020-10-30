from supporting_classes import Bus, Line
from supporting_methods import print_title1, print_title2, print_title3, create_simplified_y_bus, get_from_and_to_bus
import numpy as np
from distribution_factors_and_IMML.support import IMML_algorithm, calculate_distribution_factors
import copy
"""
Settings
"""
outage_task_1 = "2-3"
error = 1e-5
basecase_dispatch = {"1": 0.8, "2": 0, "3": 1.3, "4": 1}
basecase_congested_lines = []
"""
Input values
"""
slack_bus_number = 4
V = {"1": 1, "2": 1, "3": 1, "4": 1 }
delta = {"1": 0, "2": 0, "3": 0, "4": 0}
# Q values from assignment
Q = {"1": 0, "2": 0, "3": 0, "4": None}
# P values from project
P = {"1": -1.6, "2": -0.9, "3": -0.6, "4": None}
# From LTsolve
#P = {"1": -0.8, "2": -0.9, "3": 0.7, "4": 0}
# line data
r = {"1-2": 0.0, "1-3": 0.0, "2-3": 0.0, "3-4": 0}
x = {"1-2": 0.2, "1-3": 0.1, "2-3": 0.2, "3-4": 0.25}
gen_cost = {"1": 4, "2": 5, "3": 3, "4": 2}
trans_cap = {"1-2": 1, "1-3": 1, "2-3": 1.5, "3-4": 1}


# Create buses
bus1 = Bus(1, P["1"], Q["1"], V["1"], delta["1"], dispatch=basecase_dispatch["1"])
bus2 = Bus(2, P["2"], Q["2"], V["2"], delta["2"], dispatch=basecase_dispatch["2"])
bus3 = Bus(3, P["3"], Q["3"], V["3"], delta["3"], dispatch=basecase_dispatch["3"])
bus4 = Bus(4, P["4"], Q["4"], V["4"], delta["4"], dispatch=basecase_dispatch["4"])
buses = {1: bus1, 2: bus2, 3: bus3, 4: bus4}
# Add lines
line_12 = Line(buses[1], buses[2], r["1-2"], x["1-2"], trans_cap["1-2"])
line_13 = Line(buses[1], buses[3], r["1-3"], x["1-3"], trans_cap["1-3"])
line_23 = Line(buses[2], buses[3], r["2-3"], x["2-3"], trans_cap["2-3"])
line_34 = Line(buses[3], buses[4], r["3-4"], x["3-4"], trans_cap["3-4"])

lines = [line_12, line_13, line_23, line_34]

# Add the lines that were congested in the basecase solution
basecase_congested_lines.extend([line_13, line_34])


P_array = np.zeros([len(P)-1, 1])
system_load = 0
voltage_angles_labels = []
for index, p_spec in enumerate(P.values()):
    if p_spec:
        P_array[index] = basecase_dispatch[str(index +1)] + p_spec
        system_load -=p_spec
        voltage_angles_labels.insert(index, "\u03B4{}".format(index+1))

P_array_basecase = copy.deepcopy(P_array)

B_p_basecase = np.zeros([len(P), len(P)])
B_p_basecase = create_simplified_y_bus(B_p_basecase, lines, slack_bus_number)
_, a_dict_basecase = calculate_distribution_factors(B_p_basecase, P_array, buses, lines, slack_bus_number, printing=False)

lines_basecase = copy.deepcopy(lines)

print_title1("Task 1")

# Update
from_bus, to_bus = get_from_and_to_bus(outage_task_1)
IMML_algorithm(P_array, buses, lines, slack_bus_number, from_bus, to_bus, h_modification=0.5, printing=False)
# Removal of one line means the equivalent impedance on the remaining line is doubled
line_23.reactance *= 2
line_23.transfer_capacity /=2

# Congestion check
congested = False
congested_lines = []
print_title3("Power flow on lines after running an IMML")
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

# For convenience the distribution factors are calculated
B_p = np.zeros([len(P), len(P)])
B_p = create_simplified_y_bus(B_p, lines, slack_bus_number)
_, a_dict = calculate_distribution_factors(B_p, P_array, buses, lines, slack_bus_number, printing=False)

# The subproblem is formulated

# Create the objective function for the subproblem

print_title3("Objective function for subproblem:")
print("\nmin: Ks = PG1_up + PG2_up + PG3_up + PG4_up + PG1_down + PG2_down + PG3_down + PG4_down;")
# Create constraints based on IMML congestion analysis to be considered in the subproblem
print("\n/* Constraint(s) to be added to the subproblem formulation in addition to the load balance: */ \n")
for line in congested_lines:
    line.p_power_flow = line.p_power_flow[0][0]
    bus_numbers = "{}{}".format(line.from_bus.bus_number, line.to_bus.bus_number)
    constraint = "P{0}: {1} <= {2} PG1_up + {3} PG2_up + {4} PG3_up - {2} PG1_down - {3} PG2_down - {4} PG3_down" \
                 " <= {5};".format(bus_numbers, round(-line.transfer_capacity - line.p_power_flow, 3), round(a_dict[line.name][0][0], 3),
                                                round(a_dict[line.name][1][0], 3), round(a_dict[line.name][2][0],3),
                                                round(line.transfer_capacity - line.p_power_flow, 3))

    print(constraint)
    print("\nLoad balance: PG1_up + PG2_up + PG3_up + PG4_up - PG1_down - PG2_down - PG3_down - PG4_down = 0;\n" )

for bus in buses.values():
    print("PG{1}_up_lower  : PG{1}_up >= 0;".format(bus.bus_number, bus.bus_number))
    print("PG{0}_down_upper: PG{0}_down <= {1};".format(bus.bus_number, basecase_dispatch[str(bus.bus_number)]))


print_title1("Task 3")

# Input this objective function and these constraints to LPsolve and get new optimal solution
"""
Output from LPsolve:
PG1_up                   0.233372
PG2_up                          0
PG3_up                          0
PG4_up                          0
PG1_down                        0
PG2_down                        0
PG3_down                 0.233372
PG4_down                        0
"""
# Get the output from LPsolve
k = 0.46674446
dispatch_change ={"1": 0.233372, "2": 0, "3": -0.233372, "4": 0}
new_dispatch = {}
duals = {"1": 0, "2": 0.667445, "3": 2, "4": 2}
costs = {"1": 1, "2": 1, "3": 1, "4": 1}
for bus in buses.values():
    bus.sensitivity = duals["{}".format(bus.bus_number)] - costs["{}".format(bus.bus_number)]
# Print the results
print_title3("Subproblem LPsolve output")
print("\nOptimal objective function value, from LPsolve:", round(k, 3))
for bus_number in dispatch_change.keys():
    new_dispatch[bus_number] = basecase_dispatch[bus_number] + dispatch_change[bus_number]
    print("\nChange in dispatch for bus {}: {} pu".format(int(bus_number), round(dispatch_change[bus_number], 3)))
    print("Dispatch for bus {}: {} pu".format(int(bus_number), round(new_dispatch[bus_number], 3)))
    print("Sensitivity for bus {}: {} ".format(int(bus_number), round(buses[int(bus_number)].sensitivity, 3)))

print_title1("Task 4")

# Formulate the constraint to be added to master problem in order to make a feasible solution with the contingency

master_constraint = "Benders_cut: Ks + C1 * dKs/d\u0394PG1 * \u0394PG1 + C2 * dKs/d\u0394PG2 * \u0394PG2 + C3 * dKs/d\u0394PG3 * \u0394PG3 + C4 * dKs/d\u0394PG4" \
                    " * \u0394PG4 <= 0;"


master_constraint_values = "Benders_cut: {} ".format(round(k, 3))
rhs = -k
master_constraint_values_final = "Benders_cut: "
for index, bus in enumerate(buses.values()):
    master_constraint_values += " + {}*(PG{}-{})".format(round(bus.sensitivity,3),bus.bus_number , basecase_dispatch["{}".format(bus.bus_number)])
    if index != 0:
        master_constraint_values_final += " + "

    master_constraint_values_final += "{} PG{}".format(round(bus.sensitivity,3), bus.bus_number)
    rhs += bus.sensitivity * basecase_dispatch["{}".format(bus.bus_number)]

master_constraint_values += " <= 0;"
master_constraint_values_final += " <= {};".format(round(rhs, 3))
print_title3("Final constraint to be added to master problem: ")
print("\n")
print(master_constraint)
print(master_constraint_values)
print(master_constraint_values_final)

print_title3("Final master formulation:")

print("\nmin: Ks = {} PG1 + {} PG2 + {} PG3 + {} PG4;".format(gen_cost["1"], gen_cost["2"], gen_cost["3"], gen_cost["4"]))
# Create constraints based on IMML congestion analysis to be considered in the subproblem

# Either get the constraints from task 2 or make new with a_dict before contingency, if so just calculate before
#  removing the line in the beginning of the code and create with code below
print("\n/*Constraints: */\n")
for line in basecase_congested_lines:
    bus_numbers = "{}{}".format(line.from_bus.bus_number, line.to_bus.bus_number)
    line_constraint_basecase = "P{0}: {1} <= {2} PG1 + {3} + {4} PG2 + {5} + {6} PG3 + {7} <= {8};".format(bus_numbers,
        round(-line.transfer_capacity,2), round(a_dict_basecase[line.name][0][0], 2), round(a_dict_basecase[line.name][0][0] *
        bus1.p_spec, 3), round(a_dict_basecase[line.name][1][0], 2), round(a_dict_basecase[line.name][1][0] * bus2.p_spec, 3),
        round(a_dict_basecase[line.name][2][0], 2), round(a_dict_basecase[line.name][2][0] * bus3.p_spec, 3),
        round(line.transfer_capacity, 2))
    print(line_constraint_basecase)

print("\n{};".format(master_constraint_values_final))
print("\nLoad_balance: PG1 + PG2 + PG3 + PG4 = {}; \n".format(system_load))

for bus in buses.values():
    print("PG{1}_limit: PG{1} >= 0;".format(bus.bus_number, bus.bus_number))

print_title1("Task 5")

"""
Value of objective function: 9.3335

Actual values of the variables:
PG1                        1.0335
PG2                             0
PG3                        1.0665
PG4                             1
"""
k_final = 9.3335
bus1.dispatch = 1.0335
bus2.dispatch = 0
bus3.dispatch = 1.0665
# Update new P_array for BC dispatch
for bus in buses.values():
    if bus.bus_number != slack_bus_number:
        P_array[bus.bus_number -1] = bus.dispatch + bus.p_spec

print_title3("LP-solve output")
for bus in buses.values():
    print("Bus {} dispatch: {} pu".format(bus.bus_number, round(bus.dispatch, 4)))
print("\nThe new master formulation with the benders cut has a small increase in the objective function. In the \n"
      "basecase formulation the cost was 9.1, and it has now increased to {}. This is because the master formulation\n"
      "now has the advantage that a outage on line 2-3 should not cause any congestions. This was done by increasing\n"
      "the production on bus 1 from 0.8 pu to 1.03 pu, and decreasing production on bus 3 from 1.3 pu to 1.07 pu. This\n"
      "will be verified.".format(round(k_final, 1))  )

voltage_angles = np.linalg.solve(B_p, P_array)
voltage_angles_basecase = np.linalg.solve(B_p_basecase, P_array)

# Calculate power flow
# For contingency case
for line in lines:
    if line.from_bus.bus_number != slack_bus_number:
        from_angle = voltage_angles[line.from_bus.bus_number-1][0]
    else: from_angle = 0 #  Slack value
    if line.to_bus.bus_number  != slack_bus_number:
        to_angle = voltage_angles[line.to_bus.bus_number - 1][0]
    else: to_angle = 0  # Slack angle
    line.p_power_flow = (from_angle - to_angle)/line.reactance

# For basecase
for line in lines_basecase:
    if line.from_bus.bus_number != slack_bus_number:
        from_angle = voltage_angles_basecase[line.from_bus.bus_number - 1][0]
    else: from_angle = 0 #  Slack value
    if line.to_bus.bus_number  != slack_bus_number:
        to_angle = voltage_angles_basecase[line.to_bus.bus_number - 1][0]
    else: to_angle = 0  # Slack angle
    line.p_power_flow = (from_angle - to_angle)/line.reactance

# Congestion check
congested = False
congested_lines = []
print_title2("Power flow on lines after running benders cut")
print_title3("Contingency case on line {}".format(outage_task_1))
for line in lines:
    print("\n{} power flow: {} pu".format(line.name, round(line.p_power_flow, 4)))
    print("Transfer capacity {} pu".format(line.transfer_capacity))
    if line.transfer_capacity <= (abs(line.p_power_flow) - error):
        print(line.name)
        congested = True
        congested_lines.append(line)

congested_basecase = False
congested_lines_basecase = []
print_title3("Basecase")
for line in lines_basecase:
    print("\n{} power flow: {} pu".format(line.name, round(line.p_power_flow, 4)))
    print("Transfer capacity {} pu".format(line.transfer_capacity))
    if line.transfer_capacity <= (abs(line.p_power_flow) - error):
        congested_basecase = True
        congested_lines_basecase.append(line)

if congested or congested_basecase:
    if congested:
        print("\nCongested lines for contingency case:")
        for line in congested_lines:
            print(line.name)
    if congested_basecase:
        print("\nCongested lines for basecase:")
        for line in congested_lines_basecase:
            print(line.name)
    print("\nOne or more lines are congested. This means the benders decomposition was not run properly..")
else:
    print("\nNo lines are congested. This means the benders decomposition was done properly. The system is now secured\n"
          "for an outage on line {}. The results are hence verified.".format(outage_task_1))

