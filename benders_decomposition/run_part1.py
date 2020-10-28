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
        congested_lines.append(line.name)

if congested:
    print("\nCongested line(s):\n", congested_lines)
    print("\nSince one or more lines are congested, linear programming is required.")
else:
    print("\nNo lines are congested.")

print("Only the constraints corresponding to the congested lines are included in the OPF using LPsolve.")

# Calculated values from LP-solve for part 1 task 1
dispatch = {"1": 0.8, "2": 0, "3": 1.3, "4": 1}
k = 9.1
dispatch_duals   = {"1": 0, "2": 0, "3": 0, "4": 0} # Dual values for dispatch limits set to prod > 0
for index, bus in enumerate(buses.values()):
    bus.gen_cost = gen_cost["{}".format(bus.bus_number)]
    bus.p_gen = dispatch["{}".format(bus.bus_number)]
    bus.marginal_cost = dispatch_duals["{}".format(bus.bus_number)]

print("\nOptimal objective function value, from LPsolve:", k)
for bus_number in dispatch.keys():
    print("Dispatch for bus {} : {} pu".format(int(bus_number)+1, dispatch[bus_number]))

print_title1("Task 2")
#Task 2
#Check the marginal costs (reduced cost, dual variables) and check these against the operating cost at each bus.



"""
#The optimization part with Pyomo (Ignored because LTsolve is used)
print_title1("Optimization part")
# The set for optimization
bus_set = []
line_set = []
loads = {}
for bus in buses.values():
    bus_set.append("Bus {}".format(bus.bus_number))
    loads["Bus {}".format(bus.bus_number)] = bus.p_spec
for line in lines:
    line_set.append(line.name)
# We declare the model (ie. creating a pyo object called model)
model = pyo.ConcreteModel()
"""
"""
Sets:
    - Buses
    - Lines
"""
"""
# We declare the set for buses
model.Bus = pyo.Set(ordered=True, initialize=bus_set)
model.Line = pyo.Set(ordered=True, initialize=line_set)
""""""
Parameters:
    - distribution factors
    - loads
    - generation costs
    - transfer capacities
""""""
dist_factors = {}
line_a = {}
for line in lines:
    for index, bus in enumerate(buses.values()):
        if bus.bus_number != slack_bus_number:
            line_a["Bus {}".format(bus.bus_number)] = a_dict[line.name][index][0]
    dist_factors[line.name] = copy.deepcopy(line_a)  # <-- value må være en ny dict for flere busser
model.a = pyo.Param(model.Line, initialize=dist_factors, within=pyo.Any)
model.loads = pyo.Param(model.Bus, initialize=loads, within=pyo.Reals)
model.gen_cost = pyo.Param(model.Bus, initialize=gen_cost, within=pyo.Reals)
model.trans_cap = pyo.Param(model.Line, initialize=trans_cap, within=pyo.Reals)
""""""
Variables:
    - Variable for generation at a bus
""""""
# Declare the variable for generation at a bus
model.generation = pyo.Var(model.Bus, within=pyo.NonNegativeReals)
# Objective: Minimize cost of generation
def Objective(model):
    return (sum(model.gen_cost[n] * model.generation[n] for n in model.Bus))
model.OBJ = pyo.Objective(rule=Objective, sense=pyo.minimize)
""""""
Constraints:
    production cannot be more than inflow and stored capacity
    reservoir cannot hold more than inflow minus reservoir capacity (rest is spillage)
""""""
# Constraint for power balance
# Constraint for power balance during the first period. Production is equal to inflow - stored for next period - spilled water
def power_balance(model):
    return (sum(model.generation[n] for n in model.Bus) == sum(model.loads[n] for n in model.Bus))
model.power_balance = pyo.Constraint(rule=power_balance)
# Constraint for flow on line
def line_upper_limit(model, line):
    print(line)
    flow = 0
    for index, bus in enumerate(model.Bus):
        if index != slack_bus_number - 1:
            flow += model.a[line][bus] * (model.generation[bus] - model.loads[bus])
    return flow <= model.trans_cap[line]
model.line_upper_limit = pyo.Constraint(model.Line, rule=line_upper_limit)
def line_lower_limit(model, line):
    flow = 0
    for index, bus in enumerate(model.Bus):
        if index != slack_bus_number -1:
            flow += model.a[line][bus] * (model.generation[bus] - model.loads[bus])
    return flow >= -model.trans_cap[line]
    #return (sum(model.a[line][bus] * (model.generation[bus] - model.loads[bus]) for bus in model.Bus[:-1]) >= -model.trans_cap[line])
model.line_lower_limit = pyo.Constraint(model.Line, rule=line_lower_limit)
# Solver used (since quadratic term in objective function, gurobi must be used)
opt = SolverFactory("gurobi")
# Include dual values, the dual values of the constraints will help us find the system price
model.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT)
# Solve the problem
results = opt.solve(model, load_solutions=True)
print("Displaying results:")
# Display values
model.display()
print("Displaying dual results:")
# Display dual values
model.dual.display()
"""