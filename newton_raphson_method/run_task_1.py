﻿from supporting_classes import Bus, Line
from supporting_methods import print_title1, print_title2
from newton_raphson_method.support import run_newton_raphson, Load_Flow
"""
Input values
"""
slack_bus_number = 3
V = {"1": 1, "2": 1, "3": 1}
delta = {"1": 0, "2": 0, "3": 0}
# Specified active load data
P = {"1": -1, "2": -0.5, "3": None}
# Specified reactive load data
Q = {"1": -0.5, "2": -0.5, "3": None}
# Line data
r = {"1-2": 0.05, "1-3": 0.05, "2-3": 0.05}
x = {"1-2": 0.2, "1-3": 0.25, "2-3": 0.15}

"""
Program
"""
print_title1("Newton Raphson method iteration")

# Create buses
buses = {}
for bus_number in V:
    buses[int(bus_number)] = Bus(int(bus_number), P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])
# Add lines
line_12 = Line(buses[1], buses[2], r["1-2"], x["1-2"])
line_13 = Line(buses[1], buses[3], r["1-3"], x["1-3"])
line_23 = Line(buses[2], buses[3], r["2-3"], x["2-3"])

lines = [line_12, line_13, line_23]

# Initialize a system object (stores information about the grid)
N_R = Load_Flow(buses, slack_bus_number, lines)

# Iterate NR
run_newton_raphson(N_R)

print_title2("Iteration completed")
print("Iterations: {}".format(N_R.iteration))
# Get post analysis results
N_R.calculate_line_data()
N_R.print_line_data()
N_R.calculate_slack_values()
N_R.print_buses()
print("Total losses: P = {} pu, Q = {} pu".format(round(N_R.total_losses_p, 3), round(N_R.total_losses_q, 3)))