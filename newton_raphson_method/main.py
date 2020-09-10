﻿import Ybus as ad
from classes import NR_Method
import numpy as np
# heisann
"""
Settings:
    outage  can be set to
    0: no outage
    1: line 2-3
    2: line 2-5
    3: line 3-5

Q_limit, lim_node and lim_size decides if there is a limit, which node it applies to and the limit size
Only works for one node.
"""
q_limit = True
lim_node = 5
lim_size = 1

"""
Initial values
"""
#  Voltages for 2,3,4 and the delta are guess initial values
slack_node = 1
V = {"1": 1, "2": 1, "3": 1, "4": 1, "5": 1}
delta = {"1": 0, "2": 0, "3": 0, "4": 0, "5": 0}
# Q values from project
Q = {"1": None, "2": -0.6, "3": -0.4, "4": -0.2, "5": None}
# P values from project
P = {"1": None, "2": -0.4, "3": -0.8, "4": -0.3, "5": 0.9}

"""
Program
"""

print("\n*--- Newton Raphson method iteration ---*\n")


y_bus = ad.y_bus

# Initialize iteration counter
iter = 0

# Initialize a system object (stores information about the grid)
N_R = NR_Method(P, Q, V, delta, slack_node, y_bus)

# Iterate NS
while N_R.power_error() > 0.0001:
    print("\nIteration: {}\n".format(iter))
    N_R.calc_new_power()
    N_R.check_limit(q_limit, lim_node, lim_size)
    N_R.error_specified_vs_calculated()
    N_R.print_nodes()
    N_R.create_jacobian()
    print(N_R.jacobian)
    N_R.update_values()
    iter += 1
    if iter > 1000:
        print("No convergence")
        break

print("*--- ITERATION COMPLETED ---*")
print("Iterations: {}".format(iter))

# Get post analysis results
N_R.calculate_line_data()
N_R.print_line_data()
N_R.calculate_slack_values()
N_R.print_nodes()
print("Total losses: P={}pu, Q={}pu".format(round(N_R.total_losses_p, 5), round(N_R.total_losses_q, 5)))