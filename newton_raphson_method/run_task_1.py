from classes import Load_Flow, Bus, Line
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
q_limit = False
lim_node = 3
lim_size = 1

"""
Initial values
"""
#  Voltages for 1,2 and the delta are guessed initial values
slack_bus_number = 3
V = {"1": 1, "2": 1, "3": 1}
delta = {"1": 0, "2": 0, "3": 0}
# Q values from project
Q = {"1": -0.5, "2": -0.5, "3": None}
# P values from project
P = {"1": -1, "2": -0.5, "3": None}
# line data
r = {"1-2": 0.05, "1-3": 0.05, "2-3": 0.05}
x = {"1-2": 0.2, "1-3": 0.1 , "2-3": 0.15}

# Create buses
buses = {}
for bus_number in V:
    buses[int(bus_number)] = Bus(int(bus_number), P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])
# Add lines
line_12 = Line(buses[1], buses[2], r["1-2"], x["1-2"])
line_13 = Line(buses[1], buses[3], r["1-3"], x["1-3"])
line_23 = Line(buses[2], buses[3], r["2-3"], x["2-3"])

lines = [line_12, line_13, line_23]

"""
Program
"""

print("\n*--- Newton Raphson method iteration ---*\n")


# Initialize a system object (stores information about the grid)
N_R = Load_Flow(buses, slack_bus_number, lines)

# Iterate NR
while N_R.power_error() > 0.0001:
    N_R.iteration += 1
    N_R.reset_values()
    print("\nIteration: {}\n".format(N_R.iteration))
    N_R.calc_new_power_injections()
    N_R.check_limit(q_limit, lim_node, lim_size)
    N_R.error_specified_vs_calculated()
    N_R.print_buses()
    N_R.jacobian.create()
    N_R.update_values()
    N_R.print_matrices()
    if N_R.diverging():
        print("No convergence")
        break

print("*--- ITERATION COMPLETED ---*")
print("Iterations: {}".format(N_R.iteration))
# Get post analysis results
N_R.calculate_line_data()
N_R.print_line_data()
N_R.calculate_slack_values()
N_R.print_buses()
print("Total losses: P={}pu, Q={}pu".format(round(N_R.total_losses_p, 5), round(N_R.total_losses_q, 5)))