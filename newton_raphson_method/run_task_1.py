from classes import Load_Flow, Bus, Line
from supporting_methods import print_title1, print_title2, print_title3
"""
Input values
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
while N_R.power_error() > 0.0001:
    N_R.iteration += 1
    N_R.reset_values()
    print("\nIteration: {}\n".format(N_R.iteration))
    N_R.calc_new_power_injections()
    N_R.error_specified_vs_calculated()
    N_R.print_buses()
    N_R.jacobian.create()
    N_R.find_x_diff()
    N_R.update_values()
    N_R.print_matrices()
    if N_R.diverging():
        print("No convergence")
        break

print_title2("Iteration completed")
print("Iterations: {}".format(N_R.iteration))
# Get post analysis results
N_R.calculate_line_data()
N_R.print_line_data()
N_R.calculate_slack_values()
N_R.print_buses()
print("Total losses: P={}pu, Q={}pu".format(round(N_R.total_losses_p, 5), round(N_R.total_losses_q, 5)))