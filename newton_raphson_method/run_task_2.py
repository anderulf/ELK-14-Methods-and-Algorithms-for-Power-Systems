
from classes import Load_Flow, Bus, Line
import matplotlib.pyplot as plt

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
updating_values = 1
flat_start = 1
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
P = {"1": -0.8, "2": -0.4, "3": None}
# line data
r = {"1-2": 0.1, "1-3": 0.05, "2-3": 0.05}
x = {"1-2": 0.2, "1-3": 0.25 , "2-3": 0.15}

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

print("\n*--- Newton Raphson method iteration with load increases ---*\n")


# Initialize iteration counter
iter = 1


# Assignment 1 Task 2
# initializing vectors for plot
P_increase = []
V_vector_bus1 = []
V_vector_bus2 = []

# Initialize a system object (stores information about the grid)
convergence = 1
if flat_start:
    start = 1
    V_vector_bus1.append(V["1"])
    V_vector_bus2.append(V["2"])
    P_increase.append(-(P["1"]+P["2"]))
    P["1"] -= 0.06
    P["2"] -= 0.14
else:
    start = 0
# start is used to get the correct order of the list which are used to make plot
# If flat start, then the first elements are initialized earlier in the code.
while convergence:
    # Recreate Bus-objects
    for bus_number in V:
        buses[int(bus_number)].update_values(P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])
    N_R = Load_Flow(buses, slack_bus_number, lines)
    # Iterate NS
    while N_R.power_error() > 0.0001:
        print("\nIteration: {}\n".format(iter))
        N_R.calc_new_power_injections()
        N_R.check_limit(q_limit, lim_node, lim_size)
        N_R.error_specified_vs_calculated()
        N_R.print_buses()
        N_R.jacobian.create()
        N_R.update_values()
        N_R.print_matrices()
        iter += 1
        if N_R.diverging():
            print("No convergence")
            convergence = 0
            break

    print("*--- ITERATION COMPLETED ---*")
    print("Iterations: {}".format(iter))
    # Get post analysis results
    N_R.calculate_line_data()
    N_R.print_line_data()
    N_R.calculate_slack_values()
    N_R.print_buses()
    print("Total losses: P={}pu, Q={}pu".format(round(N_R.total_losses_p, 5), round(N_R.total_losses_q, 5)))

    if updating_values:
        V["1"] = N_R.buses_dict[1].voltage
        V["2"] = N_R.buses_dict[2].voltage
        delta["1"] = N_R.buses_dict[1].delta
        delta["2"] = N_R.buses_dict[2].delta
        Q["1"] = N_R.buses_dict[1].q_spec
        Q["2"] = N_R.buses_dict[2].q_spec
        P["1"] = N_R.buses_dict[1].p_spec
        P["2"] = N_R.buses_dict[2].p_spec

    if convergence and start:
        # Adding additional load 0.2pu 30% at Bus 1, and 70% at Bus 2.
        V_vector_bus1.append(N_R.buses_dict[1].voltage)
        V_vector_bus2.append(N_R.buses_dict[2].voltage)
        P_increase.append(-(P["1"]+P["2"]))
        P["1"] -= 0.06
        P["2"] -= 0.14

    # initializing for next iteration (while loop)

    iter = 1
    start = 1


plt.plot(P_increase,V_vector_bus1,label='V_Bus_1')
plt.plot(P_increase,V_vector_bus2, label='V_Bus_2')
plt.xlabel('Total load power drawn from the system')
plt.ylabel('Voltage [pu]')
plt.legend()
plt.show()
