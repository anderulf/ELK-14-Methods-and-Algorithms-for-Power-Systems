from classes import Load_Flow, Bus, Line
import matplotlib.pyplot as plt
from supporting_methods import run_newton_raphson, print_title1, print_title2
"""
Settings:
    flat_start can be True or False
        True: For load increases the load flow starts from the input values
        False: For load increases the load flow starts from the previous load flow
"""
flat_start = False

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
P = {"1": -0.8, "2": -0.4, "3": None}
# line data
r = {"1-2": 0.1, "1-3": 0.05, "2-3": 0.05}
x = {"1-2": 0.2, "1-3": 0.25 , "2-3": 0.15}

"""
Program
"""
print_title1("Newton Raphson method iteration with load increases and flat start = {}".format(flat_start))

# Create buses
buses = {}
for bus_number in V:
    buses[int(bus_number)] = Bus(int(bus_number), P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])
# Add lines
line_12 = Line(buses[1], buses[2], r["1-2"], x["1-2"])
line_13 = Line(buses[1], buses[3], r["1-3"], x["1-3"])
line_23 = Line(buses[2], buses[3], r["2-3"], x["2-3"])

lines = [line_12, line_13, line_23]

# Initialize iteration counter
total_iterations = 1


# Assignment 1 Task 2
# initializing vectors for plot
Q_increase = []
V_vector_bus1 = []
V_vector_bus2 = []

# Initialize a system object (stores information about the grid)
convergence = True
updating_values = True

if flat_start:
    updating_values = False
    start = 1
    V_vector_bus1.append(V["1"])
    V_vector_bus2.append(V["2"])
    Q_increase.append(-(Q["1"]+Q["2"]))
    Q["1"] -= 0.06
    Q["2"] -= 0.14
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
    total_iterations, convergence = run_newton_raphson(N_R, total_iterations=total_iterations, convergence=True)

    print_title2("Iteration completed")
    print("\nIterations: {}".format(N_R.iteration))
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
        Q_increase.append(-(Q["1"]+Q["2"]))
        Q["1"] -= 0.06
        Q["2"] -= 0.14

    # initializing for next iteration (while loop)

    N_R.iteration = 1
    start = 1

print("Total iterations: ", total_iterations)
plt.plot(Q_increase,V_vector_bus1,label='V_Bus_1')
plt.plot(Q_increase,V_vector_bus2, label='V_Bus_2')
plt.xlabel('Total Reactive power drawn from the system')
plt.ylabel('Voltage [pu]')
plt.legend()
plt.show()
