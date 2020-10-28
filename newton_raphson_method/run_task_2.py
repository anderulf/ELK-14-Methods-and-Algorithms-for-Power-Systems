import matplotlib.pyplot as plt
from classes import Bus, Line
from supporting_methods import print_title1, print_title2, print_title3
from newton_raphson_method.newton_raphson_support import run_newton_raphson, Load_Flow
"""
Settings:

flat_start can be True or False
    True: For load increases the load flow starts from the input values
    False: For load increases the load flow starts from the previous load flow
"""
flat_start = False
# Additional active load 0.2pu 30% at Bus 1, and 70% at Bus 2.
bus_1_load_increase = 0.06
bus_2_load_increase = 0.14

"""
Input values
"""
#  Voltages for 1,2 and the delta are guessed initial values
slack_bus_number = 3
V = {"1": 1, "2": 1, "3": 1}
delta = {"1": 0, "2": 0, "3": 0}
# Specified active load data
P = {"1": -0.8, "2": -0.4, "3": None}
# Specified reactive load data
Q = {"1": -0.5, "2": -0.5, "3": None}
# line data
r = {"1-2": 0.1, "1-3": 0.05, "2-3": 0.05}
x = {"1-2": 0.2, "1-3": 0.25 , "2-3": 0.15}

"""
Program
"""
print_title1("Task 3")
print_title2("NR iteration with load increases. Flat start={}".format(flat_start))

# Create buses
buses = {}
for bus_number in V:
    buses[int(bus_number)] = Bus(int(bus_number), P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])
# Add lines
line_12 = Line(buses[1], buses[2], r["1-2"], x["1-2"])
line_13 = Line(buses[1], buses[3], r["1-3"], x["1-3"])
line_23 = Line(buses[2], buses[3], r["2-3"], x["2-3"])

lines = [line_12, line_13, line_23]

total_iterations = 1
# Assignment 1 Task 2
# initializing vectors for plot
P_increase = []
V_vector_bus1 = []
V_vector_bus2 = []

# Initialize a system object (stores information about the grid)
convergence = True

if flat_start:
    start = True
    V_vector_bus1.append(V["1"])
    V_vector_bus2.append(V["2"])
    P_increase.append(-(P["1"]+P["2"]))
    P["1"] -= bus_1_load_increase
    P["2"] -= bus_2_load_increase
else:
    start = False
# start is used to get the correct order of the list which are used to make plot
# If flat start, then the first elements are initialized earlier in the code.
while convergence:
    # Recreate Bus-objects
    for bus_number in V:
        buses[int(bus_number)].update_values(P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])
    print_title2("Current load: Bus 1 = {} pu, bus 2 = {} pu".format(round(P["1"], 3), round(P["2"], 3)))
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

    if not flat_start:
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
        P["1"] -= bus_1_load_increase
        P["2"] -= bus_2_load_increase

    # initializing for next iteration (while loop)

    N_R.iteration = 1
    start = True

# After divergence
print_title2("Analysis completed")
print("\nTotal iterations: ", total_iterations)
print("Setting for analysis: flat start = {}".format(flat_start))
plt.plot(P_increase,V_vector_bus1,label='V_Bus_1')
plt.plot(P_increase,V_vector_bus2, label='V_Bus_2')
plt.xlabel('Total load power drawn from the system')
plt.ylabel('Voltage [pu]')
plt.legend()
plt.show()