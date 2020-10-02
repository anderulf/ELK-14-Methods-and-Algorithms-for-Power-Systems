from classes import Bus, Line, Continuation, Load_Flow

"""
Continuation settings
"""
max_voltage_step = 0.05
max_load_step = 1 #0.3 = S? step er hvertfall alltid li 0.3 hos oss da vi forenkler beregning av step

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
# Continuation load increase parameters
beta = {"1": 0.5, "2": 0.5, "3": None}
alpha = {"1": 0, "2": 0, "3": None}

# line data
r = {"1-2": 0.05, "1-3": 0.05, "2-3": 0.05}
x = {"1-2": 0.2, "1-3": 0.1 , "2-3": 0.15}

# Create buses
buses = {}
for bus_number in V:
    buses[int(bus_number)] = Bus(int(bus_number), P[bus_number], Q[bus_number], V[bus_number], delta[bus_number], beta[bus_number], alpha[bus_number])
# Add lines
line_12 = Line(buses[1], buses[2], r["1-2"], x["1-2"])
line_13 = Line(buses[1], buses[3], r["1-3"], x["1-3"])
line_23 = Line(buses[2], buses[3], r["2-3"], x["2-3"])

lines = [line_12, line_13, line_23]

continuation = Continuation(buses, slack_bus_number, lines)
continuation.initialize(max_voltage_step, max_load_step)

# Calculate initial load flow (Point A)

while continuation.power_error() > 0.0001:
    continuation.iteration += 1
    continuation.reset_values()
    print("\nIteration: {}\n".format(continuation.iteration))
    continuation.calc_new_power_injections()
    continuation.error_specified_vs_calculated()
    continuation.print_buses()
    continuation.jacobian.create()
    continuation.update_values()
    continuation.print_matrices()
    if continuation.diverging():
        print("No convergence")
        break

# Start predictions and corrections
#2.
continuation.initialize_predictor_phase()
continuation.reset_values()
continuation.find_x_diff()
continuation.increment_values()
continuation.print_matrices()

#3.
continuation.step=0.3;
continuation.initialize_corrector_phase("load")
while continuation.power_error() > 0.0001:
    continuation.iteration += 1
    continuation.reset_values()
    print("\nIteration: {}\n".format(continuation.iteration))
    continuation.calc_new_power_injections()
    continuation.error_specified_vs_calculated()
    continuation.print_buses()
    continuation.jacobian.create()
    continuation.update_values()
    continuation.print_matrices()
    if continuation.diverging():
        print("No convergence")
        break



continuation.initialize_predictor_phase()




