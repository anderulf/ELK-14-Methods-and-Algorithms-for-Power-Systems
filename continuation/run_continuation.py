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
beta = {"1": 0.3, "2": 0.7, "3": None}
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

print("Task 1.")
while continuation.power_error() > 0.0001:
    continuation.iteration += 1
    continuation.reset_values()
    print("\nIteration: {}\n".format(continuation.iteration))
    continuation.calc_new_power_injections()
    continuation.error_specified_vs_calculated()
    continuation.jacobian.create()
    continuation.find_x_diff()
    continuation.update_values()
    if continuation.diverging():
        print("No convergence")
        break

# Start predictions and corrections
#2.
print("Task 2.")
print("\n*--- Predictor phase ---*\n")
continuation.initialize_predictor_phase()
continuation.reset_values()
continuation.find_x_diff()

#3.
print("Task 3.")
continuation.step = 0.3;
continuation.update_continuation_values()

print("\n*--- Corrector phase ---*\n")
continuation.iteration = 0
continuation_parameter = "load"
continuation.initialize_corrector_phase(continuation_parameter)
continuation.error_specified_vs_calculated()
while continuation.power_error() > 0.0001:
    continuation.iteration += 1
    continuation.reset_values()
    print("\nIteration: {}\n".format(continuation.iteration))
    continuation.calc_new_power_injections()
    continuation.error_specified_vs_calculated()
    continuation.jacobian.reset_original_matrix()
    continuation.jacobian.create()
    continuation.jacobian.continuation_expand(continuation_parameter, continuation.buses_dict)
    continuation.find_x_diff()
    continuation.update_values()
    if continuation.diverging():
        print("No convergence")
        break
continuation.print_matrices()
#4.
print("Task 4.")
print("\n*--- Predictor phase ---*\n")
continuation.store_values()
continuation.initialize_predictor_phase()
continuation.reset_values()
continuation.find_x_diff()
continuation.update_continuation_values()
continuation.print_matrices()
constant_bus_index = continuation.constant_voltage_bus()

#5.
print("Task 5.")
print("\n*--- Corrector phase ---*\n")
continuation.iteration = 0
continuation_parameter = "voltage"
continuation.initialize_corrector_phase(continuation_parameter, constant_bus_index)
continuation.error_specified_vs_calculated()
while continuation.power_error() > 0.0001:
    continuation.iteration += 1
    continuation.reset_values()
    print("\nIteration: {}\n".format(continuation.iteration))
    continuation.calc_new_power_injections()
    continuation.error_specified_vs_calculated()
    continuation.jacobian.reset_original_matrix()
    continuation.jacobian.create()
    continuation.jacobian.continuation_expand(continuation_parameter, continuation.buses_dict, constant_bus_index)
    continuation.find_x_diff()
    print(continuation.x_diff)
    continuation.S = continuation.x_diff[-1][0]
    continuation.update_continuation_values()
    continuation.print_matrices()
    #if continuation.diverging():
    if continuation.iteration > 100:
        print("No convergence")
        break
