from classes import Bus, Line
import numpy as np
from supporting_methods import print_title1, print_title2, print_title3
from newton_raphson_method.newton_raphson_support import run_newton_raphson
from continuation.continuation_support import Continuation
"""
Continuation settings
"""
max_voltage_step = 0.05
max_load_step = 1

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
# Continuation load increase parameters
beta = {"1": 0.3, "2": 0.7, "3": None}
alpha = {"1": 0, "2": 0, "3": None}

# line data
r = {"1-2": 0.1, "1-3": 0.05, "2-3": 0.05}
x = {"1-2": 0.2, "1-3": 0.25 , "2-3": 0.15}

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
print_title1("Task 1")
run_newton_raphson(continuation, printing=False)
print_title2("Base case condition assuming a flat start")
continuation.print_matrices()
# Start predictions and corrections
#2.
print_title1("Task 2")
print_title2("Predictor phase")
continuation.initialize_predictor_phase()
continuation.reset_values()
continuation.find_x_diff()
print("\nSensitivities/Prediction vector:")
print(np.c_[continuation.correction_vector_labels, np.round(continuation.x_diff, 4)])

#3.
print_title1("Task 3")
continuation.step = 0.3; #1 hos Fosso
continuation.update_continuation_values()
continuation.print_matrices()

continuation.iteration = 0
continuation_parameter = "load"
print_title2("Corrector phase for {} continuation parameter".format(continuation_parameter))
continuation.initialize_corrector_phase(continuation_parameter)
continuation.error_specified_vs_calculated() #needs this one in order to reset the power_error
while continuation.power_error() > 0.0001:
    continuation.iteration += 1
    continuation.reset_values()
    print_title3("Iteration: {}".format(continuation.iteration))
    continuation.calc_new_power_injections()
    continuation.error_specified_vs_calculated()
    continuation.jacobian.reset_original_matrix()
    continuation.jacobian.create()
    continuation.jacobian.continuation_expand(continuation_parameter, continuation.buses_dict)
    continuation.find_x_diff()
    continuation.update_values()
    continuation.print_matrices()
    if continuation.diverging():
        print_title3("No convergence")
        break
continuation.print_matrices()
#4.
print("")
print_title1("Task 4")
print_title2("Predictor phase")
continuation.store_values()
continuation.initialize_predictor_phase()
continuation.reset_values()
continuation.find_x_diff()
continuation.update_continuation_values()
continuation.print_matrices()
constant_bus_index = continuation.constant_voltage_bus()

#5.
print_title1("Task 5")
continuation.iteration = 0
continuation_parameter = "voltage"
print_title2("Corrector phase for {} continuation parameter".format(continuation_parameter))
continuation.step = 1 # Increase correction step so that the iterations converges faster
continuation.initialize_corrector_phase(continuation_parameter, constant_bus_index)
continuation.error_specified_vs_calculated() #needs this one in order to rset the power_error
while continuation.power_error() > 0.0001:
    continuation.iteration += 1
    continuation.reset_values()
    print_title3("Iteration: {}".format(continuation.iteration))
    continuation.calc_new_power_injections()
    continuation.error_specified_vs_calculated()
    continuation.jacobian.reset_original_matrix()
    continuation.jacobian.create()
    continuation.jacobian.continuation_expand(continuation_parameter, continuation.buses_dict, constant_bus_index)
    continuation.find_x_diff()
    continuation.S = continuation.x_diff[-1][0]
    continuation.update_continuation_values()
    continuation.print_matrices()
    #if continuation.diverging():
    if continuation.diverging():
        print_title3("No convergence")
        break


print("\nThis implementation shows the different steps in the continuation load flow method. A fully integrated \n"
      "algorithm is able to make these steps automatically.")