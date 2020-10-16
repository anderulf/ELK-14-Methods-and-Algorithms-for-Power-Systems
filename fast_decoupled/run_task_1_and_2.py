from classes import Bus, Line, Fast_Decoupled
from fast_decoupled.fast_decoupled_methods import run_primal_method, run_dual_method, run_standard_method
import numpy as np
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

#Husk å oppdater continuation til run_fast_decoupled!!!
fast_dec = Fast_Decoupled(buses, slack_bus_number, lines)
fast_dec.set_up_matrices()

# Calculate initial load flow (Point A)
print("")
print("Task 1.")
while fast_dec.power_error() > 0.0001:
    fast_dec.iteration += 1
    fast_dec.reset_values()
#    print("\nIteration: {}\n".format(continuation.iteration))
    fast_dec.calc_new_power_injections()
    fast_dec.error_specified_vs_calculated()
    fast_dec.jacobian.create()
    fast_dec.find_x_diff()
    fast_dec.update_values()
    if fast_dec.diverging():
        print("No convergence")
        break
print("")
print("Base case condition assuming a flat start")
fast_dec.print_matrices()

#2. #Viktig at vi deler opp, oppdaterer P_calc og Q_calc underveis med oppdatert theta og v.
# Ikke løs samtidig (for store system er det mye mindre effektivt å løse samtidig)
print("\n",150*"#")
print("Task 2.")
print("Primal Fast Decoupled Power Flow")

#Reset X-vector to flat start
for bus_number in V:
    buses[int(bus_number)].update_values(P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])
fast_dec = Fast_Decoupled(buses, slack_bus_number, lines)

phase = "Primal"
#Initialize primal jacobian
fast_dec.set_up_matrices(phase)
run_primal_method(fast_dec, printing=True)
fast_dec.print_final_solution(phase)

print("\n", 100*"*")
print("Dual Fast Decoupled Power Flow")
# Reset to flat start
for bus_number in V:
    buses[int(bus_number)].update_values(P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])
fast_dec = Fast_Decoupled(buses, slack_bus_number, lines)

phase = "Dual"
# Initialize primal jacobian (phase)
fast_dec.set_up_matrices(phase)
iteration = 1
run_dual_method(fast_dec, printing=True)
fast_dec.print_final_solution(phase)

print("\n", 100*"*")
print("Standard Decoupled Power Flow")
# Reset to flat start
for bus_number in V:
    buses[int(bus_number)].update_values(P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])
fast_dec = Fast_Decoupled(buses, slack_bus_number, lines)

phase = "Standard"
# Initialize primal jacobian (phase)
fast_dec.set_up_matrices(phase)
run_standard_method(fast_dec, printing=True)
fast_dec.print_final_solution(phase)