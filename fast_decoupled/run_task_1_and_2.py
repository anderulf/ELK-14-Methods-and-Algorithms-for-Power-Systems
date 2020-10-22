from classes import Bus, Line, Fast_Decoupled
from fast_decoupled.fast_decoupled_methods import run_primal_method, run_dual_method, run_standard_method
from supporting_methods import print_title1, print_title2
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
r = {"1-2": 0.1, "1-3": 0.05, "2-3": 0.05}
x = {"1-2": 0.2, "1-3": 0.2 , "2-3": 0.15}

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
print_title1("Task 1")
while fast_dec.power_error() > 0.0001:
    fast_dec.iteration += 1
    fast_dec.reset_values()
#    print("\nIteration: {}\n".format(continuation.iteration))
    fast_dec.calc_new_power_injections()
    fast_dec.error_specified_vs_calculated()
    fast_dec.jacobian.create()
    fast_dec.find_x_diff()
    fast_dec.update_values()
    fast_dec.print_matrices()
    if fast_dec.diverging():
        print("No convergence")
        break

NR_iterations = fast_dec.iteration
print("")
print("Base case condition assuming a flat start")
fast_dec.print_matrices()

#2.
print_title1("Task 2")
print_title2("Primal Fast Decoupled Power Flow")

#Reset X-vector to flat start
for bus_number in V:
    buses[int(bus_number)].update_values(P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])

fast_dec = Fast_Decoupled(buses, slack_bus_number, lines)
#reset calc values
fast_dec.calc_new_power_injections()

phase = "Primal"
#Initialize primal jacobian
fast_dec.set_up_matrices(phase)
fast_dec.print_matrices()
primal_iterations = run_primal_method(fast_dec, printing=True)
fast_dec.print_final_solution(phase)

print_title2("Dual Fast Decoupled Power Flow")
# Reset to flat start
for bus_number in V:
    buses[int(bus_number)].update_values(P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])
fast_dec = Fast_Decoupled(buses, slack_bus_number, lines)
#reset calc values
fast_dec.calc_new_power_injections()

phase = "Dual"
# Initialize primal jacobian (phase)
fast_dec.set_up_matrices(phase)
fast_dec.print_matrices()
dual_iterations = run_dual_method(fast_dec, printing=True)
fast_dec.print_final_solution(phase)

print_title2("Standard Decoupled Power Flow")
# Reset to flat start
for bus_number in V:
    buses[int(bus_number)].update_values(P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])
fast_dec = Fast_Decoupled(buses, slack_bus_number, lines)
#reset calc values
fast_dec.calc_new_power_injections()

phase = "Standard"
# Initialize primal jacobian (phase)
fast_dec.set_up_matrices(phase)
standard_iterations = run_standard_method(fast_dec, printing=True)
fast_dec.print_final_solution(phase)
print("\nNewton raphson iterations: ", NR_iterations)
print("Primal iterations: ", primal_iterations)
print("Dual iterations: ", dual_iterations)
print("Standard iterations: ", standard_iterations)