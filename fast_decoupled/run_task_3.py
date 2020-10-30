from supporting_classes import Bus, Line
from fast_decoupled.support import run_primal_method, run_dual_method, Fast_Decoupled
from supporting_methods import print_title1, print_title2
"""
Initial values
"""
#  Voltages and angles are set to flat start
slack_bus_number = 3
V = {"1": 1, "2": 1, "3": 1}
delta = {"1": 0, "2": 0, "3": 0}
# Specified active load data
P = {"1": -1, "2": -0.5, "3": None}
# Specified reactive load data
Q = {"1": -0.5, "2": -0.5, "3": None}
# Line data
x = {"1-2": 0.2, "1-3": 0.2 , "2-3": 0.15}

# Task 3
# Set Rij = Xij
r = x

# Create buses
buses = {}
for bus_number in V:
    buses[int(bus_number)] = Bus(int(bus_number), P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])
# Add lines
line_12 = Line(buses[1], buses[2], r["1-2"], x["1-2"])
line_13 = Line(buses[1], buses[3], r["1-3"], x["1-3"])
line_23 = Line(buses[2], buses[3], r["2-3"], x["2-3"])

lines = [line_12, line_13, line_23]

fast_dec = Fast_Decoupled(buses, slack_bus_number, lines)

print_title1("Task 3")

print_title2("Primal Fast Decoupled Power Flow with R = X")
phase = "Primal"

#Initialize primal jacobian
fast_dec.set_up_matrices(phase)
fast_dec.print_matrices()
primal_iterations = run_primal_method(fast_dec, printing=True)

print_title2("Dual Fast Decoupled Power Flow with R = X")
# Reset to flat start
for bus_number in V:
    buses[int(bus_number)].update_values(P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])

fast_dec = Fast_Decoupled(buses, slack_bus_number, lines)
# Reset calculated injections
fast_dec.calc_new_power_injections()

phase = "Dual"
fast_dec.set_up_matrices(phase)
fast_dec.print_matrices()

dual_iterations = run_dual_method(fast_dec, printing=True)

print("\nPrimal method iterations: ", primal_iterations)
print("Dual method iterations: ", dual_iterations)
print("\nDual method converges faster for Rij = Xij in larger systems. In this small 3 bus system the difference is minor.")
print("The fast decoupled power flow method generally does not perform well when R/X is high.")