from supporting_classes import Bus, Line
from fast_decoupled.support import run_primal_method, run_dual_method, Fast_Decoupled
from supporting_methods import print_title1, print_title2
"""
Settings
"""

load_increase_at_bus_1 = 0.5

"""
Initial values
"""
#  Voltages for 1,2 and the delta are guessed initial values
slack_bus_number = 3
V = {"1": 1, "2": 1, "3": 1}
delta = {"1": 0, "2": 0, "3": 0}
# Specified active load data
P = {"1": -1, "2": -0.5, "3": None}
# Specified reactive load data
Q = {"1": -0.5, "2": -0.5, "3": None}
# line data
r = {"1-2": 0.05, "1-3": 0.05, "2-3": 0.05}
x = {"1-2": 0.2, "1-3": 0.1 , "2-3": 0.15}

# Load is increased on the bus
P["1"] -= load_increase_at_bus_1

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

# Reset calculated injections
fast_dec.calc_new_power_injections()

fast_dec.set_up_matrices()

#4.
print_title1("Task 4")

print_title2("Primal Decoupled Power Flow with {} pu load on bus 1".format(abs(P["1"])))

for bus_number in V:
    buses[int(bus_number)] = Bus(int(bus_number), P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])

# # Create lines with new data
line_12 = Line(buses[1], buses[2], r["1-2"], x["1-2"])
line_13 = Line(buses[1], buses[3], r["1-3"], x["1-3"])
line_23 = Line(buses[2], buses[3], r["2-3"], x["2-3"])

lines = [line_12, line_13, line_23]

# Recreate Fast Decoupled object
fast_dec = Fast_Decoupled(buses, slack_bus_number, lines)
phase = "Primal"
fast_dec.set_up_matrices(phase)

primal_iterations = run_primal_method(fast_dec, printing=True)

print_title2("Dual Decoupled Power Flow with {} pu load on bus 1".format(abs(P["1"])))

# Reset to flat start
for bus_number in V:
    buses[int(bus_number)].update_values(P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])

# Recreate object
fast_dec = Fast_Decoupled(buses, slack_bus_number, lines)
phase = "Dual"
fast_dec.set_up_matrices(phase)
dual_iterations = run_dual_method(fast_dec, printing=True)

print("\nPrimal method iterations: ", primal_iterations)
print("Dual method iterations: ", primal_iterations)
print("\nPrimal and dual methods uses same amount of iterations but this is not likely for bigger systems. In bigger \n"
      "systems the dual method is preferred.")
