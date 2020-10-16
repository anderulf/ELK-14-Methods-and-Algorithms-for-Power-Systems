from classes import Bus, Line, Fast_Decoupled
from fast_decoupled.fast_decoupled_methods import run_primal_method, run_dual_method
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
x = {"1-2": 0.2, "1-3": 0.1 , "2-3": 0.15}
# 3
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
fast_dec.set_up_matrices()


#3.
print("\n", 150*"#")
print("Task 3.")

# Reset buses
# Create buses
buses = {}
for bus_number in V:
    buses[int(bus_number)] = Bus(int(bus_number), P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])

# # Create lines with new data
line_12 = Line(buses[1], buses[2], r["1-2"], x["1-2"])
line_13 = Line(buses[1], buses[3], r["1-3"], x["1-3"])
line_23 = Line(buses[2], buses[3], r["2-3"], x["2-3"])

lines = [line_12, line_13, line_23]

# Recreate object
fast_dec = Fast_Decoupled(buses, slack_bus_number, lines)
print("Primal Fast Decoupled Power Flow with R = X")
phase = "Primal"

#Initialize primal jacobian
fast_dec.set_up_matrices(phase)
primal_iterations = run_primal_method(fast_dec, printing=True)

print("\n", 100*"*")
print("Dual Fast Decoupled Power Flow with R = X")
# Reset to flat start
for bus_number in V:
    buses[int(bus_number)].update_values(P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])
fast_dec = Fast_Decoupled(buses, slack_bus_number, lines)

phase = "Dual"
# Initialize primal jacobian (phase)
fast_dec.set_up_matrices(phase)

dual_iterations = run_dual_method(fast_dec, printing=True)
print("\nPrimal method iterations: ", primal_iterations)
print("Dual method iterations: ", dual_iterations)
print("\nDual method converges faster for Rij = Xij in bigger systems. In this small 3 bus system the difference is minor")