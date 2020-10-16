from classes import Bus, Line, Fast_Decoupled
from fast_decoupled.fast_decoupled_methods import run_primal_method, run_dual_method
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
iteration = 1
while fast_dec.power_error() > 0.0001:
    # Calculate active power injections
    fast_dec.calculate_P_injections()
    # Calculate the mismatches
    fast_dec.calculate_fast_decoupled_mismatches("P")
    # Calculate the corrections for angles
    theta_correction = np.linalg.solve(fast_dec.B_p, fast_dec.mismatch.get_P())
    # Update theta i X-vector, theta_new = theta_correction + theta_old
    fast_dec.update_fast_decoupled_voltage_or_angle(angles=theta_correction)
    # Calculate reactive power injections based on new angles
    fast_dec.calculate_Q_injections()
    fast_dec.calculate_fast_decoupled_mismatches("Q")
    # Voltage_correction = B_dp.invers*Q_mismatch(v, theta_updated) (samme som over)
    voltage_correction = np.linalg.solve(fast_dec.B_dp, fast_dec.mismatch.get_Q())
    # Update voltage i X-vector, V_new = V_correction + v_old
    fast_dec.update_fast_decoupled_voltage_or_angle(voltages=voltage_correction)
    # Printing
    print("\nIteration ", iteration)
    fast_dec.print_data(theta_correction, voltage_correction)
    if fast_dec.diverging():
        print("No convergence")
        break
    iteration += 1
print(fast_dec.buses_dict[1].p_calc)
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
while fast_dec.power_error() > 0.0001:
    # Calculate reactive power injections
    fast_dec.calculate_Q_injections()
    fast_dec.calculate_fast_decoupled_mismatches("Q")
    # Voltage_correction = B_dp.invers*Q_mismatch(v, theta_updated) (samme som over)
    voltage_correction = np.linalg.solve(fast_dec.B_dp, fast_dec.mismatch.get_Q())
    # Update voltage i X-vector, V_new = V_correction + v_old
    fast_dec.update_fast_decoupled_voltage_or_angle(voltages=voltage_correction)

    # Calculate active power injections based on new voltage corrections
    fast_dec.calculate_P_injections()
    # Calculate the mismatches
    fast_dec.calculate_fast_decoupled_mismatches("P")
    # Calculate the corrections for angles
    theta_correction = np.linalg.solve(fast_dec.B_p, fast_dec.mismatch.get_P())
    # Update theta in X-vector, theta_new = theta_correction + theta_old
    fast_dec.update_fast_decoupled_voltage_or_angle(angles=theta_correction)
    # Printing
    print("\nIteriation ", iteration)
    fast_dec.print_data(theta_correction, voltage_correction)
    if fast_dec.diverging():
        print("No convergence")
        break
    iteration += 1

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
iteration = 1
while fast_dec.power_error() > 0.0001:
    # Calculate active power injections
    fast_dec.calculate_P_injections()
    # Calculate the mismatches
    fast_dec.calculate_fast_decoupled_mismatches("P")
    # Calculate the corrections for angles
    theta_correction = np.linalg.solve(fast_dec.B_p, fast_dec.mismatch.get_P())
    # Update theta i X-vector, theta_new = theta_correction + theta_old
    fast_dec.update_fast_decoupled_voltage_or_angle(angles=theta_correction)
    # Calculate reactive power injections based on new angles
    fast_dec.calculate_Q_injections()
    fast_dec.calculate_fast_decoupled_mismatches("Q")
    # Voltage_correction = B_dp.invers*Q_mismatch(v, theta_updated) (samme som over)
    voltage_correction = np.linalg.solve(fast_dec.B_dp, fast_dec.mismatch.get_Q())
    # Update voltage i X-vector, V_new = V_correction + v_old
    fast_dec.update_fast_decoupled_voltage_or_angle(voltages=voltage_correction)
    # Printing
    print("\nIteration ", iteration)
    fast_dec.print_data(theta_correction, voltage_correction)
    if fast_dec.diverging():
        print("No convergence")
        break
    iteration += 1

fast_dec.print_final_solution(phase)

#3.
print("\n", 150*"#")
print("Task 3.")

# Reset buses
# Create buses
buses = {}
for bus_number in V:
    buses[int(bus_number)] = Bus(int(bus_number), P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])
# Set Rij = Xij
r = x
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
iteration = 1
while fast_dec.power_error() > 0.0001:
    # Calculate active power injections
    fast_dec.calculate_P_injections()
    # Calculate the mismatches
    fast_dec.calculate_fast_decoupled_mismatches("P")
    # Calculate the corrections for angles
    theta_correction = np.linalg.solve(fast_dec.B_p, fast_dec.mismatch.get_P())
    # Update theta i X-vector, theta_new = theta_correction + theta_old
    fast_dec.update_fast_decoupled_voltage_or_angle(angles=theta_correction)
    # Calculate reactive power injections based on new angles
    fast_dec.calculate_Q_injections()
    fast_dec.calculate_fast_decoupled_mismatches("Q")
    # Voltage_correction = B_dp.invers*Q_mismatch(v, theta_updated) (samme som over)
    voltage_correction = np.linalg.solve(fast_dec.B_dp, fast_dec.mismatch.get_Q())
    # Update voltage i X-vector, V_new = V_correction + v_old
    fast_dec.update_fast_decoupled_voltage_or_angle(voltages=voltage_correction)
    if fast_dec.diverging():
        print("No convergence")
        break
    iteration += 1

print("Primal iterations: ", iteration)

print("\n")
print("Dual Fast Decoupled Power Flow with R = X")
# Reset to flat start
for bus_number in V:
    buses[int(bus_number)].update_values(P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])
fast_dec = Fast_Decoupled(buses, slack_bus_number, lines)

phase = "Dual"
# Initialize primal jacobian (phase)
fast_dec.set_up_matrices(phase)
iteration = 1
while fast_dec.power_error() > 0.0001:
    # Calculate reactive power injections
    fast_dec.calculate_Q_injections()
    fast_dec.calculate_fast_decoupled_mismatches("Q")
    # Voltage_correction = B_dp.invers*Q_mismatch(v, theta_updated) (samme som over)
    voltage_correction = np.linalg.solve(fast_dec.B_dp, fast_dec.mismatch.get_Q())
    # Update voltage i X-vector, V_new = V_correction + v_old
    fast_dec.update_fast_decoupled_voltage_or_angle(voltages=voltage_correction)

    # Calculate active power injections based on new voltage corrections
    fast_dec.calculate_P_injections()
    # Calculate the mismatches
    fast_dec.calculate_fast_decoupled_mismatches("P")
    # Calculate the corrections for angles
    theta_correction = np.linalg.solve(fast_dec.B_p, fast_dec.mismatch.get_P())
    # Update theta in X-vector, theta_new = theta_correction + theta_old
    fast_dec.update_fast_decoupled_voltage_or_angle(angles=theta_correction)
    if fast_dec.diverging():
        print("No convergence")
        break
    iteration += 1

print("Dual iterations: ", iteration)
print("\nDual method converges faster for Rij = Xij in bigger systems. In this small 3 bus system the difference is minor")


#Sett Rij=Xij og da vil man se at antall iterasjoner øker litt, ikke heng deg opp i resultatet da vi bare har et tre-bus system.
#For store system skal dual konvergere raskere. Holder å sjekke performance på antall iterasjoner.
#Tipper vi bruker 4-5 iterasjoner, kommer ann på hvor strengt kravet vårt er

#4.
print("\n", 150*"#")
print("Task 4.")

# Reset resistances back to their initial value
r = {"1-2": 0.05, "1-3": 0.05, "2-3": 0.05}

P["1"] -= 2.1

print("Primal Decoupled Power Flow with {}pu load on bus 1".format(abs(P["1"])))

for bus_number in V:
    buses[int(bus_number)] = Bus(int(bus_number), P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])

# # Create lines with new data
line_12 = Line(buses[1], buses[2], r["1-2"], x["1-2"])
line_13 = Line(buses[1], buses[3], r["1-3"], x["1-3"])
line_23 = Line(buses[2], buses[3], r["2-3"], x["2-3"])

lines = [line_12, line_13, line_23]

# Recreate object
fast_dec = Fast_Decoupled(buses, slack_bus_number, lines)
phase = "Primal"
fast_dec.set_up_matrices(phase)

primal_iterations = run_primal_method(fast_dec)
print("Primal iterations: ", primal_iterations)

print("\n", 100*"*")
print("Dual Decoupled Power Flow with {}pu load on bus 1".format(abs(P["1"])))

# Reset to flat start
for bus_number in V:
    buses[int(bus_number)].update_values(P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])

# Recreate object
fast_dec = Fast_Decoupled(buses, slack_bus_number, lines)
phase = "Dual"
fast_dec.set_up_matrices(phase)
dual_iterations = run_dual_method(fast_dec)
print("Dual iterations: ", primal_iterations)

#Øk last på bus 1 gradvis, og test performance ut fra antall iterasjoner, men igjen ikke heng deg opp i resultatet
# da vi bare har et tre-bus system.
#Tipper vi bruker 4-5 iterasjoner
#Tror vi skal se at dual konvergerer raskere når vi nærmer oss høy belastning for da vil vi ha en sterkere kopling mellom
#P og V og Q og theta








