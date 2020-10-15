from classes import Bus, Line, Fast_Decoupled
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
print("")
print("Task 2.")
print("Primal Fast Decoupled Power Flow")

#Reset X-vector to flat start
for bus_number in V:
    buses[int(bus_number)].update_values(P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])
fast_dec = Fast_Decoupled(buses, slack_bus_number, lines)

phase = "Primal"
#Initialize primal jacobian
fast_dec.set_up_matrices(phase)

while fast_dec.power_error() > 0.001:
    # Calculate active power injections
    fast_dec.calculate_P_injections()
    # Calculate the mismatches
    fast_dec.calculate_fast_decoupled_mismatches()
    # Calculate the corrections for angles
    theta_correction = np.linalg.solve(fast_dec.B_p, fast_dec.mismatch.get_P())
    # Update the angles based on the corrections
    #update theta i X-vector, theta_new = theta_correction + theta_old

    # Calculate reactive power injections based on new angles

    #voltage_correction = B_dp.invers*Q_mismatch(v, theta_updated) (samme som over)
    #update voltage i X-vector, V_new = V_correction + v_old

    if fast_dec.diverging():
        print("No convergence")
        break

print("Dual Fast Decoupled Power Flow")
# Reset X-vector til flat start
for bus_number in V:
    buses[int(bus_number)].update_values(P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])
fast_dec = Fast_Decoupled(buses, slack_bus_number, lines)

phase = "Dual"
# Initialize primal jacobian (phase)
fast_dec.set_up_matrices(phase)

while fast_dec.power_error() > 0.0001:
    # Samme som primal men motsatt retning
    # Q_calc(basert på forrige iterasjon oppdaterte x-verdier)
    # voltage_correction = B_dp.invers*Q_mismatch(v, theta)
    # update voltage i X-vector, V_new = V_correction + v_old

    # p_calc (basert på oppdatert v og forrige theta)
    # theta correction = B_p.invers*P_mismatch(v_updated, theta)
    # update theta i X-vector, theta_new = theta_correction + theta_old

    if fast_dec.diverging():
        print("No convergence")
        break

print("Standard Decoupled Power Flow") #Standard er lik primal med en enklere jacobian (helt uavhengig av r_ij)
# Reset X-vector til flat start
for bus_number in V:
    buses[int(bus_number)].update_values(P[bus_number], Q[bus_number], V[bus_number], delta[bus_number])
fast_dec = Fast_Decoupled(buses, slack_bus_number, lines)

phase = "Standard"
# Initialize primal jacobian (phase)
fast_dec.set_up_matrices(phase)

while fast_dec.power_error() > 0.0001:
    # p_calc
    # theta correction = B_p.invers*P_mismatch(v, theta)
    # update theta i X-vector, theta_new = theta_correction + theta_old

    # Q_calc(basert på oppdatert theta og forrige spenning)
    # voltage_correction = B_dp.invers*Q_mismatch(v, theta_updated)
    # update voltage i X-vector, V_new = V_correction + v_old

    if fast_dec.diverging():
        print("No convergence")
        break
#3.
print("")
print("Task 3.")
#Sett Rij=Xij og da vil man se at antall iterasjoner øker litt, ikke heng deg opp i resultatet da vi bare har et tre-bus system.
#For store system skal dual konvergere raskere. Holder å sjekke performance på antall iterasjoner.
#Tipper vi bruker 4-5 iterasjoner, kommer ann på hvor strengt kravet vårt er

#4.
print("")
print("Task 4.")
#Øk last på bus 1 gradvis, og test performance ut fra antall iterasjoner, men igjen ikke heng deg opp i resultatet
# da vi bare har et tre-bus system.
#Tipper vi bruker 4-5 iterasjoner
#Tror vi skal se at dual konvergerer raskere når vi nærmer oss høy belastning for da vil vi ha en sterkere kopling mellom
#P og V og Q og theta








