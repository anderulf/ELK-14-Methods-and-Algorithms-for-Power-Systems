from classes import Bus, Line
from Distribution_factors_and_IMML.IMML_algorithm import IMML_algorithm

"""
Settings
"""
outage_task_2 = "1-2"
outage_task_3 = "1-3"
"""
Input values
"""
slack_bus_number = 4
V = {"1": 1, "2": 1, "3": 1, "4": 1 }
delta = {"1": 0, "2": 0, "3": 0, "4": 0}
# Q values from assignment
Q = {"1": 0, "2": 0, "3": 0, "4": None}
# P values from project
P = {"1": -1.25, "2": -0.4, "3": -0.6, "4": None}
# line data
r = {"1-2": 0.0, "1-3": 0.0, "2-3": 0.0, "3-4": 0}
x = {"1-2": 0.2, "1-3": 0.1, "2-3": 0.25, "3-4": 0.25}

# Create buses
buses = {}
for bus_number in V:
    buses[int(bus_number)] = Bus(int(bus_number), P[bus_number], Q[bus_number], V[bus_number], delta[bus_number] )
# Add lines
line_12 = Line(buses[1], buses[2], r["1-2"], x["1-2"])
line_13 = Line(buses[1], buses[3], r["1-3"], x["1-3"])
line_23 = Line(buses[2], buses[3], r["2-3"], x["2-3"])
line_34 =Line(buses[3], buses[4], r["3-4"], x["3-4"])

lines = [line_12, line_13, line_23, line_34]

"""
Program
"""
#Part 2
print("\nTask 1.")
print("\nIMML is a fast and general technique for simulating impacts of modifying a network topolgy, such as outages and compensations. "
      "The method is based on creating three sub-matrices, H, delta_h and M, which will be constant during the entire process. "
      "H describes the original network topology, delta_h reflects what to me manipulated and M reflects which lines/buses to be affected by the "
      "manipulation. The results/solutions after applying an IMML are the state of the new voltage angles as well as the new power flow for the given modification\n")

print("\nTask 2.")
#Find voltage angles and the power flow when the line 1-2 is disconnected by using the IMML
from_bus, to_bus = IMML_algorithm(P, buses, lines, slack_bus_number, outage_task_2)
for line in lines:
    if from_bus == line.from_bus.bus_number and to_bus == line.to_bus.bus_number:
        pass
    else:
        print("P on", line.name, ":", (line.from_bus.delta - line.to_bus.delta)/line.reactance)

print("\nTask 3.")

from_bus, to_bus  = IMML_algorithm(P, buses, lines, slack_bus_number, outage_task_3, h_modification=0.5)
# Removal of one line means the equivalent impedance on the remaining line is doubled
line_13.reactance *= 2
for line in lines:
    print("P on", line.name, ":", (line.from_bus.delta - line.to_bus.delta)/line.reactance)


