from classes import NR_Method
import numpy as np

error_limit = 0.001
q_limit = True
limit_bus = 5
limit_size = 1

# Creating ybus
z12 = complex(0.05, 0.25)
z23 = complex(0.02, 0.15)
z25 = complex(0.025, 0.1)
z34 = complex(0.01, 0.1)
z35 = complex(0.02, 0.1)
# Creating off-diagonal elements of ybus
y12 = 1/z12
y23 = 1/z23
y25 = 1/z25
y34 = 1/z34
y35 = 1/z35
# Creating diagonal elements of ybus
y11 = y12
y22 = y12 + y23 + y25
y33 = y34 + y23 + y35
y44 = y34
y55 = y25 + y35
# Creating ybus
ybus = np.matrix([[y11, -y12, 0, 0, 0], [-y12, y22, -y23, 0, -y25], [0, -y23, y33, -y34, -y35], [0, 0, -y34, y44, 0],
                [0, -y25, -y35, 0, y55]])

slack_bus_number = 1
V = {"1": 1, "2": 1, "3": 1, "4": 1, "5": 1}
delta = {"1": 0, "2": 0, "3": 0, "4": 0, "5": 0}
# Q values from project
Q = {"1": None, "2": -0.6, "3": -0.4, "4": -0.2, "5": None}
# P values from project
P = {"1": None, "2": -0.4, "3": -0.8, "4": -0.3, "5": 0.9}

N_R = NR_Method(P, Q, V, delta, ybus, slack_bus_number)

while N_R.maximum_power_injection_error() > error_limit:
    N_R.create_jacobian()
    N_R.calculate_new_power_injections()
    N_R.calculate_error_specified_vs_calculated()
    N_R.print_post_results(True)
    N_R.update_values()
    N_R.iteration += 1
    print("Error is now: {}".format(round(N_R.maximum_power_injection_error(),4)))
    if N_R.iteration > 3:
        print("The system did not converge")
        break

N_R.calculate_line_data()
N_R.calculate_slack_values()
N_R.print_line_data()

