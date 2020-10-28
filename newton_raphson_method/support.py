import numpy as np
import cmath as ma
from supporting_methods import print_title3, polar_to_rectangular
from supporting_classes import Jacobian, Mismatch

# Numpy printing options
np.set_printoptions(suppress=True)  # suppress scientific notations

def run_newton_raphson(N_R, printing=True, total_iterations=None, convergence=False):
    """
    Runs the general newton raphson method
    N_R: Load_Flow, Contination or Fast_Decoupled object
    set printing=False to suppress printing
    Set input total_iterations if external iterations should be counted
    Set convergence=True to track convergence
    """
    if convergence:
        converging = True
    while N_R.power_error() > 0.0001:
        N_R.iteration += 1
        if total_iterations:
            total_iterations += 1
        N_R.reset_values()
        N_R.calc_new_power_injections()
        N_R.error_specified_vs_calculated()
        N_R.jacobian.create()
        N_R.find_x_diff()
        N_R.update_values()
        if printing:
            print_title3("Iteration: {}".format(N_R.iteration))
            N_R.print_matrices()
        if N_R.diverging():
            if printing:
                print_title3("No convergence")
                converging = False
            break
    if total_iterations and convergence:
        return total_iterations, converging
    elif convergence and not total_iterations:
        return converging
    elif total_iterations and not convergence:
        return total_iterations

class Load_Flow:
    def __init__(self, buses, slack_bus_number, lines):
        """
        Initializing the class. p_dict and q_dict should be lists of dictionary types holding key equal to bus number 1, 2, .. and values equal
        their rated values in pu. If the rated active or reactive power is not given set value to None.

        The object calculates the number of PQ, PV, slack buses based on the input data

        For voltage and delta set value to flat start value or initial value

        The input data should have the same length (all buses except slack)

        self.buses_dict: dictionary holding key: bus number value: bus object
            Dictionary type is used because it gives better control over indexing. Now indexing can be done per bus number (key), and not from
            the position in the list

        The line objects in lines contains bus objects, which means that changing the buses will also change the lines.
        Hence updating the lines is not necessary if the buses are updated.

        continuation_flag is used to give further functionality for the continuation power flow method. It is set to
        False by default which means it does not need to be inputed. Setting it to True will activate continuation
        power flow functions which is necessary if doing this approach.
        """
        self.iteration = 0 # Is incremented to one on the first iteration
        self.lines = lines
        self.buses_dict = buses
        self.slack_bus_number = slack_bus_number
        self.create_y_bus()
        self.n_pq = 0
        self.n_pv = 0
        self.n_slack = 0
        self.n = 0
        self.m = 0
        self.calculate_n_values()
        self.jacobian = Jacobian(self.n_pq, self.n, self.m, buses, self.y_bus)
        self.total_losses_p = 0
        self.total_losses_q = 0
        self.x_new = np.zeros([self.m, 1])
        self.x_diff = np.zeros([self.m, 1])
        self.x_old = np.zeros([self.m, 1])
        self.x_vector_labels = []
        self.mismatch_vector_labels = []
        self.net_injections_vector = np.zeros([2 * self.n_pq + 2 * self.n_pv + 2 * self.n_slack, 1]) #Adding 2 for each bus (Pi, Qi)
        self.net_injections_vector_labels = []
        self.correction_vector_labels = []
        self.create_label_vectors()
        self.mismatch = Mismatch(self.m, self.n, self.mismatch_vector_labels)
        self.net_losses_p = 0
        self.net_losses_q = 0
        self.error_history = [] # A list of all load flow errors (maximum of delta_p and delta_q) for every iteration taken
        # Variables for continuation load flow
        self.max_voltage_step = None
        self.max_load_step = None
        self.old_buses_dict = None
        self.continuation_parameter = None
        # Variables for fast decoupled power flow
        self.P_mismatches = np.zeros([self.n, 1])
        self.Q_mismatches = np.zeros([self.n_pq, 1])

    def calculate_n_values(self):
        """
        Calculate number of different bus types
        """
        for bus_number in self.buses_dict: #bus_number is the key of each element in the dictionary
            if self.buses_dict[bus_number].bus_type == "PQ":
                self.n_pq += 1
            elif self.buses_dict[bus_number].bus_type == "PV":
                self.n_pv += 1
            else:
                self.n_slack += 1
        if self.n_slack > 1 or self.n_slack == 0:
            print("WARNING INVALID SYSTEM: system has {} slack buses".format(self.n_slack))
        else:
            self.m = 2 * self.n_pq + self.n_pv
            self.n = self.n_pq + self.n_pv

    def create_y_bus(self):
        """
        Creates the y_bus from the line data
        First the off diagonal elements are added by using the line resistance and reactance
        Then the diagonal elements are found by summing the rows of each diagonal element
        """
        self.y_bus = np.zeros([len(self.buses_dict), len(self.buses_dict)], dtype=complex)
        for line in self.lines:
            self.y_bus[line.from_bus.bus_number -1, line.to_bus.bus_number -1] = -1/complex(line.resistance, line.reactance)
            self.y_bus[line.to_bus.bus_number - 1, line.from_bus.bus_number - 1] = self.y_bus[line.from_bus.bus_number -1, line.to_bus.bus_number -1]
        # Get the sum of the rows
        diagonal_elements = np.sum(self.y_bus, axis = 1) # axis 1 meaning the we sum each colomn along the rows
        for i, Y_ii in enumerate(diagonal_elements):
            self.y_bus[i,i] = -Y_ii # subracting because the off diagonal elements are negative (--=+)

    def calc_new_power_injections(self):
        """
        Calculate power values based on voltage and delta for all buses except slack.
        The values are saved in the object

        Additionally add these values to the net injections vector
        """
        buses = self.buses_dict
        for number, i in enumerate(buses):
            buses[i].p_calc = 0  # Resets the value of p_calc/q_calc so the loop works
            buses[i].q_calc = 0
            # Skip slack bus
            if i == self.slack_bus_number:
                for j in buses:
                    # Adding the Q's from the lines. Note that Ybus is offset with -1 because Python uses 0-indexing and the buses are indexed from 1
                    buses[i].p_calc += abs(self.y_bus[i - 1, j - 1]) * buses[i].voltage * buses[j].voltage * np.cos(
                        buses[i].delta - buses[j].delta - ma.phase(self.y_bus[i - 1, j - 1]))
                    buses[i].q_calc += abs(self.y_bus[i - 1, j - 1]) * buses[i].voltage * buses[j].voltage * np.sin(
                        buses[i].delta - buses[j].delta - ma.phase(self.y_bus[i - 1, j - 1]))
                # Add values to net injection vector
                self.net_injections_vector[i-1] = round(buses[i].p_calc, 3)
                self.net_injections_vector[i + self.n] = round(buses[i].q_calc,3)
                self.net_losses_p += round(buses[i].p_calc, 3)
                self.net_losses_q += round(buses[i].q_calc,3)


            else:
                for j in buses:
                    try:
                        # Adding the Q's from the lines. Note that Ybus is offset with -1 because Python uses 0-indexing and the buses are indexed from 1
                        buses[i].p_calc += abs(self.y_bus[i - 1, j - 1]) * buses[i].voltage * buses[j].voltage * np.cos(
                            buses[i].delta - buses[j].delta - ma.phase(self.y_bus[i - 1, j - 1]))
                        buses[i].q_calc += abs(self.y_bus[i - 1, j - 1]) * buses[i].voltage * buses[j].voltage * np.sin(
                            buses[i].delta - buses[j].delta - ma.phase(self.y_bus[i - 1, j - 1]))
                    except Exception as e:
                        print(e)
                # Add values to net injection vector
                self.net_injections_vector[i-1] = round(buses[i].p_calc, 3)
                self.net_injections_vector[i + self.n] = round(buses[i].q_calc,3)
                self.net_losses_p += round(buses[i].p_calc, 3)
                self.net_losses_q += round(buses[i].q_calc, 3)

    def error_specified_vs_calculated(self):
        """
        Finds all the error terms and store in object
        """
        for bus in self.buses_dict: #bus is equal to the key for each element in the dict.
            if bus == self.slack_bus_number:
                pass
            else:
                # Ignore all buses which doesn't have a specified value (value == None)
                if self.buses_dict[bus].p_spec:
                    self.buses_dict[bus].delta_p = self.buses_dict[bus].p_spec - self.buses_dict[bus].p_calc
                if self.buses_dict[bus].q_spec:
                    self.buses_dict[bus].delta_q = self.buses_dict[bus].q_spec - self.buses_dict[bus].q_calc

    def power_error(self):
        """
        Adds all the error terms to a list and returns the maximum value of the list
        """
        #print("\nPower error:")
        error_list = []
        buses = self.buses_dict
        for i in buses:
            # Ignore all buses which doesn't have a specified value (value == None)
            if not buses[i].p_spec == None:
                error_list.append(abs(buses[i].delta_p))
            if not buses[i].q_spec == None:
                error_list.append(abs(buses[i].delta_q))
        self.error_history.append(max(error_list))
        return max(error_list)

    def diverging(self):
        """
        Determines if the load flow is diverging by looking at the last three errors. If they are increasing or not
        changing the method determines that the system is diverging. Does not compare before the third iteration
        """
        if len(self.error_history) < 3:
            # Continue if less than three iterations has been run
            return False
        else:
            # Check if error increased or was unchanged the last two iterations
            if (self.error_history[-1] >= self.error_history[-2]) and (self.error_history[-2] >= self.error_history[-3]):
                return True
            elif self.iteration > 7:
                # Stop if too many iterations were run
                return True
            else: return False

    def find_x_diff(self):
        for i in range(self.n):
            self.x_old[i, 0] = self.buses_dict[i + 1].delta
            self.mismatch.vector[i, 0] = self.buses_dict[i + 1].delta_p
        for i in range(self.n_pq):
            self.mismatch.vector[i + self.n, 0] = self.buses_dict[i + 1].delta_q
            self.x_old[i + self.n, 0] = self.buses_dict[i + 1].voltage
        self.x_diff = np.linalg.solve(self.jacobian.matrix, self.mismatch.vector)

    def update_values(self):
        """
        Update the the system values by calculating the next step in the NS method

        This method uses a default offset of 1 for the buses
        """
        self.x_new = self.x_old + self.x_diff
        for i in range(self.n):
            self.buses_dict[i + 1].delta = self.x_new[i, 0]
        for i in range(self.n_pq):
            self.buses_dict[i + 1].voltage = self.x_new[i + self.n, 0]

    def calculate_line_data(self):
        """
        Calculate the line data for all the line objects

        Note that writing to line updates self.lines because line and lines[i] is the same object ie. points to same
        location in memory

        """
        for line in self.lines:
            v_from = polar_to_rectangular(line.from_bus.voltage, line.from_bus.delta) # v_i
            v_to = polar_to_rectangular(line.to_bus.voltage, line.to_bus.delta) # v_j
            line.from_current = self.y_bus[line.to_bus.bus_number -1, line.from_bus.bus_number -1] * (v_to - v_from) # i_ij
            line.to_current = self.y_bus[line.from_bus.bus_number - 1, line.to_bus.bus_number - 1] * (v_from - v_to) # i_ji
            apparent_loss = v_from * line.from_current.conjugate() + v_to * line.to_current.conjugate() # v_i * I_ij + v_j * I_ji
            line.p_loss = apparent_loss.real
            line.q_loss = apparent_loss.imag
            line.real_power_flow = (-v_from * line.from_current).real
            line.reactive_power_flow = (-v_from * line.from_current).imag
            self.total_losses_p += line.p_loss
            self.total_losses_q += line.q_loss

    def calculate_slack_values(self):
        """
        Calculate the slack bus values based on the NS iteration

        This function is currently not in use
        """
        for i in self.buses_dict:
            # Skip slack
            if i == self.slack_bus_number:
                pass
            else:
                # Append injections for each bus. Using negative values because slack should cover loads but are covered by other generators
                self.buses_dict[self.slack_bus_number].p_calc += -self.buses_dict[i].p_calc
                self.buses_dict[self.slack_bus_number].q_calc += -self.buses_dict[i].q_calc
        # Add losses
        self.buses_dict[self.slack_bus_number].p_calc += self.total_losses_p
        self.buses_dict[self.slack_bus_number].q_calc += self.total_losses_q
        self.net_injections_vector[self.slack_bus_number - 1] = round(self.buses_dict[self.slack_bus_number].p_calc, 3)
        self.net_injections_vector[self.slack_bus_number + self.n] = round(self.buses_dict[self.slack_bus_number].q_calc, 3)

    def reset_values(self):
        """
        Reset values before new iteration
        """
        self.total_losses_p = 0
        self.total_losses_q = 0

        self.net_losses_p = 0
        self.net_losses_q = 0

    def create_label_vectors(self):
        for i in self.buses_dict:
            self.net_injections_vector_labels.insert(i-1, "P" + str(i))
            self.net_injections_vector_labels.insert(i - 1 + self.n_pq + self.n_pv + self.n_slack, "Q" + str(i))

            if i == self.slack_bus_number:
                pass
            elif self.buses_dict[i].bus_type == "PV":
                self.mismatch_vector_labels.insert(i-1, "\u0394P" + str(i))
                self.correction_vector_labels.insert(i-1, "\u0394\u03B4" + str(i)) # \u0394 = Delta (greek), \u03B4 = delta (greek)
                self.x_vector_labels.insert(i-1, "\u03B4" + str(i)) # \u03B4 = delta (greek)
            else:
                self.mismatch_vector_labels.insert(i-1, "\u0394P" + str(i))
                self.mismatch_vector_labels.insert(i -1 + self.n_pq + self.n_pv, "\u0394Q" + str(i))

                self.correction_vector_labels.insert(i - 1, "\u0394\u03B4" + str(i))
                self.correction_vector_labels.insert(i -1 + self.n_pq + self.n_pv, "\u0394V" + str(i))

                self.x_vector_labels.insert(i-1, "\u03B4" + str(i))
                self.x_vector_labels.insert(i -1 + self.n_pq + self.n_pv, "V" + str(i))

    def print_buses(self):
        """
        Prints the data for all the bus. Mostly needed for debugging
        """
        print(
            "System has {0} slack bus, {1} PQ-bus(es) and {2} PV-bus(es). The jacobian matrix has dimension: {3}x{3}".format(
                self.n_slack, self.n_pq, self.n_pv, self.m))
        for bus_number in self.buses_dict:
            self.buses_dict[bus_number].print_data(self.slack_bus_number)

    def print_line_data(self):
        """
        Print line values (losses, flows, current)
        """
        for line in self.lines:
            print("Line {}-{} has I={}, P_flow={}, Q_flow={}, P_loss={} and Q_loss={}".format(line.from_bus.bus_number, line.to_bus.bus_number, round(line.from_current,3), round(line.real_power_flow,3), round(line.reactive_power_flow,3), round(line.p_loss,3), round(line.q_loss,3)))

    def print_matrices(self):
        print("\nJacobi matrix:")
        print(np.round(self.jacobian.matrix, 4))
        print("\nNet injections")
        print(np.c_[self.net_injections_vector_labels, np.round(self.net_injections_vector, 4)])
        print("\nMismatches")
        print(np.c_[self.mismatch_vector_labels, np.round(self.mismatch.vector, 4)])
        print("\nCorrection vector")
        print(np.c_[self.correction_vector_labels, np.round(self.x_diff, 4)])
        print("\nNew x vector")
        print(np.c_[self.x_vector_labels, np.round(self.x_new, 4)])