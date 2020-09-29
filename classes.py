import numpy as np
import cmath as ma
from supporting_methods import polar_to_rectangular


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

        The limit flag is used to know if the reactive power limit has been reached

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
        self.n_pd = 0
        self.n = 0
        self.m = 0
        self.calculate_n_values()
        self.jacobian = Jacobian(self.n_pq, self.n, self.m, buses, self.y_bus)
        self.total_losses_p = 0
        self.total_losses_q = 0
        self.limit_flag = 0
        self.x_new = np.zeros([self.m, 1])
        self.x_diff = np.zeros([self.m, 1])
        self.x_old = np.zeros([self.m, 1])
        self.x_vector_labels = []
        self.mismatch = Mismatch(self.m)
        self.mismatch_vector_labels = []
        self.net_injections_vector = np.zeros([2*self.n_pq + 2*self.n_pv + 2*self.n_pd, 1]) #Adding 2 for each bus (Pi, Qi)
        self.net_injections_vector_labels = []
        self.correction_vector_labels = []
        self.create_label_vectors()
        self.net_losses_p = 0
        self.net_losses_q = 0
        self.error_history = [] # A list of all load flow errors (maximum of delta_p and delta_q) for every iteration taken
        # Variables for continuation load flow
        self.max_voltage_step = None
        self.max_load_step = None
        self.old_buses_dict = None
        self.continuation_parameter = None

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
                self.n_pd += 1
        if self.n_pd > 1 or self.n_pd == 0:
            print("WARNING INVALID SYSTEM: system has {} slack buses".format(self.n_pd))
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


    def check_limit(self, q_limit, lim_bus, lim_size):
        """
        Check if a certain bus has reached the q_limit, and change the bus if it has
        """
        if q_limit:
            if self.buses_dict[lim_bus].q_calc > lim_size and not self.limit_flag:
                self.limit_flag = 1
                self.buses_dict[lim_bus].q_spec = lim_size
            if self.buses_dict[lim_bus].q_calc < lim_size and self.limit_flag:
                self.limit_flag = 0
                self.buses_dict[lim_bus].q_spec = None
            self.n_pq = 0
            self.n_pv = 0
            self.n_pd = 0
            self.calculate_n_values()

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
            self.net_injections_vector_labels.insert(i - 1 + self.n_pq + self.n_pv + self.n_pd, "Q" + str(i))

            if i == self.slack_bus_number:
                pass
            elif self.buses_dict[i].bus_type == "PV":
                self.mismatch_vector_labels.insert(i-1, "P" + str(i))
                self.correction_vector_labels.insert(i-1, "\u0394\u03B4" + str(i)) # \u0394 = Delta (greek), \u03B4 = delta (greek)
                self.x_vector_labels.insert(i-1, "\u03B4" + str(i)) # \u03B4 = delta (greek)
            else:
                self.mismatch_vector_labels.insert(i-1, "P" + str(i))
                self.mismatch_vector_labels.insert(i -1 + self.n_pq + self.n_pv, "Q" + str(i))

                self.correction_vector_labels.insert(i - 1, "\u0394\u03B4" + str(i))
                self.correction_vector_labels.insert(i -1 + self.n_pq + self.n_pv, "V" + str(i))

                self.x_vector_labels.insert(i-1, "\u03B4" + str(i))
                self.x_vector_labels.insert(i -1 + self.n_pq + self.n_pv, "V" + str(i))

    def print_buses(self):
        """
        Prints the data for all the bus. Mostly needed for debugging
        """
        print(
            "System has {0} slack bus, {1} PQ-bus(es) and {2} PV-bus(es). The jacobian matrix has dimension: {3}x{3}".format(
                self.n_pd, self.n_pq, self.n_pv, self.m))
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
        print(self.jacobian.matrix)
        print("\nNet injections")
        print(np.c_[self.net_injections_vector_labels, self.net_injections_vector])
        print("\nMismatches")
        print(np.c_[self.mismatch_vector_labels, self.mismatch.vector])
        print("\nCorrection vector")
        print(np.c_[self.correction_vector_labels, self.x_new-self.x_old])
        print("\nNew x vector")
        print(np.c_[self.x_vector_labels, self.x_new])
        print("P_total_losses {}".format(self.net_losses_p))
        print("Q_total_losses {}".format(self.net_losses_q))
        print("P_total_losses_2nd_method {}".format(self.total_losses_p))
        print("Q_total_losses_2nd_method {}".format(self.total_losses_q))

class Bus:
    """
    Object holding data for a bus

    Takes OPTIONAL inputs alpha and beta defaulting to None.

    alpha and beta should only be inputed for continuation load flow method
    """
    def __init__(self, bus_number, p_spec, q_spec, voltage, delta, beta=None, alpha=None):
        self.bus_number = bus_number
        self.p_spec = p_spec
        self.q_spec = q_spec
        self.voltage = voltage
        self.delta = delta
        self.p_calc = 0
        self.q_calc = 0
        self.delta_p = 1
        self.delta_q = 1
        self.bus_type = None
        self.classify_bus_type()
        self.beta = beta
        self.alpha = alpha

    def classify_bus_type(self):
        """
        Sets the correct bus type for the bus object
        """
        if self.p_spec and self.q_spec:
            self.bus_type = "PQ"
        elif self.p_spec and not self.q_spec:
            self.bus_type = "PV"
        else:
            self.bus_type = "PD"

    def update_values(self, p_spec, q_spec, voltage, delta):
        self.p_spec = p_spec
        self.q_spec = q_spec
        self.voltage = voltage
        self.delta = delta
        self.delta_p = 1
        self.delta_q = 1

    def print_data(self, slack_bus_number):
        """
        Print the data for a bus
        """
        if self.bus_number == slack_bus_number:
            s = " **** SLACK BUS ****"
        else:
            s = ""
        print(
            "Bus {}: P_spec = {}, Q_spec = {}, voltage = {}, delta = {} deg, P_calc = {}, Q_calc = {}, deltaP = {}, deltaQ = {}".format(
                self.bus_number, self.p_spec, self.q_spec, round(self.voltage, 4), round(self.delta * 180 / np.pi, 4),
                round(self.p_calc, 4), round(self.q_calc, 4), round(self.delta_p, 6), round(self.delta_q, 6)) + s)

class Line:
    def __init__(self, from_bus, to_bus, resistance, reactance):
        self.name = "Line {}-{}".format(from_bus.bus_number, to_bus.bus_number)
        self.from_bus = from_bus
        self.to_bus = to_bus
        self.resistance = resistance
        self.reactance = reactance
        self.impedance = complex(resistance, reactance)
        self.conductance = self.impedance.real
        self.susceptance = self.impedance.imag
        self.p_loss = 0
        self.q_loss = 0
        # to and from currents usually denoted I_ij and I_ji are not necessary equal but opposite due to shunts
        self.to_current = 0
        self.from_current = 0
        self.p_power_flow = 0
        self.q_power_flow = 0


    def __str__(self):
        return self.name

    def print_line_data(self):
        """
        Prints the line data on the format
        Line #-#: Z=R+jX, Y=G+jB

        where Z is the impedance and Y is the admittance
        """
        print(self.name + ": Z={}+j{}, Y = {}+j{}".format(round(self.resistance,2), round(self.reactance,2), round(self.conductance,2), round(self.susceptance, 2)))

class Jacobian:
    """
    The jacobian class holds the matrix, its dimensions and methods for creating and altering the matrix for use
    with the Newton Raphson method and continuation power flow
    """
    def __init__(self,n_pq, n, m, buses_dict, y_bus):
        self.n_pq = n_pq
        self.n = n
        self.m = m
        self.dimensions = [m, m]
        self.cols = m
        self.rows = m
        self.buses_dict = buses_dict
        self.matrix = np.zeros([m, m])
        self.y_bus = y_bus
        self.create()

    def create(self):
        """
        Calculate and return the jacobian matrix for the current iteration.
        The full matrix is constructed from the submatrix parts; J1, J2, J3, J4.
        n = n_pq + n_pv?

        Read:
        The offset is set because a element in the matrix, ie. 0,0 should hold d P_2/d delta_2 if bus 1 is slack. Thus,
        in this example i and j must be offset with 1 because slack bus is bus 3.

        """
        buses = self.buses_dict
        #i_offset = 1
        #j_offset = 0
        for i in range(self.n):
            # J1 of Jacobian
            #if i == self.slack_bus_number-1:
            #    i_offset = 2
            for j in range(self.n):
                #if j == self.slack_bus_number-1:
                #    j_offset = 2
                if i == j:
                    self.matrix[i, j] = -buses[i + 1].q_calc - self.y_bus[i, j].imag * buses[i + 1].voltage * buses[i + 1].voltage
                else:
                    self.matrix[i, j] = abs(buses[i + 1].voltage) * abs(buses[j + 1].voltage) * (self.y_bus[i, j].real * np.sin(buses[i + 1].delta - buses[j + 1].delta) -self.y_bus[i, j].imag * np.cos(buses[i + 1].delta - buses[j + 1].delta))
            #j_offset = 1
            # J2 of Jacobian
            for j in range(self.n_pq):
                #if j == self.slack_bus_number-1:
                #    j_offset = 2
                if i == j:
                    self.matrix[i, j + self.n] = buses[i + 1].p_calc / abs(buses[i + 1].voltage) + \
                                                   self.y_bus[i, i].real * abs(buses[i + 1].voltage)
                else:
                    self.matrix[i, j + self.n] = abs(buses[i + 1].voltage) * (self.y_bus[i, j].real * np.cos(buses[i + 1].delta - buses[j + 1].delta) + self.y_bus[i, j].imag * np.sin(buses[i + 1].delta - buses[j + 1].delta))
            #j_offset = 1
        #i_offset = 1
        #j_offset = 1
        for i in range(self.n_pq):
            # J3 of Jacobian
            #if i == self.slack_bus_number-1:
            #    i_offset = 2
            for j in range(self.n):
                #if j == self.slack_bus_number-1:
                #    j_offset = 2
                if i == j:
                    self.matrix[i + self.n, j] = buses[i + 1].p_calc - self.y_bus[i, i].real * buses[i + 1].voltage * buses[i + 1].voltage
                else:
                    self.matrix[i + self.n, j] = -abs(buses[i + 1].voltage) * abs(buses[j + 1].voltage) * (self.y_bus[i, j].real * np.cos(buses[i + 1].delta - buses[j + 1].delta) + self.y_bus[i, j].imag * np.sin(buses[i + 1].delta - buses[j + 1].delta))
            #j_offset = 1
            for j in range(self.n_pq):
                # J4 of Jacobian
                #if j == self.slack_bus_number-1:
                #    j_offset = 2
                if i == j:
                    self.matrix[i + self.n, j + self.n] = buses[i + 1].q_calc / abs(buses[i + 1].voltage) - self.y_bus[i, i].imag * abs(buses[i + 1].voltage)
                else:
                    self.matrix[i + self.n, j + self.n] = abs(buses[i + 1].voltage) * (self.y_bus[i, j].real * np.sin(buses[i + 1].delta - buses[j + 1].delta) - self.y_bus[i, j].imag * np.cos(buses[i + 1].delta - buses[j + 1].delta))
            #j_offset = 1

    def continuation_expand(self, parameter, buses, constant_voltage_bus_number=None):
        """
        This method is used for Continium Power Flow where the jacobian matrix is expanded with one row and one column

        Parameter should be a string type input with value either "voltage" or "load". "load" input is also used in
        predictor phase.

        buses is the buses dictionary from the load flow class
        """
        self.reset_original_matrix()
        self.cols = self.m + 1
        self.rows = self.m + 1
        new_row = [0] * (self.m +1)
        new_col = [0] * self.m
        for bus in buses.values():
            if bus.bus_type == "PD":
                pass # skip slack bus
            else:
                new_col[bus.bus_number-1] = [bus.beta]
                new_col[bus.bus_number-1+self.n] = [bus.alpha]
        # Add row element which is based on phase
        if parameter == "load":
            new_row[-1] = 1 # index -1 is the last element ie. the diagonal element
        elif parameter == "voltage":
            new_row[constant_voltage_bus_number] = 1 # index -2 is the second last element ie. the last voltage value (in 3 bus system bus 2)
        else:
            print("Error: The parameter \"{}\" in Jacobian.continium_expand is not a valid input.".format(parameter))
            return
        self.matrix = np.hstack([self.matrix, new_col]) # Add the new column
        self.matrix = np.vstack([self.matrix, new_row]) # Add the new row

    def reset_original_matrix(self):
        """
        Remove the last column and row from the altered jacobian matrix to get the original jacobian matrix

        in np.delete obj is the row or column to delete
        """
        if self.cols > self.m and self.rows > self.m:
            # Delete last row
            self.matrix = np.delete(self.matrix, obj=-1, axis=0)  # obj=-1 is the last element, axis=0 means row
            # Delete last col
            self.matrix = np.delete(self.matrix, obj=-1, axis=1)  # obj=-1 is the last element, axis=1 means column
            self.cols = self.m
            self.rows = self.m
        else:
            return

class Mismatch:
    """
    A class for the mismatch vector in load flow methods.

    The class has continium method extensions which makes for a cleaner implementation
    """
    def __init__(self, m):
        self.m = m
        self.rows = m
        self.vector = np.zeros([self.m, 1])

    def continuation_expansion(self, phase):
        """
        expands the vector for continium

        inputs phase which is either "predictor" or "correction"
        """
        self.reset_original_vector()
        self.rows = self.m + 1
        if phase == "predictor":
            self.vector = np.vstack([self.vector, 1])
        elif phase == "correction":
            self.vector = np.vstack([self.vector, 0])
        else:
            self.rows = self.m
            print("Error: phase {} is not defined in mismatch.continium_expansion".format(phase))

    def reset_original_vector(self):
        """
        Removes the added row if it was added
        Do nothing if no rows have been added
        """
        if self.rows > self.m:
            self.vector = np.delete(self.vector, obj=-1, axis=0)  # obj=-1 is the last element, axis=0 means row
            self.rows = self.m
        else: return

class Continuation(Load_Flow):
    """
    The continuation class is a subclass of the Load_Flow class which means Load_Flow (it's super) is
    available for it to use. Ie. inside Continuation you can use self.jacobian, even though jacobian is
    a member variable of the Load_Flow class. Basically everything in Load_Flow is available for
    Continuation, but Continuation is not available for Load_Flow.

    It is implemented as a subclass to have access to these variables and methods, while getting some
    own variables and methods used for the continuation load flow method.

    Note that the class does not have a __init__ because it is done in the super
    """
    def initialize(self, max_voltage_step, max_load_step):
        """
        Initialize some values
        """
        self.max_voltage_step = max_voltage_step
        self.max_load_step = max_load_step
        self.old_buses_dict = self.buses_dict # Stores the last step to be able to reverse
        self.old_mismatch = self.mismatch
        self.continuation_parameter = None
        self.phase = None
        self.step = 1
        self.S = 1

    def initialize_predictor_phase(self):
        """
        Set up the predictor phase

        expand jacobian matrix and mismatch vector
        """
        self.phase = "predictor"
        self.jacobian.continuation_expand("load", self.buses_dict)
        self.mismatch.continuation_expansion(self.phase)

    def initialize_corrector_phase(self, parameter):
        """
        Set up the corrector phase

        input parameter should be "load" or "voltage"
        """
        self.phase = "corrector"
        self.jacobian.continuation_expand(parameter, self.buses_dict)
        self.mismatch.continuation_expansion(self.phase)

    def determine_continuation_parameter(self):
        """
        Determine if NR load flow converges at the current step

        Computational heavy way of determining this. Could perhaps allow higher error limit? Other ways of determining?
        """
        self.jacobian.reset_original_matrix()
        self.mismatch.reset_original_vector()
        convergence = True
        while self.power_error() > 0.0001:
            self.iteration += 1
            self.reset_values()
            self.calc_new_power_injections()
            self.error_specified_vs_calculated()
            self.jacobian.create()
            self.update_values()
            if self.diverging():
                convergence = False
                break
        if convergence:
            self.continuation_parameter = "load"
            # It converged so save the new solution to avoid doing this step again
            self.store_values()
        else:
            self.continuation_flag = "voltage"
            # It did not converge, so don't save the new solution
            self.reverse_step()

        return self.continuation_parameter

    def increment_values(self):
        """
        Increase x and P,Q
        """
        for bus in self.buses_dict.values():
            if bus.bus_number == self.slack_bus_number:
                pass
            else:
                self.bus.p_spec -= self.step * self.bus.beta * self.S
                self.bus.q_spec -= self.step * self.bus.alpha * self.S
                if self.phase == "predictor":
                    self.bus.delta +=  self.x_diff[bus.bus_number-1] * self.step
                    self.bus.voltage += self.x_diff[bus.bus_number - 1 + self.n] * self.step
                else:
                    self.bus.delta += self.x_diff[bus.bus_number - 1]
                    self.bus.voltage += self.x_diff[bus.bus_number - 1 + self.n]

    def constant_voltage_bus(self):
        """
        Return the highest voltage change
        """
        voltage_change_list = []
        for bus_number in self.buses_dict:
            max_voltage_temp = self.buses_dict[bus_number].voltage - self.old_buses_dict[bus_number].voltage
            if max_voltage_temp > max_voltage:
                max_voltage = max_voltage_temp
                max_voltage_bus_number = bus_number;
            else:
                pass
        return max_voltage_bus_number

    def store_values(self):
        """
        If a step was successful store the new values in the old variable
        :return:
        """
        self.old_buses_dict = self.buses_dict
        self.old_mismatch = self.mismatch
        # some other values? mismatch etc.

    def reverse_step(self):
        """
        If a step was unsuccessful revert the values to the last step
        """
        self.buses_dict = self.old_buses_dict
        self.mismatch = self.old_mismatch
        # some other values? mismatch etc.
