import numpy as np
import cmath as ma


class NR_Method:
    def __init__(self, p_dict, q_dict, voltage_dict, delta_dict, slack_bus_number, y_bus, lines):
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
        """
        self.lines = lines
        self.buses_dict = {}
        self.fill_buses_dict(p_dict, q_dict, voltage_dict, delta_dict)
        self.slack_bus_number = slack_bus_number
        self.y_bus = y_bus
        self.n_pq = 0
        self.n_pv = 0
        self.n_pd = 0
        self.n = 0
        self.m = 0
        self.calculate_n_values()
        self.jacobian = np.zeros([self.m, self.m])
        self.loss_matrix_p = None
        self.loss_matrix_q = None
        self.power_flow_matrix = None
        self.current_matrix = None
        self.total_losses_p = 0
        self.total_losses_q = 0
        self.limit_flag = 0
        self.x_new = np.zeros([self.m, 1])
        self.x_old = np.zeros([self.m, 1])
        self.x_vector_labels = []
        self.diff_b = np.zeros([self.m, 1])
        self.mismatch_vector_labels = []
        self.net_injections_vector = np.zeros([2*self.n_pq + 2*self.n_pv + 2*self.n_pd, 1]) #Adding 2 for each bus (Pi, Qi)
        self.net_injections_vector_labels = []
        self.correction_vector_labels = []
        self.create_label_vectors()


    def fill_buses_dict(self, p_dict, q_dict, voltage_dict, delta_dict):
        """
        Initialize buses_dict based on input data
        """
        for bus_number in p_dict:
            self.buses_dict[int(bus_number)] = Bus(bus_number, p_dict[bus_number], q_dict[bus_number], voltage_dict[bus_number],
                                                 delta_dict[bus_number])

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

    def calc_new_power_injections(self):
        """
        Calculate power values based on voltage and delta for all buses except slack.
        The values are saved in the object

        Additionally add these values to the net injections vector
        """
        buses = self.buses_dict
        self.sum_real_power_injections = 0
        self.sum_ractive_power_injections = 0
        for number, i in enumerate(buses):
            buses[i].p_calc = 0  # Resets the value of p_calc/q_calc so the loop works
            buses[i].q_calc = 0
            # Skip slack bus
            if i == self.slack_bus_number:
                pass
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
        return max(error_list)

    def create_jacobian(self):
        """
        Calculate and return the jacobian matrix for the current iteration.
        The full matrix is constructed from the submatrix parts; J1, J2, J3, J4.
        n = n_pq + n_pv?

        Read:
        The offset is set because a element in the matrix, ie. 0,0 should hold d P_2/d delta_2 if bus 1 is slack. Thus,
        in this example i and j must be offset with 2.

        """
        buses = self.buses_dict
        self.jacobian = np.zeros([self.m, self.m])
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
                    self.jacobian[i, j] = -buses[i + 1].q_calc - self.y_bus[i, j].imag * buses[
                        i + 1].voltage * buses[i + 1].voltage
                else:
                    self.jacobian[i, j] = abs(buses[i + 1].voltage) * abs(buses[j + 1].voltage) * (
                            self.y_bus[i, j].real * np.sin(buses[i + 1].delta - buses[j + 1].delta) -
                            self.y_bus[i, j].imag * np.cos(
                        buses[i + 1].delta - buses[j + 1].delta))
            #j_offset = 1
            # J2 of Jacobian
            for j in range(self.n_pq):
                #if j == self.slack_bus_number-1:
                #    j_offset = 2
                if i == j:
                    self.jacobian[i, j + self.n] = buses[i + 1].p_calc / abs(buses[i + 1].voltage) + \
                                                   self.y_bus[i, i].real * abs(buses[i + 1].voltage)
                else:
                    self.jacobian[i, j + self.n] = abs(buses[i + 1].voltage) * (self.y_bus[i, j].real * np.cos(
                            buses[i + 1].delta - buses[j + 1].delta) + self.y_bus[i, j].imag * np.sin(
                            buses[i + 1].delta - buses[j + 1].delta))
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
                    self.jacobian[i + self.n, j] = buses[i + 1].p_calc - self.y_bus[i, i].real * buses[
                        i + 1].voltage * buses[i + 1].voltage
                else:
                    self.jacobian[i + self.n, j] = -abs(buses[i + 1].voltage) * abs(
                        buses[j + 1].voltage) * (self.y_bus[i, j].real * np.cos(
                        buses[i + 1].delta - buses[j + 1].delta) + self.y_bus[i, j].imag * np.sin(
                        buses[i + 1].delta - buses[j + 1].delta))
            #j_offset = 1
            for j in range(self.n_pq):
                # J4 of Jacobian
                #if j == self.slack_bus_number-1:
                #    j_offset = 2
                if i == j:
                    self.jacobian[i + self.n, j + self.n] = buses[i + 1].q_calc / abs(
                        buses[i + 1].voltage) - self.y_bus[i, i].imag * abs(buses[i + 1].voltage)
                else:
                    self.jacobian[i + self.n, j + self.n] = abs(buses[i + 1].voltage) * (
                                self.y_bus[i, j].real * np.sin(
                            buses[i + 1].delta - buses[j + 1].delta) - self.y_bus[i, j].imag * np.cos(
                            buses[i + 1].delta - buses[j + 1].delta))
            #j_offset = 1

    def update_values(self):
        """
        Update the the system values by calculating the next step in the NS method

        This method uses a default offset of 1 for the buses
        """
        for i in range(self.n):
            self.x_old[i, 0] = self.buses_dict[i + 1].delta
            self.diff_b[i, 0] = self.buses_dict[i + 1].delta_p
        for i in range(self.n_pq):
            self.diff_b[i + self.n, 0] = self.buses_dict[i + 1].delta_q
            self.x_old[i + self.n, 0] = self.buses_dict[i + 1].voltage
        self.x_new = self.x_old + np.linalg.solve(self.jacobian, self.diff_b)
        for i in range(self.n):
            self.buses_dict[i + 1].delta = self.x_new[i, 0]
        for i in range(self.n_pq):
            self.buses_dict[i + 1].voltage = self.x_new[i + self.n, 0]

    def calculate_line_data(self):
        """
        Calculate the line data for all the line objects
        """
        for line in self.lines:
            v_from = polar_to_rectangular(line.from_bus.voltage, line.to_bus.delta)
            v_to = polar_to_rectangular(line.to_bus.voltage, line.to_bus.delta)
            line.from_current = self.y_bus[line.to_bus.bus_number -1, line.from_bus.bus_number -1] * (v_to - v_from)
            line.to_current = self.y_bus[line.from_bus.bus_number - 1, line.to_bus.bus_number - 1] * (v_from - v_to)
            apparent_loss = v_from * line.from_current.conjugate() + v_to * line.to_current.conjugate() # v_i * I_ij + v_j * I_ji
            line.p_loss = abs(apparent_loss.real)
            line.q_loss = apparent_loss.imag # shunts can supply reactive power and increase losses
            line.real_power_flow = (-v_from * line.from_current).real
            line.reactive_power_flow = (-v_from * line.from_current).imag
            self.total_losses_p += line.p_loss
            self.total_losses_q += line.q_loss

    def calculate_slack_values(self):
        """
        Calculate the slack bus values based on the NS iteration
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

        #prints when debugging
        print(*self.net_injections_vector_labels, sep = '\n')
        print(self.mismatch_vector_labels)
        print(self.correction_vector_labels)
        print(self.x_vector_labels)





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
            print("Line {}-{} has I={}, P_flow={}, Q_flow={}, P_loss={} and Q_loss={}".format(line.from_bus.bus_number, line.to_bus.bus_number, line.from_current, line.real_power_flow, line.reactive_power_flow, line.p_loss, line.q_loss))

    def print_matrices(self):
        print("\nJacobi matrix:")
        print(self.jacobian)
        print("\nNet injections")
        print(np.c_[self.net_injections_vector_labels, self.net_injections_vector])
        print("\nMismatches")
        print(np.c_[self.mismatch_vector_labels, self.diff_b])
        print("\nCorrection vector")
        print(np.c_[self.correction_vector_labels, self.x_new-self.x_old])
        print("\nNew x vector")
        print(np.c_[self.x_vector_labels, self.x_new])

class Bus:
    """
    Object holding data for a bus
    """

    def __init__(self, bus_number, p_spec, q_spec, voltage, delta):
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
        self.name = "Line {}-{}".format(to_bus.bus_number, to_bus.bus_number)
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
        self.real_power_flow = 0
        self.reactive_power_flow = 0


    def __str__(self):
        return self.name

    def print_line_data(self):
        """
        Prints the line data on the format
        Line #-#: Z=R+jX, Y=G+jB

        where Z is the impedance and Y is the admittance
        """
        print(self.name + ": Z={}+j{}, Y = {}+j{}".format(round(self.resistance,2), round(self.reactance,2), round(self.conductance,2), round(self.susceptance, 2)))

"""
Auxillary methods
"""
def rectangular_to_polar(complex_number):
    r = np.sqrt(complex_number.real * complex_number.real + complex_number.imag * complex_number.imag).real
    if complex_number.real < 0 and complex_number.imag > 0:  # second quadrant
        angle = np.pi - np.arctan(abs(complex_number.imag / complex_number.real))
    elif complex_number.real < 0 and complex_number.imag < 0:  # third quadrant
        angle = -(np.pi - np.arctan(abs(complex_number.imag / complex_number.real)))
    elif complex_number.real > 0 and complex_number.imag < 0:  # fourth quadrant
        angle = -np.arctan(abs(complex_number.imag / complex_number.real))
    else:  # first quadrant
        angle = np.arctan(abs(complex_number.imag / complex_number.real))
    return r, angle


def polar_to_rectangular(abs, angle):
    a = abs * np.cos(angle)
    b = abs * np.sin(angle)
    return complex(a, b)


def complex_angle(complex_number):
    r, angle = rectangular_to_polar(complex_number)
    return angle


def complex_radius(complex_number):
    r, angle = rectangular_to_polar(complex_number)
    return r