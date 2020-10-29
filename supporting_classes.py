import numpy as np
# Numpy printing options
np.set_printoptions(suppress=True)  # suppress scientific notations

class Bus:
    """
    Object holding data for a bus

    Takes OPTIONAL inputs alpha and beta defaulting to None.

    alpha and beta should only be inputed for continuation load flow method
    """
    def __init__(self, bus_number, p_spec, q_spec, voltage, delta, beta=None, alpha=None, dispatch=None):
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
        # Variables for continuation
        self.beta = beta
        self.alpha = alpha
        # Variables for benders decomposition
        self.dispatch = dispatch
        self.p_gen = None
        self.gen_cost = None
        self.sensitivity = None

    def classify_bus_type(self):
        """
        Sets the correct bus type for the bus object
        """
        if self.p_spec and self.q_spec:
            self.bus_type = "PQ"
        elif self.p_spec and not self.q_spec:
            self.bus_type = "PV"
        else:
            self.bus_type = "Slack"

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
        if self.p_spec:
            p_spec = round(self.p_spec, 4)
        else:
            p_spec = None
        if self.q_spec:
            q_spec = round(self.q_spec, 4)
        else:
            q_spec = None
        print(
            "Bus {}: P_spec = {}pu, Q_spec = {}pu, voltage = {}pu, delta = {}\u00B0, P_calc = {}pu, Q_calc = {}pu, deltaP = {}pu, deltaQ = {}pu".format(
                self.bus_number, p_spec, q_spec, round(self.voltage, 4), round(self.delta * 180 / np.pi, 4),
                f"{self.p_calc:.4f}", f"{self.q_calc:.4f}", f"{self.delta_p:.4f}", f"{self.delta_q:.4f}") + s)

class Line:
    def __init__(self, from_bus, to_bus, resistance, reactance, transfer_capacity=None):
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
        self.transfer_capacity = transfer_capacity


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

    def create(self, fast_decoupled=False):
        """
        Calculate and return the jacobian matrix for the current iteration.
        The full matrix is constructed from the submatrix parts; J1, J2, J3, J4.
        n = n_pq + n_pv?

        Read:
        The offset is set because a element in the matrix, ie. 0,0 should hold d P_2/d delta_2 if bus 1 is slack. Thus,
        in this example i and j must be offset with 1 because slack bus is bus 3.

        fast_decoupled should be set to True if the fast decoupled power flow method is used. It skips J2 and J3 jacobi
        submatrices.
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
            if not fast_decoupled:
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
            if not fast_decoupled:
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

    def continuation_expand(self, parameter, buses, constant_voltage_bus_index=None):
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
            if bus.bus_type == "Slack":
                pass # skip slack bus
            else:
                new_col[bus.bus_number-1] = [bus.beta]
                new_col[bus.bus_number-1+self.n] = [bus.alpha]
        # Add row element which is based on phase
        if parameter == "load":
            new_row[-1] = 1 # index -1 is the last element ie. the diagonal element
        elif parameter == "voltage":
            new_row[constant_voltage_bus_index] = 1 # index -2 is the second last element ie. the last voltage value (in 3 bus system bus 2)
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

    The class has continuation method extensions which makes for a cleaner implementation

    Also supports extracting only active or reactive parts for decoupled power flow implemenations
    """
    def __init__(self, m, n, labels):
        self.m = m
        self.n = n
        self.rows = m
        self.vector = np.zeros([self.m, 1])
        self.labels = labels

    def continuation_expansion(self, phase):
        """
        expands the vector for continium

        inputs phase which is either "predictor" or "correction"
        """
        self.reset_original_vector()
        self.rows = self.m + 1
        if phase == "predictor":
            self.vector = np.vstack([self.vector, 1])
        elif phase == "corrector":
            self.vector = np.vstack([self.vector, 0])
        else:
            self.rows = self.m
            print("Error: phase {} is not defined in mismatch.continuation_expansion".format(phase))

    def reset_original_vector(self):
        """
        Removes the added row if it was added
        Do nothing if no rows have been added
        """
        if self.rows > self.m:
            self.vector = np.delete(self.vector, obj=-1, axis=0)  # obj=-1 is the last element, axis=0 means row
            self.rows = self.m
        else: return

    def get_P(self):
        """
        The active mismatches only for fast decoupled power flow
        """
        return self.vector[:self.n]

    def get_Q(self):
        """
        The reactive mismatches only for fast decoupled power flow
        """
        return self.vector[self.n:]

    def get_P_label(self):
        """
        The active mismatch labels only for fast decoupled power flow
        """
        return self.labels[:self.n]

    def get_Q_label(self):
        """
        The reactive mismatch labels only for fast decoupled power flow
        """
        return self.labels[self.n:]
