import cmath as math
import numpy as np

from supporting_methods import polar_to_rectangular

class NR_Method:
    """
    The newton raphson method class contains all necessary values and operations associated with the method
    """
    def __init__(self, real_power_dict, reactive_power_dict, voltage_dict, angle_dict, ybus, slack_bus_number):
        self.buses_dict = {}
        self.ybus = ybus
        self.slack_bus_number = slack_bus_number
        self.iteration = 1
        self.fill_buses_dict(real_power_dict, reactive_power_dict, voltage_dict, angle_dict)
        self.n_pq = 0
        self.n_pv = 0
        self.n_pd = 0
        self.calculate_number_of_bus_types()
        self.m = 2 * self.n_pq + self.n_pv
        self.n = self.n_pq + self.n_pv
        self.jacobian = np.zeros([self.m, self.m])
        self.loss_matrix_real_power = None
        self.loss_matrix_reactive_power = None
        self.power_flow_matrix = None
        self.current_matrix = None
        self.total_losses_real_power = 0
        self.total_losses_reactive_power = 0
        self.limit_flag = False

    def fill_buses_dict(self, p_dict, q_dict, voltage_dict, angle_dict):
        """
        Initialize buses_dict based on input data
        """
        for bus_number in p_dict:
            self.buses_dict[int(bus_number)] = Bus(voltage_dict[bus_number],
                                                   angle_dict[bus_number], p_dict[bus_number], q_dict[bus_number], bus_number, self.slack_bus_number)

    def calculate_number_of_bus_types(self):
        """
        Calculates the number of PQ, PD and PV buses. Also calculates the dimensions of the jacobian matrix mxn
        """
        for bus in self.buses_dict:
            if self.buses_dict[bus].real_power_specified and self.buses_dict[bus].reactive_power_specified:
                self.n_pq += 1
            elif self.buses_dict[bus].real_power_specified and not self.buses_dict[bus].reactive_power_specified:
                self.n_pv += 1
            else:
                self.n_pd += 1
        if self.n_pd > 1 or self.n_pd == 0:
            print("Warning: The system is not valid. It has {} slack buses".format(self.n_pd))

    def create_jacobian(self):
        """
        Calculate and return the jacobian matrix for the current iteration.
        The full matrix is constructed from the submatrix parts; J1, J2, J3, J4.
        n = n_pq + n_pv?

        Read:
        The offset is set because a element in the matrix, ie. 0,0 should hold d P_2/d angle_2 if bus 1 is slack. Thus,
        in this example i and j must be offset with 2.
        """
        buses = self.buses_dict
        i_offset = 2
        j_offset = 2
        for i in range(self.n):
            # J1 of Jacobian
            for j in range(self.n):
                if i == j:
                    self.jacobian[i, j] = -buses[i+i_offset].reactive_power_calculated - self.ybus[i+1, i+1].imag * buses[i+i_offset].voltage * buses[i+i_offset].voltage
                else:
                    self.jacobian[i, j] = abs(buses[i+i_offset].voltage) * abs(buses[j+j_offset].voltage) * (
                                self.ybus[i+1, j+1].real * np.sin(buses[i+i_offset].angle - buses[j+j_offset].angle) - self.ybus[i+1, j+1].imag * np.cos(
                        buses[i+i_offset].angle - buses[j+j_offset].angle))
            # J2 of Jacobian
            for j in range(self.n_pq):
                if i == j:
                    self.jacobian[i, j+self.n] = buses[i+i_offset].real_power_calculated / abs(buses[i+i_offset].voltage) + self.ybus[i+1, i+1].real * abs(buses[i+i_offset].voltage)
                else:
                    self.jacobian[i, j+self.n] = abs(buses[i+i_offset].voltage) * (self.ybus[i+1, j+1].real * np.cos(buses[i+i_offset].angle - buses[j+j_offset].angle) + self.ybus[i+1, j+1].imag * np.sin(
                        buses[i+i_offset].angle - buses[j+j_offset].angle))

        for i in range(self.n_pq):
            # J3 of Jacobian
            for j in range(self.n):
                if i == j:
                    self.jacobian[i+self.n, j] = buses[i+i_offset].real_power_calculated - self.ybus[i+1, i+1].real * buses[i+i_offset].voltage * buses[i+i_offset].voltage
                else:
                    self.jacobian[i+self.n, j] = -abs(buses[i+i_offset].voltage) * abs(buses[j+j_offset].voltage) * (self.ybus[i+1, j+1].real * np.cos(buses[i+i_offset].angle - buses[j+j_offset].angle) + self.ybus[i+1, j+1].imag * np.sin(buses[i+i_offset].angle - buses[j+j_offset].angle))

            for j in range(self.n_pq):
                # J4 of Jacobian
                if i == j:
                    self.jacobian[i+self.n, j + self.n] = buses[i+i_offset].reactive_power_calculated / abs(buses[i+i_offset].voltage) - self.ybus[i+1, i+1].imag * abs(buses[i+i_offset].voltage)
                else:
                    self.jacobian[i+self.n, j + self.n] = abs(buses[i+i_offset].voltage) * (self.ybus[i+1, j+1].real * np.sin(buses[i+i_offset].angle - buses[j+j_offset].angle) - self.ybus[i+1, j+1].imag * np.cos(buses[i+i_offset].angle - buses[j+j_offset].angle))


    def calculate_new_power_injections(self):
        """
        Calculate power values based on voltage and delta for all buses except slack.
        The values are saved in the object
        """
        buses = self.buses_dict
        for i in buses:
            buses[i].p_calc = 0  # Resets the value of p_calc/q_calc so the loop works
            buses[i].q_calc = 0
            # Skip slack bus
            if i == self.slack_bus_number:
                pass
            else:
                for j in buses:
                    try:
                        # Adding the Q's from the lines. Node that Ybus is offset with -1 because Python uses 0-indexing and the buses are indexed from 1
                        # Active power
                        buses[i].real_power_calculated += abs(self.ybus[i - 1, j - 1]) * buses[i].voltage * buses[j].voltage * np.cos(
                            buses[i].angle - buses[j].angle - math.phase(self.ybus[i - 1, j - 1]))
                        # Reactive power
                        buses[i].reactive_power_calculated += abs(self.ybus[i - 1, j - 1]) * buses[i].voltage * buses[j].voltage * np.sin(
                            buses[i].angle - buses[j].angle - math.phase(self.ybus[i - 1, j - 1]))
                    except Exception as e:
                        print(e)

    def check_limit(self, q_limit, lim_node, lim_size):
        """
        Check if a certain bus has reached the q_limit, and change the bus if it has
        """
        if q_limit:
            if self.buses_dict[lim_node].reactive_power_calculated > lim_size and not self.limit_flag:
                self.limit_flag = True
                self.buses_dict[lim_node].reactive_power_specified = lim_size
            if self.buses_dict[lim_node].reactive_power_calculated < lim_size and self.limit_flag:
                self.limit_flag = False
                self.buses_dict[lim_node].reactive_power_specified = None
            self.n_pq = 0
            self.n_pv = 0
            self.n_pd = 0
            self.calculate_number_of_bus_types()

    def calculate_error_specified_vs_calculated(self):
        """
        Finds all the error terms and store in object
        """
        for bus_number in self.buses_dict:
            if bus_number == self.slack_bus_number:
                pass
            else:
                # Ignore all buses which doesn't have a specified value (value == None)
                if self.buses_dict[bus_number].real_power_specified:
                    self.buses_dict[bus_number].real_power_delta = self.buses_dict[bus_number].real_power_specified-self.buses_dict[bus_number].real_power_calculated
                if self.buses_dict[bus_number].reactive_power_specified:
                    self.buses_dict[bus_number].reactive_power_delta = self.buses_dict[bus_number].reactive_power_specified - self.buses_dict[bus_number].reactive_power_calculated

    def update_values(self):
        """
        Update the the system values by calculating the next step in the NS method
        """
        x_old = np.zeros([self.m, 1])
        b_matrix = np.zeros([self.m, 1])
        for i in range(self.n):
            x_old[i, 0] = self.buses_dict[i+2].angle
            b_matrix[i, 0] = self.buses_dict[i + 2].real_power_delta
        for i in range(self.n_pq):
            b_matrix[i+self.n, 0] = self.buses_dict[i+2].reactive_power_delta
            x_old[i+self.n, 0] = self.buses_dict[i+2].voltage
        #inv_jac = np.linalg.inv(self.jacobian)
        #x_new = x_old + np.matmul(inv_jac, b_matrix)
        x_new = x_old + np.linalg.solve(self.jacobian, b_matrix) # Faster method compared to inverting
        for i in range(self.n):
            self.buses_dict[i + 2].angle = x_new[i, 0]
        for i in range(self.n_pq):
            self.buses_dict[i + 2].voltage = x_new[i+self.n, 0]

    def maximum_power_injection_error(self):
        """
        Adds all the error terms to a list and returns the maximum value of the list
        """
        error_list = []
        buses = self.buses_dict
        for i in buses:
            # Ignore all buses which doesn't have a specified value (value == None)
            if buses[i].real_power_specified:
                error_list.append(abs(buses[i].real_power_delta))
            if buses[i].reactive_power_specified:
                error_list.append(abs(buses[i].reactive_power_delta))
        return max(error_list)

    def calculate_line_data(self):
        """
        Calculate the losses for all the top right values (skipping diagonals) and store in seperate matrices for active and reactive losses
        """
        buses = self.buses_dict
        rows, cols = self.ybus.shape
        self.loss_matrix_real_power = np.zeros([rows, cols])
        self.loss_matrix_reactive_power = np.zeros([rows, cols])
        self.power_flow_matrix = np.zeros([rows, cols],dtype=np.complex_)
        self.current_matrix = np.zeros([rows, cols],dtype=np.complex_)
        for i in range(rows - 1):
            for j in range(rows - 1, i, -1):
                # Ignore non-existing lines
                if self.ybus[i, j]:
                    # Get rectangular values
                    v_i = polar_to_rectangular(buses[i+1].voltage, buses[i+1].angle)
                    v_j = polar_to_rectangular(buses[j+1].voltage, buses[j+1].angle)
                    # Calculate losses
                    current_ij = self.ybus[i,j] * (v_i - v_j)
                    current_ji = self.ybus[i,j] * (v_j - v_i)
                    self.power_flow_matrix[i,j] = -v_i * current_ij.conjugate() # Negative values because power injections are defined into branch, but this is not intuitive for flow between branches
                    self.power_flow_matrix[j,i] = -v_j * current_ji.conjugate() # Negative values because power injections are defined into branch, but this is not intuitive for flow between branches
                    loss = v_i * current_ij.conjugate() + v_j * current_ji.conjugate()
                    self.current_matrix[i,j] = current_ij
                    self.current_matrix[j,i] = current_ji
                    self.loss_matrix_real_power[i,j] = abs(loss.real)
                    self.loss_matrix_reactive_power[i,j] = abs(loss.imag)
                    # Copy to other side of diagonal just in case
                    self.loss_matrix_real_power[j,i] = self.loss_matrix_real_power[i,j]
                    self.loss_matrix_reactive_power[j,i] = self.loss_matrix_reactive_power[i,j]
                    # Save total losses
                    self.total_losses_real_power += self.loss_matrix_real_power[i,j]
                    self.total_losses_reactive_power += self.loss_matrix_reactive_power[i,j]

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
                self.buses_dict[self.slack_bus_number].real_power_calculated += -self.buses_dict[i].real_power_calculated
                self.buses_dict[self.slack_bus_number].reactive_power_calculated += -self.buses_dict[i].reactive_power_calculated
        # Add losses
        self.buses_dict[self.slack_bus_number].real_power_calculated += self.total_losses_real_power
        self.buses_dict[self.slack_bus_number].reactive_power_calculated += self.total_losses_reactive_power

    def print_line_data(self):
        """
        Print line values (losses, flows, current)
        """
        rows, cols = self.loss_matrix_real_power.shape
        print("Power losses in lines: ")
        for i in range(rows - 1):
            for j in range(rows - 1, i, -1):
                if self.loss_matrix_real_power[i, j]:
                    print("Line {}-{} has I={}, P_flow={}, Q_flow={}, P_loss={} and Q_loss={}".format(i + 1, j + 1,round(self.current_matrix[i, j], 4),round(self.power_flow_matrix[i, j].real,4), round(self.power_flow_matrix[i, j].imag, 4), round(self.loss_matrix_real_power[i, j], 4), round(self.loss_matrix_reactive_power[i, j], 4)))

    def print_post_results(self, jacobian=False):
        """
        Print the relevant data after each iteration
        """
        print(75*"*")
        print("Iteration number: {}\n".format(self.iteration))
        print(
            "The system has {0} slack bus, {1} PQ-bus(es) and {2} PV-bus(es). The jacobian matrix has dimension: {3}x{3}".format(self.n_pd, self.n_pq, self.n_pv, self.m))
        if jacobian:
            print("Jacobian matrix: ")
            print(self.jacobian)
        print("\nBus data:")
        for bus_number in self.buses_dict:
            print(self.buses_dict[bus_number])
        print("\nSystem data:")
        # print something
        print(75*"*")

class Bus:
    """
    The bus class contains all the information of a specific bus
    """
    def __init__(self, voltage, angle, real_power_specified, reactive_power_specified, bus_number, slack_bus_number):
        self.voltage = voltage
        self.angle = angle
        self.real_power_specified = real_power_specified
        self.reactive_power_specified = reactive_power_specified
        self.real_power_calculated = 0
        self.reactive_power_calculated = 0
        self.real_power_delta = 1
        self.reactive_power_delta = 1
        self.bus_number = int(bus_number)
        if self.bus_number == slack_bus_number:
            self.is_slack_bus = True
        else:
            self.is_slack_bus = False

    def __str__(self):
        """
        Overwrite the string type for this class so that a class object can be printed with print-method.
        """
        if self.is_slack_bus:
            s = 10*"*" + " SLACK BUS (bus {}) ".format(self.bus_number) + 10 * "*"
        else:
            s = "Bus {}: voltage = {} pu, angle = {}, P_calc = {}, Q_calc = {}, delta P  = {}, delta Q = {}".format(
                round(self.bus_number,3), round(self.voltage, 3), round(self.angle, 3), round(self.real_power_calculated, 3), round(self.reactive_power_calculated, 3),
                round(self.real_power_delta, 3), round(self.reactive_power_delta, 3))
        return s



