
import numpy as np

class Power_Network:
    """
    The power network class contains all the information about the network, and Bus objects for each bus in the network.
    """
    def __init__(self, real_power_dict, reactive_power_dict, voltage_dict, delta_dict, slack_bus_number):
        self.buses_dict = {}
        self.fill_buses_dict(real_power_dict, reactive_power_dict, voltage_dict, delta_dict)
        self.slack_bus_number = slack_bus_number
        self.n_pq = 0
        self.n_pv = 0
        self.n_pd = 0
        self.n = 0
        self.m = 0
        self.calculate_number_of_bus_types()
        self.loss_matrix_real_power = None
        self.loss_matrix_reactive_power = None
        self.power_flow_matrix = None
        self.current_matrix = None
        self.total_losses_real_power = 0
        self.total_losses_reactive_power = 0
        self.limit_flag = 0

    def fill_buses_dict(self, p_dict, q_dict, voltage_dict, delta_dict):
        """
        Initialize buses_dict based on input data
        """
        for bus_number in p_dict:
            self.buses_dict[int(bus_number)] = Bus(voltage_dict[bus_number],
                                                  delta_dict[bus_number], p_dict[bus_number], q_dict[bus_number], bus_number, slack_bus_number)

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

class NR_Method(Power_Network):
    """
    The newton raphson method class contains all necessary values and operations associated with the method
    It is a subclass of the Power_Network class
    """
    def __init__(self):
        pass

class Jacobian(Power_Network):
    """
    The jacobian class contains the jacobian matrix builder and all values associated with the jacobian matrix
    It is a subclass of the NR_Method class
    """

    def __init__(self):
        self.m = 2 * self.n_pq + self.n_pv
        self.n = self.n_pq + self.n_pv
        self.jacobian = np.zeros([self.m, self.m])

    def create_jacobian(self, y_bus):
        """
        Calculate and return the jacobian matrix for the current iteration.
        The full matrix is constructed from the submatrix parts; J1, J2, J3, J4.
        n = n_pq + n_pv

        Read:
        The offset is set because a element in the matrix, ie. 0,0 should hold P_2/d delta_2 if bus 1 is slack. Thus,
        if bus 1 is slack bus then offset with 2. Ie. one index offset due to the inherit zero-indexing of python, and
        another for the slack bus which is not included in the matrix. With this simple method the algorithm must be
        altered if the slack bus is not bus 1.
        """
        buses = self.buses_dict
        self.jacobian = np.zeros([self.m, self.m])
        i_offset = 2
        j_offset = 2
        for i in range(self.n):
            # Generalized slack bus method proposal
            #if i == self.slack_bus_number:
            #    i_offset = 2
            # J1 of Jacobian
            for j in range(self.n):
                #if j == self.slack_bus_number:
                #    j_offset = 2
                # Diagonal element
                if i == j:
                    self.jacobian[i, j] = -buses[i + i_offset].q_calc - y_bus[i + 1, i + 1].imag * buses[
                        i + i_offset].voltage * buses[i + i_offset].voltage
                # Non-diagonal element
                else:
                    self.jacobian[i, j] = abs(buses[i + i_offset].voltage) * abs(buses[j + j_offset].voltage) * (
                            y_bus[i + 1, j + 1].real * np.sin(buses[i + i_offset].delta - buses[j + j_offset].delta) -
                            y_bus[i + 1, j + 1].imag * np.cos(
                        buses[i + i_offset].delta - buses[j + j_offset].delta))
            #j_offset = 1
            # J2 of Jacobian
            for j in range(self.n_pq):
                # if j == self.slack_bus_number:
                #    j_offset += 1
                # Diagonal element
                if i == j:
                    self.jacobian[i, j + self.n] = buses[i + i_offset].p_calc / abs(buses[i + i_offset].voltage) + \
                                                   y_bus[i + 1, i + 1].real * abs(buses[i + i_offset].voltage)
                # Non-diagonal element
                else:
                    self.jacobian[i, j + self.n] = abs(buses[i + i_offset].voltage) * (
                                y_bus[i + 1, j + 1].real * np.cos(
                        buses[i + i_offset].delta - buses[j + j_offset].delta) + y_bus[i + 1, j + 1].imag * np.sin(
                        buses[i + i_offset].delta - buses[j + j_offset].delta))
        #i_offset = 1
        #j_offset = 1
        for i in range(self.n_pq):
            # J3 of Jacobian
            for j in range(self.n):
                # if j == self.slack_bus_number:
                #    j_offset = 2
                # Diagonal element
                if i == j:
                    self.jacobian[i + self.n, j] = buses[i + i_offset].p_calc - y_bus[i + 1, i + 1].real * buses[
                        i + i_offset].voltage * buses[i + i_offset].voltage
                # Non-diagonal element
                else:
                    self.jacobian[i + self.n, j] = -abs(buses[i + i_offset].voltage) * abs(
                        buses[j + j_offset].voltage) * (y_bus[i + 1, j + 1].real * np.cos(
                        buses[i + i_offset].delta - buses[j + j_offset].delta) + y_bus[i + 1, j + 1].imag * np.sin(
                        buses[i + i_offset].delta - buses[j + j_offset].delta))
            #j_offset = 1
            # J4 of Jacobian
            for j in range(self.n_pq):
                # if j == self.slack_bus_number:
                #    j_offset = 2
                # Diagonal element
                if i == j:
                    self.jacobian[i + self.n, j + self.n] = buses[i + i_offset].q_calc / abs(
                        buses[i + i_offset].voltage) - y_bus[i + 1, i + 1].imag * abs(buses[i + i_offset].voltage)
                # Non-diagonal element
                else:
                    self.jacobian[i + self.n, j + self.n] = abs(buses[i + i_offset].voltage) * (
                                y_bus[i + 1, j + 1].real * np.sin(
                        buses[i + i_offset].delta - buses[j + j_offset].delta) - y_bus[i + 1, j + 1].imag * np.cos(
                        buses[i + i_offset].delta - buses[j + j_offset].delta))


class Bus(Power_Network):
    """
    The bus class contains all the information of a specific bus
    It is a subclass of the Power_Network class
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
        self.bus_number = bus_number
        if bus_number == slack_bus_number:
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
            s = "Bus {}: voltage = {} pu, angle =  {}, P_calc = {}, Q_calc = {}, delta P  = {}, delta Q = {}".format(
                self.bus_number, self.voltage, self.angle, self.real_power_calculated, self.reactive_power_calculated,
                self.real_power_delta, self.reactive_power_delta)
        return s
