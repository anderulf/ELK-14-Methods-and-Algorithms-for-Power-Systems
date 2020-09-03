
import numpy as np

class Power_Network:
    """
    The power network class contains all the information about the network, and Bus objects for each bus in the network.
    """
    def __init__(self, real_power_dict, reactive_power_dict, voltage_dict, angle_dict, slack_bus_number):
        self.slack_bus_number = slack_bus_number
        self.buses_dict = {}
        self.fill_buses_dict(real_power_dict, reactive_power_dict, voltage_dict, angle_dict)
        self.n_pq = 0
        self.n_pv = 0
        self.n_pd = 0
        self.n = 0
        self.m = 0
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
        self.limit_flag = 0

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


    def create_jacobian(self, y_bus):
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
                    self.jacobian[i, j] = -buses[i+i_offset].reactive_power_calculated - y_bus[i+1, i+1].imag * buses[i+i_offset].voltage * buses[i+i_offset].voltage
                else:
                    self.jacobian[i, j] = abs(buses[i+i_offset].voltage) * abs(buses[j+j_offset].voltage) * (
                                y_bus[i+1, j+1].real * np.sin(buses[i+i_offset].angle - buses[j+j_offset].angle) - y_bus[i+1, j+1].imag * np.cos(
                        buses[i+i_offset].angle - buses[j+j_offset].angle))
            # J2 of Jacobian
            for j in range(self.n_pq):
                if i == j:
                    self.jacobian[i, j+self.n] = buses[i+i_offset].real_power_calculated / abs(buses[i+i_offset].voltage) + y_bus[i+1, i+1].real * abs(buses[i+i_offset].voltage)
                else:
                    self.jacobian[i, j+self.n] = abs(buses[i+i_offset].voltage) * (y_bus[i+1, j+1].real * np.cos(buses[i+i_offset].angle - buses[j+j_offset].angle) + y_bus[i+1, j+1].imag * np.sin(
                        buses[i+i_offset].angle - buses[j+j_offset].angle))

        for i in range(self.n_pq):
            # J3 of Jacobian
            for j in range(self.n):
                if i == j:
                    self.jacobian[i+self.n, j] = buses[i+i_offset].real_power_calculated - y_bus[i+1, i+1].real * buses[i+i_offset].voltage * buses[i+i_offset].voltage
                else:
                    self.jacobian[i+self.n, j] = -abs(buses[i+i_offset].voltage) * abs(buses[j+j_offset].voltage) * (y_bus[i+1, j+1].real * np.cos(buses[i+i_offset].angle - buses[j+j_offset].angle) + y_bus[i+1, j+1].imag * np.sin(buses[i+i_offset].angle - buses[j+j_offset].angle))

            for j in range(self.n_pq):
                # J4 of Jacobian
                if i == j:
                    self.jacobian[i+self.n, j + self.n] = buses[i+i_offset].reactive_power_calculated / abs(buses[i+i_offset].voltage) - y_bus[i+1, i+1].imag * abs(buses[i+i_offset].voltage)
                else:
                    self.jacobian[i+self.n, j + self.n] = abs(buses[i+i_offset].voltage) * (y_bus[i+1, j+1].real * np.sin(buses[i+i_offset].angle - buses[j+j_offset].angle) - y_bus[i+1, j+1].imag * np.cos(buses[i+i_offset].angle - buses[j+j_offset].angle))


class NR_Method:
    """
    The newton raphson method class contains all necessary values and operations associated with the method
    """
    def __init__(self):
        pass


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
        self.real_power_delta = None
        self.reactive_power_delta = None
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
                self.bus_number, self.voltage, self.angle, self.real_power_calculated, self.reactive_power_calculated,
                self.real_power_delta, self.reactive_power_delta)
        return s


