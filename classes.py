from newton_raphson_method.nr_classes import NR_Method

class Power_Network:
    """
    The power network class contains all the information about the network, and Bus objects for each bus in the network.
    """
    def __init__(self, p_dict, q_dict, voltage_dict, delta_dict, slack_bus_number):
        pass


class Jacobian(NR_Method):
    """
    The jacobian class contains the jacobian matrix builder and all values associated with the jacobian matrix
    It is a subclass of the NR_Method class
    """

    def __init__(self):
        pass

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