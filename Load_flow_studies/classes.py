# One general class for the NR method

# One class for the busses

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
            s = " /*/*/*/ SLACK BUS (bus {}) /*/*/*/".format(self.bus_number)
        else:
            s = "Bus {}: voltage = {} pu, angle =  {}, P_calc = {}, Q_calc = {}, delta P  = {}, delta Q = {}".format(
                self.bus_number, self.voltage, self.angle, self.real_power_calculated, self.reactive_power_calculated,
                self.real_power_delta, self.reactive_power_delta)
        return s