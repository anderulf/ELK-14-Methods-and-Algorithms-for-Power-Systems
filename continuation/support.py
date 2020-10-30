import numpy as np
import copy
from newton_raphson_method.support import Load_Flow, run_newton_raphson

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

        max voltage step is not implemented in the assignment. It is used in continuation as a upper limit on how much
        a voltage can change through a prediction step
        """
        self.max_voltage_step = max_voltage_step
        self.max_load_step = max_load_step
        self.old_buses_dict = copy.deepcopy(self.buses_dict) # Stores the last step to be able to reverse
        self.old_mismatch = copy.deepcopy(self.mismatch)
        self.continuation_parameter = None
        self.phase = None
        self.step = max_load_step
        self.S = 1

    def initialize_predictor_phase(self):
        """
        Set up the predictor phase

        expand jacobian matrix and mismatch vector
        """
        self.phase = "predictor"
        self.jacobian.continuation_expand("load", self.buses_dict)
        self.mismatch.continuation_expansion(self.phase)
        self.expand_other_vectors()

    def initialize_corrector_phase(self, parameter, constant_voltage_bus=None):
        """
        Set up the corrector phase

        input parameter should be "load" or "voltage"
        """
        self.phase = "corrector"
        self.jacobian.continuation_expand(parameter, self.buses_dict, constant_voltage_bus)
        self.mismatch.continuation_expansion(self.phase)

    def determine_continuation_parameter(self):
        """
        Determine if NR load flow converges at the current step

        Computational heavy way of determining this. Could perhaps allow higher error limit? Other ways of determining?

        This method is used for the continuation method but not implemented in the ELK-14 assignment
        """
        self.jacobian.reset_original_matrix()
        self.mismatch.reset_original_vector()
        convergence = run_newton_raphson(N_R=self, printing=False, convergence=True)
        if convergence:
            self.continuation_parameter = "load"
            # It converged so save the new solution to avoid doing this step again
            self.store_values()
        else:
            self.continuation_flag = "voltage"
            # It did not converge, so don't save the new solution
            self.reverse_step()
        return self.continuation_parameter

    def update_continuation_values(self, parameter = None):
        """
        Increase x and P,Q
        """
        for bus in self.buses_dict.values():
            if bus.bus_number == self.slack_bus_number:
                pass
            else:
                if parameter == "load":
                    pass
                else: #parameter = "voltage eller none og dermed predictor phase"
                    bus.p_spec -= self.step * bus.beta * self.S
                    bus.q_spec -= self.step * bus.alpha * self.S
                if self.phase == "predictor":
                    bus.delta +=  self.x_diff[bus.bus_number-1][0] * self.step
                    bus.voltage += self.x_diff[bus.bus_number - 1 + self.n][0] * self.step
                    self.x_new[bus.bus_number - 1] = bus.delta
                    self.x_new[bus.bus_number -1 + self.n] = bus.voltage
                else: #phase = corrector
                    bus.delta += self.x_diff[bus.bus_number - 1][0]
                    bus.voltage += self.x_diff[bus.bus_number - 1 + self.n][0]
                    self.x_new[bus.bus_number - 1] = bus.delta
                    self.x_new[bus.bus_number - 1 + self.n] = bus.voltage

    def constant_voltage_bus(self):
        """
        Return the index in the jacobian matrix for the highest voltage change
        The returned value is offset with the number of pq and pv buses corresponding to the number of angle elements before
        the voltage elements in the jacobian matrix.

        This method is used for the continuation method but not implemented in the ELK-14 assignment
        """
        max_voltage_bus_index = 0
        max_voltage = 0
        index = 0
        for bus_number in self.buses_dict:
            if bus_number == self.slack_bus_number:
                pass
            else:
                max_voltage_temp = self.old_buses_dict[bus_number].voltage - self.buses_dict[bus_number].voltage
                if max_voltage_temp > max_voltage:
                    max_voltage = max_voltage_temp
                    max_voltage_bus_index = index;
                else:
                    pass
                index += 1
        return self.n + max_voltage_bus_index

    def store_values(self):
        """
        If a step was successful store the new values in the old variable
        """
        self.old_buses_dict = copy.deepcopy(self.buses_dict)
        self.old_mismatch = copy.deepcopy(self.mismatch)

    def reverse_step(self):
        """
        If a step was unsuccessful revert the values to the last step
        """
        self.buses_dict = self.old_buses_dict
        self.mismatch = self.old_mismatch

    def expand_other_vectors(self):
        """
        Expands the x vectors (x_old and x_new) with one extra row

        Expands the load flow label vectors for the continuation power flow method.
        """
        if len(self.mismatch_vector_labels) > self.m:
            return
        else:
            self.mismatch_vector_labels.append("S")
            self.x_vector_labels.append("S")
            self.correction_vector_labels.append("\u0394S")
            self.x_old = np.vstack([self.x_old, 0])
            self.x_new = np.vstack([self.x_new, 0])
            self.x_diff = np.vstack([self.x_diff, 0])
