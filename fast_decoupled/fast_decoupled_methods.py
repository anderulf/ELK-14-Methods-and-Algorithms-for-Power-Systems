import numpy as np
import cmath as ma
from supporting_methods import print_title3, create_simplified_y_bus
from newton_raphson_method.newton_raphson_support import Load_Flow

# Numpy printing options
np.set_printoptions(suppress=True)  # suppress scientific notations

def run_primal_method(fast_dec, printing=False):
    """
    Inputs a Fast_Decoupled object fast_dec
    returns the number of iterations
    """
    fast_dec.iteration = 0
    # Calculate initial active power injections
    fast_dec.calculate_P_injections()
    # Calculate the initial mismatches
    fast_dec.calculate_fast_decoupled_mismatches("P")
    while fast_dec.power_error() > 0.0001:
        fast_dec.iteration += 1
        # Calculate the corrections for angles
        theta_correction = np.linalg.solve(fast_dec.B_p, fast_dec.mismatch.get_P())
        # Update theta i X-vector, theta_new = theta_correction + theta_old
        fast_dec.update_fast_decoupled_voltage_or_angle(angles=theta_correction)
        # Calculate reactive power injections based on new angles
        fast_dec.calculate_Q_injections()
        fast_dec.calculate_fast_decoupled_mismatches("Q")
        # Voltage_correction = B_dp.invers*Q_mismatch(v, theta_updated) (samme som over)
        voltage_correction = np.linalg.solve(fast_dec.B_dp, fast_dec.mismatch.get_Q())
        # Update voltage i X-vector, V_new = V_correction + v_old
        fast_dec.update_fast_decoupled_voltage_or_angle(voltages=voltage_correction)

        # Calculate updated active power injections and reactive power injections
        fast_dec.calculate_P_injections()
        fast_dec.calculate_Q_injections()
        # Calculate the mismatches
        fast_dec.calculate_fast_decoupled_mismatches("P")
        fast_dec.calculate_fast_decoupled_mismatches("Q")

        if printing:
            print_title3("Iteration {}".format(fast_dec.iteration))
            fast_dec.print_data(theta_correction, voltage_correction)
        if fast_dec.diverging():
            print_title3("No convergence")
            break
    return fast_dec.iteration

def run_dual_method(fast_dec, printing=False):
    """
    Inputs a Fast_Decoupled object fast_dec
    Returns the number of iterations
    """
    fast_dec.iteration = 0
    # Calculate initial reactive power injections and mismatches
    fast_dec.calculate_Q_injections()
    fast_dec.calculate_fast_decoupled_mismatches("Q")
    while fast_dec.power_error() > 0.0001:
        fast_dec.iteration += 1
        # Voltage_correction = B_dp.invers*Q_mismatch(v, theta_updated) (samme som over)
        voltage_correction = np.linalg.solve(fast_dec.B_dp, fast_dec.mismatch.get_Q())
        # Update voltage i X-vector, V_new = V_correction + v_old
        fast_dec.update_fast_decoupled_voltage_or_angle(voltages=voltage_correction)

        # Calculate active power injections based on new voltage corrections
        fast_dec.calculate_P_injections()
        # Calculate the mismatches
        fast_dec.calculate_fast_decoupled_mismatches("P")
        # Calculate the corrections for angles
        theta_correction = np.linalg.solve(fast_dec.B_p, fast_dec.mismatch.get_P())
        # Update theta in X-vector, theta_new = theta_correction + theta_old
        fast_dec.update_fast_decoupled_voltage_or_angle(angles=theta_correction)

        # Calculate updated active power injections and reactive power injections
        fast_dec.calculate_P_injections()
        fast_dec.calculate_Q_injections()
        # Calculate the mismatches
        fast_dec.calculate_fast_decoupled_mismatches("P")
        fast_dec.calculate_fast_decoupled_mismatches("Q")

        if printing:
            print_title3("Iteration {}".format(fast_dec.iteration))
            fast_dec.print_data(theta_correction, voltage_correction)
        if fast_dec.diverging():
            print_title3("No convergence")
            break
    return fast_dec.iteration

def run_standard_method(fast_dec, printing=False):  # Standard is implemented in the same way as the primal algorithm
    fast_dec.iteration = 0
    # Calculate initial active power injections
    fast_dec.calculate_P_injections()
    # Calculate initial active power mismatches
    fast_dec.calculate_fast_decoupled_mismatches("P")
    while fast_dec.power_error() > 0.0001:
        fast_dec.iteration += 1
        # Calculate the corrections for angles
        theta_correction = np.linalg.solve(fast_dec.B_p, fast_dec.mismatch.get_P())
        # Update theta i X-vector, theta_new = theta_correction + theta_old
        fast_dec.update_fast_decoupled_voltage_or_angle(angles=theta_correction)
        # Calculate reactive power injections based on new angles
        fast_dec.calculate_Q_injections()
        fast_dec.calculate_fast_decoupled_mismatches("Q")
        # Voltage_correction = B_dp.invers*Q_mismatch(v, theta_updated) (samme som over)
        voltage_correction = np.linalg.solve(fast_dec.B_dp, fast_dec.mismatch.get_Q())
        # Update voltage i X-vector, V_new = V_correction + v_old
        fast_dec.update_fast_decoupled_voltage_or_angle(voltages=voltage_correction)
        # Calculate updated active power injections and reactive power injections
        fast_dec.calculate_P_injections()
        fast_dec.calculate_Q_injections()
        # Calculate the mismatches
        fast_dec.calculate_fast_decoupled_mismatches("P")
        fast_dec.calculate_fast_decoupled_mismatches("Q")

        if printing:
            print_title3("Iteration {}".format(fast_dec.iteration))
            fast_dec.print_data(theta_correction, voltage_correction)
        if fast_dec.diverging():
            print_title3("No convergence")
            break
    return fast_dec.iteration


class Fast_Decoupled(Load_Flow):
    """
    Subclass of Load Flow class

    From the jobian matrix we define the following submatrices:
    H = J1
    N = J2
    M = J3
    L = J4
    """
    def set_up_matrices(self, phase = None):
        """
        Set up matrices needed for fast decoupled power flow with inputed phase
        """
        self.jacobian.create(fast_decoupled=True)
        # Store submatrices of jacobian
        self.H = self.jacobian.matrix[0:self.n, 0:self.n]
        self.L = self.jacobian.matrix[self.n:, self.n:]
        # Create zero matrices for correction matrices building
        self.N_zeros = np.zeros([self.n_pq, self.n_pq])
        self.M_zeros = np.zeros([self.n_pq, self.n])
        if phase == "Primal":
            self.B_p = self.H
        elif phase == "Dual":
            self.B_dp = self.L
        self.create_modified_jacobians(phase)

    def create_modified_jacobians(self, phase):
        if phase == "Primal":
            ## B_dp
            self.B_dp = np.zeros([len(self.buses_dict), len(self.buses_dict)], dtype=float)
            self.B_dp = create_simplified_y_bus(self.B_dp, self.lines, self.slack_bus_number)
        elif phase == "Dual":
            ## B_p
            self.B_p = np.zeros([len(self.buses_dict), len(self.buses_dict)], dtype=float)
            self.B_p = create_simplified_y_bus(self.B_p, self.lines, self.slack_bus_number)
        elif phase == "Standard":
            ## B_p
            self.B_p = np.zeros([len(self.buses_dict), len(self.buses_dict)], dtype=float)
            self.B_p = create_simplified_y_bus(self.B_p, self.lines, self.slack_bus_number)
            self.B_dp = self.B_p

    def calculate_P_injections(self):
        """
        Calculates the new active power injections based on previous voltage and angles
        """
        buses = self.buses_dict
        for number, i in enumerate(buses):
            buses[i].p_calc = 0  # Resets the value of p_calc so the loop works
            # Skip slack bus
            if i == self.slack_bus_number:
                pass
            else:
                for j in buses:
                    # Adding the Q's from the lines. Note that Ybus is offset with -1 because Python uses 0-indexing and the buses are indexed from 1
                    buses[i].p_calc += abs(self.y_bus[i - 1, j - 1]) * buses[i].voltage * buses[j].voltage * np.cos(
                        buses[i].delta - buses[j].delta - ma.phase(self.y_bus[i - 1, j - 1]))

    def calculate_Q_injections(self):
        """
        Calculates the new reactive power injections based on previous voltage and angles
        """
        buses = self.buses_dict
        for number, i in enumerate(buses):
            buses[i].q_calc = 0 # Resets the value of q_calc so the loop works
            # Skip slack bus
            if i == self.slack_bus_number:
                pass
            else:
                for j in buses:
                    # Adding the Q's from the lines. Note that Ybus is offset with -1 because Python uses 0-indexing and the buses are indexed from 1
                    buses[i].q_calc += abs(self.y_bus[i - 1, j - 1]) * buses[i].voltage * buses[j].voltage * np.sin(
                        buses[i].delta - buses[j].delta - ma.phase(self.y_bus[i - 1, j - 1]))

    def calculate_fast_decoupled_mismatches(self, parameter):
        """
        Calculate the mismatches and store in mismatch vector for either P or Q
        """
        if parameter == "P":
            for i in range(len(self.buses_dict)):
                if self.buses_dict[i + 1].bus_type == "PQ" or self.buses_dict[i + 1].bus_type == "PV":
                    self.buses_dict[i + 1].delta_p = self.buses_dict[i+1].p_spec - self.buses_dict[i + 1].p_calc
                    self.mismatch.vector[i, 0] = self.buses_dict[i + 1].delta_p
        if parameter == "Q":
            for i in range(len(self.buses_dict)):
                if self.buses_dict[i + 1].bus_type == "PQ":
                    self.buses_dict[i + 1].delta_q = self.buses_dict[i+1].q_spec - self.buses_dict[i + 1].q_calc
                    self.mismatch.vector[i + self.n, 0] = self.buses_dict[i + 1].delta_q

    def update_fast_decoupled_voltage_or_angle(self, voltages=None, angles=None):
        """
        Updates either voltages or angles based on which are inputed

        Skips the slack bus
        """
        number = 0
        for bus in self.buses_dict.values():
            if bus.bus_number != self.slack_bus_number:
                if voltages is not None:
                    bus.voltage += voltages[number][0]
                elif angles is not None:
                    bus.delta += angles[number][0]
                number += 1
            else: # slack
                pass

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
            elif self.iteration >= 30:
                # Stop if too many iterations were run
                return True
            else: return False

    def print_data(self, theta_correction, voltage_correction):
        """
        Overwritten from Load Flow class to match the interesting parts of the fast decoupled power flow
        """
        print("\nJacobian: B':\n", np.round(self.B_p, 4))
        print("\nJacobian: B\":\n", np.round(self.B_dp, 4))
        print("\nP injections:")
        for bus in self.buses_dict.values():
            if bus.bus_number == self.slack_bus_number:
                pass
            else:
                print("P{} = {}".format(bus.bus_number, round(bus.p_calc, 3)))
        print("\nQ injections:")
        for bus in self.buses_dict.values():
            if bus.bus_number == self.slack_bus_number:
                pass
            else:
                print("Q{} = {}".format(bus.bus_number, round(bus.q_calc, 3)))
        print("\nP Mismatches:")
        print(np.c_[self.mismatch.get_P_label(), np.round(self.mismatch.get_P(), 4)])
        print("\nQ Mismatches:")
        print(np.c_[self.mismatch.get_Q_label(), np.round(self.mismatch.get_Q(), 4)])
        print("\nAngle corrections:")
        print(np.c_[self.correction_vector_labels[:self.n], np.round(theta_correction, 4)])
        print("\nVoltage corrections:")
        print(np.c_[self.correction_vector_labels[self.n:], np.round(voltage_correction, 4)])

    def print_final_solution(self, phase):
        """
        Prints data for the finished solution after iterations has completed
        """
        print("\nFinal solution for {} method: ".format(phase))
        for bus in self.buses_dict.values():
            bus.print_data(slack_bus_number=self.slack_bus_number)