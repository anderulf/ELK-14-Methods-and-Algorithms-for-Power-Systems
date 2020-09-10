import numpy as np
import cmath as ma


class NR_Method:
    def __init__(self, p_dict, q_dict, voltage_dict, delta_dict, slack_node, y_bus):
        """
        Initializing the class. p_dict and q_dict should be lists of dictionary types holding key equal to node number 1, 2, .. and values equal
        their rated values in pu. If the rated active or reactive power is not given set value to None.

        The object calculates the number of PQ, PV, slack buses based on the input data

        For voltage and delta set value to flat start value or initial value

        The input data should have the same length (all nodes except slack)

        self.nodes_dict: dictionary holding key: node number value: node object
            Dictionary type is used because it gives better control over indexing. Now indexing can be done per node number (key), and not from
            the position in the list

        The limit flag is used to know if the reactive power limit has been reached
        """
        self.nodes_dict = {}
        self.fill_nodes_dict(p_dict, q_dict, voltage_dict, delta_dict)
        self.slack_node = slack_node
        self.y_bus = y_bus
        self.n_pq = 0
        self.n_pv = 0
        self.n_pd = 0
        self.n = 0
        self.m = 0
        self.get_n_values()
        self.jacobian = np.zeros([self.m, self.m])
        self.loss_matrix_p = None
        self.loss_matrix_q = None
        self.power_flow_matrix = None
        self.current_matrix = None
        self.total_losses_p = 0
        self.total_losses_q = 0
        self.limit_flag = 0

    def fill_nodes_dict(self, p_dict, q_dict, voltage_dict, delta_dict):
        """
        Initialize nodes_dict based on input data
        """
        for node_num in p_dict:
            self.nodes_dict[int(node_num)] = Node(p_dict[node_num], q_dict[node_num], voltage_dict[node_num],
                                                  delta_dict[node_num])

    def get_n_values(self):
        """
        Calculate number of different bus types
        """
        for node in self.nodes_dict:
            if self.nodes_dict[node].p_spec and self.nodes_dict[node].q_spec:
                self.n_pq += 1
            elif self.nodes_dict[node].p_spec and not self.nodes_dict[node].q_spec:
                self.n_pv += 1
            else:
                self.n_pd += 1
        if self.n_pd > 1 or self.n_pd == 0:
            print("WARNING INVALID SYSTEM: system has {} slack buses".format(self.n_pd))
        else:
            self.m = 2 * self.n_pq + self.n_pv
            self.n = self.n_pq + self.n_pv

    def print_nodes(self):
        """
        Prints the data for all the node. Mostly needed for debugging
        """
        print(
            "System has {0} slack bus, {1} PQ-bus(es) and {2} PV-bus(es). The jacobian matrix has dimension: {3}x{3}".format(
                self.n_pd, self.n_pq, self.n_pv, self.m))
        for node_num in self.nodes_dict:
            self.nodes_dict[node_num].print_data(node_num, self.slack_node)

    def calc_new_power(self):
        """
        Calculate power values based on voltage and delta for all nodes except slack.
        The values are saved in the object
        """
        nodes = self.nodes_dict
        for i in nodes:
            nodes[i].p_calc = 0  # Resets the value of p_calc/q_calc so the loop works
            nodes[i].q_calc = 0
            # Skip slack bus
            if i == self.slack_node:
                pass
            else:
                for j in nodes:
                    try:
                        # Adding the Q's from the lines. Note that Ybus is offset with -1 because Python uses 0-indexing and the nodes are indexed from 1
                        nodes[i].p_calc += abs(self.y_bus[i - 1, j - 1]) * nodes[i].voltage * nodes[j].voltage * np.cos(
                            nodes[i].delta - nodes[j].delta - ma.phase(self.y_bus[i - 1, j - 1]))
                        nodes[i].q_calc += abs(self.y_bus[i - 1, j - 1]) * nodes[i].voltage * nodes[j].voltage * np.sin(
                            nodes[i].delta - nodes[j].delta - ma.phase(self.y_bus[i - 1, j - 1]))
                    except Exception as e:
                        print(e)
                if i == 1:
                    print("Q_calc,1")
                    print(abs(self.y_bus[0, 0]))
                    print(nodes[1].delta - nodes[1].delta - ma.phase(self.y_bus[0, 0]))
                    print("+")
                    print(abs(self.y_bus[0, 1]))
                    print(nodes[1].delta - nodes[2].delta - ma.phase(self.y_bus[0, 1]))
                    print("+")
                    print(abs(self.y_bus[0, 2]))
                    print(nodes[1].delta - nodes[3].delta - ma.phase(self.y_bus[0, 2]))
                    print("=")
                    print(nodes[1].q_calc)


    def check_limit(self, q_limit, lim_node, lim_size):
        """
        Check if a certain bus has reached the q_limit, and change the bus if it has
        """
        if q_limit:
            if self.nodes_dict[lim_node].q_calc > lim_size and not self.limit_flag:
                self.limit_flag = 1
                self.nodes_dict[lim_node].q_spec = lim_size
            if self.nodes_dict[lim_node].q_calc < lim_size and self.limit_flag:
                self.limit_flag = 0
                self.nodes_dict[lim_node].q_spec = None
            self.n_pq = 0
            self.n_pv = 0
            self.n_pd = 0
            self.get_n_values()

    def error_specified_vs_calculated(self):
        """
        Finds all the error terms and store in object
        """
        for node in self.nodes_dict:
            if node == self.slack_node:
                pass
            else:
                # Ignore all nodes which doesn't have a specified value (value == None)
                if self.nodes_dict[node].p_spec:
                    self.nodes_dict[node].delta_p = self.nodes_dict[node].p_spec - self.nodes_dict[node].p_calc
                if self.nodes_dict[node].q_spec:
                    self.nodes_dict[node].delta_q = self.nodes_dict[node].q_spec - self.nodes_dict[node].q_calc

    def power_error(self):
        """
        Adds all the error terms to a list and returns the maximum value of the list
        """
        error_list = []
        nodes = self.nodes_dict
        for i in nodes:
            # Ignore all nodes which doesn't have a specified value (value == None)
            if not nodes[i].p_spec == None:
                error_list.append(abs(nodes[i].delta_p))
            if not nodes[i].q_spec == None:
                error_list.append(abs(nodes[i].delta_q))
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
        nodes = self.nodes_dict
        self.jacobian = np.zeros([self.m, self.m])
        #i_offset = 1
        #j_offset = 0
        for i in range(self.n):
            # J1 of Jacobian
            #if i == self.slack_node-1:
            #    i_offset = 2
            for j in range(self.n):
                print("Jacobian J1 i,j = {},{}".format(i, j))
                #if j == self.slack_node-1:
                #    j_offset = 2
                if i == j:
                    self.jacobian[i, j] = -nodes[i + 1].q_calc - self.y_bus[i, j].imag * nodes[
                        i + 1].voltage * nodes[i + 1].voltage
                    if i == 0 and j == 0:
                        print("Jacobian 1,1 is:")
                        print(-nodes[1].q_calc)
                        print("-")
                        print(self.y_bus[0,0].imag)
                        print("*")
                        print(nodes[1].voltage)
                        print("*")
                        print(nodes[1].voltage)
                        print("=")
                        print(-nodes[1].q_calc - self.y_bus[0,0].imag * nodes[1].voltage * nodes[1].voltage)
                else:
                    self.jacobian[i, j] = abs(nodes[i + 1].voltage) * abs(nodes[j + 1].voltage) * (
                            self.y_bus[i, j].real * np.sin(nodes[i + 1].delta - nodes[j + 1].delta) -
                            self.y_bus[i, j].imag * np.cos(
                        nodes[i + 1].delta - nodes[j + 1].delta))
            #j_offset = 1
            # J2 of Jacobian
            for j in range(self.n_pq):
                #if j == self.slack_node-1:
                #    j_offset = 2
                if i == j:
                    self.jacobian[i, j + self.n] = nodes[i + 1].p_calc / abs(nodes[i + 1].voltage) + \
                                                   self.y_bus[i, i].real * abs(nodes[i + 1].voltage)
                else:
                    self.jacobian[i, j + self.n] = abs(nodes[i + 1].voltage) * (self.y_bus[i, j].real * np.cos(
                            nodes[i + 1].delta - nodes[j + 1].delta) + self.y_bus[i, j].imag * np.sin(
                            nodes[i + 1].delta - nodes[j + 1].delta))
            #j_offset = 1
        #i_offset = 1
        #j_offset = 1
        for i in range(self.n_pq):
            # J3 of Jacobian
            #if i == self.slack_node-1:
            #    i_offset = 2
            for j in range(self.n):
                #if j == self.slack_node-1:
                #    j_offset = 2
                if i == j:
                    self.jacobian[i + self.n, j] = nodes[i + 1].p_calc - self.y_bus[i, i].real * nodes[
                        i + 1].voltage * nodes[i + 1].voltage
                else:
                    self.jacobian[i + self.n, j] = -abs(nodes[i + 1].voltage) * abs(
                        nodes[j + 1].voltage) * (self.y_bus[i, j].real * np.cos(
                        nodes[i + 1].delta - nodes[j + 1].delta) + self.y_bus[i, j].imag * np.sin(
                        nodes[i + 1].delta - nodes[j + 1].delta))
            #j_offset = 1
            for j in range(self.n_pq):
                # J4 of Jacobian
                #if j == self.slack_node-1:
                #    j_offset = 2
                if i == j:
                    self.jacobian[i + self.n, j + self.n] = nodes[i + 1].q_calc / abs(
                        nodes[i + 1].voltage) - self.y_bus[i, i].imag * abs(nodes[i + 1].voltage)
                else:
                    self.jacobian[i + self.n, j + self.n] = abs(nodes[i + 1].voltage) * (
                                self.y_bus[i, j].real * np.sin(
                            nodes[i + 1].delta - nodes[j + 1].delta) - self.y_bus[i, j].imag * np.cos(
                            nodes[i + 1].delta - nodes[j + 1].delta))
            #j_offset = 1

    def update_values(self):
        """
        Update the the system values by calculating the next step in the NS method
        """
        x_new = np.zeros([self.m, 1])
        x_old = np.zeros([self.m, 1])
        diff_b = np.zeros([self.m, 1])
        for i in range(self.n):
            x_old[i, 0] = self.nodes_dict[i + 2].delta
            diff_b[i, 0] = self.nodes_dict[i + 2].delta_p
        for i in range(self.n_pq):
            diff_b[i + self.n, 0] = self.nodes_dict[i + 2].delta_q
            x_old[i + self.n, 0] = self.nodes_dict[i + 2].voltage
        print("delta_x")
        print(np.linalg.solve(self.jacobian, diff_b))
        x_new = x_old + np.linalg.solve(self.jacobian, diff_b)
        for i in range(self.n):
            self.nodes_dict[i + 2].delta = x_new[i, 0]
        for i in range(self.n_pq):
            self.nodes_dict[i + 2].voltage = x_new[i + self.n, 0]

    def calculate_line_data(self):
        """
        Calculate the losses for all the top right values (skipping diagonals) and store in seperate matrices for active and reactive losses
        """
        nodes = self.nodes_dict
        rows, cols = self.y_bus.shape
        self.loss_matrix_p = np.zeros([rows, cols])
        self.loss_matrix_q = np.zeros([rows, cols])
        self.power_flow_matrix = np.zeros([rows, cols], dtype=np.complex_)
        self.current_matrix = np.zeros([rows, cols], dtype=np.complex_)
        for i in range(rows - 1):
            for j in range(rows - 1, i, -1):
                # Ignore non-existing lines
                if self.y_bus[i, j]:
                    # Get rectangular values
                    v_i = polar_to_rectangular(nodes[i + 1].voltage, nodes[i + 1].delta)
                    v_j = polar_to_rectangular(nodes[j + 1].voltage, nodes[j + 1].delta)
                    # Calculate losses
                    current_ij = self.y_bus[i, j] * (v_i - v_j)
                    current_ji = self.y_bus[i, j] * (v_j - v_i)
                    self.power_flow_matrix[
                        i, j] = -v_i * current_ij.conjugate()  # Negative values because power injections are defined into branch, but this is not intuitive for flow between branches
                    self.power_flow_matrix[
                        j, i] = -v_j * current_ji.conjugate()  # Negative values because power injections are defined into branch, but this is not intuitive for flow between branches
                    loss = v_i * current_ij.conjugate() + v_j * current_ji.conjugate()
                    self.current_matrix[i, j] = current_ij
                    self.current_matrix[j, i] = current_ji
                    self.loss_matrix_p[i, j] = abs(loss.real)
                    self.loss_matrix_q[i, j] = abs(loss.imag)
                    # Copy to other side of diagonal just in case
                    self.loss_matrix_p[j, i] = self.loss_matrix_p[i, j]
                    self.loss_matrix_q[j, i] = self.loss_matrix_q[i, j]
                    # Save total losses
                    self.total_losses_p += self.loss_matrix_p[i, j]
                    self.total_losses_q += self.loss_matrix_q[i, j]

    def print_line_data(self):
        """
        Print line values (losses, flows, current)
        """
        rows, cols = self.loss_matrix_p.shape
        print("Power losses in lines: ")
        for i in range(rows - 1):
            for j in range(rows - 1, i, -1):
                if self.loss_matrix_p[i, j]:
                    print("Line {}-{} has I={}, P_flow={}, Q_flow={}, P_loss={} and Q_loss={}".format(i + 1, j + 1,
                                                                                                      round(
                                                                                                          self.current_matrix[
                                                                                                              i, j], 4),
                                                                                                      round(
                                                                                                          self.power_flow_matrix[
                                                                                                              i, j].real,
                                                                                                          4), round(
                            self.power_flow_matrix[i, j].imag, 4), round(self.loss_matrix_p[i, j], 4), round(
                            self.loss_matrix_q[i, j], 4)))

    def calculate_slack_values(self):
        """
        Calculate the slack bus values based on the NS iteration
        """
        for i in self.nodes_dict:
            # Skip slack
            if i == self.slack_node:
                pass
            else:
                # Append injections for each bus. Using negative values because slack should cover loads but are covered by other generators
                self.nodes_dict[self.slack_node].p_calc += -self.nodes_dict[i].p_calc
                self.nodes_dict[self.slack_node].q_calc += -self.nodes_dict[i].q_calc
        # Add losses
        self.nodes_dict[self.slack_node].p_calc += self.total_losses_p
        self.nodes_dict[self.slack_node].q_calc += self.total_losses_q


class Node:
    """
    Object holding data for a node
    """

    def __init__(self, p_spec, q_spec, voltage, delta):
        self.p_spec = p_spec
        self.q_spec = q_spec
        self.voltage = voltage
        self.delta = delta
        self.p_calc = 0
        self.q_calc = 0
        self.delta_p = 1
        self.delta_q = 1

    def print_data(self, node_num, slack_node):
        """
        Print the data for a node
        """
        if node_num == slack_node:
            s = " **SLACK NODE**"
        else:
            s = ""
        print(
            "Node {}: P_spec = {}, Q_spec = {}, voltage = {}, delta = {} deg, P_calc = {}, Q_calc = {}, deltaP = {}, deltaQ = {}".format(
                node_num, self.p_spec, self.q_spec, round(self.voltage, 4), round(self.delta * 180 / np.pi, 4),
                round(self.p_calc, 4), round(self.q_calc, 4), round(self.delta_p, 6), round(self.delta_q, 6)) + s)


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