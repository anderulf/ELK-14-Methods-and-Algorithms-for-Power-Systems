import numpy as np

def run_primal_method(fast_dec, printing=False):
    """
    Inputs a Fast_Decoupled object fast_dec
    returns the number of iterations
    """
    iteration = 0
    while fast_dec.power_error() > 0.0001:
        iteration += 1
        # Calculate active power injections
        fast_dec.calculate_P_injections()
        # Calculate the mismatches
        fast_dec.calculate_fast_decoupled_mismatches("P")
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
        if printing:
            print("\nIteriation ", iteration)
            fast_dec.print_data(theta_correction, voltage_correction)
        if fast_dec.diverging():
            print("No convergence")
            break
    return iteration

def run_dual_method(fast_dec, printing=False):
    """
    Inputs a Fast_Decoupled object fast_dec
    Returns the number of iterations
    """
    iteration = 0
    while fast_dec.power_error() > 0.0001:
        iteration += 1
        # Calculate reactive power injections
        fast_dec.calculate_Q_injections()
        fast_dec.calculate_fast_decoupled_mismatches("Q")
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
        if printing:
            print("\nIteriation ", iteration)
            fast_dec.print_data(theta_correction, voltage_correction)
        if fast_dec.diverging():
            print("No convergence")
            break
    return iteration

def run_standard_method(fast_dec, printing=False):
    iteration = 0
    while fast_dec.power_error() > 0.0001:
        iteration += 1
        # Calculate active power injections
        fast_dec.calculate_P_injections()
        # Calculate the mismatches
        fast_dec.calculate_fast_decoupled_mismatches("P")
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
        if printing:
            print("\nIteration ", iteration)
            fast_dec.print_data(theta_correction, voltage_correction)
        if fast_dec.diverging():
            print("No convergence")
            break
    return iteration
