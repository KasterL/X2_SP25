# region imports
from scipy.integrate import solve_ivp, quad
import numpy as np
import matplotlib.pyplot as plt
from math import sin
# endregion

# region class definitions
class circuit():
    def __init__(self, R=10, L=20, C=0.05, A=20, w=20, p=0):
        """
        Initializes the RLC circuit parameters with default values.
        :param R: Resistance in ohms
        :param L:Inductance in henreys
        :param C:Capacitance in Farads
        :param A:Amplitude of the input voltage source
        :param w:Angular frequency for the input voltage
        :param p:Phase of the input voltage
        :param t:Time array for the simulation
        :param x:Solution array containing system state variables
        """
        #region attributes

        self.R = R
        self.L = L
        self.C = C
        self.A = A
        self.w = w
        self.p = p
        self.t = None  # Stores time values after simulation
        self.X = None  # Stores system state solutions after simulation
        # endregion

    # region methods
    def ode_system(self, t, X):
        """
        this is the odeSystem callback I'm using for solve_ivp()
        :param X: the current values of the state variables
        :param t: the current time
        :return: list of derivatives of state variables
        """
        i1, i2, vc = X  # State variables: current in loop 1, current in loop 2, capacitor voltage
        v_t = self.A * sin(self.w * t + self.p)  # Input voltage v(t)

        # Kirchhoff’s Law for the two loops
        di1_dt = (v_t - self.R * (i1 - i2)) / self.L  # Left loop
        di2_dt = (v_t - self.R * (i1 - i2)) / self.L - i2 / (self.R * self.C)  # Right loop
        dvc_dt = -i2 / self.C  # Capacitor voltage equation

        return [di1_dt, di2_dt, dvc_dt]

    def simulate(self, t=10, pts=500):
        """
        For simulating transient behavior of circuit.
        :param: time over which to carry out the simulation in seconds
        :param pts: number of points in the simulation
        :return: nothing, just store I
        """
        time = np.linspace(0, t, pts)  # Time array
        X0 = [0, 0, 0]  # Initial conditions: i1(0) = 0, i2(0) = 0, vc(0) = 0

        # Solve the ODE system using `solve_ivp` with RK45 method
        solution = solve_ivp(self.ode_system, (0, t), X0, t_eval=time, method='RK45')

        # Store results
        self.t = solution.t
        self.X = solution.y

    def doPlot(self):
        """
        Re-written on 4/21/2022 to adapt to plotting on GUI if ax is not None
        :param args: contains ((R, list of time values, and results of solve_ivp))
        :param ax: Matplotlib axis object
        :return: none
        """
        fig, ax = plt.subplots()  # New figure and primary axis

        # Extract solution variables
        i1 = self.X[0]  # Current in loop 1
        i2 = self.X[1]  # Current in loop 2
        vc = self.X[2]  # Capacitor voltage

        # Plot currents on primary y-axis
        ax.plot(self.t, i1, 'k-', label=r'$i_1(t)$')  # Solid line for i1
        ax.plot(self.t, i2, 'k--', label=r'$i_2(t)$')  # Dashed line for i2

        # Configure x and y-axis labels and limits
        ax.set_xlim(0, 10)  # Expected x-axis range
        ax.set_ylim(-0.06, 0.10)  # Expected y-axis range
        ax.set_xticks(np.arange(0,11,1))
        ax.set_xlabel("t (s)")
        ax.set_ylabel(r"$i_1, i_2$ (A)")
        ax.grid(True, linestyle='--', linewidth=0.5)  # Dashed grid lines

        # Create secondary y-axis for voltage across the capacitor
        ax2 = ax.twinx()
        ax2.plot(self.t, vc, 'k:', label=r'$v_C(t)$')  # Dotted line for vc
        ax2.set_ylim(-0.5, 0.1)  # Expected y-axis range for voltage
        ax2.set_ylabel(r"$v_C(t)$ (V)")

        # Add legends
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left')

        plt.show()  # Display the plot
    # endregion

# endregion
# region function definitions

def main():
    """
    For solving problem 2 on exam
    :return: none
    """
    goAgain = True
    default_params = {'R': 10, 'L': 20, 'C': 0.05, 'A': 20, 'w': 20, 'p': 0}  # create a circuit object with default values

    while goAgain:
        # Prompt user for circuit parameters with default values
        print("Press Enter to reuse previous values or enter new ones.")
        R = input(f"Enter resistance R (Ohms) [{default_params['R']}]: ")
        L = input(f"Enter inductance L (Henreys) [{default_params['L']}]: ")
        C = input(f"Enter capacitance C (Farads) [{default_params['C']}]: ")
        A = input(f"Enter voltage amplitude A [{default_params['A']}]: ")
        w = input(f"Enter voltage frequency w [{default_params['w']}]: ")
        p = input(f"Enter voltage phase p [{default_params['p']}]: ")

        # Convert inputs to float values or keep defaults
        R = float(R) if R else default_params['R']
        L = float(L) if L else default_params['L']
        C = float(C) if C else default_params['C']
        A = float(A) if A else default_params['A']
        w = float(w) if w else default_params['w']
        p = float(p) if p else default_params['p']

        # Update default values for next iteration
        default_params = {'R': R, 'L': L, 'C': C, 'A': A, 'w': w, 'p': p}

        # Run simulation with user-provided parameters
        CircuitObj = circuit(R, L, C, A, w, p)
        CircuitObj.simulate(t=10, pts=500)
        CircuitObj.doPlot()

        # Prompt user to run another simulation with new parameters
        choice = input("Run the simulation again with new parameters? (y/n): ")
        if choice.lower() != 'y':
            goAgain = False  # Exit loop if user selects 'n'

# endregion

# region function calls
if __name__ == "__main__":
    main()
# endregion
