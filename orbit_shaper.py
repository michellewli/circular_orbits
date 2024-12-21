import numpy as np
import matplotlib.pyplot as plt
from fractions import Fraction


# Function to plot an orbit given the frequencies
def plot_orbit(omega_r, omega_phi, r0=1, time_duration=10, num_points=1000):
    """
    Plot an orbit based on the radial frequency (omega_r) and angular frequency (omega_phi).

    Parameters:
        omega_r (float): Radial frequency
        omega_phi (float): Angular frequency
        r0 (float): Initial radial displacement
        time_duration (float): Total time for the simulation
        num_points (int): Number of points for the plot
    """
    # Check if the ratio is rational
    ratio = omega_phi / omega_r
    if Fraction(ratio).denominator == 1:
        print(f"The orbit is closed (rational ratio: {ratio}).")
        # Calculate least common multiple (LCM) of the periods
        T_r = 2 * np.pi / omega_r
        T_phi = 2 * np.pi / omega_phi
        time_duration = np.lcm(int(T_r * 1000), int(T_phi * 1000)) / 1000  # LCM in seconds
        print(f"Adjusted time duration for closed orbit: {time_duration:.2f}")
    else:
        print(f"The orbit is not closed (irrational ratio: {ratio}).")

    # Time array
    t = np.linspace(0, time_duration, num_points)

    # Radial displacement (r) as a sinusoidal function of time
    r = r0 + 0.1 * np.sin(omega_r * t)  # r oscillates around r0

    # Angular displacement (phi) as a function of time
    phi = omega_phi * t

    # Parametric equations for x and y
    x = r * np.cos(phi)
    y = r * np.sin(phi)

    # Plot the orbit
    plt.figure(figsize=(6, 6))
    plt.plot(x, y, label=f"$\omega_r = {omega_r}$, $\omega_\phi = {omega_phi}$")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Orbit Plot")
    plt.axis("equal")
    plt.legend()
    plt.grid()
    plt.show()


# Example usage
omega_r = 1  # Radial frequency
omega_phi = 2.5  # Angular frequency
plot_orbit(omega_r, omega_phi)
