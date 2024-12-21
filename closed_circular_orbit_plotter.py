import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from math import isclose

# Define U_eff'(x) = 0 to solve for circular orbits
def u_eff_derivative(x, beta):
    return (1 / x**2) * np.exp(-x) + (1 / x) * np.exp(-x) - (2 * beta) / x**3

# Define the second derivative of U_eff
def u_eff_second_derivative(x, beta, E0=1):
    term1 = (-2 / x**3) * np.exp(-x)
    term2 = (-1 / x**2) * np.exp(-x)
    term3 = (2 / x**2) * np.exp(-x)
    term4 = (1 / x) * np.exp(-x)
    term5 = (-6 * beta) / x**4
    return E0 * (term1 + term2 + term3 + term4 + term5)

# Calculate omega_r^2 (radial frequency squared)
def omega_r_squared(x, beta, E0=1, m=1):
    return (1 / m) * u_eff_second_derivative(x, beta, E0)

# Calculate omega_phi^2 (angular frequency squared)
def omega_phi_squared(x, beta, E0=1, m=1):
    return (2 * E0 * beta) / (m * x**4)

# Check if a number is approximately rational
def is_rational(value, tol=1e-3):
    for denominator in range(1, 100):  # Check fractions with small denominators
        if isclose(value, round(value * denominator) / denominator, abs_tol=tol):
            return True
    return False

# Prompt user for beta and x values
while True:
    try:
        beta_user = float(input("Enter the value of beta: "))
        x_user = float(input("Enter the value of x: "))
        break
    except ValueError:
        print("Invalid input. Please enter numerical values for beta and x.")

# Perform calculations
omega_r2_user = omega_r_squared(x_user, beta_user)
omega_phi2_user = omega_phi_squared(x_user, beta_user)

# Print the results
print(f"\nFor beta = {beta_user} and x = {x_user}:")

if omega_r2_user > 0 and omega_phi2_user > 0:
    omega_r_user = np.sqrt(omega_r2_user)
    omega_phi_user = np.sqrt(omega_phi2_user)
    ratio_user = omega_phi_user / omega_r_user

    print(f"Radial frequency (omega_r): {omega_r_user:.5f}")
    print(f"Angular frequency (omega_phi): {omega_phi_user:.5f}")
    print(f"Frequency ratio (omega_phi / omega_r): {ratio_user:.5f}")

    # Determine if orbit is circular
    if is_rational(ratio_user):
        print("The orbit is circular and closed (rational ratio).")
    else:
        print("The orbit is circular but not closed (irrational ratio).")
else:
    print("The orbit is not circular (invalid frequencies).")

# Always print the frequencies regardless of conditions
print("\nAdditional Information:")
print(f"omega_r^2 = {omega_r2_user:.5f}")
print(f"omega_phi^2 = {omega_phi2_user:.5f}")

# Plotting frequencies as a function of beta
beta_values = np.linspace(0.01, 0.42, 100)
omega_r_values = []
omega_phi_values = []
ratios = []

for beta in beta_values:
    try:
        # Solve for x
        result = root_scalar(u_eff_derivative, args=(beta,), bracket=[0.5, 5], method='brentq')
        x = result.root
        # Calculate frequencies
        omega_r2 = omega_r_squared(x, beta)
        omega_phi2 = omega_phi_squared(x, beta)
        if omega_r2 > 0 and omega_phi2 > 0:
            omega_r = np.sqrt(omega_r2)
            omega_phi = np.sqrt(omega_phi2)
            omega_r_values.append(omega_r)
            omega_phi_values.append(omega_phi)
            ratios.append(omega_phi / omega_r)
        else:
            omega_r_values.append(np.nan)
            omega_phi_values.append(np.nan)
            ratios.append(np.nan)
    except ValueError:
        omega_r_values.append(np.nan)
        omega_phi_values.append(np.nan)
        ratios.append(np.nan)


# Plot the ratio
plt.figure(figsize=(10, 6))
plt.plot(beta_values, ratios, label=r'$\omega_\phi / \omega_r$', color='green')
plt.axhline(y=1, color='red', linestyle='--', label='Rational = 1')
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\omega_\phi / \omega_r$')
plt.title(r'Frequency Ratio $\omega_\phi / \omega_r$ vs $\beta$')
plt.legend()
plt.grid()
plt.show()
