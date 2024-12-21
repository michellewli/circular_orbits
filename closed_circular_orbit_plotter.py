import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from math import isclose

# Define U_eff'(x) = 0 to solve for x in terms of beta
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

# Calculate omega_r and omega_phi
def omega_r_squared(x, beta, E0=1, m=1):
    return (1 / m) * u_eff_second_derivative(x, beta, E0)

def omega_phi_squared(x, beta, E0=1, m=1):
    return (2 * E0 * beta) / (m * x**4)

# Determine if a number is approximately rational
def is_rational(value, tol=1e-3):
    for denominator in range(1, 100):  # Check fractions with small denominators
        if isclose(value, round(value * denominator) / denominator, abs_tol=tol):
            return True
    return False

# Prompt user for beta and x values
while True:
    try:
        beta_user = float(input("Enter the value of beta for a single calculation: "))
        x_user = float(input("Enter the value of x for a single calculation: "))
        break
    except ValueError:
        print("Invalid input. Please enter numerical values for beta and x.")

# Single calculation for user-provided beta and x
omega_r2_user = omega_r_squared(x_user, beta_user)
omega_phi2_user = omega_phi_squared(x_user, beta_user)

if omega_r2_user > 0 and omega_phi2_user > 0:
    ratio_user = np.sqrt(omega_phi2_user / omega_r2_user)
    print(f"For beta = {beta_user} and x = {x_user}:")
    print(f"omega_r^2 = {omega_r2_user}")
    print(f"omega_phi^2 = {omega_phi2_user}")
    print(f"omega_phi / omega_r = {ratio_user}")
    if is_rational(ratio_user):
        print("The orbit is circular and closed (rational ratio).")
    else:
        print("The orbit is circular but not closed (irrational ratio).")
else:
    print(f"The values result in an open orbit for beta = {beta_user} and x = {x_user}.")

# Plotting for a range of beta values
beta_values = np.linspace(0.01, 0.42, 100)
ratios = []

for beta in beta_values:
    # Solve for x
    try:
        result = root_scalar(u_eff_derivative, args=(beta,), bracket=[0.5, 5], method='brentq')
        x = result.root
        # Calculate omega_r and omega_phi
        omega_r2 = omega_r_squared(x, beta)
        omega_phi2 = omega_phi_squared(x, beta)
        if omega_r2 > 0 and omega_phi2 > 0:
            ratios.append(np.sqrt(omega_phi2 / omega_r2))
        else:
            ratios.append(np.nan)
    except ValueError:
        ratios.append(np.nan)

# Plot the results
plt.plot(beta_values, ratios, label=r'$\omega_\phi / \omega_r$')
plt.axhline(y=1, color='r', linestyle='--', label='Rational = 1')
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\omega_\phi / \omega_r$')
plt.title(r'$\omega_\phi / \omega_r$ vs $\beta$')
plt.legend()
plt.grid()
plt.show()

