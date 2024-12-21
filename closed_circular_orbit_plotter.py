import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

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

# Calculate omega_r and omega_theta
def omega_r_squared(x, beta, E0=1, m=1):
    return (1 / m) * u_eff_second_derivative(x, beta, E0)

def omega_theta_squared(x, beta, L=1, m=1):
    return (L**2) / (m**2 * x**4)

# Define beta range and compute values
beta_values = np.linspace(0.01, 0.42, 100)
ratios = []

for beta in beta_values:
    # Solve for x
    try:
        result = root_scalar(u_eff_derivative, args=(beta,), bracket=[0.5, 5], method='brentq')
        x = result.root
        # Calculate omega_r and omega_theta
        omega_r2 = omega_r_squared(x, beta)
        omega_theta2 = omega_theta_squared(x, beta)
        # Append the ratio
        if omega_r2 > 0 and omega_theta2 > 0:
            ratios.append(np.sqrt(omega_theta2 / omega_r2))
        else:
            ratios.append(np.nan)
    except ValueError:
        ratios.append(np.nan)

# Plot the results
plt.plot(beta_values, ratios, label=r'$\omega_\theta / \omega_r$')
plt.axhline(y=1, color='r', linestyle='--', label='Rational = 1')
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\omega_\theta / \omega_r$')
plt.title(r'$\omega_\theta / \omega_r$ vs $\beta$')
plt.legend()
plt.grid()
plt.show()
