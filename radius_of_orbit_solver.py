import numpy as np
from scipy.optimize import fsolve

# Define the equation to solve
def equation(x, beta):
    return 2*x + 1 - 2*beta*np.exp(x)

# Given beta value
beta_value = float(input("Beta value: "))

# Initial guesses to find multiple roots
initial_guesses = [-1.0, 0.0, 1.0, 2.0]

# Solve for x using fsolve
solutions = [fsolve(equation, guess, args=(beta_value))[0] for guess in initial_guesses]

# Remove duplicates and sort solutions
unique_solutions = sorted(set(np.round(solutions, 6)))

# Display the results
print(f"The solutions for x when beta = {beta_value} are: {unique_solutions}")
