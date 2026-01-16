import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""Old reference formulas used in testing. Will keep here for archiving.
# Head equation used for comparison
c1 = -0.5
c2 = 14317.3

y = lambda x: c1*np.exp(-2*x) + c2*x*np.exp(-2*x) + 0.25*x + 0.5

c1 = 1.41122
c2 = 0.057534

y = lambda x: c1*np.exp(-4*x) + c2*np.exp(x) - ((3/4)*x**2 + (17/8)*x + (47/32))

"""

# Name of the solver output file
solver_file = "solution_output.dat"

dataframe = pd.read_csv(solver_file, sep="\t")
x_0 = dataframe["x"].iat[0]
x_1 = dataframe["x"].iat[-1]

x = np.linspace(x_0, x_1, 200)
#y_ref = y(x)


# Plot and output
plt.figure(figsize = (15, 10))
plt.plot(dataframe["x"], dataframe["y"], "-ok", label="Solver")
# plt.plot(x, y_ref, label="Reference") # Used for any reference solution; originally used for the deprecated functions above.
plt.legend()
plt.grid()
plt.savefig("results.png")
plt.show()






