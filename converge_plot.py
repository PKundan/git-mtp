import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("convergePlot.csv", delimiter=",");
plt.semilogy(df.iter, df.RSS, '.r', label=r"RSS for $\rho$")
plt.xlabel(f"iterations"); plt.ylabel(f"RSS")
plt.title("convergence plot")
# plt.grid()
plt.legend()
plt.show()


