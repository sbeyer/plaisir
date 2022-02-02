import pandas as pd
import numpy as np

results = pd.read_csv('results.csv')

feasible = results[results.time < np.inf]
print("== feasible ==")
print(feasible.describe())

optimal = results[results.optimal]
print("== optimal ==")
print(optimal.describe())
