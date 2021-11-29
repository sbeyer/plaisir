#!/usr/bin/env python3
#
# Run with arguments results/.../out_*.txt

import enum
import numpy as np
import os
import pandas as pd
import sys

results = pd.read_csv("results.csv")
results.fillna(np.inf, inplace=True)

_, *args = sys.argv

for filepath in args:
    filedir, filename = os.path.split(filepath)
    _, commit = os.path.split(filedir)

    basename, extension = os.path.splitext(filename)
    if extension != ".txt":
        continue

    prefix, instance = basename.split("_", 1)
    if prefix != "out":
        continue

    try:
        with open(filepath, "r") as infile:
            solution_data = infile.readlines()
    except Exception as err:
        print(f"Failed to read file {filepath}: {err}")
        continue

    class State(enum.Enum):
        Ignore = 1
        Days = 2
        CostA = 3
        CostB = 4
        CostTotal = 5
        CPU = 6
        Time = 7

    state = State.Ignore
    total_cost = 0
    total_time = 0
    best_cost = np.inf
    best_time = np.inf
    optimum = False
    best_changed = False
    for line in solution_data:
        line = line.rstrip()
        if state == State.Ignore:
            if line == "# Final solution":
                optimum = True
            elif line == "Day 1":
                state = State.Days
        elif state == State.Days:
            if line[:4] != "Day " and line[:6] != "Route ":
                state = State.CostA
        elif state == State.CostA:
            state = State.CostB
        elif state == State.CostB:
            state = State.CostTotal
        elif state == State.CostTotal:
            total_cost = float(line)
            if total_cost < best_cost:
                best_changed = True
                best_cost = total_cost
            state = State.CPU
        elif state == State.CPU:
            state = state.Time
        elif state == State.Time:
            total_time = float(line)
            state = State.Ignore
            if best_changed:
                best_time = total_time
                best_changed = False

    row = results.loc[results.instance == instance].squeeze()

    if best_cost != row.bestsol or best_time != row.time or optimum != row.optimal:
        print(f"NEW: {instance}\t{best_cost}\t{commit}\t{best_time}\t{optimum}")
        print(
            f"OLD: {instance}\t{row.bestsol}\t{row.commit}\t{row.time}\t{row.optimal}"
        )

    if best_cost < row.bestsol or (best_cost == row.bestsol and best_time < row.time):
        print(f" `-> IMPROVED!")
        results.loc[
            results.instance == instance, ("commit", "bestsol", "time", "optimal")
        ] = (commit, best_cost, best_time, optimum)

    print()

results.to_csv("results2.csv", index=False, float_format="%.14g")
