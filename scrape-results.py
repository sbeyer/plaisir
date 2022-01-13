#!/usr/bin/env python3
#
# Run with arguments results/.../out_*.txt

import enum
import numpy as np
import os
import pandas as pd
import sys

bks = pd.read_csv("reference.csv")
results = pd.read_csv("results.csv")

_, *args = sys.argv


def compute_score(instance, bestsol):
    refsol_series = bks.loc[bks.instance == instance]
    if refsol_series.empty:
        return np.nan

    refsol = refsol_series.squeeze().bestsol
    return min(10, 100 * (bestsol / refsol - 1.0))


def get_solution_from_data(data):
    class State(enum.Enum):
        Ignore = 1
        Days = 2
        CostA = 3
        CostB = 4
        CostTotal = 5
        CPU = 6
        Time = 7

    output_list = []
    feasible_output = ""
    state = State.Ignore
    total_cost = 0
    total_time = 0
    best_cost = np.inf
    best_time = np.inf
    optimum = False
    best_changed = False
    for line in solution_data:
        line = line.rstrip()

        if output_list:
            output_list.append(line)

        if state == State.Ignore:
            if line == "# Final solution":
                optimum = True
            elif line == "Day 1":
                output_list.append(line)
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
                feasible_output = "\n".join(output_list)
            output_list = []

    return {
        "bestsol": best_cost,
        "time": best_time,
        "optimal": optimum,
        "output": feasible_output,
    }


def save_output(instance, content):
    if content:
        try:
            latest_solution_file = f"results_best/out_{instance}.txt"
            with open(latest_solution_file, "w") as outfile:
                outfile.writelines(content)
        except Exception as err:
            print(f"Failed to write solution to file {latest_solution_file}: {err}")


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

    solution = get_solution_from_data(solution_data)

    row = results.loc[results.instance == instance].squeeze()

    update = False
    if row.empty:
        print(
            f"NEW: {instance}\t{solution['bestsol']}\t{commit}\t{solution['time']}\t{solution['optimal']}"
        )
        print(" `-> NEWLY INSERTED!")
        print()
        results.loc[
            len(results.index), ("instance", "commit", "bestsol", "time", "optimal")
        ] = (
            instance,
            commit,
            solution["bestsol"],
            solution["time"],
            solution["optimal"],
        )

        save_output(instance, solution["output"])
    elif (
        solution["bestsol"] != row.bestsol
        or solution["time"] != row.time
        or solution["optimal"] != row.optimal
    ):
        print(
            f"NEW: {instance}\t{solution['bestsol']}\t{commit}\t{solution['time']}\t{solution['optimal']}"
        )
        print(
            f"OLD: {instance}\t{row.bestsol}\t{row.commit}\t{row.time}\t{row.optimal}"
        )

        if solution["bestsol"] <= row.bestsol:
            if solution["bestsol"] == row.bestsol and solution["time"] > row.time:
                print(" `-> SOLUTION KEPT ALTHOUGH WORSE TIME!")
            elif row.optimal and not solution["optimal"]:
                print(
                    " `-> SOLUTION KEPT ALTHOUGH WE LOST PROOF OF OPTIMALITY! (We keep info of optimality)"
                )
                solution["optimal"] = True
            else:
                print(" `-> IMPROVED!")

            score = compute_score(instance, solution['bestsol'])
            results.loc[
                results.instance == instance,
                ("commit", "bestsol", "time", "optimal", "score"),
            ] = (
                commit,
                solution["bestsol"],
                solution["time"],
                solution["optimal"],
                score,
            )

            save_output(instance, solution["output"])

        print()

results.to_csv("results2.csv", index=False, float_format="%.14g")
