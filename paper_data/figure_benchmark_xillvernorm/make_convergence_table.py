#!/usr/bin/env python3
"""Extract the five final outer-iteration diagnostics; no model runs."""

import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
KEYS = ["temperature", "radiation", "local_balance", "integrated_balance", "boundary_flux"]
THRESHOLDS = [100 * math.expm1(0.003 * math.log(10)), 0.3, 0.3, 1.0, 1.0]


def percentages(row):
    # max|log10(T_new/T_old)| -> max|T_new-T_old|/min(T_new,T_old).
    return [100 * math.expm1(row[1] * math.log(10)),
            100 * row[2], 100 * row[3], 100 * row[4], 100 * abs(row[5] - 1)]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest = json.loads((HERE / "benchmark_inputs.json").read_text())
    records = []
    for run in manifest["runs"]:
        directory = ROOT / "results" / run["DAO_hash"]
        history = directory / "thermal_history.dat"
        status_path = directory / "thermal_status.json"
        status = json.loads(status_path.read_text())
        lines = history.read_text().splitlines()
        assert "max_relative_L1_dJ" in lines[0], "Expected angle-averaged run"
        rows = [[float(x) for x in line.split()]
                for line in lines if line.strip() and not line.startswith("#")]
        final = rows[-1]
        iteration = int(final[0])
        assert status["state"] == "converged"
        assert status["iteration"] == iteration == run["DAO_iteration"]
        last_three = rows[-3:]
        assert [int(row[0]) for row in last_three] == list(range(iteration - 2, iteration + 1))
        for row in last_three:
            assert all(math.isfinite(x) for x in row)
            assert all(value < limit for value, limit in zip(percentages(row), THRESHOLDS))

        budget = directory / f"thermal_budget_iter{iteration}.dat"
        budget_lines = budget.read_text().splitlines()
        cells = [[float(x) for x in line.split()]
                 for line in budget_lines if line.strip() and not line.startswith("#")]
        footer = budget_lines[-1].split()[1:]
        totals = dict(zip(footer[0::2], map(float, footer[1::2])))
        # Independently check the three energy diagnostics against the depth budget.
        columns = budget_lines[2].removeprefix("# ").split()
        relative_index, residual_index, width_index = [columns.index(k) for k in ["relative_R", "R", "dr"]]
        checks = [max(abs(cell[relative_index]) for cell in cells),
                  sum(abs(cell[residual_index]) * cell[width_index] for cell in cells) / totals["Fin"],
                  abs(totals["closure"])]
        for actual, expected in zip(checks, [final[3], final[4], abs(final[5] - 1)]):
            assert math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-13)
        spectrum = ROOT / run["source_spectrum"]
        assert sha256(spectrum) == sha256(HERE / run["DAO_file"]) == run["source_sha256"]
        records.append({
            "log_xi": run["log_xi"], "DAO_hash": run["DAO_hash"],
            "final_iteration": iteration,
            "values_percent": dict(zip(KEYS, percentages(final))),
            "last_three_iterations": [int(row[0]) for row in last_three],
            "last_three_history_rows": last_three,
            "last_three_all_five_pass": True,
            "history_columns": lines[0].removeprefix("# ").split(),
            "source_sha256": {str(p.relative_to(ROOT)): sha256(p)
                              for p in [history, status_path, budget, spectrum]},
        })

    output = {
        "units": "percent",
        "temperature_conversion": "100 * (10**max_dlogT - 1); denominator min(T_new,T_old)",
        "radiation_measure": "maximum cell-wise energy-integrated relative L1 change in J",
        "local_balance": "max_d abs(R_d)/(absorption + emission + abs(net Compton exchange))",
        "integrated_balance": "sum_d abs(R_d)*dr_d / Fin",
        "boundary_flux": "abs(Fout/Fin - 1); both slab faces included",
        "energy_diagnostics_stage": "after damping and the full-slab transfer solve",
        "strict_upper_limits_percent": dict(zip(KEYS, THRESHOLDS)),
        "runs": records,
    }
    (HERE / "convergence_values.json").write_text(json.dumps(output, indent=2) + "\n")

    table = r"""% Requires booktabs, graphicx, and xcolor; all table values are percentages.
\begin{table}[t]
\color{red}
\centering
\setlength{\tabcolsep}{3pt}
\renewcommand{\arraystretch}{1.15}
\caption{\textcolor{red}{Final values of the five outer-iteration
convergence measures. All values and limits are in per cent;
each measure must lie below its limit. The final iterations
are 63, 35, and 33 for $\log\xi=1$, 2, and 3, respectively.
All five criteria passed in the last three consecutive
outer iterations of each run.}}
\label{tab:benchmark_convergence}
\resizebox{\columnwidth}{!}{%
\begin{tabular}{@{}lrrrr@{}}
\toprule
 & & \multicolumn{3}{c}{$\log\xi$} \\
\cmidrule(l){3-5}
Measure & Limit & 1 & 2 & 3 \\
\midrule
"""
    labels = ["Temperature change", "Radiation-field change", "Local balance",
              "Integrated residual", "Boundary flux mismatch"]
    limits = ["0.6932", "0.3", "0.3", "1.0", "1.0"]
    for key, label, limit in zip(KEYS, labels, limits):
        entries = " & ".join(f"{record['values_percent'][key]:.4f}" for record in records)
        table += f"{label} & {limit} & {entries}" + r" \\" + "\n"
    table += r"""\bottomrule
\end{tabular}%
}
\end{table}
"""
    (HERE / "convergence_table.tex").write_text(table)
    for record in records:
        print(record["DAO_hash"], record["final_iteration"], record["values_percent"])
    print("Verified: all five criteria in the last three iterations; final energy budgets; spectrum hashes.")


if __name__ == "__main__":
    main()
