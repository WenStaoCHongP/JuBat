import os
import csv
import numpy as np

try:
    import pybamm
except ImportError as exc:
    raise SystemExit(
        "pybamm is not installed. Install it with: pip install pybamm"
    ) from exc


def run_spme_curve(c_rate, t_end_s, out_path):
    model = pybamm.lithium_ion.SPMe({"thermal": "lumped"})
    params = pybamm.ParameterValues("Chen2020")

    experiment = pybamm.Experiment(
        [f"Discharge at {c_rate}C for {t_end_s} seconds"],
        period="10 seconds",
    )

    sim = pybamm.Simulation(model, parameter_values=params, experiment=experiment)
    sol = sim.solve()

    t = sol["Time [s]"].entries
    v = sol["Terminal voltage [V]"].entries
    cap = sol["Discharge capacity [A.h]"].entries
    temp = sol["X-averaged cell temperature [K]"].entries

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time [s]", "capacity [A.h]", "voltage [V]", "temperature [K]"])
        for i in range(len(t)):
            writer.writerow([f"{t[i]:.6g}", f"{cap[i]:.6g}", f"{v[i]:.6g}", f"{temp[i]:.6g}"])


def main():
    cases = [
        (0.3, 12000.0),
        (0.7, 4800.0),
        (1.0, 3600.0),
        (2.0, 1800.0),
        (5.0, 720.0),
    ]

    base_dir = os.path.join(os.path.dirname(__file__), "src", "data")
    for c_rate, t_end in cases:
        out_file = os.path.join(base_dir, f"pybamm_SPMe_LGM50_{c_rate:.1f}C.csv")
        run_spme_curve(c_rate, t_end, out_file)
        print(f"Saved: {out_file}")


if __name__ == "__main__":
    main()
