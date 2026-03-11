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

    def get_series(candidates):
        """Try multiple PyBaMM variable names and return entries or NaNs."""
        for name in candidates:
            try:
                return sol[name].entries
            except Exception:
                continue
        return np.full_like(t, np.nan, dtype=float)

    # Heat-source breakdown (naming differs slightly across PyBaMM versions)
    q_total_vol = get_series([
        "X-averaged total heating [W.m-3]",
        "X-averaged total heating [W.m-3]",
    ])
    q_ohmic_vol = get_series([
        "X-averaged Ohmic heating [W.m-3]",
        "X-averaged ohmic heating [W.m-3]",
    ])
    q_irrev_vol = get_series([
        "X-averaged irreversible electrochemical heating [W.m-3]",
        "X-averaged irreversible heating [W.m-3]",
    ])
    q_rev_vol = get_series([
        "X-averaged reversible heating [W.m-3]",
        "X-averaged reversible electrochemical heating [W.m-3]",
    ])

    p_total = get_series([
        "Total heating [W]",
        "Volume-averaged total heating [W]",
    ])
    p_ohmic = get_series([
        "Ohmic heating [W]",
        "Volume-averaged Ohmic heating [W]",
    ])
    p_irrev = get_series([
        "Irreversible electrochemical heating [W]",
        "Irreversible heating [W]",
    ])
    p_rev = get_series([
        "Reversible heating [W]",
        "Reversible electrochemical heating [W]",
    ])

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "time [s]",
            "capacity [A.h]",
            "voltage [V]",
            "temperature [K]",
            "q_total [W.m-3]",
            "q_ohmic [W.m-3]",
            "q_irreversible [W.m-3]",
            "q_reversible [W.m-3]",
            "Q_total [W]",
            "Q_ohmic [W]",
            "Q_irreversible [W]",
            "Q_reversible [W]",
        ])
        for i in range(len(t)):
            writer.writerow([
                f"{t[i]:.6g}",
                f"{cap[i]:.6g}",
                f"{v[i]:.6g}",
                f"{temp[i]:.6g}",
                f"{q_total_vol[i]:.6g}",
                f"{q_ohmic_vol[i]:.6g}",
                f"{q_irrev_vol[i]:.6g}",
                f"{q_rev_vol[i]:.6g}",
                f"{p_total[i]:.6g}",
                f"{p_ohmic[i]:.6g}",
                f"{p_irrev[i]:.6g}",
                f"{p_rev[i]:.6g}",
            ])


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
