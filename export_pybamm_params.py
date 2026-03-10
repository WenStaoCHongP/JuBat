import datetime
import json
import os
import sys


def try_load_param_set():
    import pybamm

    candidates = [
        "ORegan2022",
    ]

    last_err = None
    for name in candidates:
        try:
            pv = pybamm.ParameterValues(name)
            return name, pv
        except Exception as exc:  # pylint: disable=broad-except
            last_err = exc

    raise RuntimeError(
        "Failed to load any parameter set. Last error: %r" % (last_err,)
    )


def dump_params_md(out_path, param_set_name, pv):
    # Override geometry for 21700 (diameter 21 mm, height 70 mm)
    pv.update(
        {
            "Cell diameter [m]": 0.021,
            "Cell height [m]": 0.07,
        },
        check_already_exists=False,
    )
    params = pv.copy()

    # Sort keys for stable output
    keys = sorted(params.keys())

    lines = []
    lines.append("# PyBaMM Parameters Export (LGM50T)")
    lines.append("")
    lines.append("- parameter_set: %s" % param_set_name)
    lines.append("- exported_at: %s" % datetime.datetime.now().isoformat(timespec="seconds"))
    lines.append("- geometry_override: 21700 (diameter 0.021 m, height 0.07 m)")
    lines.append("")
    lines.append("## Sources")
    lines.append("")
    lines.append("The parameter set in Table 8 has been made available in the PyBaMM software package.")
    lines.append("Further information can be found at https://www.pybamm.org/.")
    lines.append("The data repository with parameter values for the electrode solid-phase diffusivity,")
    lines.append("entropic term, open circuit voltages, exchange current density, electronic conductivity,")
    lines.append("specific heat capacity, and thermal conductivity can be found at Zenodo under DOI")
    lines.append("10.5281/zenodo.5171874.")
    lines.append("The data repository containing validation data for cells tested under different operating")
    lines.append("conditions can be found at Zenodo under DOI 10.5281/zenodo.4864437.")
    citations = params.get("citations", [])
    if citations:
        lines.append("")
        lines.append("### PyBaMM citations field")
        for citation in citations:
            lines.append("- %s" % citation)
    lines.append("")
    lines.append("## Full Parameter List")
    lines.append("")
    lines.append("| Key | Value |")
    lines.append("| --- | --- |")
    for key in keys:
        val = params[key]
        try:
            if hasattr(val, "__call__"):
                val_str = "<function>"
            else:
                val_str = repr(val)
        except Exception:  # pylint: disable=broad-except
            val_str = "<unprintable>"
        lines.append("| %s | %s |" % (key, val_str))

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def main():
    try:
        import pybamm  # noqa: F401
    except Exception as exc:  # pylint: disable=broad-except
        print("PyBaMM import failed:", exc)
        sys.exit(1)

    param_set_name, pv = try_load_param_set()

    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = os.path.join(repo_root, "output")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "pybamm_params_LGM50T.md")

    dump_params_md(out_path, param_set_name, pv)
    print("Exported parameters to", out_path)


if __name__ == "__main__":
    main()
