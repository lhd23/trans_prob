"""Run either the transitions or peak workflow from settings in configs/."""

import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
from settings import load_config, run_directory


def run(command):
    print("+", " ".join(map(str, command)), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workflow", choices=("transitions", "peak"))
    parser.add_argument("--config", type=Path, help="use an alternative configuration")
    parser.add_argument("--plot-only", action="store_true", help="reuse the saved simulation")
    parser.add_argument("--linear-y", action="store_true", help="linear transition histograms")
    args = parser.parse_args()
    if args.linear_y and args.workflow != "transitions":
        parser.error("--linear-y applies to the transitions workflow")
    config_path = (args.config or ROOT / "configs" / f"{args.workflow}.json").resolve()
    config = load_config(config_path)
    directory = run_directory(config, args.workflow)
    if not args.plot_only:
        run([sys.executable, f"src/simulate_{args.workflow}.py", "--config", str(config_path)])
    if args.workflow == "transitions":
        # Reuse an existing analysis, including any deliberately changed binning.
        if not (directory / "analysis.npz").exists():
            run([sys.executable, "src/analyse.py", str(directory)])
        command = [sys.executable, "plots/plot_transitions.py", str(directory)]
        if args.linear_y:
            command.append("--linear-y")
        run(command)
    else:
        run([sys.executable, "plots/plot_peak.py", str(directory)])


if __name__ == "__main__":
    main()
