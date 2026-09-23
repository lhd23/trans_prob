"""Draw histograms and cumulants from the same analysed transition ensemble."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from initial_conditions import linear_growth_at


def probability_band(edges, density, mean):
    """Bounds enclosing 0.34 probability on each side of the mean.

    A sufficiently skewed distribution may not have 0.34 on each side.
    Return None in that case rather than silently changing the definition.
    """
    cumulative = np.concatenate(([0.0], np.cumsum(density * np.diff(edges))))
    at_mean = np.interp(mean, edges, cumulative)
    if at_mean < 0.34 or cumulative[-1] - at_mean < 0.34:
        return None
    return np.interp([at_mean - 0.34, at_mean + 0.34], cumulative, edges)


def clipped_histogram(edges, density, lower, upper):
    vertices = []
    for left, right, value in zip(edges[:-1], edges[1:], density):
        left, right = max(left, lower), min(right, upper)
        if right > left:
            vertices.extend(((left, value), (right, value)))
    return np.asarray(vertices)


def plot(directory, linear_y=False):
    metadata = json.loads((directory / "metadata.json").read_text())
    if metadata["workflow"] != "transitions":
        raise ValueError("this figure needs transition data")
    config = metadata["configuration"]
    length = config["box_length_mpc_over_h"]
    with np.load(directory / "analysis.npz", allow_pickle=False) as archive:
        data = {key: archive[key] for key in archive.files}
    times, targets = data["times"], data["targets"]
    if targets.size != 4:
        raise ValueError("the histogram layout has four initial-separation rows")
    indices = []
    for time in (0.1, 0.3, 0.6, 1.0):
        matches = np.flatnonzero(np.isclose(times, time, rtol=0, atol=1e-14))
        if matches.size != 1:
            raise ValueError(f"missing histogram time {time}")
        indices.append(int(matches[0]))
    tables = ROOT / "output/tables" / directory.name
    tables.mkdir(parents=True, exist_ok=True)
    edges = data["plot_edges"]
    centers = (edges[:-1] + edges[1:]) / 2
    density = data["density"][:, indices]
    normalizer = density[-1, 0].max() if linear_y else 1.0
    if normalizer <= 0:
        raise ValueError("the reference histogram has zero height")
    with np.load(directory / "spectrum.npz", allow_pickle=False) as spectrum:
        wavenumbers, power = spectrum["wavenumbers"], spectrum["power_1d_z0"]
    initial_growth = linear_growth_at(config["initial_scale_factor"],
                                     config["omega_matter"], config["omega_lambda"])[0]
    growth = [linear_growth_at(times[i], config["omega_matter"], config["omega_lambda"])[0]
              for i in indices]
    smooth_changes = np.linspace(edges[0], edges[-1], 1601)
    bands = []
    band_flags = []
    for row, target in enumerate(targets):
        with (tables / f"trajectories_{row + 1}.dat").open("w") as stream:
            for slot, trajectory in enumerate(data["trajectories"][row]):
                pair_index = data["trajectory_pair_indices"][row, slot]
                if pair_index < 0:
                    continue
                stream.write(f"# pair index within this target: {pair_index}\n")
                np.savetxt(stream, np.column_stack((data["trajectory_times"], trajectory * length)))
                stream.write("\n")
        # Gaussian relative-displacement approximation evaluated from the
        # initialized Fourier spectrum, without fitting simulation moments.
        relative_variance = 4 / length * np.sum(
            power / wavenumbers**2 * (1 - np.cos(wavenumbers * target * length))) / length**2
        for column, snapshot in enumerate(indices):
            mean = data["cumulants"][row, snapshot, 0] * target
            np.savetxt(tables / f"density_{row + 1}_{column + 1}.dat",
                       np.column_stack((centers * length, density[row, column] / normalizer,
                                        np.full(centers.size, mean * length))),
                       header="r_minus_q_mpc_over_h density mean_mpc_over_h")
            variance = (growth[column] - initial_growth)**2 * relative_variance
            model = np.exp(-smooth_changes**2 / (2 * variance)) / np.sqrt(2 * np.pi * variance)
            np.savetxt(tables / f"model_{row + 1}_{column + 1}.dat",
                       np.column_stack((smooth_changes * length, model / normalizer)))
            band = probability_band(edges, density[row, column], mean)
            band_flags.append(int(band is not None))
            bands.append(None if band is None else (band * length).tolist())
            if band is not None:
                vertices = clipped_histogram(edges, density[row, column], *band)
                vertices *= np.array([length, 1 / normalizer])
                np.savetxt(tables / f"fill_{row + 1}_{column + 1}.dat", vertices)
            else:
                print(f"q={target * length:g}, a={times[snapshot]:g}: cannot enclose 34% on each side; shading omitted")
    (tables / "probability_bands.json").write_text(json.dumps({
        "bounds_mpc_over_h_row_major": bands, "probability_on_each_side_of_mean": 0.34,
        "linear_density_normalizer": float(normalizer)}, indent=2) + "\n")

    # All signed measured cumulants are retained; the logarithmic figure
    # shows their magnitudes and omits the four extra histogram-only times.
    np.savetxt(tables / "cumulants.dat", np.column_stack(
        (times, data["cumulants"][0], data["cumulants"][-1])),
        header="a kappa1_first kappa2_first kappa3_first kappa1_last kappa2_last kappa3_last")
    pdf_directory = ROOT / "output/pdf"
    pdf_directory.mkdir(parents=True, exist_ok=True)
    suffix = "_linear" if linear_y else ""
    histogram_path = pdf_directory / f"{directory.name}_histograms{suffix}.pdf"
    shared = (f"TABLES={json.dumps(str(tables))};"
              "array Q[4] = [" + ",".join(f"{q * length:g}" for q in targets) + "];")
    variables = (shared + f"OUTPUT={json.dumps(str(histogram_path))};LINEAR_Y={int(linear_y)};"
                 f"DENSITY_MAX={max(1.05, 1.05 * density.max() / normalizer) if linear_y else 500:g};"
                 "array BANDS[16] = [" + ",".join(map(str, band_flags)) + "]")
    subprocess.run(["gnuplot", "-e", variables, str(ROOT / "plots/histograms.gp")], check=True)
    cumulant_path = pdf_directory / f"{directory.name}_cumulants.pdf"
    variables = shared + f"OUTPUT={json.dumps(str(cumulant_path))};AMIN={times[0]:g};AMAX={times[-1]:g}"
    subprocess.run(["gnuplot", "-e", variables, str(ROOT / "plots/cumulants.gp")], check=True)
    print(f"wrote {histogram_path}\nwrote {cumulant_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--linear-y", action="store_true")
    args = parser.parse_args()
    plot(args.directory.resolve(), args.linear_y)
