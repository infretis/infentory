#!/usr/bin/env python3
"""
Standalone script to compute and plot histograms and free energies from infretis data.

Computes WHAM-weighted histograms and per-ensemble data, then generates plots.
Excludes WHAM vs 0+ comparison plot.

Usage:
    python standalone_histograms.py --toml infretis.toml --data infretis_data.txt --trajdir load
    python standalone_histograms.py --help
"""
import argparse
import os
import gzip
from pathlib import Path
from typing import Optional, List, Dict, Tuple
import numpy as np
import matplotlib.pyplot as plt

try:
    import tomli
except ImportError:
    print("Error: tomli package required. Install with: pip install tomli")
    exit(1)

# Try to import scienceplots for publication-ready figures
try:
    import scienceplots  # type: ignore
    SCIENCEPLOTS_AVAILABLE = True
except Exception:
    SCIENCEPLOTS_AVAILABLE = False

# Default color cycle
COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]


# ============================================================================
# INFTOOLS HELPER FUNCTIONS (extracted for standalone use)
# ============================================================================

def read_toml(toml_path):
    """Read TOML configuration file."""
    toml = Path(toml_path)
    if toml.exists():
        with open(toml, "rb") as rfile:
            config = tomli.load(rfile)
        return config
    else:
        return False


def data_reader(inp):
    """Read infretis_data.txt type of file."""
    paths = []

    # Handle .gz files
    if inp[-3:] == ".gz":
        oopen = gzip.open
        readmode = "rt"
    else:
        oopen = open
        readmode = "r"

    with oopen(inp, readmode) as read:
        ensl = 0

        for line in read:
            if line.startswith("#"):
                continue

            rip = line.rstrip().split()
            # Get num ensembles from first non-comment line
            if not ensl:
                ensl = int(len(rip[3:]) / 2)

            # Get line data
            pn, len0, max_op = rip[:3]
            path = {"pn": pn, "len": len0, "max_op": max_op}
            path["cols"] = {}
            f0l, w0l = rip[3:ensl + 3], rip[3 + ensl:2 * ensl + 3]

            # Skip if no weights
            if set(f0l) == set(w0l) == set(("----",)):
                continue

            # Store only the weights based on col
            for col, (f0, w0) in enumerate(zip(f0l, w0l)):
                if '----' in (f0, w0):
                    continue
                path["cols"][col] = [f0, w0]

            paths.append(path)
        return paths


def get_WHAMfactors(matrix, lambda_interfaces, i0plus, Q, lm1):
    """Compute WHAM factors for path weighting."""
    imax = 2  # index where lambda_max is stored
    intfQ = lambda_interfaces[:-1]  # interfaces for determining Q-index
    numC = len(intfQ)
    WHAMfactors = []

    for x in matrix:
        lmax = x[imax]
        indexQ = max(
            (i for i, val in enumerate(intfQ) if val < lmax), default=-1
        )
        if indexQ == -1:
            if lmax == intfQ[0]:
                indexQ == 0
            elif lm1 is not None:
                indexQ == 0
            else:
                print("Error: lambda_max is lower or equal to all TIS interfaces")
                print("data line=", x)
                print("lmax=", lmax)
                print("interfaces except last: ", intfQ)
                exit()
        Qmax = Q[indexQ]
        sumC = sum(x[i0plus : i0plus + numC])
        Chi_X = Qmax * sumC
        WHAMfactors.append(Chi_X)
    return WHAMfactors


def extract(trajfile, xcol, ycol=None):
    """Extract trajectory data from order.txt file."""
    traj = np.loadtxt(trajfile)
    data = traj[1:-1, xcol]  # remove first and last frames
    if ycol is not None:
        data = np.vstack((data, traj[1:-1, ycol]))
    return data


def update_histogram(data, factor, histogram, Minx, Miny, dx, dy):
    """Update histogram with weighted data."""
    if Miny is not None and dy is not None:
        x = data[0]
        y = data[1]
        ix = ((x - Minx) / dx).astype(int)
        iy = ((y - Miny) / dy).astype(int)
        np.add.at(histogram, (ix, iy), factor)
    else:
        x = data if data.ndim == 1 else data[:, 0]
        mask = (x >= Minx) & (x <= Minx + len(histogram) * dx)
        ix = ((x[mask] - Minx) / dx).astype(int)
        ix[ix == len(histogram)] = len(histogram) - 1
        np.add.at(histogram, ix, factor)

    return histogram


# ============================================================================
# HISTOGRAM COMPUTATION
# ============================================================================

def compute_all_histograms(
    toml: str = "infretis.toml",
    data: str = "infretis_data.txt",
    trajdir: str = "load",
    outdir: str = "histograms",
    nskip: int = 0,
    dlambda: Optional[float] = None,
    nbins: Optional[int] = None,
    lmin: Optional[float] = None,
    lmax: Optional[float] = None,
    xcol: int = 1,
    dt: Optional[float] = None,
    subcycles: Optional[int] = None,
    normalize: str = "none",
    lm1: bool = True,
):
    """
    Compute histograms and free energies for an infretis simulation.

    Parameters
    ----------
    toml : str
        Path to infretis.toml configuration file.
    data : str
        Path to infretis_data.txt file.
    trajdir : str
        Directory containing path folders.
    outdir : str
        Output directory for CSV files.
    nskip : int
        Skip first N paths.
    dlambda : float, optional
        Bin width (mutually exclusive with nbins).
    nbins : int, optional
        Number of bins (mutually exclusive with dlambda).
    lmin : float, optional
        Minimum order parameter.
    lmax : float, optional
        Maximum order parameter.
    xcol : int
        Order parameter column in order.txt.
    dt : float, optional
        Time step (auto-detected from toml if not set).
    subcycles : int, optional
        MD subcycles per frame (auto-detected from toml if not set).
    normalize : str
        Normalization: none, time, density, or probability.
    lm1 : bool
        Use lambda_-1 correction for WHAM.

    Returns
    -------
    dict
        Info about computed histograms.
    """
    # Create output directory
    outdir_path = Path(outdir)
    try:
        outdir_path.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        print(f"Warning: Could not create output directory {outdir}: {e}")
        print("Using current directory instead...")
        outdir = "."
        outdir_path = Path(".")

    print(f"Computing histograms for simulation in {os.path.dirname(toml) or '.'}")
    print(f"Output directory: {outdir}")

    # Get interfaces from toml
    toml_dict = read_toml(toml)
    interfaces = [float(i) for i in toml_dict["simulation"]["interfaces"]]
    nintf = len(interfaces)
    nplus_ens = nintf - 1

    # Get dt and subcycles from toml if not provided
    if dt is None:
        dt = float(toml_dict.get("engine", {}).get("timestep", 1.0))
        print(f"Auto-detected dt={dt} from toml")
    if subcycles is None:
        subcycles = int(toml_dict.get("engine", {}).get("subcycles", 1))
        print(f"Auto-detected subcycles={subcycles} from toml")

    # Configuration
    i0plus, i0min = 4, 3
    lamres = 0.005
    lambdaA, lambdaB = interfaces[0], interfaces[-1]

    # Check for lambda_-1
    lambda_m1 = None
    if lm1 and "tis_set" in toml_dict.get("simulation", {}):
        lambda_m1 = toml_dict["simulation"]["tis_set"].get("lambda_minus_one")
        if lambda_m1 is not None:
            print(f"Using lambda_-1 = {lambda_m1}")

    # Determine binning
    if dlambda is None and nbins is None:
        nbins = 100
    if dlambda is not None and nbins is not None:
        raise ValueError("Specify either -dlambda or -nbins, not both.")

    bin_lmin = lmin if lmin is not None else lambdaA - 0.3 * (lambdaB - lambdaA)
    bin_lmax = lmax if lmax is not None else lambdaB + 0.1 * (lambdaB - lambdaA)

    if nbins is not None:
        actual_nbins = nbins
        actual_dlambda = (bin_lmax - bin_lmin) / nbins
    else:
        if dlambda is None or dlambda <= 0:
            raise ValueError("Specify a positive -dlambda when -nbins is not provided")
        actual_dlambda = dlambda
        actual_nbins = int(np.floor((bin_lmax - bin_lmin) / actual_dlambda))
        if actual_nbins <= 0:
            raise ValueError("Computed number of bins is non-positive")
        bin_lmax = bin_lmin + actual_nbins * actual_dlambda

    print(f"Binning: {actual_nbins} bins, width={actual_dlambda:.6f}")
    print(f"  Range: [{bin_lmin:.6f}, {bin_lmax:.6f}]")

    # Create histo_stuff dict
    histo_stuff = {
        "nbx": actual_nbins,
        "minx": bin_lmin,
        "maxx": bin_lmax,
        "xcol": xcol,
        "nby": None,
        "miny": None,
        "maxy": None,
        "ycol": None,
    }

    # Read and process data matrix
    print("\n--- Reading data matrix ---")
    print(f"  Reading from: {data}")
    with open(data) as f:
        matrix = [
            [float(x) if x != "----" else 0.0 for x in line.strip().split()]
            for line in f if not line.startswith("#")
        ][nskip:]
    print(f"  Loaded {len(matrix)} paths (skipped first {nskip})")

    # Initialize eta and v_alpha
    print("\n--- Initializing WHAM computation ---")
    eta = [0.0] * nplus_ens
    lambda_values = [i * lamres for i in range(round(lambdaA / lamres), round(lambdaB / lamres) + 1)]
    v_alpha = [0.0] * len(lambda_values)
    v_alpha[0] = 1.0
    print(f"  Interfaces: {[f'{intf:.6f}' for intf in interfaces]}")
    print(f"  n_ensembles: {nintf}, n_plus_ensembles: {nplus_ens}")

    # Unweight matrix with HA-weights
    print("  Unweighting matrix with HA-weights...")
    sumPxy = [0.0] * nintf
    sumPxy_afterw = [0.0] * nintf

    for x in matrix:
        for y in range(nintf):
            y1, y2 = i0min + y, i0min + y + nintf
            P_xy = x[y1]
            sumPxy[y] += P_xy
            x[y1] = P_xy / x[y2] if x[y2] > 0 else 0.0
            sumPxy_afterw[y] += x[y1]

    # Normalize by average inverse HA-weight
    for y in range(nintf):
        if sumPxy[y] > 0:
            y1 = i0min + y
            AvinvwHA = sumPxy_afterw[y] / sumPxy[y]
            for x in matrix:
                x[y1] /= AvinvwHA

    # Compute eta and v_alpha from paths
    for x in matrix:
        lambdamax = x[2]
        for i in range(nplus_ens):
            Cxy = x[i0plus + i]
            eta[i] += Cxy

            lambda_i = interfaces[i]
            alpha_max = int(np.floor((lambdamax - lambdaA) / lamres))
            alpha_min = round((lambda_i - lambdaA) / lamres)
            if alpha_max > len(v_alpha) - 1:
                alpha_max = len(v_alpha) - 1
            alpha_min += 1
            for alpha in range(alpha_min, alpha_max + 1):
                v_alpha[alpha] += Cxy

    # Compute Q factors for WHAM
    def WHAM_PQ(npe, interf, res, eta, v_alpha):
        P, Q, invQ = [0.0] * npe, [0.0] * npe, [0.0] * npe
        P[0], invQ[0] = 1.0, eta[0]
        if invQ[0] == 0:
            return P, Q
        Q[0] = 1 / invQ[0]

        for i in range(1, npe):
            alpha = round((interf[i] - interf[0]) / res)
            if alpha >= len(v_alpha):
                alpha = len(v_alpha) - 1
            P[i] = v_alpha[alpha] * Q[i - 1]
            if P[i] == 0:
                return P, Q
            invQ[i] = invQ[i - 1] + eta[i] / P[i]
            Q[i] = 1 / invQ[i]
        return P, Q

    Pi0_wham, Q = WHAM_PQ(nplus_ens, interfaces, lamres, eta, v_alpha)

    print("\n--- Computing WHAM factors ---")
    print(f"  eta (sampling per ensemble): {[f'{e:.2f}' for e in eta]}")
    print(f"  Q factors: {[f'{q:.6f}' for q in Q]}")
    print(f"  Pi0_wham (crossing probs): {[f'{p:.6e}' for p in Pi0_wham]}")

    WHAMfactors = get_WHAMfactors(matrix, interfaces, i0plus, Q, lambda_m1 if lm1 else None)
    sum_wham = sum(WHAMfactors)
    WHAMfactors = [w / sum_wham for w in WHAMfactors] if sum_wham > 0 else WHAMfactors

    # Get [0-] ensemble weights
    WHAMfactors_0min = [x[i0min] for x in matrix]
    sum_0min = sum(WHAMfactors_0min)
    WHAMfactors_0min = [w / sum_0min for w in WHAMfactors_0min] if sum_0min > 0 else WHAMfactors_0min

    trajlabels = [int(x[0]) for x in matrix]
    dt_frame = dt * subcycles

    print(f"Processed {len(trajlabels)} paths")
    print(f"  [i+]: {sum(1 for w in WHAMfactors if w > 0)} paths with non-zero weight")
    print(f"  [0-]: {sum(1 for w in WHAMfactors_0min if w > 0)} paths with non-zero weight")
    print(f"  dt_frame = dt * subcycles = {dt} * {subcycles} = {dt_frame}")
    print(f"  Normalization mode: {normalize}")

    # Path data cache
    path_cache = {}

    def get_filtered_trajectory(path_nr, max_x=None):
        """Load trajectory and optionally filter by max_x."""
        cache_key = (path_nr, max_x)
        if cache_key in path_cache:
            return path_cache[cache_key]

        trajfile = os.path.join(trajdir, str(path_nr), "order.txt")
        if not os.path.exists(trajfile):
            path_cache[cache_key] = np.array([])
            return path_cache[cache_key]

        try:
            traj_data = np.asarray(extract(trajfile, xcol))
            if max_x is not None:
                traj_data = traj_data[traj_data < max_x]
            path_cache[cache_key] = traj_data
            return traj_data
        except Exception as e:
            print(f"Warning: Could not process {trajfile}: {e}")
            path_cache[cache_key] = np.array([])
            return path_cache[cache_key]

    # Helper to compute histogram with given weights
    def compute_histogram(weights, label, do_normalize=True, max_x_allowed=None):
        Nbinsx = histo_stuff["nbx"]
        Minx, Maxx = histo_stuff["minx"], histo_stuff["maxx"]
        dx = (Maxx - Minx) / Nbinsx
        edges = np.linspace(Minx, Maxx, Nbinsx + 1)
        xval = 0.5 * (edges[:-1] + edges[1:])
        histogram = np.zeros(Nbinsx)

        print(f"  {label}:")
        print(f"    Bin edges: [{Minx:.6f}, {Maxx:.6f}], Nbins={Nbinsx}, dx={dx:.6g}")
        if max_x_allowed is not None:
            print(f"    Filtering: keeping only frames with x < {max_x_allowed:.6f}")

        total_frames = 0.0
        npaths_weight = 0.0
        frames_filtered_out = 0.0

        for path_nr, weight in zip(trajlabels, weights):
            if weight == 0:
                continue

            traj_data = get_filtered_trajectory(path_nr, max_x_allowed)
            if len(traj_data) == 0:
                continue

            if max_x_allowed is not None:
                original_data = get_filtered_trajectory(path_nr, None)
                frames_filtered_out += weight * (len(original_data) - len(traj_data))

            traj_data = traj_data[(traj_data >= Minx) & (traj_data < Maxx)]
            if len(traj_data) == 0:
                continue

            histogram = update_histogram(traj_data, weight, histogram, Minx, None, dx, None)
            total_frames += weight * len(traj_data)
            npaths_weight += weight

        # Mask forbidden bins
        if max_x_allowed is not None:
            mask = xval >= max_x_allowed
            n_masked_bins = np.sum(mask)
            weight_before_mask = np.sum(histogram[mask])
            histogram[mask] = 0.0
            if n_masked_bins > 0:
                print(f"    Masked {n_masked_bins} bins (x >= {max_x_allowed:.6f}), removed weight={weight_before_mask:.6g}")
            if frames_filtered_out > 0:
                print(f"    Filtered out {frames_filtered_out:.6g} weighted frames")

        # Normalize histogram
        npaths = npaths_weight if npaths_weight > 0 else 1.0
        total_weight_sum = np.sum(histogram)

        print(f"    Effective paths: {npaths_weight:.6g}, weighted frames: {total_frames:.6g}")
        print(f"    Total histogram weight: {total_weight_sum:.6g}")

        if do_normalize and total_weight_sum > 0:
            if normalize == "time":
                histogram = histogram / total_weight_sum * total_frames * dt_frame / npaths / dx
                print(f"    Time normalization applied")
            elif normalize == "density":
                time_ens = sum(
                    weight * max(0, len(get_filtered_trajectory(pn, max_x_allowed)) - 2)
                    for pn, weight in zip(trajlabels, weights) if weight > 0
                ) * dt_frame

                if time_ens > 0:
                    histogram = histogram / total_weight_sum * total_frames * dt_frame / npaths / dx / time_ens
                    print(f"    Density normalization: time_ens={time_ens:.6g}")
                else:
                    histogram = histogram / total_weight_sum / dx
                    print(f"    Density normalization: time_ens=0")
                integral = np.sum(histogram * dx)
                print(f"    Integral={integral:.6g}")
            elif normalize == "probability":
                histogram = histogram / total_weight_sum
                print(f"    Probability normalization: sum={np.sum(histogram):.6g}")

        return xval, histogram

    # Helper to compute free energy
    def compute_free_energy(histogram):
        max_value = np.max(histogram)
        if max_value > 0:
            prob = histogram / max_value
        else:
            prob = histogram.copy()

        with np.errstate(divide="ignore"):
            fe = -np.log(prob)
        fe[np.isinf(fe)] = np.nan

        return prob, fe

    # Determine column header
    hist_header = {
        "time": "avg_time_per_path",
        "density": "probability_density",
        "probability": "probability",
    }.get(normalize, "weighted_counts")

    # Helper to save files
    def save_histogram_files(prefix, xval, histogram, prob, fe):
        hist_file = os.path.join(outdir, f"{prefix}_histogram.csv")
        fe_file = os.path.join(outdir, f"{prefix}_free_energy.csv")
        np.savetxt(hist_file, np.c_[xval, histogram, prob],
                   delimiter=",", header=f"order_parameter,{hist_header},probability", comments="")
        np.savetxt(fe_file, np.c_[xval, fe],
                   delimiter=",", header="order_parameter,free_energy_kBT", comments="")
        print(f"  Saved {prefix}_histogram.csv and {prefix}_free_energy.csv")

    # Compute WHAM [i+] histogram
    print("\n--- Computing WHAM [i+] histogram ---")
    xval, hist_wham_plus = compute_histogram(WHAMfactors, "WHAM [i+]")
    prob_wham_plus, fe_wham_plus = compute_free_energy(hist_wham_plus)
    save_histogram_files("wham_plus", xval, hist_wham_plus, prob_wham_plus, fe_wham_plus)

    # Compute [0-] ensemble histogram
    print("\n--- Computing [0-] ensemble histogram ---")
    xval, hist_0min = compute_histogram(WHAMfactors_0min, "[0-]", max_x_allowed=interfaces[0])
    prob_0min, fe_0min = compute_free_energy(hist_0min)
    save_histogram_files("ens_0min", xval, hist_0min, prob_0min, fe_0min)

    # Compute per-ensemble histograms
    print("\n--- Computing per-ensemble histograms (without WHAM) ---")

    paths = data_reader(data)[nskip:]
    ens_paths = {i: [] for i in range(nintf)}
    ens_weights = {i: [] for i in range(nintf)}
    print(f"  Building ensemble path lists from {len(paths)} paths...")

    for path in paths:
        pn = int(path["pn"])
        for col, (f0, w0) in path["cols"].items():
            f0, w0 = float(f0), float(w0)
            if w0 > 0:
                ens_paths[col].append(pn)
                ens_weights[col].append(f0 / w0)

    # Normalize weights per ensemble
    for i in range(nintf):
        w = np.array(ens_weights[i])
        if len(w) > 0 and np.sum(w) > 0:
            ens_weights[i] = w / np.sum(w)
            ens_label = f"[{i}-]" if i == 0 else f"[{i-1}+]"
            print(f"  Ensemble {i} ({ens_label}): {len(w)} paths, sum(weights)={np.sum(ens_weights[i]):.6f}")

    all_ens_hists = []
    all_ens_fes = []
    headers_hist = []
    headers_fe = []

    # Compute histogram for each ensemble
    def compute_ensemble_histogram(ens_idx, path_nrs, weights):
        ens_label = f"[{ens_idx}-]" if ens_idx == 0 else f"[{ens_idx-1}+]"
        max_x = interfaces[0] if ens_idx == 0 else None

        xval, histogram = compute_histogram(weights, f"Ensemble {ens_idx} ({ens_label})", do_normalize=True, max_x_allowed=max_x)
        prob, fe = compute_free_energy(histogram)

        # Save files
        np.savetxt(
            os.path.join(outdir, f"ens_{ens_idx:03d}_histogram.csv"),
            np.c_[xval, histogram, prob],
            delimiter=",",
            header=f"order_parameter,{hist_header}_{ens_label},probability_{ens_label}",
            comments="",
        )
        np.savetxt(
            os.path.join(outdir, f"ens_{ens_idx:03d}_free_energy.csv"),
            np.c_[xval, fe],
            delimiter=",",
            header=f"order_parameter,free_energy_{ens_label}_kBT",
            comments="",
        )

        return histogram, fe, ens_label

    for ens_idx in range(nintf):
        if len(ens_paths[ens_idx]) == 0:
            ens_label = f"[{ens_idx}-]" if ens_idx == 0 else f"[{ens_idx-1}+]"
            print(f"  Ensemble {ens_idx} ({ens_label}): no paths")
            continue

        # Temporarily swap trajlabels
        orig_trajlabels = trajlabels
        trajlabels = ens_paths[ens_idx]

        hist, fe, label = compute_ensemble_histogram(ens_idx, ens_paths[ens_idx], ens_weights[ens_idx])

        trajlabels = orig_trajlabels

        all_ens_hists.append(hist)
        all_ens_fes.append(fe)
        headers_hist.append(label)
        headers_fe.append(label)

    # Save combined per-ensemble data
    if all_ens_hists:
        np.savetxt(
            os.path.join(outdir, "all_ensembles_histogram.csv"),
            np.column_stack([xval] + all_ens_hists),
            delimiter=",", header="order_parameter," + ",".join(headers_hist), comments="",
        )
        np.savetxt(
            os.path.join(outdir, "all_ensembles_free_energy.csv"),
            np.column_stack([xval] + all_ens_fes),
            delimiter=",", header="order_parameter," + ",".join(headers_fe), comments="",
        )
        print(f"  Saved combined data with {len(all_ens_hists)} ensembles")

    print("\n--- Histogram computation done ---")
    print(f"All output files saved to: {outdir}")
    return {"outdir": outdir, "interfaces": interfaces, "n_paths": len(trajlabels)}


# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

def load_csv(filepath: str) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """Load a CSV file with header."""
    with open(filepath, "r") as f:
        header_line = f.readline().strip()
        headers = header_line.split(",")

    data = np.loadtxt(filepath, delimiter=",", skiprows=1)
    if data.ndim == 1:
        data = data.reshape(-1, 1)

    x = data[:, 0]
    y = data[:, 1:] if data.shape[1] > 1 else data[:, 0:1]

    return x, y, headers[1:] if len(headers) > 1 else headers


def detect_normalization(headers: List[str]) -> str:
    """Detect normalization type from CSV headers."""
    header_str = ",".join(headers).lower()
    if "time" in header_str:
        return "time"
    elif "probability_density" in header_str:
        return "density"
    else:
        return "none"


def load_histogram_data(datadir: str) -> Dict:
    """Load all histogram data from a directory."""
    datadir = Path(datadir)
    result = {"datadir": str(datadir), "normalization": "none"}

    # Load WHAM [i+] ensembles data
    wham_plus_hist_file = datadir / "wham_plus_histogram.csv"
    wham_plus_fe_file = datadir / "wham_plus_free_energy.csv"

    if wham_plus_hist_file.exists():
        x, y, headers = load_csv(wham_plus_hist_file)
        result["wham_plus_x"] = x
        result["wham_plus_hist"] = y[:, 0] if y.ndim > 1 else y.flatten()
        if y.shape[1] > 1:
            result["wham_plus_prob"] = y[:, 1]
        result["normalization"] = detect_normalization(headers)
    if wham_plus_fe_file.exists():
        x, y, _ = load_csv(wham_plus_fe_file)
        result["wham_plus_fe"] = y[:, 0] if y.ndim > 1 else y.flatten()

    # Load [0-] ensemble data
    ens_0min_hist_file = datadir / "ens_0min_histogram.csv"
    ens_0min_fe_file = datadir / "ens_0min_free_energy.csv"

    if ens_0min_hist_file.exists():
        x, y, headers = load_csv(ens_0min_hist_file)
        result["ens_0min_x"] = x
        result["ens_0min_hist"] = y[:, 0] if y.ndim > 1 else y.flatten()
        if y.shape[1] > 1:
            result["ens_0min_prob"] = y[:, 1]
    if ens_0min_fe_file.exists():
        x, y, _ = load_csv(ens_0min_fe_file)
        result["ens_0min_fe"] = y[:, 0] if y.ndim > 1 else y.flatten()

    # Load all-ensembles data
    all_hist_file = datadir / "all_ensembles_histogram.csv"
    all_fe_file = datadir / "all_ensembles_free_energy.csv"

    if all_hist_file.exists():
        x, y, headers = load_csv(all_hist_file)
        result["ens_x"] = x
        result["ens_hist"] = y
        result["ens_labels"] = headers

    if all_fe_file.exists():
        x, y, headers = load_csv(all_fe_file)
        result["ens_fe"] = y
        result["ens_fe_labels"] = headers

    return result


def plot_wham_results(
    data: Dict,
    interfaces: Optional[List[float]] = None,
    title: str = "",
    figsize: Tuple[float, float] = (14, 5),
    save_path: Optional[str] = None,
    show: bool = True,
    log_scale: bool = False,
) -> plt.Figure:
    """Plot WHAM histogram and free energy side by side."""
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    has_new_format = "wham_plus_x" in data and "ens_0min_x" in data

    # Determine y-axis label
    normalization = data.get("normalization", "none")
    if normalization == "time":
        ylabel_hist = "Time per dlambda per path"
        title_hist = "WHAM Histogram (Time Normalized)"
    elif normalization == "density":
        ylabel_hist = "Probability density"
        title_hist = "WHAM Histogram (Density Normalized)"
    else:
        ylabel_hist = "Probability P(λ)"
        title_hist = "WHAM Histogram (Probability)"

    # Plot histogram
    ax = axes[0]
    if has_new_format:
        x_0min = data.get("ens_0min_x", np.array([]))
        hist_0min = data.get("ens_0min_hist", np.array([]))
        x_plus = data.get("wham_plus_x", np.array([]))
        hist_plus = data.get("wham_plus_hist", np.array([]))

        if len(x_0min) > 0 and len(hist_0min) > 0:
            ax.plot(x_0min, hist_0min, color="tab:green", lw=2, marker="o", markersize=3, label="[0-] ensemble")
        if len(x_plus) > 0 and len(hist_plus) > 0:
            ax.plot(x_plus, hist_plus, color="gold", lw=2, marker="o", markersize=3, label="WHAM [i+] ensembles")
        ax.legend(loc="best")

    ax.set_xlabel("Order parameter λ")
    ax.set_ylabel(ylabel_hist)
    if log_scale:
        ax.set_yscale("log")
    ax.set_title(title_hist)
    ax.grid(True, alpha=0.3)

    if interfaces:
        for intf in interfaces:
            ax.axvline(intf, color="k", alpha=0.3, ls="--")

    # Plot free energy
    ax = axes[1]
    if has_new_format:
        x_0min = data.get("ens_0min_x", np.array([]))
        fe_0min = data.get("ens_0min_fe", np.array([]))
        x_plus = data.get("wham_plus_x", np.array([]))
        fe_plus = data.get("wham_plus_fe", np.array([]))

        if len(x_0min) > 0 and len(fe_0min) > 0:
            ax.plot(x_0min, fe_0min, color="tab:green", lw=2, marker="o", markersize=3, label="[0-] ensemble")
        if len(x_plus) > 0 and len(fe_plus) > 0:
            ax.plot(x_plus, fe_plus, color="gold", lw=2, marker="o", markersize=3, label="WHAM [i+] ensembles")
        ax.legend(loc="best")

    ax.set_xlabel("Order parameter λ")
    ax.set_ylabel("Free energy F(λ) [kBT]")
    ax.set_title("Conditional free energy")
    ax.grid(True, alpha=0.3)

    if interfaces:
        for intf in interfaces:
            ax.axvline(intf, color="k", alpha=0.3, ls="--")

    if title:
        fig.suptitle(title, fontsize=14)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")

    if show:
        plt.show()

    return fig


def plot_ensemble_histograms(
    data: Dict,
    interfaces: Optional[List[float]] = None,
    title: str = "",
    figsize: Tuple[float, float] = (14, 5),
    save_path: Optional[str] = None,
    show: bool = True,
    log_scale: bool = False,
) -> plt.Figure:
    """Plot per-ensemble histograms and free energies."""
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    x = data.get("ens_x", np.array([]))
    hist = data.get("ens_hist", np.array([]))
    fe = data.get("ens_fe", np.array([]))
    labels = data.get("ens_labels", [])

    n_ens = hist.shape[1] if hist.ndim > 1 else 1

    # Determine y-axis label
    normalization = data.get("normalization", "none")
    if normalization == "time":
        ylabel_hist = "Time per bin width"
        title_hist = "Per-Ensemble Histograms (Time Normalized)"
    elif normalization == "density":
        ylabel_hist = "Probability density"
        title_hist = "Per-Ensemble Histograms (Density Normalized)"
    else:
        ylabel_hist = "Probability P(λ)"
        title_hist = "Per-Ensemble Histograms"

    # Plot histograms
    ax = axes[0]
    for i in range(n_ens):
        h = hist[:, i] if hist.ndim > 1 else hist
        label = labels[i] if i < len(labels) else f"Ens {i}"
        color = COLORS[i % len(COLORS)]
        ax.plot(x, h, color=color, lw=2, marker="o", markersize=2, label=label, alpha=0.8)

    ax.set_xlabel("Order parameter λ")
    ax.set_ylabel(ylabel_hist)
    if log_scale:
        ax.set_yscale("log")
    ax.set_title(title_hist)
    ax.legend(loc="best", fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)

    if interfaces:
        for intf in interfaces:
            ax.axvline(intf, color="k", alpha=0.2, ls="--")

    # Plot free energies
    ax = axes[1]
    for i in range(n_ens):
        f = fe[:, i] if fe.ndim > 1 else fe
        label = labels[i] if i < len(labels) else f"Ens {i}"
        color = COLORS[i % len(COLORS)]
        ax.plot(x, f, color=color, lw=1.5, marker="o", markersize=2, label=label, alpha=0.8)

    ax.set_xlabel("Order parameter λ")
    ax.set_ylabel("Free Energy F(λ) [kBT]")
    ax.set_title("Per-Ensemble Free Energies")
    ax.legend(loc="best", fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)

    if interfaces:
        for intf in interfaces:
            ax.axvline(intf, color="k", alpha=0.2, ls="--")

    if title:
        fig.suptitle(title, fontsize=14)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")

    if show:
        plt.show()

    return fig


# ============================================================================
# MAIN
# ============================================================================

def main():
    """Main function for standalone execution."""
    parser = argparse.ArgumentParser(
        description="Compute and plot histograms and free energies from infretis data"
    )

    # Computation arguments
    parser.add_argument("--toml", default="infretis.toml", help="Path to infretis.toml")
    parser.add_argument("--data", default="infretis_data.txt", help="Path to infretis_data.txt")
    parser.add_argument("--trajdir", default="load", help="Directory with path folders")
    parser.add_argument("--outdir", default="histograms", help="Output directory for CSV files")
    parser.add_argument("--nskip", type=int, default=0, help="Skip first N paths")
    parser.add_argument("--dlambda", type=float, help="Bin width (mutually exclusive with --nbins)")
    parser.add_argument("--nbins", type=int, help="Number of bins (mutually exclusive with --dlambda)")
    parser.add_argument("--lmin", type=float, help="Minimum order parameter")
    parser.add_argument("--lmax", type=float, help="Maximum order parameter")
    parser.add_argument("--xcol", type=int, default=1, help="Order parameter column in order.txt")
    parser.add_argument("--dt", type=float, help="Time step (auto-detected if not set)")
    parser.add_argument("--subcycles", type=int, help="MD subcycles per frame (auto-detected if not set)")
    parser.add_argument("--normalize", default="none", choices=["none", "time", "density", "probability"],
                        help="Normalization mode")

    # Plotting arguments
    parser.add_argument("--plot", action="store_true", help="Generate plots after computation")
    parser.add_argument("--plot-only", action="store_true", help="Skip computation, only plot existing data")
    parser.add_argument("--interfaces", help='Comma-separated interface values. Use: --interfaces=-0.1,0.0,0.1 or --interfaces "-0.1,0.0,0.1"')
    parser.add_argument("--title", default="", help="Plot title")
    parser.add_argument("--save-plots", help="Directory to save plot figures")
    parser.add_argument("--no-show", action="store_true", help="Don't display figures interactively")
    parser.add_argument("--log", action="store_true", help="Use log scale for histograms")

    args = parser.parse_args()

    # Parse interfaces
    intf_list = None
    if args.interfaces:
        intf_list = [float(x.strip()) for x in args.interfaces.split(",")]

    # Compute histograms (unless --plot-only)
    if not args.plot_only:
        print("=" * 70)
        print("COMPUTING HISTOGRAMS")
        print("=" * 70)
        compute_all_histograms(
            toml=args.toml,
            data=args.data,
            trajdir=args.trajdir,
            outdir=args.outdir,
            nskip=args.nskip,
            dlambda=args.dlambda,
            nbins=args.nbins,
            lmin=args.lmin,
            lmax=args.lmax,
            xcol=args.xcol,
            dt=args.dt,
            subcycles=args.subcycles,
            normalize=args.normalize,
            lm1=True,
        )

    # Generate plots (if --plot or --plot-only)
    if args.plot or args.plot_only:
        print("\n" + "=" * 70)
        print("GENERATING PLOTS")
        print("=" * 70)

        print(f"Loading data from {args.outdir}...")
        data = load_histogram_data(args.outdir)

        # Build save paths
        if args.save_plots:
            wham_save = os.path.join(args.save_plots, "wham_histogram.png")
            ens_save = os.path.join(args.save_plots, "per_ensemble_histogram.png")
            os.makedirs(args.save_plots, exist_ok=True)
        else:
            wham_save = None
            ens_save = None

        # Plot WHAM results
        print("\nPlotting WHAM results...")
        plot_wham_results(
            data,
            intf_list,
            args.title or "WHAM Results",
            save_path=wham_save,
            show=not args.no_show,
            log_scale=args.log
        )

        # Plot ensemble results
        if "ens_x" in data:
            print("Plotting per-ensemble results...")
            plot_ensemble_histograms(
                data,
                intf_list,
                args.title or "Per-Ensemble Results",
                save_path=ens_save,
                show=not args.no_show,
                log_scale=args.log
            )
        else:
            print("No per-ensemble data found, skipping ensemble plots")

        print("\nPlots complete.")


if __name__ == "__main__":
    main()
