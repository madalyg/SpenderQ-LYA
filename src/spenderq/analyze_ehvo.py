"""
This program creates ML-reconstructed unabsorbed continuua for SDSS quasar spectra
using a modified version of SpenderQ (located in SpenderQ.py). On start-up it asks
"Analyze multi-epoch spectra using ratios?":

    Y  Reconstruct both observations of each quasar, then write epoch-vs-epoch ratios,
       ratio and difference plots, and recon_norm_ratios/ text files (multi_epoch/).
    N  Reconstruct every observation independently and save its flux / reconstruction
       to recon_norm_individual/, without dividing by another epoch.

A CSV row with only one observation is still reconstructed either way.

It then asks "Convert output to rest frame (Y) or keep in observed frame (N)?". Y writes the
wavelength column of the recon_norm_ratios/ or recon_norm_individual/ txt files as
wavelength / (1 + z); N keeps observed wavelengths. Plots are unaffected.

@author: Madaly
"""

from pathlib import Path
import sys

import torch

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from spenderq.multi_epoch.compare import compare_epochs
from spenderq.multi_epoch.ratios import recon_norm_ratio_paths
from spenderq.plotting.epoch_plots import masked_forest_pixels, plot_flux_over_continuum, plot_reconstruction
from spenderq.reconstruction.continuum import reconstruct_epoch
from spenderq.reconstruction.outputs import save_absorption_csv, save_recon_norm_txt
from spenderq.reconstruction.spectrum_io import get_mjd, quasar_dir, read_quasar_csv
from spenderq.spenderq import SpenderQ

########################################################################################################################################

CASE_DIR = Path("dat/Dense_Lya")
#CASE_DIR = Path("dat/Variability/NV/")
CSV_PATH = CASE_DIR / "J2318_list_fits.csv"
#CSV_PATH = CASE_DIR / "NV_variability_list_fits.csv"
SHOW_PLOT = False # Set to True to show all plots

########################################################################################################################################

MODEL_NAME = "qso.dr1.hiz"
ANALYSIS_SUBDIR = "spenderq_analysis"  # per-quasar plots (and *_ratio_norm.txt in ratio mode)
RECON_NORM_RATIOS_DIR = CASE_DIR / "recon_norm_ratios"  # Y: case-level epoch-vs-epoch (flux/recon) ratios
RECON_NORM_INDIVIDUAL_DIR = CASE_DIR / "recon_norm_individual"  # N: (flux/recon) spectra, not divided by another epoch


def ask_use_ratios():
    """Prompt until the user answers Y or N; True means run the multi-epoch ratio analysis."""
    while True:
        answer = input("Analyze multi-epoch spectra using ratios? (Y/N): ").strip().upper()
        if answer in ("Y", "N"):
            return answer == "Y"
        print("Please enter Y or N.")


def ask_rest_frame():
    """Prompt until the user answers Y or N; True means save the txt outputs in the rest frame."""
    while True:
        answer = input("Convert output to rest frame (Y) or keep in observed frame (N)? ").strip().upper()
        if answer in ("Y", "N"):
            return answer == "Y"
        print("Please enter Y or N.")


def frame_name(rest_frame):
    return "rest" if rest_frame else "observed"


def mark_complete(marker, rest_frame):
    """Touch marker, recording the wavelength frame its txt outputs were written in."""
    marker.write_text(frame_name(rest_frame))


def is_up_to_date(marker, inputs, outputs, rest_frame=False):
    """True when marker is newer than every input, every output exists, and the outputs
    were written in the requested frame (an empty, older marker means observed).

    Delete the marker to force a reprocess.
    """
    return (
        marker.exists()
        and (marker.read_text().strip() or "observed") == frame_name(rest_frame)
        and all(p.exists() and marker.stat().st_mtime >= p.stat().st_mtime for p in inputs)
        and all(p.exists() for p in outputs)
    )


def existing_fits(qso_dir, observations):
    """Paths of the CSV observations that are present in qso_dir."""
    paths = []
    for filename in observations:
        path = qso_dir / filename
        if path.exists():
            paths.append(path)
        else:
            print(f"  [skip] FITS not found: {path}")
    return paths


def reconstruct_and_plot(fits_path, z_qso, model, analysis_dir):
    """Reconstruct one epoch and save the plots made in both modes."""
    epoch = reconstruct_epoch(fits_path, z_qso, model)
    rest_png = plot_reconstruction(epoch, analysis_dir, observed=False, show=SHOW_PLOT)
    observed_png = plot_reconstruction(epoch, analysis_dir, observed=True, show=SHOW_PLOT)
    ratio_png = plot_flux_over_continuum(epoch, analysis_dir, show=SHOW_PLOT)
    n_masked = int(masked_forest_pixels(epoch).sum())
    print(f"{fits_path.name}: normalization={epoch.norm:.4g}")
    print(f"Saved {rest_png} ({n_masked} Ly-alpha forest pixels masked)")
    print(f"Saved {observed_png}")
    print(f"Saved {ratio_png}")
    return epoch


def print_quasar_header(quasar, qso_dir):
    print(f"\n{'='*60}")
    print(f"Quasar: {quasar['name']}  z={quasar['redshift']}")
    print(f"Dir:    {qso_dir}")


def run_quasar_ratios(quasar, model, rest_frame=False):
    """Y: reconstruct both epochs and write the epoch-vs-epoch outputs.

    Returns the recon_norm_ratios txt paths (empty if fewer than two spectra are available).
    """
    z_qso = quasar["redshift"]
    qso_dir = CASE_DIR / quasar_dir(quasar["name"])
    analysis_dir = qso_dir / ANALYSIS_SUBDIR
    analysis_dir.mkdir(parents=True, exist_ok=True)

    requested = [qso_dir / filename for filename in quasar["observations"]]
    done = None
    if len(requested) == 2:
        mjd_a, mjd_b = get_mjd(requested[0]), get_mjd(requested[1])
        # One marker per pair so quasars with 3+ epochs run every requested pair.
        done = analysis_dir / f".complete_{mjd_a}_{mjd_b}"
        ratio_txts = recon_norm_ratio_paths(RECON_NORM_RATIOS_DIR, *requested)
        if is_up_to_date(done, requested, ratio_txts, rest_frame):
            print(f"\n{'='*60}\n  [skip] {qso_dir.name} {mjd_a} vs {mjd_b} — delete {done.name} to rerun\n")
            return ratio_txts

    print_quasar_header(quasar, qso_dir)
    epochs = []
    for fits_path in existing_fits(qso_dir, quasar["observations"]):
        epochs.append(reconstruct_and_plot(fits_path, z_qso, model, analysis_dir))

    if len(epochs) < 2:
        print(f"  [skip] ratios need 2 spectra, got {len(epochs)} — skipping cross-file outputs")
        return []

    ratio_txts = compare_epochs(epochs[0], epochs[1], analysis_dir, RECON_NORM_RATIOS_DIR, show=SHOW_PLOT, rest_frame=rest_frame)
    mark_complete(done, rest_frame)
    return ratio_txts


def run_quasar_individual(quasar, model, rest_frame=False):
    """N: reconstruct every listed epoch on its own and save its flux / recon.

    Returns the recon_norm_individual txt paths.
    """
    z_qso = quasar["redshift"]
    qso_dir = CASE_DIR / quasar_dir(quasar["name"])
    analysis_dir = qso_dir / ANALYSIS_SUBDIR
    analysis_dir.mkdir(parents=True, exist_ok=True)

    print_quasar_header(quasar, qso_dir)
    saved = []
    for fits_path in existing_fits(qso_dir, quasar["observations"]):
        txt = RECON_NORM_INDIVIDUAL_DIR / f"{fits_path.stem}_over_recon.txt"
        done = analysis_dir / f".complete_individual_{fits_path.stem}"
        if is_up_to_date(done, [fits_path], [txt], rest_frame):
            print(f"  [skip] {fits_path.name} — delete {done.name} to rerun")
        else:
            epoch = reconstruct_and_plot(fits_path, z_qso, model, analysis_dir)
            print(f"Saved {save_recon_norm_txt(epoch, txt, rest_frame)}")
            mark_complete(done, rest_frame)
        saved.append(txt)
    return saved


if __name__ == "__main__":
    torch.set_grad_enabled(False)  # We are not training the model, only running inference (aka running the model forward)

    use_ratios = ask_use_ratios()
    rest_frame = ask_rest_frame()
    if use_ratios:
        run_quasar, out_dir, csv_suffix = run_quasar_ratios, RECON_NORM_RATIOS_DIR, "recon_ratios_norm"
    else:
        run_quasar, out_dir, csv_suffix = run_quasar_individual, RECON_NORM_INDIVIDUAL_DIR, "recon_norm_individual"
    out_dir.mkdir(parents=True, exist_ok=True)

    quasars = read_quasar_csv(CSV_PATH)
    print(f"Found {len(quasars)} quasar(s) in {CSV_PATH}")
    model = SpenderQ(MODEL_NAME)

    summary_entries = []
    for quasar in quasars:
        if not quasar["observations"]:
            print(f"  [skip] {quasar['name']}: no observation filenames in {CSV_PATH.name}")
            continue
        for txt in run_quasar(quasar, model, rest_frame):
            summary_entries.append((txt, quasar["redshift"]))

    if summary_entries:
        # Column layout is the input format of the downstream absorption search (absorption.py).
        save_absorption_csv(CASE_DIR / f"{CASE_DIR.resolve().name}_{csv_suffix}.csv", summary_entries)
