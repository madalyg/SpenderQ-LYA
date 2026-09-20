# SpenderQ-LYA

Batch pipeline for finding **extremely high-velocity outflow (EHVO)** absorption in SDSS quasar spectra, including troughs from outflows at v > 20% the speed of light that are buried deep in the Lyman-Alpha (Ly-α) forest.

Only two such forest EHVOs have been identified. The absorption is mixed with the Ly-α forest, so they cannot be identified from a single spectrum by eye. This repo uses [SpenderQ](https://github.com/changhoonhahn/SpenderQ) to reconstruct the unabsorbed continuum, then compares **two epochs of the same quasar**.

**First principles.** Outflow absorption varies between epochs due to the fast-moving nature of the gases in outflows. Lyα forest lines do not: they stay relatively constant over the timescales in which we observe quasars. Dividing two spectra of the same quasar therefore isolates the variable (outflow) absorption and cancels the forest lines. Each spectrum is then divided by its SpenderQ reconstruction, which normalizes it if the underlying continuum is the same in both epochs. Dividing those continuum-normalized spectra by each other amplifies the variable absorption between epochs.

It is built to run on large datasets: one CSV of quasars and two FITS filenames per object, for any SDSS data release that uses the usual `spec-PLATE-MJD-FIBER.fits` naming.

## Pipeline

1. **CSV** — one row per quasar, two observations.
2. **`init_directories.py`** — one folder per quasar; downloads the two SDSS spectra.
3. **`analyze_ehvo.py`** — SpenderQ reconstruction for each spectrum, analysis plots, and a direct comparison of the two epochs.

Optional:

- **`analyze_ehvo_remove_lya.py`** — same analysis after interpolating sharp Lyα forest spikes (`remove_lya.py`).
- **`summarize_continuum_ratios.py`** — median continuum ratios in rest-frame bands across a case directory.
- **`convert_fits.py`** — rest-frame text spectra → observed-frame FITS (one-off helper).

## Input CSV

Columns:

```
Quasar Name, Observation 1, Observation 2, Redshift
```

`Observation 1` / `Observation 2` are SDSS spectrum filenames, e.g. `spec-0948-52428-0370.fits`.

## Setup

```bash
pip install -e .
pip install matplotlib scipy
```

Trained SpenderQ weights (`qso.dr1.hiz*.pt`) are under `src/spenderq/dat/`. The original SpenderQ training and DESI/PICCA scripts in `bin/` and `doc/` are unchanged upstream code.

## Run

From the repo root:

```bash
python src/spenderq/init_directories.py path/to/list.csv path/to/case_dir
```

Set `CASE_DIR` and `CSV_PATH` at the top of `analyze_ehvo.py`, then:

```bash
python src/spenderq/analyze_ehvo.py
```

`init_directories` tries SDSS DR16 SAS (eBOSS, then SDSS-I/II). Already-downloaded FITS files are skipped.

## Outputs (per quasar)

Under `case_dir/<quasar>/`:

| Path | Contents |
|------|----------|
| `spec-*.fits` | Downloaded spectra |
| `spenderq_analysis/` | Reconstruction plots (rest + observed), continuum-ratio and flux/continuum-ratio panels, `*_ratio_norm.txt` |
| `../recon_norm_ratios/` | `(flux/cont)_A / (flux/cont)_B` text files, named by MJD |
| `../flux_over_recon/` | Per-epoch `flux / continuum` text files |

A case-level summary CSV of the recon-ratio files is written next to the case directory.

## Layout

```
src/spenderq/
  init_directories.py      # CSV → folders + SDSS downloads
  analyze_ehvo.py          # reconstructions + epoch comparison
  analyze_ehvo_remove_lya.py
  remove_lya.py            # median-filter Lyα spike cleaning
  summarize_continuum_ratios.py
  continuum_regions.py     # rest-frame band masks
  convert_fits.py
  spenderq.py              # SpenderQ model wrapper (upstream)
  lyalpha.py, desi_qso.py, util.py
```

## SpenderQ

Continuum reconstruction is [SpenderQ](https://github.com/changhoonhahn/SpenderQ) (Hahn et al.): a data-driven spectrum autoencoder that masks absorption and reconstructs the intrinsic quasar spectrum. This repository adds SDSS data ingestion, per-quasar directory layout, and detailed epoch-to-epoch EHVO analysis on top of that model.
