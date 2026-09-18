"""Read SpenderQ continuum ratio spectra (*_ratio_norm.txt from analyze_ehvo) and compute
median ratios in standard rest-frame bands for each exposure pair. Prints a batch summary
and writes continuum_ratio_summary.csv when run on a case or quasar directory."""

import argparse
import csv
from collections.abc import Iterator
from pathlib import Path

import numpy as np

from spenderq.continuum_regions import REGIONS

RATIO_LO, RATIO_HI = 0.0, 5.0


def load_redshifts(directory: Path) -> dict[str, float]:
    """Load quasar redshifts from *list_fits.csv under directory or its parents."""
    for root in [directory, *directory.parents]:
        for csv_path in sorted(root.glob("*list_fits.csv")):
            zs = {}
            with open(csv_path, newline="", encoding="utf-8") as fh:
                for row in csv.DictReader(fh):
                    name = row.get("Quasar Name", "").strip()
                    z = float(row["Redshift"])
                    zs[name] = z
                    short = name.split("+")[0].split(".")[0]
                    zs[short] = z
                    if name.startswith("J") and "+" in name:
                        zs[name.split("+")[0]] = z
            if zs:
                return zs
    return {}


def quasar_key(path: Path) -> str:
    """Quasar folder name for a ratio_norm.txt under .../<quasar>/spenderq_analysis/."""
    if path.parent.name == "spenderq_analysis":
        return path.parent.parent.name
    return path.parent.name


def unique_ratio_files(directory: Path) -> list[Path]:
    seen: set[tuple[str, str]] = set()
    files: list[Path] = []
    for path in sorted(directory.rglob("*_ratio_norm.txt")):
        stem = path.stem.replace("_ratio_norm", "")
        if "_vs_" not in stem:
            continue
        a, b = stem.split("_vs_", 1)
        key = tuple(sorted([a, b]))
        if key in seen:
            continue
        seen.add(key)
        files.append(path)
    return files


def resolve_z(quasar: str, redshifts: dict[str, float], override: float | None) -> float:
    if override is not None:
        return override
    if quasar in redshifts:
        return redshifts[quasar]
    for key, z in redshifts.items():
        if quasar in key or key in quasar:
            return z
    raise ValueError(f"No redshift for {quasar!r}; pass --redshift or add *list_fits.csv")


def summarize_pair(path: Path, z: float) -> list[dict]:
    data = np.loadtxt(path)
    wave_obs, ratio = data[:, 0], data[:, 1]
    wave_rest = wave_obs / (1.0 + z)
    good = np.isfinite(ratio) & (ratio > RATIO_LO) & (ratio < RATIO_HI)

    rows = []
    for region in REGIONS:
        mask = region.mask_fn(wave_rest, wave_obs, good)
        if mask.sum() < 5:
            continue
        r = ratio[mask]
        rows.append(
            {
                "file": path.name,
                "region_id": region.region_id,
                "region": region.label,
                "n_pix": int(mask.sum()),
                "median": float(np.median(r)),
                "p16": float(np.percentile(r, 16)),
                "p84": float(np.percentile(r, 84)),
                "min": float(r.min()),
                "max": float(r.max()),
            }
        )
    return rows


def iter_region_groups(pair_rows: list[dict]) -> Iterator[tuple[str, str, list[dict]]]:
    for region in REGIONS:
        sub = [r for r in pair_rows if r["region_id"] == region.region_id]
        if sub:
            yield region.region_id, region.label, sub


def aggregate_pair_medians(sub: list[dict]) -> dict[str, float | int]:
    meds = np.fromiter((r["median"] for r in sub), dtype=float)
    return {
        "n_pairs": len(meds),
        "med_of_medians": float(np.median(meds)),
        "pair_min": float(np.min(meds)),
        "pair_max": float(np.max(meds)),
        "p16_pairs": float(np.percentile(meds, 16)),
        "p84_pairs": float(np.percentile(meds, 84)),
    }


def print_summary(pair_rows: list[dict]) -> None:
    print(f"{'Region':<38} {'N':>3} {'Median':>8} {'Pair range':>14} {'16–84%':>14}")
    print("-" * 82)
    for _, region_label, sub in iter_region_groups(pair_rows):
        agg = aggregate_pair_medians(sub)
        print(
            f"{region_label:<38} {agg['n_pairs']:3d} {agg['med_of_medians']:8.3f} "
            f"{agg['pair_min']:6.3f}–{agg['pair_max']:5.3f}   "
            f"{agg['p16_pairs']:6.3f}–{agg['p84_pairs']:5.3f}"
        )


def write_csv(outpath: Path, pair_rows: list[dict]) -> None:
    with open(outpath, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["quasar", "z", "file", "region_id", "region", "n_pix", "median", "p16", "p84", "min", "max"])
        for row in pair_rows:
            w.writerow(
                [
                    row["quasar"],
                    row["z"],
                    row["file"],
                    row["region_id"],
                    row["region"],
                    row["n_pix"],
                    f"{row['median']:.6f}",
                    f"{row['p16']:.6f}",
                    f"{row['p84']:.6f}",
                    f"{row['min']:.6f}",
                    f"{row['max']:.6f}",
                ]
            )
        w.writerow([])
        w.writerow(["region_id", "region", "n_pairs", "med_of_medians", "pair_min", "pair_max", "p16_pairs", "p84_pairs"])
        for region_id, region_label, sub in iter_region_groups(pair_rows):
            agg = aggregate_pair_medians(sub)
            w.writerow(
                [
                    region_id,
                    region_label,
                    agg["n_pairs"],
                    f"{agg['med_of_medians']:.6f}",
                    f"{agg['pair_min']:.6f}",
                    f"{agg['pair_max']:.6f}",
                    f"{agg['p16_pairs']:.6f}",
                    f"{agg['p84_pairs']:.6f}",
                ]
            )
    print(f"Wrote {outpath}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="Case dir or quasar spenderq_analysis dir")
    parser.add_argument("-z", "--redshift", type=float, default=None, help="Redshift if a single quasar dir")
    parser.add_argument("-o", "--output", type=Path, default=None, help="Output CSV path")
    args = parser.parse_args()

    directory = args.directory.resolve()
    if not directory.is_dir():
        raise SystemExit(f"Not a directory: {directory}")

    files = unique_ratio_files(directory)
    if not files:
        raise SystemExit(f"No *_ratio_norm.txt files under {directory}")

    redshifts = load_redshifts(directory)
    pair_rows: list[dict] = []
    for path in files:
        quasar = quasar_key(path)
        z = resolve_z(quasar, redshifts, args.redshift)
        for row in summarize_pair(path, z):
            pair_rows.append({"quasar": quasar, "z": z, **row})

    print(f"Found {len(files)} unique pair(s) under {directory}\n")
    print_summary(pair_rows)

    outpath = args.output or (directory / "continuum_ratio_summary.csv")
    write_csv(outpath, pair_rows)


if __name__ == "__main__":
    main()
