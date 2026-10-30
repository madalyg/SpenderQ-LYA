"""Success criteria for SpenderQ continuum reconstructions.

Each check returns one row per continuum_regions band with the median and 16–84%
spread, so per-quasar runs and batch summaries report the same numbers.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from spenderq.continuum_regions import REGIONS

RATIO_LO, RATIO_HI = 0.0, 5.0
MIN_PIX = 5


def safe_ratio(num, den):
    """num / den, NaN where either is non-finite or den is zero."""
    return np.divide(
        num, den,
        out=np.full_like(num, np.nan, dtype=float),
        where=np.isfinite(num) & np.isfinite(den) & (den != 0),
    )


def region_ratio_stats(wave_obs, ratio, z):
    """Median, 16th and 84th percentile of a ratio spectrum in each band of REGIONS."""
    wave_rest = wave_obs / (1.0 + z)
    good = np.isfinite(ratio) & (ratio > RATIO_LO) & (ratio < RATIO_HI)
    rows = []
    for region in REGIONS:
        mask = region.mask_fn(wave_rest, wave_obs, good)
        if mask.sum() < MIN_PIX:
            continue
        r = ratio[mask]
        rows.append({
            "region_id": region.region_id,
            "region": region.label,
            "n_pix": int(mask.sum()),
            "median": float(np.median(r)),
            "p16": float(np.percentile(r, 16)),
            "p84": float(np.percentile(r, 84)),
            "min": float(r.min()),
            "max": float(r.max()),
        })
    return rows


def epoch_continuum_stability(wave_obs, cont_a, cont_b, z):
    """Ratio of two epochs' independent reconstructions; a median near 1 means both
    epochs recovered the same underlying continuum."""
    return region_ratio_stats(wave_obs, safe_ratio(cont_a, cont_b), z)


def _region_label(label):
    """Break a band name before its wavelength range so the cell shows both lines."""
    if " (" in label:
        name, rest = label.split(" (", 1)
        return f"{name}\n({rest}"
    return label


def save_stats_pdf(rows, title, outpath):
    """Write one success-criteria table to a PDF in the quasar analysis folder."""
    outpath = Path(outpath)
    columns = ["Region", "n", "median", "16%", "84%"]
    if rows:
        cells = [
            [
                _region_label(row["region"]),
                f"{row['n_pix']}",
                f"{row['median']:.3f}",
                f"{row['p16']:.3f}",
                f"{row['p84']:.3f}",
            ]
            for row in rows
        ]
    else:
        cells = [["no band with\nenough pixels", "—", "—", "—", "—"]]

    n_rows = len(cells) + 1
    fig_h = 0.52 * n_rows + 0.7
    fig = plt.figure(figsize=(8.5, max(3.2, fig_h)))
    ax = fig.add_axes([0.05, 0.16, 0.90, 0.74])
    ax.axis("off")
    fig.suptitle(title, fontsize=11, x=0.05, ha="left")
    table = ax.table(
        cellText=cells,
        colLabels=columns,
        loc="upper center",
        cellLoc="center",
        colWidths=[0.46, 0.12, 0.14, 0.14, 0.14],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    # Matplotlib sizes rows for one line and ignores breaks inside a cell.
    row_height = table._approx_text_height() * 2.6
    for (r, c), cell in table.get_celld().items():
        cell.set_height(row_height)
        text = cell.get_text()
        if c == 0:
            text.set_ha("left")
        if r == 0:
            cell.set_facecolor("#e6e6e6")
            text.set_weight("bold")
    fig.text(
        0.05,
        0.05,
        "n is the number of wavelength pixels in the band.\n"
        "16% and 84% are the 16th and 84th percentiles of the continuum ratio.",
        fontsize=9,
        ha="left",
        va="bottom",
    )
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, format="pdf", bbox_inches="tight")
    plt.close(fig)
    return outpath
