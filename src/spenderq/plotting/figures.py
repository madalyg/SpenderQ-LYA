"""Figure style, spectral-line markers, and saving for every SpenderQ-LYA plot."""

import matplotlib as mpl
import matplotlib.pyplot as plt

from spenderq.spectral_lines import ION_LINES, LYA

mpl.rcParams.update({"axes.labelsize": 14, "axes.titlesize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 11, "figure.titlesize": 14})


def mark_lines(ax, z_qso=0.0, labeled=tuple(ION_LINES), label_y=None):
    """Draw Lyα and the ION_LINES doublets, labeling the ions in labeled.

    Pass z_qso=0 for a rest-frame axis. Labels sit at label_y, or the top of the axis.
    """
    scale = 1.0 + z_qso
    ax.axvline(LYA * scale, color="tab:purple", ls="--", lw=0.9, label="Ly-alpha")
    for waves in ION_LINES.values():
        for wave in waves:
            ax.axvline(wave * scale, color="tab:red", ls=":", lw=0.9)
    y = ax.get_ylim()[1] if label_y is None else label_y
    for ion in labeled:
        ax.text(ION_LINES[ion][0] * scale, y, ion, color="tab:red", fontsize=14, va="top")


def finish_figure(fig, outpath, show=False):
    """Save fig at 200 dpi, then show it or close it."""
    fig.savefig(outpath, dpi=200)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return outpath
