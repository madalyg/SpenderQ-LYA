"""Rest-frame line wavelengths and plot windows shared by the individual-epoch and
multi-epoch analyses."""

LYA = 1215.67
LYB = 1026.00
ION_LINES = {
    "N V": (1238.82, 1242.80),
    "Si IV": (1393.755, 1402.770),
    "C IV": (1548.1950, 1550.770),
}
RECON_PLOT_XLIM_REST = (1000.0, 1700.0)


def lya_region_observed(z_qso):
    """Return observed-frame bounds for the Lyα forest region (Lyβ–Lyα rest)."""
    return LYB * (1.0 + z_qso), LYA * (1.0 + z_qso)
