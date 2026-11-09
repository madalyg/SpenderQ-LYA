"""Every two-epoch output for one quasar pair."""

from spenderq.multi_epoch.plots import (
    plot_continuum_ratio_panels,
    plot_recon_norm_difference,
    plot_recon_norm_ratio,
)
from spenderq.multi_epoch.ratios import (
    continuum_ratios,
    recon_norm_differences,
    recon_norm_ratios,
    save_continuum_ratio_txt,
    save_recon_norm_ratio_txt,
)
from spenderq.multi_epoch.success_criteria import epoch_continuum_stability, save_stats_pdf


def compare_epochs(a, b, analysis_dir, ratios_dir, show=False, rest_frame=False):
    """Write continuum-stability and (flux/recon) ratio outputs for epochs a and b.

    Plots and *_ratio_norm.txt go to analysis_dir; returns the recon_norm_ratios txt paths
    written to ratios_dir (rest-frame wavelengths if rest_frame, otherwise observed).
    """
    cont_ab, cont_ba = continuum_ratios(a, b)
    for path in save_continuum_ratio_txt(a, b, cont_ab, cont_ba, analysis_dir):
        print(f"Saved {path}")

    stability_pdf = analysis_dir / f"{a.path.stem}_vs_{b.path.stem}_epoch_continuum_stability.pdf"
    save_stats_pdf(
        epoch_continuum_stability(a.wave, a.continuum, b.continuum, a.z_qso),
        f"Epoch continuum stability ({a.path.stem} / {b.path.stem})",
        stability_pdf,
    )
    print(f"Saved {stability_pdf}")
    print(f"Saved {plot_continuum_ratio_panels(a, b, cont_ab, cont_ba, analysis_dir, show)}")

    ratio_ab, ratio_ba = recon_norm_ratios(a, b)
    print(f"Saved {plot_recon_norm_ratio(a, b, ratio_ab, ratio_ba, analysis_dir, show)}")
    diff_ab, diff_ba = recon_norm_differences(a, b)
    print(f"Saved {plot_recon_norm_difference(a, b, diff_ab, diff_ba, analysis_dir, show)}")

    ratio_txts = save_recon_norm_ratio_txt(a, b, ratio_ab, ratio_ba, ratios_dir, rest_frame)
    for path in ratio_txts:
        print(f"Saved {path}")
    return ratio_txts
