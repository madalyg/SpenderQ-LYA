"""SpenderQ continuum reconstruction of one epoch on SpenderQ's observed wavelength grid."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch

from spenderq import util as U
from spenderq.desi_qso import DESI
from spenderq.reconstruction.spectrum_io import load_sdss_like_spectrum


@dataclass
class EpochReconstruction:
    """One spectrum and its SpenderQ continuum, all on SpenderQ's observed grid."""

    path: Path
    z_qso: float
    norm: float                  # flux scale divided out before reconstruction
    wave: np.ndarray             # observed wavelength [Å]
    flux: np.ndarray
    ivar: np.ndarray
    continuum: np.ndarray
    recon_norm: np.ndarray       # flux / continuum
    recon_norm_err: np.ndarray
    weight_clipped: np.ndarray   # ivar after SpenderQ zeroed Lyα-forest pixels


def normalize_spectrum(wave_obs, flux, ivar, z_qso):
    """Normalize near a clean continuum window."""
    wave_rest = wave_obs / (1.0 + z_qso)
    # Anchor bands in mostly line-free continuum (1450 A main, 1700 A backup)
    windows = [(1445.0, 1455.0), (1700.0, 1705.0)]

    norm = np.nan
    for lo, hi in windows:
        use = (wave_rest > lo) & (wave_rest < hi) & np.isfinite(flux) & (ivar > 0)
        if np.count_nonzero(use) > 5:
            norm = np.nanmedian(flux[use])
            break

    if not np.isfinite(norm) or norm == 0:
        use = np.isfinite(flux) & (ivar > 0)
        norm = np.nanmedian(flux[use])

    # Dividing flux by norm scales ivar by 1/(norm)^2 so S/N is unchanged.
    return flux / norm, ivar * norm**2, norm


def resample_to_spenderq_grid(wave_obs, flux, ivar):
    """Interpolate a single spectrum onto SpenderQ's DESI wavelength grid."""
    # SpenderQ was trained on fixed observed wavelength bins (DESI). SDSS must land on the same grid
    wave_grid = np.asarray(DESI._wave_obs, dtype=float)
    good = np.isfinite(wave_obs) & np.isfinite(flux) & np.isfinite(ivar)
    good &= ivar > 0 # Check values are finite and that ivar > 0

    order = np.argsort(wave_obs[good])
    wave_in = wave_obs[good][order]
    flux_in = flux[good][order]
    ivar_in = ivar[good][order]

    flux_grid = np.interp(wave_grid, wave_in, flux_in, left=np.nan, right=np.nan)
    ivar_grid = np.interp(wave_grid, wave_in, ivar_in, left=0.0, right=0.0)

    # Outside SDSS coverage, flux/weight = 0 so SpenderQ ignores those pixels.
    bad = ~np.isfinite(flux_grid)
    flux_grid[bad] = 0.0
    ivar_grid[bad] = 0.0

    return wave_grid, flux_grid.astype(np.float32), ivar_grid.astype(np.float32)


def rescale_continuum_to_flux(wave_obs, flux, continuum, ivar, z_qso):
    """Rescale continuum so its median matches the observed flux at the normalization window.

    Removes the small shift between SpenderQ's normalization and the SDSS normalization,
    setting the continuum-normalized ratio to ~1 at the normalization window
    (1445-1455 Å restframe).
    """
    wave_rest = wave_obs / (1.0 + z_qso)
    windows = [(1445.0, 1455.0), (1700.0, 1705.0), (1275.0, 1290.0)]
    scale = np.nan
    for lo, hi in windows:
        m = (wave_rest > lo) & (wave_rest < hi) & np.isfinite(flux) & np.isfinite(continuum) & (ivar > 0) & (continuum > 0)
        if m.sum() > 5:
            scale = np.nanmedian(flux[m]) / np.nanmedian(continuum[m])
            break
    if not np.isfinite(scale) or scale == 0:
        scale = 1.0
    return continuum * scale, scale


def rebin_reconstruction_to_observed(wave_rest_model, recon_rest, wave_obs, z_qso):
    # Model output is in restframe. Rebin to this spectrum's observed wavelengths before comparing to data
    wave_model_obs = wave_rest_model * (1.0 + z_qso)
    try:
        return U.trapz_rebin(wave_model_obs, recon_rest, xnew=wave_obs)
    except ValueError:
        return np.interp(wave_obs, wave_model_obs, recon_rest, left=np.nan, right=np.nan)


def run_spenderq(model, wave_obs, flux, ivar, z_qso):
    """Return the observed-frame continuum and the ivar left after Lyα clipping."""
    spec = torch.from_numpy(flux[None, :])
    weight = torch.from_numpy(ivar[None, :])
    z = torch.tensor([z_qso], dtype=torch.float32)

    weight_after_lya_clip = weight.clone()
    # eval() iteratively zeroes weights on Lya forest pixels so the recon fits the continuum
    _, recon = model.eval(spec, weight_after_lya_clip, z)

    continuum = rebin_reconstruction_to_observed(
        model.wave_recon(), np.asarray(recon[0]), wave_obs, z_qso
    )
    continuum, _ = rescale_continuum_to_flux(wave_obs, flux, continuum, ivar, z_qso)
    return continuum, np.asarray(weight_after_lya_clip[0])


def reconstruct_epoch(path, z_qso, model):
    """Load one FITS spectrum, reconstruct its continuum with model (a SpenderQ instance),
    and divide the flux by that continuum."""
    wave_obs, flux, ivar = load_sdss_like_spectrum(path, z_qso)
    flux, ivar, norm = normalize_spectrum(wave_obs, flux, ivar, z_qso)
    wave, flux, ivar = resample_to_spenderq_grid(wave_obs, flux, ivar)
    continuum, weight_clipped = run_spenderq(model, wave, flux, ivar, z_qso)

    good = np.isfinite(continuum) & (continuum != 0) & (ivar > 0)
    recon_norm = np.full_like(flux, np.nan, dtype=float)
    recon_norm_err = np.full_like(flux, np.nan, dtype=float)
    recon_norm[good] = flux[good] / continuum[good]
    recon_norm_err[good] = 1.0 / np.sqrt(ivar[good]) / np.abs(continuum[good])

    return EpochReconstruction(
        path=Path(path),
        z_qso=z_qso,
        norm=norm,
        wave=wave,
        flux=flux,
        ivar=ivar,
        continuum=continuum,
        recon_norm=recon_norm,
        recon_norm_err=recon_norm_err,
        weight_clipped=weight_clipped,
    )
