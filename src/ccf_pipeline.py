"""
src/ccf_pipeline.py
-------------------
Author: Isha Chakraborty (University of Bologna & INAF Catania)

Reusable functions for noiseless high-resolution exoplanet spectroscopy simulations:
1. Doppler-shifting planetary templates across orbital phases
2. Injecting signals into blaze-modulated synthetic stellar echelle spectra
3. Single-pass cubic spline continuum + blaze normalization
4. Computing 2D phase-resolved Cross-Correlation Function (CCF) matrices
"""

import numpy as np
from scipy.interpolate import interp1d, UnivariateSpline

C_KMS = 299792.458  # Speed of light in km/s


def doppler_shift_wavelength(wl_um: np.ndarray, rv_kms: float) -> np.ndarray:
    """Shift wavelength grid by radial velocity rv_kms (relativistic first-order)."""
    return wl_um * (1.0 + rv_kms / C_KMS)


def simulate_blaze_echelle_order(wl_um: np.ndarray) -> np.ndarray:
    """Create a synthetic stellar continuum modulated by a parabolic echelle blaze profile."""
    wl_mid = 0.5 * (wl_um[0] + wl_um[-1])
    half_span = 0.5 * (wl_um[-1] - wl_um[0])
    x_norm = (wl_um - wl_mid) / half_span
    blaze_profile = 1.0 - 0.38 * (x_norm ** 2) - 0.04 * x_norm
    return blaze_profile


def single_pass_spline_normalization(
    wl_um: np.ndarray,
    raw_flux: np.ndarray,
    n_knots: int = 12,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Fit a single-pass cubic spline directly to a blaze-uncorrected echelle order
    to simultaneously remove the stellar continuum and instrumental blaze function.
    """
    knot_edges = np.linspace(wl_um[0], wl_um[-1], n_knots + 1)
    knot_x = np.zeros(n_knots)
    knot_y = np.zeros(n_knots)

    for k in range(n_knots):
        mask = (wl_um >= knot_edges[k]) & (wl_um <= knot_edges[k + 1])
        knot_x[k] = np.median(wl_um[mask])
        knot_y[k] = np.percentile(raw_flux[mask], 90)

    spline = UnivariateSpline(knot_x, knot_y, k=3, s=0.0)
    continuum_fit = spline(wl_um)
    normalized_flux = raw_flux / continuum_fit
    return normalized_flux, continuum_fit


def compute_noiseless_2d_ccf(
    wl_obs: np.ndarray,
    flux_time_series: np.ndarray,
    wl_template: np.ndarray,
    template_absorption: np.ndarray,
    rv_grid_kms: np.ndarray,
) -> np.ndarray:
    """Compute the 2D Cross-Correlation Function matrix CCF(phase, rv) in noiseless condition."""
    n_phases = flux_time_series.shape[0]
    n_rv = rv_grid_kms.size
    ccf_matrix = np.zeros((n_phases, n_rv), dtype=float)

    residuals = 1.0 - flux_time_series
    residuals -= np.mean(residuals, axis=1, keepdims=True)

    interp_template = interp1d(
        wl_template,
        template_absorption - np.mean(template_absorption),
        kind="linear",
        bounds_error=False,
        fill_value=0.0,
    )

    for r_idx, rv in enumerate(rv_grid_kms):
        shifted_template = interp_template(wl_obs / (1.0 + rv / C_KMS))
        norm_factor = np.sqrt(np.sum(shifted_template ** 2)) + 1e-12
        ccf_matrix[:, r_idx] = (residuals @ shifted_template) / norm_factor

    return ccf_matrix
