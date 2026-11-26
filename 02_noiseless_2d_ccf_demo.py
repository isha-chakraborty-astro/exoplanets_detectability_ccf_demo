"""
02_noiseless_2d_ccf_demo.py
---------------------------
Author: Isha Chakraborty (University of Bologna & INAF Catania)

Runs an end-to-end noiseless 2D Cross-Correlation Function (CCF) simulation:
1. Generates a high-resolution planetary transmission template
2. Injects Doppler-shifted planetary signals into a blaze-modulated stellar order
3. Removes the stellar continuum and blaze in a single spline pass
4. Computes and plots the 2D phase-resolved CCF matrix and rest-frame co-added CCF
"""

import numpy as np
import matplotlib.pyplot as plt
from importlib import import_module
from src.ccf_pipeline import (
    C_KMS,
    simulate_blaze_echelle_order,
    single_pass_spline_normalization,
    compute_noiseless_2d_ccf,
)

prt_mod = import_module("01_prt_model_spectrum")


def run_demo(kp_kms: float = 150.0, vsys_kms: float = -15.0) -> None:
    wl_um, norm_trans = prt_mod.generate_prt_or_synthetic_spectrum()
    template_abs = 1.0 - norm_trans

    phases = np.linspace(-0.045, 0.045, 25)
    planet_rvs = kp_kms * np.sin(2.0 * np.pi * phases) + vsys_kms
    blaze_order = simulate_blaze_echelle_order(wl_um)

    n_phases = phases.size
    normalized_series = np.zeros((n_phases, wl_um.size))

    for p_idx, rv_p in enumerate(planet_rvs):
        shifted_abs = np.interp(wl_um / (1.0 + rv_p / C_KMS), wl_um, template_abs, left=0.0, right=0.0)
        raw_order_flux = blaze_order * (1.0 - shifted_abs)
        norm_flux, _ = single_pass_spline_normalization(wl_um, raw_order_flux)
        normalized_series[p_idx] = norm_flux

    rv_grid = np.linspace(-250.0, 250.0, 251)
    ccf_2d = compute_noiseless_2d_ccf(wl_um, normalized_series, wl_um, template_abs, rv_grid)

    coadded_ccf = np.zeros_like(rv_grid)
    for p_idx, phi in enumerate(phases):
        v_orb = kp_kms * np.sin(2.0 * np.pi * phi)
        coadded_ccf += np.interp(rv_grid + v_orb, rv_grid, ccf_2d[p_idx], left=0.0, right=0.0)
    coadded_ccf /= coadded_ccf.max()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8), dpi=150)
    im = ax1.imshow(
        ccf_2d,
        aspect="auto",
        origin="lower",
        extent=[rv_grid[0], rv_grid[-1], phases[0], phases[-1]],
        cmap="viridis",
    )
    ax1.axvline(vsys_kms, color="#38bdf8", ls="--", lw=1.2, label=f"V_sys = {vsys_kms} km/s")
    ax1.set_xlabel("Radial Velocity [km/s]")
    ax1.set_ylabel("Orbital Phase [φ]")
    ax1.set_title("Noiseless 2D Phase-Resolved CCF Matrix")
    ax1.legend(loc="upper left")
    fig.colorbar(im, ax=ax1, label="CCF Amplitude")

    ax2.plot(rv_grid, coadded_ccf, color="#0284c7", lw=2.0, label="Planet Rest-Frame Co-Added CCF")
    ax2.axvline(vsys_kms, color="#f59e0b", ls="--", lw=1.4, label=f"Injected V_sys ({vsys_kms} km/s)")
    ax2.set_xlabel("Rest-Frame Radial Velocity [km/s]")
    ax2.set_ylabel("Normalized CCF")
    ax2.set_title(f"Co-Added Rest-Frame Detection (K_p = {kp_kms} km/s)")
    ax2.legend()

    plt.tight_layout()
    plt.savefig("noiseless_2d_ccf_summary.png")
    print("Saved noiseless_2d_ccf_summary.png successfully.")


if __name__ == "__main__":
    run_demo()
