"""
01_prt_model_spectrum.py
------------------------
Author: Isha Chakraborty (University of Bologna & INAF Catania)

Demonstrates how to generate a high-resolution (lbl mode, R ~ 100,000)
exoplanet transmission spectrum in the 2.30-2.32 um K-band window using
petitRADTRANS, with an automatic analytical fallback so visitors can run
this script even without local opacity tables installed.
"""

import numpy as np
import matplotlib.pyplot as plt


def generate_prt_or_synthetic_spectrum(
    wl_min_um: float = 2.300,
    wl_max_um: float = 2.320,
    temp_iso_k: float = 1400.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (wavelength_um, normalized_transmission_depth) in [2.30, 2.32] um."""
    try:
        from petitRADTRANS import Radtrans
        from petitRADTRANS import nat_cst as nc

        atmosphere = Radtrans(
            line_species=["H2O_main_iso", "CO_main_iso"],
            rayleigh_species=["H2", "He"],
            continuum_opacities=["H2-H2", "H2-He"],
            wlen_bords_micron=[wl_min_um, wl_max_um],
            mode="lbl",
            lbl_opacity_sampling=1,
        )

        pressures = np.logspace(-6, 2, 100)
        atmosphere.setup_opa_structure(pressures)

        temperature = temp_iso_k * np.ones_like(pressures)
        mass_fractions = {
            "H2": 0.74 * np.ones_like(pressures),
            "He": 0.24 * np.ones_like(pressures),
            "H2O_main_iso": 1e-3 * np.ones_like(pressures),
            "CO_main_iso": 5e-4 * np.ones_like(pressures),
        }
        mmw = 2.33 * np.ones_like(pressures)
        r_pl = 1.2 * nc.r_jup_mean
        r_star = 1.0 * nc.r_sun
        gravity = 1e3  # cm/s^2
        p0 = 0.01

        atmosphere.calc_transm(
            temp=temperature,
            abunds=mass_fractions,
            gravity=gravity,
            mmw=mmw,
            R_pl=r_pl,
            P0_bar=p0,
        )
        wl_um = nc.c / atmosphere.freq / 1e-4
        transit_depth = (atmosphere.transm_rad / r_star) ** 2
        norm_transmission = 1.0 - (transit_depth - np.percentile(transit_depth, 10))
        return wl_um, norm_transmission

    except Exception:
        wl_um = np.linspace(wl_min_um, wl_max_um, 4000)
        rng = np.random.default_rng(42)
        line_centers = np.sort(rng.uniform(wl_min_um + 0.0005, wl_max_um - 0.0005, 45))
        line_weights = rng.uniform(0.4, 1.0, size=line_centers.size)
        sigma_um = 2.31 * (2.5 / 299792.458)

        absorption = np.zeros_like(wl_um)
        for center, weight in zip(line_centers, line_weights):
            absorption += 0.0018 * weight * np.exp(-0.5 * ((wl_um - center) / sigma_um) ** 2)

        return wl_um, 1.0 - absorption


if __name__ == "__main__":
    wl, trans = generate_prt_or_synthetic_spectrum()
    plt.figure(figsize=(10, 4), dpi=150)
    plt.plot(wl, trans, color="#0284c7", lw=1.0, label="Normalized Planetary Transmission T(λ)")
    plt.xlabel("Wavelength [µm]")
    plt.ylabel("Relative Transmission")
    plt.title("Synthetic High-Resolution Exoplanet Transmission Template (R = 100,000)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("model_spectrum_preview.png")
    print("Saved model_spectrum_preview.png successfully.")
