# High-Resolution Exoplanet Transmission Spectroscopy & Noiseless 2D CCF Pipeline

**Author:** Isha Chakraborty  
**Affiliations:** Laurea Magistrale (M.Sc.) Candidate in Astrophysics, [Alma Mater Studiorum – Università di Bologna](https://corsi.unibo.it/2cycle/Astrophysics) | Master's Thesis Researcher, [INAF – Osservatorio Astrofisico di Catania](https://www.oact.inaf.it/en/) (Supervisor: Dr. Gaetano Scandariato)

### [Launch the Live Interactive 2D CCF & Spectrum Simulator](https://isha-chakraborty-astro.github.io/exoplanets_detectability_ccf_demo/)

---

## Overview

This repository demonstrates the numerical and spectroscopic methodology behind my Master's thesis research at **INAF – Osservatorio Astrofisico di Catania**, evaluating exoplanet atmospheric detectability for the upcoming **[ANDES (ArmazoNes high Dispersion Echelle Spectrograph)](https://elt.eso.org/instrument/ANDES/)** instrument ($R \sim 100,000$) on ESO's **[Extremely Large Telescope (ELT)](https://elt.eso.org/)**.

> **Note on Data Policy:** To respect consortium boundaries and unpublished thesis simulation grids, this repository uses **purely synthetic, noiseless analytical inputs** to showcase the pipeline architecture, single-pass spline continuum normalization, and 2D Cross-Correlation Function (CCF) analysis.

---

## Repository Structure

* [`index.html`](https://isha-chakraborty-astro.github.io/exoplanets_detectability_ccf_demo/): Interactive browser simulator for exploring K-band planetary transmission templates ($\text{H}_2\text{O}$, $\text{CO}$, $\text{CH}_4$), single-pass echelle blaze normalization, and the 2D CCF velocity trail across transit phase $\phi$.
* `01_prt_model_spectrum.py`: Forward models a high-resolution ($R = 100,000$) planetary transmission spectrum in `petitRADTRANS` (`lbl` mode), with an analytical fallback for lightweight execution without local opacity tables.
* `src/ccf_pipeline.py`: Modular Python functions for relativistic Doppler shifting, parabolic echelle blaze simulation, single-pass `scipy.interpolate.UnivariateSpline` continuum normalization, and 2D CCF matrix computation.
* `02_noiseless_2d_ccf_demo.py`: End-to-end noiseless demonstration script that produces the 2D phase-resolved CCF matrix and co-added planetary rest-frame detection peak.

---

## Quickstart

```bash
git clone [https://github.com/isha-chakraborty-astro/exoplanets_detectability_ccf_demo.git](https://github.com/isha-chakraborty-astro/exoplanets_detectability_ccf_demo.git)
cd exoplanets_detectability_ccf_demo
pip install -r requirements.txt

# 1. Generate synthetic high-resolution transmission spectrum
python 01_prt_model_spectrum.py

# 2. Run end-to-end noiseless 2D CCF pipeline and export summary plots
python 02_noiseless_2d_ccf_demo.py
