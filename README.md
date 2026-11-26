# High-Resolution Exoplanet Transmission Spectroscopy & Noiseless 2D CCF Pipeline

**Author:** Isha Chakraborty  
**Affiliations:** Laurea Magistrale (M.Sc.) Candidate in Astrophysics, [Alma Mater Studiorum – Università di Bologna](https://corsi.unibo.it/2cycle/Astrophysics) | Master's Thesis Researcher, [INAF – Osservatorio Astrofisico di Catania](https://www.oact.inaf.it/en/) (Supervisor: Dr. Gaetano Scandariato)

## Overview & Public Methodology Scope

This repository showcases the numerical and spectroscopic methodology developed for my Master's thesis research at **INAF – Osservatorio Astrofisico di Catania**, which evaluates exoplanet atmospheric detectability for the upcoming **[ANDES (ArmazoNes high Dispersion Echelle Spectrograph)](https://elt.eso.org/instrument/ANDES/)** instrument ($R \sim 100,000$) on ESO's **[Extremely Large Telescope (ELT)](https://elt.eso.org/)**.

> **Note on Data Policy:** To respect consortium boundaries and unpublished thesis simulation grids (such as cloudy vs. clear 3D detectability grids), this repository uses **purely synthetic, noiseless analytical inputs** to demonstrate the software architecture, single-pass spline continuum normalization, and 2D Cross-Correlation Function (CCF) pipeline.

### Live Interactive Browser Simulator (GitHub Pages)
Once GitHub Pages is enabled on this repository (`Settings -> Pages -> Branch: main`), the included `index.html` file serves an interactive browser simulator at:
* `https://<your-github-username>.github.io/exoplanet-ccf-demo/`

## Core Files
* `index.html`: Self-contained interactive GitHub Pages web simulator (2D CCF velocity matrix, single-pass spline continuum normalization toggle, and K-band absorber templates).
* `01_prt_model_spectrum.py`: Forward models a high-resolution ($R = 100,000$) planetary transmission spectrum in `petitRADTRANS` (`lbl` mode), with an automatic analytical fallback so the code runs anywhere without local opacity tables.
* `src/ccf_pipeline.py`: Reusable Python library for Doppler shifting, parabolic echelle blaze simulation, single-pass `scipy.interpolate.UnivariateSpline` normalization, and 2D CCF matrix computation.
* `02_noiseless_2d_ccf_demo.py`: End-to-end demonstration script that generates the noiseless 2D CCF trail and co-added planetary rest-frame detection peak.
