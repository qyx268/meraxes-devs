---
title: 'Meraxes: a semi-analytic model of galaxy formation coupled to reionization and the 21-cm signal'
tags:
  - C
  - MPI
  - astronomy
  - cosmology
  - galaxy formation
  - reionization
  - 21-cm
authors:
  # TODO: confirm the author list, its order, ORCIDs and current affiliations with every
  # contributor; JOSS requires each co-author's agreement. This draft lists the code
  # contributors and the project lead. Ordered by lines of model code added over the
  # whole git history; contributors with none, then the project lead, come last.
  - name: Simon J. Mutch
    affiliation: 2  # TODO: current affiliation
  - name: Yuxiang Qin
    # orcid: TODO (add for every author)
    corresponding: true
    affiliation: 1
  - name: Fatma Shaban
    affiliation: 3  # TODO
  - name: Emanuele M. Ventura
    affiliation: 2  # TODO
  - name: Bradley Greig
    affiliation: 2  # TODO
  - name: Gregory B. Poole
    affiliation: 3  # TODO
  - name: Shihai Zhu
    affiliation: 1
  - name: Paul M. Geil
    affiliation: 2  # TODO
  - name: Yisheng Qiu
    affiliation: 2  # TODO
  - name: Sreedhar Balu
    affiliation: 3  # TODO
  - name: James E. Davies
    affiliation: 2  # TODO
  - name: Han-Seek Kim  # TODO: listed as "Hansik Kim" in the README; confirm
    affiliation: 2  # TODO
  - name: J. Stuart B. Wyithe
    affiliation: 1
affiliations:
  - name: Research School of Astronomy and Astrophysics, Australian National University, Australia  # TODO: check
    index: 1
  - name: School of Physics, The University of Melbourne, Australia  # TODO: check per author
    index: 2
  - name: TODO
    index: 3
date: 23 October 2026
bibliography: paper.bib
---

# Summary

Meraxes is a semi-analytic model of galaxy formation built for the first billion years
of cosmic history. It follows the baryons in dark matter haloes along merger trees from
cosmological N-body simulations: gas infall and cooling, star formation, supernova and
black-hole feedback, galaxy mergers and chemical enrichment. At every snapshot the
galaxies it has formed become the sources of a semi-numerical calculation of the
intergalactic medium (IGM) on a three-dimensional grid covering the same volume. That
calculation follows the ionization and thermal state of hydrogen through cosmic dawn and
the epoch of reionization (EoR) and predicts the 21-cm signal. Its results feed back on the
galaxies as photoheating in ionized regions suppresses gas accretion onto
low-mass haloes. Galaxies and the IGM therefore evolve together and self-consistently.

Meraxes is written in C, parallelized with MPI, and released under the GPL-3.0 license.
Two companion packages complete it: DRAGONS [@dragons], a Python package for reading and
analysing its output, and sector [@sector], a library that computes galaxy spectral
energy distributions and magnitudes.

# Statement of need

Galaxies are now observed beyond redshift 10, and radio interferometers aim to measure
the 21-cm signal from the same epoch. Interpreting both requires models that connect
individual galaxies to the ionization and thermal structure of the IGM. These models
must cover volumes large enough to sample reionization, hundreds of comoving
megaparsecs across, and run fast enough to be repeated many times when fitting data.
Radiation-hydrodynamic simulations resolve the relevant physics but are too expensive
for either.

Meraxes fills this gap. Galaxy histories come from merger trees, so star formation,
feedback and enrichment reflect how each galaxy assembled. The IGM is computed with
fast excursion-set methods on a grid that is spatially matched to the galaxies. The
code is aimed at researchers in high-redshift galaxy formation and 21-cm cosmology who
need galaxy populations and IGM observables from one consistent model. A single run
produces, for example, ultraviolet, [O III] and X-ray luminosity functions alongside
21-cm power spectra and light-cones (\autoref{fig:uvlf} and \autoref{fig:lightcones}).

![Rest-frame ultraviolet (1600 Å) luminosity functions of Meraxes galaxies (black lines)
from $z = 20$ to $5$, compared with measurements from Hubble, JWST and ground-based
surveys (symbols) [@finkelstein2015; @oesch2016; @livermore2017; @atek2018;
@ishigaki2018; @bouwens2021; @finkelstein2022; @kauffmann2022; @naidu2022; @donnan2023;
@harikane2023; @leethochawalit2023; @perezgonzalez2023].\label{fig:uvlf}](uvlf.jpg){ width=100% }

<!-- TODO: name the simulation and run behind the luminosity functions. -->

![Light-cones of the 21-cm brightness temperature, $\delta T_{\rm b}$, from $z \approx 30$
to $5$ in the 210 $h^{-1}$ Mpc simulations of @balu2023. L210 uses the merger trees of the
N-body simulation itself; L210_AUG extends them to haloes at the atomic-cooling limit with
the Monte Carlo algorithm darkforest [@qiu2021]; L210_AUG_lowX and L210_AUG_highX give
galaxies 100 times less and 100 times more soft X-ray luminosity per unit star formation
rate than L210_AUG. The gas is seen in absorption (red to yellow) until X-ray heating turns
the signal to emission (blue), and the signal vanishes inside ionized regions (black). With
weak X-ray heating (lowX), the gas never goes into emission.\label{fig:lightcones}](lightcones.jpg){ width=100% }

# State of the field

Semi-numerical 21-cm codes such as 21cmFAST [@mesinger2011; @murray2020] compute the IGM very efficiently. They describe the sources with analytic relations between halo mass and star formation rather than as individual galaxies with
assembly histories.

Established semi-analytic models, including L-Galaxies [@henriques2015], SAGE
[@croton2016] and Shark [@lagos2018], model galaxies in detail but treat reionization as
a uniform, prescribed background. 

Meraxes was designed around this coupling from its first release [@mutch2016]. It has
since been extended to follow the IGM through cosmic dawn, with X-ray heating,
Lyman-$\alpha$ coupling and the spin temperature [@balu2023]. It also models Population
III star formation in minihaloes, with Lyman–Werner feedback and metal enrichment
[@ventura2024; @ventura2025]. A single run therefore predicts the galaxy population,
reionization and the full 21-cm signal.

# Software design

Meraxes advances every galaxy in the volume one snapshot at a time. Galaxy physics and
the IGM calculation use two different parallel decompositions:

- **Galaxies.** Haloes are grouped into *forests*: sets of merger trees linked by mergers
  or by sharing a friends-of-friends group. Each forest lives entirely on one MPI rank,
  so the galaxy physics needs no communication. Forests are assigned to ranks using
  per-forest statistics computed in advance, which balance the work and bound the memory
  each rank needs.
- **IGM grids.** The radiation grids are divided into slabs for parallel Fourier
  transforms with FFTW. At each snapshot, galaxies deposit their star
  formation and black-hole accretion onto the grids. The excursion-set filtering of
  21cmFAST [@mesinger2011; @murray2020] is then applied in Fourier space, and the resulting
  ionization and feedback fields are returned to the galaxies on their own ranks.

**Inputs and outputs.** Meraxes reads merger trees from gbpTrees [@poole2017] and from
VELOCIraptor and TreeFrog [@elahi2019a; @elahi2019b]. This includes augmented trees,
such as those extended to lower-mass haloes with the Monte Carlo algorithm darkforest
[@qiu2021] or super-resolved with a generative model [@brennan2025]. It also reads
density and velocity grids from the N-body simulation. Output is written with parallel
HDF5, in one file per rank plus a master file that links them. It includes galaxy
catalogues, IGM grids (ionization, gas temperature, spin temperature and 21-cm
brightness temperature), 21-cm power spectra and light-cones, and mass and luminosity
functions.

**Configuration.** Run parameters come from text files layered over a defaults file and
a per-simulation file. Optional physics that changes the galaxy data structure is
selected at compile time with CMake, so runs that don't use it pay no memory for it:

- `USE_MINI_HALOS`: molecular cooling in minihaloes, Population III star formation,
  Lyman–Werner feedback and metal enrichment of the IGM [@ventura2024; @ventura2025].
- `CALC_MAGS`: galaxy magnitudes computed during the run with sector, including dust
  attenuation [@qiu2019], in rest-frame bands or in HST and JWST filters (`USE_HST`,
  `USE_JWST`).
- `USE_STOCHASTICITY`: galaxy-to-galaxy scatter in the escape fraction of ionizing
  photons and in X-ray luminosity, and recalibration of the source properties.

**Parameter inference.** Meraxes can also be built as a shared library. In this mode it
keeps haloes and grids in memory between runs, so an external driver can re-run the
model with new parameters without reading the merger trees again. Mhysa, a companion
Python package, uses this to connect Meraxes to Python samplers and optimisers for
Bayesian parameter inference [@mutch2024]; it is not yet publicly released.

**Testing and documentation.** Unit tests built with Criterion check individual physics
routines. A regression suite runs six physics configurations on a public 10 $h^{-1}$ Mpc
test set [@testdata]. It compares summary statistics against stored references:
distribution functions, quantiles of galaxy properties, grid statistics and 21-cm power
spectra. GitHub Actions builds each compile-time configuration with the address and
undefined-behaviour sanitizers, runs both suites and checks the code formatting. The
Meraxes Guide on Read the Docs documents building and running the code, every parameter
and output field, and the equations it implements. It includes a worked example that runs
Meraxes on the public test set and reads its output.

# Research impact statement

Meraxes is the galaxy-formation engine of the Dark-ages Reionization and Galaxy
Formation Simulation (DRAGONS) programme [@poole2016; @mutch2016]. It has been used to
study:

- the suppressed growth of haloes during reionization [@qin2017a]
- the contribution of quasars to reionization [@qin2017b]
- the ionizing-photon budget of galaxies [@duffy2017]
- the thermal memory of reionization [@davies2019]
- infrared excess and dust [@qiu2019]
- high-redshift black holes [@marshall2020]
- the thermal and reionization history of a large volume [@balu2023]
- Population III star formation and its imprint on the 21-cm signal [@ventura2024; @ventura2025]
- the evolution of the ionizing escape fraction, inferred jointly from galaxy and
  reionization observations [@mutch2024]
- analogues of the galaxies JWST has found at $z \gtrsim 12$ [@qin2023]
- the robustness of reionization histories inferred from the Lyman-$\alpha$ forest
  [@qin2026]
- the impact of galaxy-to-galaxy stochasticity on reionization and the 21-cm signal
  [@zhu2026]

<!-- TODO: give the number of papers in the public ADS library of Meraxes papers
(https://ui.adsabs.harvard.edu/public-libraries/CWUcYnt3TsmG6BuOKjR0Fw) and note any
users outside the core team. -->

The code has been developed openly since 2013 by contributors at several Australian
institutions, and its history is preserved in the public repository.

# AI usage disclosure

<!-- TODO: every author should check this statement and add any other use (for example
in the documentation or earlier code). -->

Generative AI was used in preparing this submission. Claude Code (Anthropic, using the
Claude Opus 5.5 model) was used in 2026 to:

- implement the regression-test suite, the continuous-integration workflows and the
  formatting tools
- fix and review code changes in recent pull requests
- update the documentation
- draft this paper

GitHub Copilot was also used to review pull requests. However, the authors reviewed all
AI-generated code and text. Code changes were built and tested on the CI and on the
authors' own computing systems, and every bibliographic entry was checked against its
DOI record. All scientific and design decisions were made by the authors.

# Acknowledgements

The development of Meraxes has received continuous support from the Australian Research
Council (ARC), including J.S.B.W.'s Australian Laureate Fellowship (FL110100072), the ARC
Centre of Excellence for All Sky Astrophysics in 3 Dimensions (ASTRO 3D; project number
CE170100013), and Y.Q.'s Discovery Early Career Researcher Award (DECRA; DE240101129).

<!-- TODO: consider acknowledging computing (e.g. Gadi at NCI Australia and OzSTAR at
Swinburne). -->

# References
