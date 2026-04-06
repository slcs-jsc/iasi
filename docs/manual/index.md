# IASI Code Collection

The IASI Code Collection provides command-line tools for processing and analyzing observations from Eumetsat's Infrared Atmospheric Sounding Interferometer (IASI).

It combines small focused utilities for spectral extraction, perturbation analysis, map generation, noise estimation, and retrieval-oriented preprocessing. The repository also includes a runnable example workflow and local copies of the external libraries needed for typical Linux builds.

## What the project provides

- extraction of radiance spectra from Level-1C granules
- computation of brightness temperature products and perturbation fields
- conversion of perturbation products into geolocated map tables
- spectral band averaging and noise diagnostics
- preparation of retrieval inputs and an MPI-enabled retrieval executable
- utility programs inherited from the JURASSIC radiative transfer model for time and date conversions

## Documentation map

- Use [Getting Started](getting-started.md) for local build and runtime setup.
- Use [Tools](tools.md) for a summary of the executables shipped in `src/`.
- Use [Example Workflow](workflow.md) for the end-to-end processing path implemented in `projects/example/run.sh`.
- Use [Links](links.md) for repository, citation, license, and related documentation references.

## Project notes

The source tree in `src/` contains both native IASI code and a few files copied from the external JURASSIC radiative transfer model for convenience. Local third-party builds are placed in `libs/build/`, CODA definitions are shipped in `data/`, and user-facing example outputs are generated in `projects/example/`.

## Contact

We are interested in sharing the IASI Code Collection for research applications. Please contact:

Dr. Lars Hoffmann  
<l.hoffmann@fz-juelich.de>

Jülich Supercomputing Centre, Forschungszentrum Jülich
