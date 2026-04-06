# Getting Started

This page summarizes the local build and runtime setup for the IASI Code Collection on Linux.

## Repository structure

The main directories are:

- `src/` for C sources, headers, and the project `Makefile`
- `libs/` for vendored third-party source archives and the local build helper
- `data/` for CODA definition files
- `tests/data/` for sample input download helpers
- `projects/example/` for an end-to-end processing and plotting workflow
- `docs/` for MkDocs and Doxygen configuration

## Requirements

The code is built with `gcc` and `make`. The processing tools depend on:

- GSL
- netCDF
- CODA
- OpenMP
- MPI for the `retrieval` executable

If these libraries are not already available on your system, build the local copies shipped with the repository:

```bash
cd libs
./build.sh
```

## Build

Compile the tools from `src/`:

```bash
cd src
make
```

Useful build variants:

- `make STATIC=1` for static linking
- `make PROF=1` for profiling flags
- `make COV=1` for coverage instrumentation
- `make INFO=1` for compiler optimization reports

By default the build uses strict warnings and treats them as errors via `-Werror`.

## Runtime environment

Several workflows expect local runtime configuration. The example pipeline uses:

```bash
export CODA_DEFINITION=../../data/
export LD_LIBRARY_PATH=../../libs/build/lib/:$LD_LIBRARY_PATH
```

`CODA_DEFINITION` should point to the directory containing the shipped `.codadef` files. If you built local libraries in `libs/build/`, add its `lib/` directory to `LD_LIBRARY_PATH`.

## Validation

There is no dedicated unit test suite in this repository. The main practical validation path is:

1. Rebuild the affected binaries in `src/`.
2. Run the example workflow in `projects/example/`.
3. Inspect the generated tabular outputs and PNG figures.

The example downloads a sample IASI Level-1C granule through `tests/data/download.sh` if it is not already available.
