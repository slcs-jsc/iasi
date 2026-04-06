# Tools

The `src/Makefile` builds the following command-line executables.

## Core IASI tools

- `spec2tab`: extracts a single IASI spectrum into a tabular text file
- `perturbation`: computes 4.3 micron and 15 micron brightness temperature products and perturbations, then writes a NetCDF file
- `map_pert`: converts perturbation NetCDF output into geolocated map tables and can optionally recompute background, filtering, and variance fields
- `bands`: averages radiances over configurable spectral bands and writes brightness temperatures to a table
- `noise`: estimates noise statistics from radiance data and reports mean brightness temperature, NEDT, and NESR
- `extract`: prepares radiance and meteorological inputs for retrieval workflows
- `retrieval`: MPI-enabled retrieval processor for IASI products

## Utility programs

- `day2doy`: converts calendar dates to day-of-year values
- `doy2day`: converts day-of-year values to calendar dates
- `time2jsec`: converts dates to Julian seconds
- `jsec2time`: converts Julian seconds back to date/time fields

The time conversion utilities and the shared files `jurassic.c` and `jurassic.h` come from the external JURASSIC radiative transfer model:

- <https://github.com/slcs-jsc/jurassic>

## Common calling pattern

Most programs follow the same structure:

```text
<program> <ctl> <output-or-input> ...
```

The first argument is a control-file placeholder interpreted through the shared `scan_ctl(...)` infrastructure in the code. In practice, many workflows pass `-` and provide options as trailing key/value pairs.

## Selected command summaries

### `spec2tab`

Usage:

```text
spec2tab <ctl> <iasi_l1b_file> [index <track> <xtrack> | geo <lon> <lat>] <spec.tab>
```

Behavior:

- reads one IASI Level-1C granule
- selects a footprint either by track/cross-track index or nearest geolocation
- writes time, geolocation, wavenumber, brightness temperature, and radiance columns

### `perturbation`

Usage:

```text
perturbation <ctl> <out.nc> <l1b_file1> [<l1b_file2> ...]
```

Behavior:

- reads one or more IASI Level-1C granules
- derives brightness temperature products in the 4 micron and 15 micron bands
- stores geolocation, brightness temperature, perturbation, and variance fields in NetCDF

### `map_pert`

Usage:

```text
map_pert <ctl> <pert.nc> <map.tab>
```

Behavior:

- reads a perturbation NetCDF product
- selects the requested perturbation family through `PERTNAME`
- can apply background fitting, smoothing, Gaussian, Hamming, median, and variance processing options
- writes a geolocated table containing time, solar zenith angle, position, brightness temperature, perturbation, variance, and scan indices

The example workflow uses:

```bash
map_pert - pert.nc map_4mu.tab PERTNAME 4mu
map_pert - pert.nc map_15mu_high.tab PERTNAME 15mu_high
map_pert - pert.nc map_15mu_low.tab PERTNAME 15mu_low
```

### `bands`

Usage:

```text
bands <ctl> <out.tab> <l1b_file1> [<l1b_file2> ...]
```

Behavior:

- reads one or more Level-1C granules
- averages radiances within user-configured spectral intervals
- converts those band means to brightness temperatures

The spectral intervals are configured through repeated `NUMIN` and `NUMAX` control parameters.

### `noise`

Usage:

```text
noise <ctl> <iasi_l1_file> <noise.tab>
```

Behavior:

- computes blockwise brightness-temperature noise diagnostics from one Level-1C granule
- writes track index, channel index, wavenumber, mean BT, NEDT, and NESR
