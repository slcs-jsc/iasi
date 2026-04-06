# Example Workflow

The repository includes a concrete end-to-end example in `projects/example/run.sh`. It demonstrates how to move from a raw IASI Level-1C granule to derived perturbation products and quick-look figures.

## Steps

### 1. Configure the runtime environment

The script points to the local executables in `src/`, exports the CODA definitions, and adds locally built libraries to `LD_LIBRARY_PATH`.

### 2. Download a sample granule

If the test input file is missing, the script enters `tests/data/` and runs `download.sh` to fetch and unpack a sample IASI Level-1C NAT file.

### 3. Extract a spectrum

`spec2tab` is used to write one selected footprint to `spec.tab`.

Example:

```bash
./spec2tab - "$f" i 0 30 spec.tab
```

The resulting table contains time, geolocation, wavenumber, brightness temperature, and radiance.

### 4. Compute perturbation products

`perturbation` reads the granule and writes a NetCDF file named `pert.nc` containing brightness temperature and perturbation products for the configured spectral regions.

Example:

```bash
./perturbation - pert.nc "$f"
```

### 5. Produce map tables

`map_pert` extracts geolocated text tables from `pert.nc` for:

- `4mu`
- `15mu_high`
- `15mu_low`

These tables are convenient for plotting, filtering, or downstream analysis.

### 6. Generate figures

The example uses `gnuplot` to create:

- a spectral plot from `spec.tab`
- geolocated brightness temperature maps
- variance maps
- perturbation maps

## Running the example

From the repository root:

```bash
cd src
make
cd ../projects/example
./run.sh
```

## Expected outputs

The example workflow produces:

- `spec.tab`
- `pert.nc`
- `map_4mu.tab`
- `map_15mu_low.tab`
- `map_15mu_high.tab`
- multiple PNG figures in `projects/example/`

## Notes

- `gnuplot` is required for the plotting step.
- The plotting commands reference an external coastline file at `~/wrk/coast/wcl.tab`.
- Generated `.tab`, `.nc`, and image files in `projects/example/` are working artifacts, not source files.
