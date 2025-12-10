## IASI Code Collection — Example Workflow

This example demonstrates a processing workflow using the **IASI Code Collection** tools. It shows how to download a sample IASI Level-1C file, extract radiance spectra, compute brightness temperature perturbation fields, generate geolocated maps, and produce visualizations using `gnuplot`.

### What This Example Does

1. **Environment Setup**  
   The script configures the required runtime environment, including the IASI executable path, CODA definitions, and shared libraries.

2. **Download Test Input**  
   It verifies the presence of an example IASI L1C NAT file.  
   If missing, it automatically downloads the file from a public data server.
   The script terminates if the download fails.

3. **Extract Spectra**  
   Using the `spec2tab` tool, the script reads a subset of IASI spectra and stores them in tabular form (`spec.tab`).  
   These data can be inspected, plotted, or used as input to other tools.

4. **Generate Perturbations**  
   The `perturbation` tool computes 4.3 and 15 micron brightness temperature perturbations and stores them in a NetCDF file (`pert.nc`).

5. **Create Maps**  
   The `map_pert` tool converts perturbation information into geolocated map tables for different spectral regions.

6. **Plot Results**  
   Finally, `gnuplot` produces several diagnostic PNG figures:
   - The extracted spectrum
   - Brightness temperature fields
   - Variance fields
   - Perturbation amplitude maps

   These visualizations illustrate how brightness temperature perturbations vary across the Earth for selected IASI channels.

### Output Products

After running the script, you should see:

- `spec.tab` — extracted spectra  
- `pert.nc` — perturbation dataset  
- `map_4mu.tab`, `map_15mu_low.tab`, `map_15mu_high.tab` — map tables  
- Multiple PNG images showing spectra and perturbation maps
