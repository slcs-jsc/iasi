#! /bin/bash

# Setup...
iasi=~/wrk/iasi/src
export CODA_DEFINITION=../../data/
export LD_LIBRARY_PATH=../../libs/build/lib/:$LD_LIBRARY_PATH

# Download L1C test file...
f=../../tests/data/IASI_xxx_1C_M02_20080101001156Z_20080101015356Z_N_O_20080101015901Z.nat
[ -s "$f" ] || (cd ../../tests/data && ./download.sh)
[ -s "$f" ] || { echo "Download failed!"; exit 1; }

# Get spectra...
$iasi/spec2tab $f i 0 30 spec.tab

# Get perturbations...
$iasi/perturbation pert.nc $f

# Get maps...
$iasi/map_pert - pert.nc map_4mu.tab PERTNAME 4mu
$iasi/map_pert - pert.nc map_15mu_high.tab PERTNAME 15mu_high
$iasi/map_pert - pert.nc map_15mu_low.tab PERTNAME 15mu_low

# Plot...
gnuplot <<EOF
set term pngcairo truecolor font "Helvetica,28" size 1600,1200 crop lw 2

set out "plot.spec.png"
set xla "Wavenumber [1/cm]"
set mxtics
set yla "Brightness temperature [K]"
set mytics
plot "spec.tab" u 6:7 w l t ""

set out "plot.map_4mu.dc.png"
set pal def
set xra [-180:180]
set yra [-90:90]
set xtics 60
set mxtics 6
set ytics 30
set mytics 3
set cbla "Brightness temperature (8 micron) [K]"
plot "map_4mu.tab" u 3:4:5 w p pt 7 lc pal z t "", "~/wrk/coast/wcl.tab" u 1:2:(250) lt -1 t ""


set out "plot.map_4mu.bt.png"
set cbra [*:*]
set cbla "Brightness temperature (4.3 micron) [K]"
plot "map_4mu.tab" u 3:4:6 w p pt 7 lc pal z t "", "~/wrk/coast/wcl.tab" u 1:2:(250) lt -1 t ""

set out "plot.map_4mu.var.png"
set cbra [0:*]
set cbla "Brightness temperature variance (4.3 micron) [K**2]"
plot "map_4mu.tab" u 3:4:8 w p pt 7 lc pal z t "", "~/wrk/coast/wcl.tab" u 1:2:(0) lt -1 t ""

set out "plot.map_4mu.pt.png"
set pal gray
set cbra [-3:3]
set cbla "Brightness temperature perturbation (4.3 micron) [K]"
plot "map_4mu.tab" u 3:4:7 w p pt 7 lc pal z t "", "~/wrk/coast/wcl.tab" u 1:2:(0) lt -1 t ""


set out "plot.map_15mu_low.bt.png"
set pal def
set cbra [*:*]
set cbla "Brightness temperature (15 micron low) [K]"
plot "map_15mu_low.tab" u 3:4:6 w p pt 7 lc pal z t "", "~/wrk/coast/wcl.tab" u 1:2:(250) lt -1 t ""

set out "plot.map_15mu_low.var.png"
set cbra [0:*]
set cbla "Brightness temperature variance (15 micron low) [K**2]"
plot "map_15mu_low.tab" u 3:4:8 w p pt 7 lc pal z t "", "~/wrk/coast/wcl.tab" u 1:2:(0) lt -1 t ""

set out "plot.map_15mu_low.pt.png"
set pal gray
set cbra [-3:3]
set cbla "Brightness temperature perturbation (15 micron low) [K]"
plot "map_15mu_low.tab" u 3:4:7 w p pt 7 lc pal z t "", "~/wrk/coast/wcl.tab" u 1:2:(0) lt -1 t ""


set out "plot.map_15mu_high.bt.png"
set pal def
set cbra [*:*]
set cbla "Brightness temperature (15 micron high) [K]"
plot "map_15mu_high.tab" u 3:4:6 w p pt 7 lc pal z t "", "~/wrk/coast/wcl.tab" u 1:2:(250) lt -1 t ""

set out "plot.map_15mu_high.var.png"
set cbra [0:*]
set cbla "Brightness temperature variance (15 micron high) [K**2]"
plot "map_15mu_high.tab" u 3:4:8 w p pt 7 lc pal z t "", "~/wrk/coast/wcl.tab" u 1:2:(0) lt -1 t ""

set out "plot.map_15mu_high.pt.png"
set pal gray
set cbra [-3:3]
set cbla "Brightness temperature perturbation (15 micron high) [K]"
plot "map_15mu_high.tab" u 3:4:7 w p pt 7 lc pal z t "", "~/wrk/coast/wcl.tab" u 1:2:(0) lt -1 t ""

EOF
