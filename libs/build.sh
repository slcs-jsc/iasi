#! /bin/bash

# ----------------------------------------------------------------------
function infomsg {
# ----------------------------------------------------------------------
    
    echo
    echo "============================================================"
    echo "Compile: $1"
    echo "============================================================"
    echo
}

# ----------------------------------------------------------------------
# Main...
# ----------------------------------------------------------------------


# Set number of cores...
THREADS=$(cat /proc/cpuinfo | grep processor | wc -l)

# Get absolute path...
target=$(mkdir -p build && cd build && pwd)

# Prepare directories...
mkdir -p $target/src $target/bin $target/lib $target/man/man1 \
    && cp *tar.gz $target/src \
    && cd $target/src \
    && for f in $(ls *tar.gz) ; do tar xvzf $f ; done \
    || exit

# GSL...
dir=gsl-2.5
infomsg $dir
cd $target/src/$dir \
    && ./configure --prefix=$target \
    && make -j$THREADS && make check && make install && make clean \
    || exit

# netCDF...
dir=netcdf-c-4.6.2
infomsg $dir
cd $target/src/$dir \
    && ./configure --prefix=$target --enable-c-only --disable-dap --disable-netcdf-4\
    && make -j$THREADS && make check && make install \
    || exit

# coda
dir=coda-2.20
infomsg $dir
cd $target/src/$dir \
    && ./configure --prefix=$target\
    && make -j$THREADS && make check && make install \
    || exit
