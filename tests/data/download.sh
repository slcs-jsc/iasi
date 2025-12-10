#! /bin/bash

# Setup...
url=https://datapub.fz-juelich.de/slcs/iasi/projects/testdata/testdata.tar.gz
f=IASI_xxx_1C_M02_20080101001156Z_20080101015356Z_N_O_20080101015901Z.nat

# Check whether file exists...
if [ ! -s $f ] ; then
    
    # Download and extract test data...
    wget $url && tar -xzf testdata.tar.gz && rm -f testdata.tar.gz || exit

fi

# Write info...
du -h $f
