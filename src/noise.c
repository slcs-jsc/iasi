/*
  This file is part of the IASI Code Collection.
  
  the IASI Code Collections is free software: you can redistribute it
  and/or modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.
  
  The IASI Code Collection is distributed in the hope that it will be
  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with the IASI Code Collection. If not, see
  <http://www.gnu.org/licenses/>.
  
  Copyright (C) 2019-2025 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Get noise estimates from radiance data.
*/

#include "libiasi.h"

int main(
  int argc,
  char *argv[]) {

  static iasi_rad_t *iasi_rad;

  static wave_t wave;

  static FILE *out;

  static double mu, nesr, sigma;

  static int ichan, itrack, ix, iy;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <iasi_l1_file> <noise.tab>");

  /* Allocate... */
  ALLOC(iasi_rad, iasi_rad_t, 1);

  /* Read IASI data... */
  printf("Read IASI data: %s\n", argv[2]);
  iasi_read(argv[2], iasi_rad);

  /* Create file... */
  printf("Write noise data: %s\n", argv[3]);
  if (!(out = fopen(argv[3], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = track index\n"
	  "# $2 = channel index\n"
	  "# $3 = wavenumber [1/cm]\n"
	  "# $4 = mean BT [K]\n"
	  "# $5 = NEDT [K]\n" "# $6 = NESR [W/(m^2 sr cm^-1)]\n");

  /* Analyze blocks of data... */
  for (itrack = 0; itrack < iasi_rad->ntrack; itrack += 60) {

    /* Write empty line... */
    fprintf(out, "\n");

    /* Loop over channels... */
    for (ichan = 0; ichan < IASI_L1_NCHAN; ichan++) {

      /* Set wave struct... */
      wave.nx = L1_NXTRACK;
      wave.ny = 0;
      for (iy = itrack; iy < GSL_MIN(itrack + 60, iasi_rad->ntrack); iy++) {
	for (ix = 0; ix < wave.nx; ix++)
	  wave.temp[ix][wave.ny] = BRIGHT(iasi_rad->Rad[iy][ix][ichan],
					  iasi_rad->freq[ichan]);
	wave.ny++;
      }

      /* Check number of data points... */
      if (wave.ny >= 55) {

	/* Get noise... */
	noise(&wave, &mu, &sigma);

	/* Get NESR... */
	nesr = PLANCK(mu + sigma, iasi_rad->freq[ichan])
	  - PLANCK(mu, iasi_rad->freq[ichan]);

	/* Write output... */
	if (gsl_finite(sigma))
	  fprintf(out, "%d %d %.4f %g %g %g\n", itrack, ichan,
		  iasi_rad->freq[ichan], mu, sigma, nesr);
      }
    }
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(iasi_rad);

  return EXIT_SUCCESS;
}
