#include "libiasi.h"

int main(
  int argc,
  char *argv[]) {

  static iasi_rad_t *iasi_rad;

  static wave_t wave;

  static FILE *out;

  static double mu, sigma;

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
	  "# $4 = mean BT [K]\n" "# $5 = NEDT [K]\n");

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
	  wave.temp[ix][wave.ny] = brightness(iasi_rad->Rad[iy][ix][ichan],
					      iasi_rad->freq[ichan]);
	wave.ny++;
      }

      /* Check number of data points... */
      if (wave.ny >= 55) {

	/* Get noise... */
	noise(&wave, &mu, &sigma);

	/* Write output... */
	if (gsl_finite(sigma))
	  fprintf(out, "%d %d %.4f %g %g\n", itrack, ichan,
		  iasi_rad->freq[ichan], mu, sigma);
      }
    }
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(iasi_rad);

  return EXIT_SUCCESS;
}
