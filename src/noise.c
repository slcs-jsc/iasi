#include "libiasi.h"

int main(
  int argc,
  char *argv[]) {

  static iasi_rad_t *iasi_rad;

  static wave_t wave;

  static FILE *out;

  static double mu, sigma;

  static int ichan, ix, iy;

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
	  "# $1 = wavenumber [1/cm]\n" "# $2 = NESR [W/(m^2 sr cm^-1)]\n\n");

  /* Loop over channels... */
  for (ichan = 0; ichan < IASI_L1_NCHAN; ichan++) {

    /* Set wave struct... */
    wave.nx = L1_NXTRACK;
    wave.ny = iasi_rad->ntrack;
    for (ix = 0; ix < wave.nx; ix++)
      for (iy = 0; iy < wave.ny; iy++)
	wave.temp[ix][iy] = iasi_rad->Rad[iy][ix][ichan];

    /* Get noise... */
    noise(&wave, &mu, &sigma);

    /* Write output... */
    fprintf(out, "%.4f %g\n", iasi_rad->freq[ichan], sigma);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(iasi_rad);

  return EXIT_SUCCESS;
}
