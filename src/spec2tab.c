#include "libiasi.h"

int main(
  int argc,
  char *argv[]) {

  static iasi_rad_t *iasi_rad;

  FILE *out;

  double dmin = 1e100, x0[3], x1[3];

  int ichan, track = -1, track2, xtrack = -1, xtrack2;

  /* Check arguments... */
  if (argc != 6)
    ERRMSG("Give parameters: <iasi_l1b_file> "
	   "[index <track> <xtrack> | geo <lon> <lat>] <spec.tab>");

  /* Allocate... */
  ALLOC(iasi_rad, iasi_rad_t, 1);

  /* Read IASI data... */
  printf("Read IASI Level-1C data file: %s\n", argv[1]);
  iasi_read(argv[1], iasi_rad);

  /* Get indices... */
  if (argv[2][0] == 'i') {
    track = atoi(argv[3]);
    xtrack = atoi(argv[4]);
  }

  /* Find nearest footprint... */
  else {
    geo2cart(0, atof(argv[3]), atof(argv[4]), x0);
    for (track2 = 0; track2 < iasi_rad->ntrack; track2++)
      for (xtrack2 = 0; xtrack2 < L1_NXTRACK; xtrack2++) {
	geo2cart(0, iasi_rad->Longitude[track2][xtrack2],
		 iasi_rad->Latitude[track2][xtrack2], x1);
	if (DIST2(x0, x1) < dmin) {
	  dmin = DIST2(x0, x1);
	  track = track2;
	  xtrack = xtrack2;
	}
      }
    if (dmin > 2500)
      ERRMSG("Geolocation not covered by granule!");
  }

  /* Check indices... */
  if (track < 0 || track >= iasi_rad->ntrack)
    ERRMSG("Along-track index out of range!");
  if (xtrack < 0 || xtrack >= L1_NXTRACK)
    ERRMSG("Across-track index out of range!");

  /* Create file... */
  printf("Write spectrum: %s\n", argv[5]);
  if (!(out = fopen(argv[5], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time (seconds since 01-JAN-2000, 00:00 UTC)\n"
	  "# $2 = satellite longitude [deg]\n"
	  "# $3 = satellite latitude [deg]\n"
	  "# $4 = footprint longitude [deg]\n"
	  "# $5 = footprint latitude [deg]\n"
	  "# $6 = wavenumber [cm^-1]\n"
	  "# $7 = brightness temperature [K]\n"
	  "# $8 = radiance [W/(m^2 sr cm^-1)]\n\n");

  /* Write data... */
  for (ichan = 0; ichan < IASI_L1_NCHAN; ichan++)
    fprintf(out, "%.2f %g %g %g %g %g %g %g\n",
	    iasi_rad->Time[track][xtrack],
	    iasi_rad->Sat_lon[track],
	    iasi_rad->Sat_lat[track],
	    iasi_rad->Longitude[track][xtrack],
	    iasi_rad->Latitude[track][xtrack],
	    iasi_rad->freq[ichan],
	    brightness(iasi_rad->Rad[track][xtrack][ichan],
		       iasi_rad->freq[ichan]),
	    iasi_rad->Rad[track][xtrack][ichan]);

  /* Close file... */
  fclose(out);

  /* Free... */
  free(iasi_rad);

  return EXIT_SUCCESS;
}
