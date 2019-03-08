#include "libiasi.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Maximum number of bands... */
#define NB 100

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static iasi_rad_t *iasi_rad;

  static FILE *out;

  static double numin[NB], numax[NB], rad[NB];

  static int iarg, ib, ichan, n, nb, track, xtrack;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <out.tab> <l1b_file1> [<l1b_file2> ...]");

  /* Allocate... */
  ALLOC(iasi_rad, iasi_rad_t, 1);

  /* Get control parameters... */
  nb = (int) scan_ctl(argc, argv, "NB", -1, "1", NULL);
  if (nb > NB)
    ERRMSG("Too many bands!");
  for (ib = 0; ib < nb; ib++) {
    numin[ib] = scan_ctl(argc, argv, "NUMIN", ib, "", NULL);
    numax[ib] = scan_ctl(argc, argv, "NUMAX", ib, "", NULL);
  }

  /* Create file... */
  printf("Write band data: %s\n", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Loop over IASI files... */
  for (iarg = 3; iarg < argc; iarg++) {

    /* Read IASI data... */
    printf("Read IASI Level-1C data file: %s\n", argv[iarg]);
    iasi_read(argv[iarg], iasi_rad);

    /* Write header... */
    if (iarg == 3) {
      fprintf(out,
	      "# $1 = time [s]\n"
	      "# $2 = footprint longitude [deg]\n"
	      "# $3 = footprint latitude [deg]\n"
	      "# $4 = satellite altitude [km]\n"
	      "# $5 = satellite longitude [deg]\n"
	      "# $6 = satellite latitude [deg]\n");
      for (ib = 0; ib < nb; ib++)
	fprintf(out,
		"# $%d = BT(%.2f/cm...%.2f/cm) [K]\n",
		7 + ib, numin[ib], numax[ib]);
    }

    /* Loop over scans... */
    for (track = 0; track < iasi_rad->ntrack; track++) {

      /* Write output... */
      fprintf(out, "\n");

      /* Loop over footprints... */
      for (xtrack = 0; xtrack < L1_NXTRACK; xtrack++) {

	/* Write output... */
	fprintf(out, "%.2f %.4f %.4f %.3f %.4f %.4f",
		iasi_rad->Time[track][xtrack],
		iasi_rad->Longitude[track][xtrack],
		iasi_rad->Latitude[track][xtrack],
		iasi_rad->Sat_z[track],
		iasi_rad->Sat_lon[track], iasi_rad->Sat_lat[track]);

	/* Loop over bands... */
	for (ib = 0; ib < nb; ib++) {

	  /* Get mean radiance... */
	  n = 0;
	  rad[ib] = 0;
	  for (ichan = 0; ichan <= IASI_L1_NCHAN; ichan++)
	    if (iasi_rad->freq[ichan] >= numin[ib]
		&& iasi_rad->freq[ichan] <= numax[ib]
		&& gsl_finite(iasi_rad->Rad[track][xtrack][ichan])) {
	      rad[ib] += iasi_rad->Rad[track][xtrack][ichan];
	      n++;
	    }
	  if (n > 0)
	    rad[ib] /= n;
	  else
	    rad[ib] = GSL_NAN;

	  /* Convert to brightness temperature... */
	  rad[ib] = brightness(rad[ib], 0.5 * (numin[ib] + numax[ib]));

	  /* Write output... */
	  fprintf(out, " %.3f", rad[ib]);
	}

	/* Write output... */
	fprintf(out, "\n");
      }
    }
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(iasi_rad);

  return EXIT_SUCCESS;
}
