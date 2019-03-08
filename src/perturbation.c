#include "libiasi.h"

/* ------------------------------------------------------------
   Constants...
   ------------------------------------------------------------ */

/* Number of 4 micron channels: */
#define N4 152

/* Number of 15 micron channels (low altitudes): */
#define N15_LOW 21

/* Number of 15 micron channels (high altitudes): */
#define N15_HIGH 2

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Add variable defintions to netCDF file. */
void addatt(
  int ncid,
  int varid,
  const char *unit,
  const char *long_name);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static iasi_rad_t *iasi_rad;

  static pert_t *pert_4mu, *pert_15mu_low, *pert_15mu_high;

  static wave_t wave;

  static double numean, radmean, var_dh = 100.;

  static int list_4mu[N4]
    = { 6711, 6712, 6713, 6714, 6715, 6716, 6717, 6718, 6719, 6720,
    6721, 6722, 6723, 6724, 6725, 6726, 6727, 6728, 6729, 6730, 6731,
    6732, 6733, 6734, 6735, 6736, 6737, 6738, 6739, 6740, 6741, 6742,
    6743, 6744, 6745, 6746, 6747, 6748, 6749, 6750, 6751, 6752, 6753,
    6754, 6755, 6756, 6757, 6758, 6759, 6760, 6761, 6762, 6763, 6764,
    6765, 6766, 6767, 6768, 6769, 6770, 6771, 6772, 6773, 6774, 6775,
    6776, 6777, 6778, 6779, 6780, 6781, 6782, 6783, 6784, 6785, 6786,
    6787, 6788, 6789, 6790, 6791, 6792, 6793, 6794, 6795, 6796, 6797,
    6798, 6799, 6800, 6801, 6802, 6803, 6804, 6830, 6831, 6832, 6833,
    6834, 6835, 6836, 6837, 6838, 6839, 6840, 6841, 6842, 6843, 6844,
    6845, 6846, 6847, 6848, 6849, 6850, 6851, 6852, 6853, 6854, 6855,
    6856, 6857, 6858, 6859, 6860, 6861, 6862, 6863, 6864, 6865, 6866,
    6867, 6868, 6869, 6870, 6871, 6872, 6873, 6874, 6875, 6876, 6877,
    6878, 6879, 6880, 6881, 6882, 6883, 6884, 6885, 6886, 6887
  };

  static int list_15mu_low[N15_LOW]
    = { 22, 28, 34, 40, 46, 52, 58, 72, 100, 105, 112, 118, 119,
    124, 125, 130, 131, 136, 137, 143, 144
  };

  static int list_15mu_high[N15_HIGH]
  = { 91, 92 };

  static int ix, iy, dimid[2], i, n, ncid, track, track0, xtrack,
    time_varid, lon_varid, lat_varid, bt_4mu_varid, bt_4mu_pt_varid,
    bt_4mu_var_varid, bt_8mu_varid, bt_15mu_low_varid, bt_15mu_low_pt_varid,
    bt_15mu_low_var_varid, bt_15mu_high_varid, bt_15mu_high_pt_varid,
    bt_15mu_high_var_varid, iarg;

  static size_t start[2], count[2];

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <out.nc> <l1b_file1> [<l1b_file2> ...]");

  /* Allocate... */
  ALLOC(iasi_rad, iasi_rad_t, 1);
  ALLOC(pert_4mu, pert_t, 1);
  ALLOC(pert_15mu_low, pert_t, 1);
  ALLOC(pert_15mu_high, pert_t, 1);

  /* ------------------------------------------------------------
     Read HDF files...
     ------------------------------------------------------------ */

  /* Loop over HDF files... */
  for (iarg = 2; iarg < argc; iarg++) {

    /* Read IASI data... */
    printf("Read IASI Level-1C data file: %s\n", argv[iarg]);
    iasi_read(argv[iarg], iasi_rad);

    /* Save geolocation... */
    pert_4mu->ntrack += iasi_rad->ntrack;
    if (pert_4mu->ntrack > PERT_NTRACK)
      ERRMSG("Too many granules!");
    pert_4mu->nxtrack = L1_NXTRACK;
    if (pert_4mu->nxtrack > PERT_NXTRACK)
      ERRMSG("Too many tracks!");
    for (track = 0; track < iasi_rad->ntrack; track++)
      for (xtrack = 0; xtrack < L1_NXTRACK; xtrack++) {
	pert_4mu->time[track0 + track][xtrack]
	  = iasi_rad->Time[track][xtrack];
	pert_4mu->lon[track0 + track][xtrack]
	  = iasi_rad->Longitude[track][xtrack];
	pert_4mu->lat[track0 + track][xtrack]
	  = iasi_rad->Latitude[track][xtrack];
      }

    pert_15mu_low->ntrack += iasi_rad->ntrack;
    if (pert_15mu_low->ntrack > PERT_NTRACK)
      ERRMSG("Too many granules!");
    pert_15mu_low->nxtrack = L1_NXTRACK;
    if (pert_15mu_low->nxtrack > PERT_NXTRACK)
      ERRMSG("Too many tracks!");
    for (track = 0; track < iasi_rad->ntrack; track++)
      for (xtrack = 0; xtrack < L1_NXTRACK; xtrack++) {
	pert_15mu_low->time[track0 + track][xtrack]
	  = iasi_rad->Time[track][xtrack];
	pert_15mu_low->lon[track0 + track][xtrack]
	  = iasi_rad->Longitude[track][xtrack];
	pert_15mu_low->lat[track0 + track][xtrack]
	  = iasi_rad->Latitude[track][xtrack];
      }

    pert_15mu_high->ntrack += iasi_rad->ntrack;
    if (pert_15mu_high->ntrack > PERT_NTRACK)
      ERRMSG("Too many granules!");
    pert_15mu_high->nxtrack = L1_NXTRACK;
    if (pert_15mu_high->nxtrack > PERT_NXTRACK)
      ERRMSG("Too many tracks!");
    for (track = 0; track < iasi_rad->ntrack; track++)
      for (xtrack = 0; xtrack < L1_NXTRACK; xtrack++) {
	pert_15mu_high->time[track0 + track][xtrack]
	  = iasi_rad->Time[track][xtrack];
	pert_15mu_high->lon[track0 + track][xtrack]
	  = iasi_rad->Longitude[track][xtrack];
	pert_15mu_high->lat[track0 + track][xtrack]
	  = iasi_rad->Latitude[track][xtrack];
      }

    /* Get 8.1 micron brightness temperature... */
    for (track = 0; track < iasi_rad->ntrack; track++)
      for (xtrack = 0; xtrack < L1_NXTRACK; xtrack++)
	pert_4mu->dc[track0 + track][xtrack]
	  = brightness(iasi_rad->Rad[track][xtrack][2345],
		       iasi_rad->freq[2345]);

    /* Get 4.3 micron brightness temperature... */
    for (track = 0; track < iasi_rad->ntrack; track++)
      for (xtrack = 0; xtrack < L1_NXTRACK; xtrack++) {
	n = 0;
	numean = radmean = 0;
	for (i = 0; i < N4; i++)
	  if (gsl_finite(iasi_rad->Rad[track][xtrack][list_4mu[i]])) {
	    radmean += iasi_rad->Rad[track][xtrack][list_4mu[i]];
	    numean += iasi_rad->freq[list_4mu[i]];
	    n++;
	  }
	if (n > 0.9 * N4)
	  pert_4mu->bt[track0 + track][xtrack]
	    = brightness(radmean / n, numean / n);
	else
	  pert_4mu->bt[track0 + track][xtrack] = GSL_NAN;
      }

    /* Get 15 micron brightness temperature (low altitudes)... */
    for (track = 0; track < iasi_rad->ntrack; track++)
      for (xtrack = 0; xtrack < L1_NXTRACK; xtrack++) {
	n = 0;
	numean = radmean = 0;
	for (i = 0; i < N15_LOW; i++)
	  if (gsl_finite(iasi_rad->Rad[track][xtrack][list_15mu_low[i]])) {
	    radmean += iasi_rad->Rad[track][xtrack][list_15mu_low[i]];
	    numean += iasi_rad->freq[list_15mu_low[i]];
	    n++;
	  }
	if (n > 0.9 * N15_LOW)
	  pert_15mu_low->bt[track0 + track][xtrack]
	    = brightness(radmean / n, numean / n);
	else
	  pert_15mu_low->bt[track0 + track][xtrack] = GSL_NAN;
      }

    /* Get 15 micron brightness temperature (high altitudes)... */
    for (track = 0; track < iasi_rad->ntrack; track++)
      for (xtrack = 0; xtrack < L1_NXTRACK; xtrack++) {
	n = 0;
	numean = radmean = 0;
	for (i = 0; i < N15_HIGH; i++)
	  if (gsl_finite(iasi_rad->Rad[track][xtrack][list_15mu_high[i]])) {
	    radmean += iasi_rad->Rad[track][xtrack][list_15mu_high[i]];
	    numean += iasi_rad->freq[list_15mu_high[i]];
	    n++;
	  }
	if (n > 0.9 * N15_HIGH)
	  pert_15mu_high->bt[track0 + track][xtrack]
	    = brightness(radmean / n, numean / n);
	else
	  pert_15mu_high->bt[track0 + track][xtrack] = GSL_NAN;
      }

    /* Increment track counter... */
    track0 += iasi_rad->ntrack;
  }

  /* ------------------------------------------------------------
     Calculate perturbations and variances...
     ------------------------------------------------------------ */

  /* Convert to wave analysis struct... */
  pert2wave(pert_4mu, &wave,
	    0, pert_4mu->ntrack - 1, 0, pert_4mu->nxtrack - 1);

  /* Estimate background... */
  background_poly(&wave, 5, 0);

  /* Compute variance... */
  variance(&wave, var_dh);

  /* Copy data... */
  for (ix = 0; ix < wave.nx; ix++)
    for (iy = 0; iy < wave.ny; iy++) {
      pert_4mu->pt[iy][ix] = wave.pt[ix][iy];
      pert_4mu->var[iy][ix] = wave.var[ix][iy];
    }

  /* Convert to wave analysis struct... */
  pert2wave(pert_15mu_low, &wave,
	    0, pert_15mu_low->ntrack - 1, 0, pert_15mu_low->nxtrack - 1);

  /* Estimate background... */
  background_poly(&wave, 5, 0);

  /* Compute variance... */
  variance(&wave, var_dh);

  /* Copy data... */
  for (ix = 0; ix < wave.nx; ix++)
    for (iy = 0; iy < wave.ny; iy++) {
      pert_15mu_low->pt[iy][ix] = wave.pt[ix][iy];
      pert_15mu_low->var[iy][ix] = wave.var[ix][iy];
    }

  /* Convert to wave analysis struct... */
  pert2wave(pert_15mu_high, &wave,
	    0, pert_15mu_high->ntrack - 1, 0, pert_15mu_high->nxtrack - 1);

  /* Estimate background... */
  background_poly(&wave, 5, 0);

  /* Compute variance... */
  variance(&wave, var_dh);

  /* Copy data... */
  for (ix = 0; ix < wave.nx; ix++)
    for (iy = 0; iy < wave.ny; iy++) {
      pert_15mu_high->pt[iy][ix] = wave.pt[ix][iy];
      pert_15mu_high->var[iy][ix] = wave.var[ix][iy];
    }

  /* ------------------------------------------------------------
     Write to netCDF file...
     ------------------------------------------------------------ */

  /* Create netCDF file... */
  printf("Write perturbation data file: %s\n", argv[1]);
  NC(nc_create(argv[1], NC_CLOBBER, &ncid));

  /* Set dimensions... */
  NC(nc_def_dim(ncid, "NTRACK", NC_UNLIMITED, &dimid[0]));
  NC(nc_def_dim(ncid, "NXTRACK", L1_NXTRACK, &dimid[1]));

  /* Add variables... */
  NC(nc_def_var(ncid, "time", NC_DOUBLE, 2, dimid, &time_varid));
  addatt(ncid, time_varid, "s", "time (seconds since 2000-01-01T00:00Z)");
  NC(nc_def_var(ncid, "lon", NC_DOUBLE, 2, dimid, &lon_varid));
  addatt(ncid, lon_varid, "deg", "footprint longitude");
  NC(nc_def_var(ncid, "lat", NC_DOUBLE, 2, dimid, &lat_varid));
  addatt(ncid, lat_varid, "deg", "footprint latitude");

  NC(nc_def_var(ncid, "bt_8mu", NC_FLOAT, 2, dimid, &bt_8mu_varid));
  addatt(ncid, bt_8mu_varid, "K", "brightness temperature at 8.1 micron");

  NC(nc_def_var(ncid, "bt_4mu", NC_FLOAT, 2, dimid, &bt_4mu_varid));
  addatt(ncid, bt_4mu_varid, "K", "brightness temperature" " at 4.3 micron");
  NC(nc_def_var(ncid, "bt_4mu_pt", NC_FLOAT, 2, dimid, &bt_4mu_pt_varid));
  addatt(ncid, bt_4mu_pt_varid, "K", "brightness temperature perturbation"
	 " at 4.3 micron");
  NC(nc_def_var(ncid, "bt_4mu_var", NC_FLOAT, 2, dimid, &bt_4mu_var_varid));
  addatt(ncid, bt_4mu_var_varid, "K^2", "brightness temperature variance"
	 " at 4.3 micron");

  NC(nc_def_var(ncid, "bt_15mu_low", NC_FLOAT, 2, dimid, &bt_15mu_low_varid));
  addatt(ncid, bt_15mu_low_varid, "K", "brightness temperature"
	 " at 15 micron (low altitudes)");
  NC(nc_def_var(ncid, "bt_15mu_low_pt", NC_FLOAT, 2, dimid,
		&bt_15mu_low_pt_varid));
  addatt(ncid, bt_15mu_low_pt_varid, "K",
	 "brightness temperature perturbation"
	 " at 15 micron (low altitudes)");
  NC(nc_def_var
     (ncid, "bt_15mu_low_var", NC_FLOAT, 2, dimid, &bt_15mu_low_var_varid));
  addatt(ncid, bt_15mu_low_var_varid, "K^2",
	 "brightness temperature variance" " at 15 micron (low altitudes)");

  NC(nc_def_var(ncid, "bt_15mu_high", NC_FLOAT, 2, dimid,
		&bt_15mu_high_varid));
  addatt(ncid, bt_15mu_high_varid, "K", "brightness temperature"
	 " at 15 micron (high altitudes)");
  NC(nc_def_var(ncid, "bt_15mu_high_pt", NC_FLOAT, 2, dimid,
		&bt_15mu_high_pt_varid));
  addatt(ncid, bt_15mu_high_pt_varid, "K",
	 "brightness temperature perturbation"
	 " at 15 micron (high altitudes)");
  NC(nc_def_var
     (ncid, "bt_15mu_high_var", NC_FLOAT, 2, dimid, &bt_15mu_high_var_varid));
  addatt(ncid, bt_15mu_high_var_varid, "K^2",
	 "brightness temperature variance" " at 15 micron (high altitudes)");

  /* Leave define mode... */
  NC(nc_enddef(ncid));

  /* Loop over tracks... */
  for (track = 0; track < pert_4mu->ntrack; track++) {

    /* Set array sizes... */
    start[0] = (size_t) track;
    start[1] = 0;
    count[0] = 1;
    count[1] = (size_t) pert_4mu->nxtrack;

    /* Write data... */
    NC(nc_put_vara_double(ncid, time_varid, start, count,
			  pert_4mu->time[track]));
    NC(nc_put_vara_double(ncid, lon_varid, start, count,
			  pert_4mu->lon[track]));
    NC(nc_put_vara_double(ncid, lat_varid, start, count,
			  pert_4mu->lat[track]));

    NC(nc_put_vara_double(ncid, bt_8mu_varid, start, count,
			  pert_4mu->dc[track]));

    NC(nc_put_vara_double(ncid, bt_4mu_varid, start, count,
			  pert_4mu->bt[track]));
    NC(nc_put_vara_double(ncid, bt_4mu_pt_varid, start, count,
			  pert_4mu->pt[track]));
    NC(nc_put_vara_double(ncid, bt_4mu_var_varid, start, count,
			  pert_4mu->var[track]));

    NC(nc_put_vara_double(ncid, bt_15mu_low_varid, start, count,
			  pert_15mu_low->bt[track]));
    NC(nc_put_vara_double(ncid, bt_15mu_low_pt_varid, start, count,
			  pert_15mu_low->pt[track]));
    NC(nc_put_vara_double(ncid, bt_15mu_low_var_varid, start, count,
			  pert_15mu_low->var[track]));

    NC(nc_put_vara_double(ncid, bt_15mu_high_varid, start, count,
			  pert_15mu_high->bt[track]));
    NC(nc_put_vara_double(ncid, bt_15mu_high_pt_varid, start, count,
			  pert_15mu_high->pt[track]));
    NC(nc_put_vara_double(ncid, bt_15mu_high_var_varid, start, count,
			  pert_15mu_high->var[track]));
  }

  /* Close file... */
  NC(nc_close(ncid));

  /* Free... */
  free(iasi_rad);
  free(pert_4mu);
  free(pert_15mu_low);
  free(pert_15mu_high);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void addatt(
  int ncid,
  int varid,
  const char *unit,
  const char *long_name) {

  /* Set long name... */
  NC(nc_put_att_text(ncid, varid, "long_name", strlen(long_name), long_name));

  /* Set units... */
  NC(nc_put_att_text(ncid, varid, "units", strlen(unit), unit));
}
