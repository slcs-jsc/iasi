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
  IASI Code Collection library definitions.
*/

#include "libiasi.h"

/*****************************************************************************/

void add_var(
  int ncid,
  const char *varname,
  const char *unit,
  const char *longname,
  int type,
  int dimid[],
  int *varid,
  int ndims) {

  /* Check if variable exists... */
  if (nc_inq_varid(ncid, varname, varid) != NC_NOERR) {

    /* Define variable... */
    NC(nc_def_var(ncid, varname, type, ndims, dimid, varid));

    /* Set long name... */
    NC(nc_put_att_text
       (ncid, *varid, "long_name", strlen(longname), longname));

    /* Set units... */
    NC(nc_put_att_text(ncid, *varid, "units", strlen(unit), unit));
  }
}

/*****************************************************************************/

void background_poly_help(
  double *xx,
  double *yy,
  int n,
  int dim) {

  gsl_multifit_linear_workspace *work;
  gsl_matrix *cov, *X;
  gsl_vector *c, *x, *y;

  double chisq, xx2[WX > WY ? WX : WY], yy2[WX > WY ? WX : WY];

  size_t i, i2, n2 = 0;

  /* Check for nan... */
  for (i = 0; i < (size_t) n; i++)
    if (gsl_finite(yy[i])) {
      xx2[n2] = xx[i];
      yy2[n2] = yy[i];
      n2++;
    }
  if ((int) n2 < dim || n2 < 0.9 * n) {
    for (i = 0; i < (size_t) n; i++)
      yy[i] = GSL_NAN;
    return;
  }

  /* Allocate... */
  work = gsl_multifit_linear_alloc((size_t) n2, (size_t) dim);
  cov = gsl_matrix_alloc((size_t) dim, (size_t) dim);
  X = gsl_matrix_alloc((size_t) n2, (size_t) dim);
  c = gsl_vector_alloc((size_t) dim);
  x = gsl_vector_alloc((size_t) n2);
  y = gsl_vector_alloc((size_t) n2);

  /* Compute polynomial fit... */
  for (i = 0; i < (size_t) n2; i++) {
    gsl_vector_set(x, i, xx2[i]);
    gsl_vector_set(y, i, yy2[i]);
    for (i2 = 0; i2 < (size_t) dim; i2++)
      gsl_matrix_set(X, i, i2, pow(gsl_vector_get(x, i), (double) i2));
  }
  gsl_multifit_linear(X, y, c, cov, &chisq, work);
  for (i = 0; i < (size_t) n; i++)
    yy[i] = gsl_poly_eval(c->data, (int) dim, xx[i]);

  /* Free... */
  gsl_multifit_linear_free(work);
  gsl_matrix_free(cov);
  gsl_matrix_free(X);
  gsl_vector_free(c);
  gsl_vector_free(x);
  gsl_vector_free(y);
}

/*****************************************************************************/

void background_poly(
  wave_t *wave,
  int dim_x,
  int dim_y) {

  double help[WX], x[WX], x2[WY], y[WX], y2[WY];

  int ix, iy;

  /* Copy temperatures to background... */
  for (ix = 0; ix < wave->nx; ix++)
    for (iy = 0; iy < wave->ny; iy++) {
      wave->bg[ix][iy] = wave->temp[ix][iy];
      wave->pt[ix][iy] = 0;
    }

  /* Check parameters... */
  if (dim_x <= 0 && dim_y <= 0)
    return;

  /* Compute fit in x-direction... */
  if (dim_x > 0)
    for (iy = 0; iy < wave->ny; iy++) {
      for (ix = 0; ix <= 53; ix++) {
	x[ix] = (double) ix;
	y[ix] = wave->bg[ix][iy];
      }
      background_poly_help(x, y, 54, dim_x);
      for (ix = 0; ix <= 29; ix++)
	help[ix] = y[ix];

      for (ix = 6; ix <= 59; ix++) {
	x[ix - 6] = (double) ix;
	y[ix - 6] = wave->bg[ix][iy];
      }
      background_poly_help(x, y, 54, dim_x);
      for (ix = 30; ix <= 59; ix++)
	help[ix] = y[ix - 6];

      for (ix = 0; ix < wave->nx; ix++)
	wave->bg[ix][iy] = help[ix];
    }

  /* Compute fit in y-direction... */
  if (dim_y > 0)
    for (ix = 0; ix < wave->nx; ix++) {
      for (iy = 0; iy < wave->ny; iy++) {
	x2[iy] = (int) iy;
	y2[iy] = wave->bg[ix][iy];
      }
      background_poly_help(x2, y2, wave->ny, dim_y);
      for (iy = 0; iy < wave->ny; iy++)
	wave->bg[ix][iy] = y2[iy];
    }

  /* Recompute perturbations... */
  for (ix = 0; ix < wave->nx; ix++)
    for (iy = 0; iy < wave->ny; iy++)
      wave->pt[ix][iy] = wave->temp[ix][iy] - wave->bg[ix][iy];
}

/*****************************************************************************/

void iasi_read(
  char *filename,
  iasi_rad_t *iasi_rad) {

  const char *product_class;

  coda_product *pf;

  coda_cursor cursor;

  iasi_raw_t *iasi_raw;

  int i, j, w, tr1, tr2, tr1_lpm, tr1_rpm, tr2_lpm, tr2_rpm,
    ichan, mdr_i, num_dims = 1;

  long dim[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

  short int IDefScaleSondNbScale, IDefScaleSondNsfirst[10],
    IDefScaleSondNslast[10], IDefScaleSondScaleFactor[10];

  float sc, scaling[IASI_L1_NCHAN];

  /* Initialize CODA... */
  coda_init();

  /* Allocate... */
  ALLOC(iasi_raw, iasi_raw_t, 1);

  /* Open IASI file... */
  CODA(coda_open(filename, &pf));
  CODA(coda_get_product_class(pf, &product_class));
  CODA(coda_cursor_set_product(&cursor, pf));

  /* Get scaling parameters... */
  CODA(coda_cursor_goto_record_field_by_name(&cursor, "GIADR_ScaleFactors"));

  CODA(coda_cursor_goto_record_field_by_name
       (&cursor, "IDefScaleSondNbScale"));
  CODA(coda_cursor_read_int16(&cursor, &IDefScaleSondNbScale));
  CODA(coda_cursor_goto_parent(&cursor));

  CODA(coda_cursor_goto_record_field_by_name
       (&cursor, "IDefScaleSondNsfirst"));
  CODA(coda_cursor_read_int16_array
       (&cursor, IDefScaleSondNsfirst, coda_array_ordering_c));
  CODA(coda_cursor_goto_parent(&cursor));

  CODA(coda_cursor_goto_record_field_by_name(&cursor, "IDefScaleSondNslast"));
  CODA(coda_cursor_read_int16_array
       (&cursor, IDefScaleSondNslast, coda_array_ordering_c));
  CODA(coda_cursor_goto_parent(&cursor));

  CODA(coda_cursor_goto_record_field_by_name
       (&cursor, "IDefScaleSondScaleFactor"));
  CODA(coda_cursor_read_int16_array
       (&cursor, IDefScaleSondScaleFactor, coda_array_ordering_c));

  /* Compute scaling factors... */
  for (ichan = 0; ichan < IASI_L1_NCHAN; ichan++)
    scaling[ichan] = GSL_NAN;
  for (i = 0; i < IDefScaleSondNbScale; i++) {
    sc = (float) pow(10.0, -IDefScaleSondScaleFactor[i]);
    for (ichan = IDefScaleSondNsfirst[i] - 1;
	 ichan < IDefScaleSondNslast[i]; ichan++) {
      w = ichan - IASI_IDefNsfirst1b + 1;
      if (w >= 0 && w < IASI_L1_NCHAN)
	scaling[w] = sc;
    }
  }

  /* Get number of tracks in record... */
  CODA(coda_cursor_goto_root(&cursor));
  CODA(coda_cursor_goto_record_field_by_name(&cursor, "MDR"));
  CODA(coda_cursor_get_array_dim(&cursor, &num_dims, dim));
  iasi_raw->ntrack = dim[0];

  /* Read tracks one by one... */
  for (mdr_i = 0; mdr_i < iasi_raw->ntrack; mdr_i++) {

    /* Reset cursor position... */
    CODA(coda_cursor_goto_root(&cursor));

    /* Move cursor to radiation data... */
    CODA(coda_cursor_goto_record_field_by_name(&cursor, "MDR"));
    CODA(coda_cursor_goto_array_element_by_index(&cursor, mdr_i));
    CODA(coda_cursor_goto_record_field_by_name(&cursor, "MDR"));
    CODA(coda_cursor_goto_record_field_by_name(&cursor, "GS1cSpect"));
    CODA(coda_cursor_read_int16_array
	 (&cursor, &iasi_raw->Radiation[mdr_i][0][0][0],
	  coda_array_ordering_c));

    /* Read time... */
    CODA(coda_cursor_goto_parent(&cursor));
    CODA(coda_cursor_goto_record_field_by_name(&cursor, "OnboardUTC"));
    CODA(coda_cursor_read_double_array
	 (&cursor, &iasi_raw->Time[mdr_i][0], coda_array_ordering_c));

    /* Read coordinates... */
    CODA(coda_cursor_goto_parent(&cursor));
    CODA(coda_cursor_goto_record_field_by_name(&cursor, "GGeoSondLoc"));
    CODA(coda_cursor_read_double_array
	 (&cursor, &iasi_raw->Loc[mdr_i][0][0][0], coda_array_ordering_c));

    /* Read satellite altitude... */
    CODA(coda_cursor_goto_parent(&cursor));
    CODA(coda_cursor_goto_record_field_by_name(&cursor,
					       "EARTH_SATELLITE_DISTANCE"));
    CODA(coda_cursor_read_uint32(&cursor, &iasi_raw->Sat_z[mdr_i]));

    /* Read spectral range... */
    iasi_raw->IDefSpectDWn1b[mdr_i] = IASI_IDefSpectDWn1b / 100.0;

    CODA(coda_cursor_goto_parent(&cursor));
    CODA(coda_cursor_goto_record_field_by_name(&cursor, "IDefNsfirst1b"));
    CODA(coda_cursor_read_int32(&cursor, &iasi_raw->IDefNsfirst1b[mdr_i]));
    if (iasi_raw->IDefNsfirst1b[mdr_i] != IASI_IDefNsfirst1b)
      ERRMSG("Unexpected value for IDefNsfirst1b!");

    CODA(coda_cursor_goto_parent(&cursor));
    CODA(coda_cursor_goto_record_field_by_name(&cursor, "IDefNslast1b"));
    CODA(coda_cursor_read_int32(&cursor, &iasi_raw->IDefNslast1b[mdr_i]));
    if (iasi_raw->IDefNslast1b[mdr_i] != IASI_IDefNslast1b)
      ERRMSG("Unexpected value for IDefNslast1b!");

    /* Compute wavenumber... */
    if (mdr_i == 0)
      for (i = 0; i < IASI_L1_NCHAN; i++)
	iasi_raw->Wavenumber[i] =
	  iasi_raw->IDefSpectDWn1b[mdr_i] *
	  (float) (iasi_raw->IDefNsfirst1b[mdr_i] + i - 1);
  }

  /* Close file... */
  CODA(coda_close(pf));

  /* Finalize CODA... */
  coda_done();

  /* Set number of tracks... */
  iasi_rad->ntrack = (int) (iasi_raw->ntrack * 2);

  /* Copy wavenumbers... */
  for (ichan = 0; ichan < IASI_L1_NCHAN; ichan++)
    iasi_rad->freq[ichan] = iasi_raw->Wavenumber[ichan];

  /* Copy footprint data... */
  for (mdr_i = 0; mdr_i < iasi_raw->ntrack; mdr_i++) {
    tr1 = mdr_i * 2;
    tr2 = mdr_i * 2 + 1;
    tr1_lpm = 3;
    tr1_rpm = 0;
    tr2_lpm = 2;
    tr2_rpm = 1;

    /* Copy time (2x2 matrix has same measurement time)...  */
    for (i = 0; i < IASI_NXTRACK; i++) {
      iasi_rad->Time[tr1][i * 2] = iasi_raw->Time[mdr_i][i];
      iasi_rad->Time[tr1][i * 2 + 1] = iasi_raw->Time[mdr_i][i];
      iasi_rad->Time[tr2][i * 2] = iasi_raw->Time[mdr_i][i];
      iasi_rad->Time[tr2][i * 2 + 1] = iasi_raw->Time[mdr_i][i];
    }

    /* Copy location... */
    for (i = 0; i < IASI_NXTRACK; i++) {
      iasi_rad->Longitude[tr1][i * 2] = iasi_raw->Loc[mdr_i][i][tr1_lpm][0];
      iasi_rad->Longitude[tr1][i * 2 + 1] =
	iasi_raw->Loc[mdr_i][i][tr1_rpm][0];
      iasi_rad->Latitude[tr1][i * 2] = iasi_raw->Loc[mdr_i][i][tr1_lpm][1];
      iasi_rad->Latitude[tr1][i * 2 + 1] =
	iasi_raw->Loc[mdr_i][i][tr1_rpm][1];

      iasi_rad->Longitude[tr2][i * 2] = iasi_raw->Loc[mdr_i][i][tr2_lpm][0];
      iasi_rad->Longitude[tr2][i * 2 + 1] =
	iasi_raw->Loc[mdr_i][i][tr2_rpm][0];
      iasi_rad->Latitude[tr2][i * 2] = iasi_raw->Loc[mdr_i][i][tr2_lpm][1];
      iasi_rad->Latitude[tr2][i * 2 + 1] =
	iasi_raw->Loc[mdr_i][i][tr2_rpm][1];
    }

    /* Copy satellite location (we only have one height value)... */
    iasi_rad->Sat_lon[tr1] = iasi_rad->Longitude[tr1][28];
    iasi_rad->Sat_lat[tr1] = iasi_rad->Latitude[tr1][28];
    iasi_rad->Sat_lon[tr2] = iasi_rad->Longitude[tr2][28];
    iasi_rad->Sat_lat[tr2] = iasi_rad->Latitude[tr2][28];
    iasi_rad->Sat_z[tr1] =
      iasi_raw->Sat_z[mdr_i] / 1000.0 - wgs84(iasi_rad->Sat_lat[tr1]);
    iasi_rad->Sat_z[tr2] =
      iasi_raw->Sat_z[mdr_i] / 1000.0 - wgs84(iasi_rad->Sat_lat[tr2]);

    /* Copy radiation data... */
    for (i = 0; i < IASI_NXTRACK; i++) {
      for (ichan = 0; ichan < IASI_L1_NCHAN; ichan++) {
	sc = scaling[ichan] * 100.0f;
	iasi_rad->Rad[tr1][i * 2][ichan] =
	  iasi_raw->Radiation[mdr_i][i][tr1_lpm][ichan] * sc;
	iasi_rad->Rad[tr1][i * 2 + 1][ichan] =
	  iasi_raw->Radiation[mdr_i][i][tr1_rpm][ichan] * sc;
	iasi_rad->Rad[tr2][i * 2][ichan] =
	  iasi_raw->Radiation[mdr_i][i][tr2_lpm][ichan] * sc;
	iasi_rad->Rad[tr2][i * 2 + 1][ichan] =
	  iasi_raw->Radiation[mdr_i][i][tr2_rpm][ichan] * sc;
      }
    }
  }

  /* Check radiance data... */
  for (i = 0; i < iasi_rad->ntrack; i++)
    for (j = 0; j < L1_NXTRACK; j++)
      if (iasi_rad->Rad[i][j][6753] > iasi_rad->Rad[i][j][6757]
	  || iasi_rad->Rad[i][j][6753] < 0)
	for (ichan = 0; ichan < IASI_L1_NCHAN; ichan++)
	  iasi_rad->Rad[i][j][ichan] = GSL_NAN;

  /* Free... */
  free(iasi_raw);
}

/*****************************************************************************/

void noise(
  wave_t *wave,
  double *mu,
  double *sig) {

  int ix, ix2, iy, iy2, n = 0, okay;

  /* Init... */
  *mu = 0;
  *sig = 0;

  /* Estimate noise (Immerkaer, 1996)... */
  for (ix = 1; ix < wave->nx - 1; ix++)
    for (iy = 1; iy < wave->ny - 1; iy++) {

      /* Check data... */
      okay = 1;
      for (ix2 = ix - 1; ix2 <= ix + 1; ix2++)
	for (iy2 = iy - 1; iy2 <= iy + 1; iy2++)
	  if (!gsl_finite(wave->temp[ix2][iy2]))
	    okay = 0;
      if (!okay)
	continue;

      /* Get mean noise... */
      n++;
      *mu += wave->temp[ix][iy];
      *sig += gsl_pow_2(+4. / 6. * wave->temp[ix][iy]
			- 2. / 6. * (wave->temp[ix - 1][iy]
				     + wave->temp[ix + 1][iy]
				     + wave->temp[ix][iy - 1]
				     + wave->temp[ix][iy + 1])
			+ 1. / 6. * (wave->temp[ix - 1][iy - 1]
				     + wave->temp[ix + 1][iy - 1]
				     + wave->temp[ix - 1][iy + 1]
				     + wave->temp[ix + 1][iy + 1]));
    }

  /* Normalize... */
  *mu /= (double) n;
  *sig = sqrt(*sig / (double) n);
}

/*****************************************************************************/

void pert2wave(
  pert_t *pert,
  wave_t *wave,
  int track0,
  int track1,
  int xtrack0,
  int xtrack1) {

  double x0[3], x1[3];

  int itrack, ixtrack;

  /* Check ranges... */
  track0 = GSL_MIN(GSL_MAX(track0, 0), pert->ntrack - 1);
  track1 = GSL_MIN(GSL_MAX(track1, 0), pert->ntrack - 1);
  xtrack0 = GSL_MIN(GSL_MAX(xtrack0, 0), pert->nxtrack - 1);
  xtrack1 = GSL_MIN(GSL_MAX(xtrack1, 0), pert->nxtrack - 1);

  /* Set size... */
  wave->nx = xtrack1 - xtrack0 + 1;
  if (wave->nx > WX)
    ERRMSG("Too many across-track values!");
  wave->ny = track1 - track0 + 1;
  if (wave->ny > WY)
    ERRMSG("Too many along-track values!");

  /* Loop over footprints... */
  for (itrack = track0; itrack <= track1; itrack++)
    for (ixtrack = xtrack0; ixtrack <= xtrack1; ixtrack++) {

      /* Get distances... */
      if (itrack == track0) {
	wave->x[0] = 0;
	if (ixtrack > xtrack0) {
	  geo2cart(0, pert->lon[itrack][ixtrack - 1],
		   pert->lat[itrack][ixtrack - 1], x0);
	  geo2cart(0, pert->lon[itrack][ixtrack],
		   pert->lat[itrack][ixtrack], x1);
	  wave->x[ixtrack - xtrack0] =
	    wave->x[ixtrack - xtrack0 - 1] + DIST(x0, x1);
	}
      }
      if (ixtrack == xtrack0) {
	wave->y[0] = 0;
	if (itrack > track0) {
	  geo2cart(0, pert->lon[itrack - 1][ixtrack],
		   pert->lat[itrack - 1][ixtrack], x0);
	  geo2cart(0, pert->lon[itrack][ixtrack],
		   pert->lat[itrack][ixtrack], x1);
	  wave->y[itrack - track0] =
	    wave->y[itrack - track0 - 1] + DIST(x0, x1);
	}
      }

      /* Save geolocation... */
      wave->time = pert->time[(track0 + track1) / 2][(xtrack0 + xtrack1) / 2];
      wave->z = 0;
      wave->lon[ixtrack - xtrack0][itrack - track0] =
	pert->lon[itrack][ixtrack];
      wave->lat[ixtrack - xtrack0][itrack - track0] =
	pert->lat[itrack][ixtrack];

      /* Save temperature data... */
      wave->temp[ixtrack - xtrack0][itrack - track0]
	= pert->bt[itrack][ixtrack];
      wave->bg[ixtrack - xtrack0][itrack - track0]
	= pert->bt[itrack][ixtrack] - pert->pt[itrack][ixtrack];
      wave->pt[ixtrack - xtrack0][itrack - track0]
	= pert->pt[itrack][ixtrack];
      wave->var[ixtrack - xtrack0][itrack - track0]
	= pert->var[itrack][ixtrack];
    }
}

/*****************************************************************************/

void variance(
  wave_t *wave,
  double dh) {

  double dh2, mu, help;

  int dx, dy, ix, ix2, iy, iy2, n;

  /* Check parameters... */
  if (dh <= 0)
    return;

  /* Compute squared radius... */
  dh2 = gsl_pow_2(dh);

  /* Get sampling distances... */
  dx =
    (int) (dh / fabs(wave->x[wave->nx - 1] - wave->x[0]) * (wave->nx - 1.0) +
	   1);
  dy =
    (int) (dh / fabs(wave->y[wave->ny - 1] - wave->y[0]) * (wave->ny - 1.0) +
	   1);

  /* Loop over data points... */
  for (ix = 0; ix < wave->nx; ix++)
    for (iy = 0; iy < wave->ny; iy++) {

      /* Init... */
      mu = help = 0;
      n = 0;

      /* Get data... */
      for (ix2 = GSL_MAX(ix - dx, 0); ix2 <= GSL_MIN(ix + dx, wave->nx - 1);
	   ix2++)
	for (iy2 = GSL_MAX(iy - dy, 0); iy2 <= GSL_MIN(iy + dy, wave->ny - 1);
	     iy2++)
	  if ((gsl_pow_2(wave->x[ix] - wave->x[ix2])
	       + gsl_pow_2(wave->y[iy] - wave->y[iy2])) <= dh2)
	    if (gsl_finite(wave->pt[ix2][iy2])) {
	      mu += wave->pt[ix2][iy2];
	      help += gsl_pow_2(wave->pt[ix2][iy2]);
	      n++;
	    }

      /* Compute local variance... */
      if (n > 1)
	wave->var[ix][iy] = help / n - gsl_pow_2(mu / n);
      else
	wave->var[ix][iy] = GSL_NAN;
    }
}

/*****************************************************************************/

double wgs84(
  double lat) {

  const double a = 6378.1370, b = 6356.7523;

  double cphi, sphi;

  cphi = cos(lat * M_PI / 180.);
  sphi = sin(lat * M_PI / 180.);

  return sqrt((gsl_pow_2(a * a * cphi) + gsl_pow_2(b * b * sphi))
	      / (gsl_pow_2(a * cphi) + gsl_pow_2(b * sphi)));
}

/*****************************************************************************/

void write_l1(
  char *filename,
  iasi_l1_t *l1) {

  int dimid[10], ncid, time_id, lon_id, lat_id,
    sat_z_id, sat_lon_id, sat_lat_id, nu_id, rad_id;

  /* Open or create netCDF file... */
  printf("Write IASI Level-1 file: %s\n", filename);
  if (nc_open(filename, NC_WRITE, &ncid) != NC_NOERR) {
    NC(nc_create(filename, NC_CLOBBER, &ncid));
  } else {
    NC(nc_redef(ncid));
  }

  /* Set dimensions... */
  if (nc_inq_dimid(ncid, "L1_NTRACK", &dimid[0]) != NC_NOERR)
    NC(nc_def_dim(ncid, "L1_NTRACK", l1->ntrack, &dimid[0]));
  if (nc_inq_dimid(ncid, "L1_NXTRACK", &dimid[1]) != NC_NOERR)
    NC(nc_def_dim(ncid, "L1_NXTRACK", L1_NXTRACK, &dimid[1]));
  if (nc_inq_dimid(ncid, "L1_NCHAN", &dimid[2]) != NC_NOERR)
    NC(nc_def_dim(ncid, "L1_NCHAN", L1_NCHAN, &dimid[2]));

  /* Add variables... */
  add_var(ncid, "l1_time", "s", "time (seconds since 2000-01-01T00:00Z)",
	  NC_DOUBLE, dimid, &time_id, 2);
  add_var(ncid, "l1_lon", "deg", "longitude", NC_DOUBLE, dimid, &lon_id, 2);
  add_var(ncid, "l1_lat", "deg", "latitude", NC_DOUBLE, dimid, &lat_id, 2);
  add_var(ncid, "l1_sat_z", "km", "satellite altitude",
	  NC_DOUBLE, dimid, &sat_z_id, 1);
  add_var(ncid, "l1_sat_lon", "deg", "(estimated) satellite longitude",
	  NC_DOUBLE, dimid, &sat_lon_id, 1);
  add_var(ncid, "l1_sat_lat", "deg", "(estimated) satellite latitude",
	  NC_DOUBLE, dimid, &sat_lat_id, 1);
  add_var(ncid, "l1_nu", "cm^-1", "channel wavenumber",
	  NC_DOUBLE, &dimid[2], &nu_id, 1);
  add_var(ncid, "l1_rad", "W/(m^2 sr cm^-1)", "channel radiance",
	  NC_FLOAT, dimid, &rad_id, 3);

  /* Leave define mode... */
  NC(nc_enddef(ncid));

  /* Write data... */
  NC(nc_put_var_double(ncid, time_id, l1->time[0]));
  NC(nc_put_var_double(ncid, lon_id, l1->lon[0]));
  NC(nc_put_var_double(ncid, lat_id, l1->lat[0]));
  NC(nc_put_var_double(ncid, sat_z_id, l1->sat_z));
  NC(nc_put_var_double(ncid, sat_lon_id, l1->sat_lon));
  NC(nc_put_var_double(ncid, sat_lat_id, l1->sat_lat));
  NC(nc_put_var_double(ncid, nu_id, l1->nu));
  NC(nc_put_var_float(ncid, rad_id, &l1->rad[0][0][0]));

  /* Close file... */
  NC(nc_close(ncid));
}

/*****************************************************************************/

void write_l2(
  char *filename,
  iasi_l2_t *l2) {

  int dimid[10], ncid, time_id, z_id, lon_id, lat_id, p_id, t_id;

  /* Create netCDF file... */
  printf("Write IASI Level-2 file: %s\n", filename);
  if (nc_open(filename, NC_WRITE, &ncid) != NC_NOERR) {
    NC(nc_create(filename, NC_CLOBBER, &ncid));
  } else {
    NC(nc_redef(ncid));
  }

  /* Set dimensions... */
  if (nc_inq_dimid(ncid, "L2_NTRACK", &dimid[0]) != NC_NOERR)
    NC(nc_def_dim(ncid, "L2_NTRACK", l2->ntrack, &dimid[0]));
  if (nc_inq_dimid(ncid, "L2_NXTRACK", &dimid[1]) != NC_NOERR)
    NC(nc_def_dim(ncid, "L2_NXTRACK", L2_NXTRACK, &dimid[1]));
  if (nc_inq_dimid(ncid, "L2_NLAY", &dimid[2]) != NC_NOERR)
    NC(nc_def_dim(ncid, "L2_NLAY", L2_NLAY, &dimid[2]));

  /* Add variables... */
  add_var(ncid, "l2_time", "s", "time (seconds since 2000-01-01T00:00Z)",
	  NC_DOUBLE, dimid, &time_id, 2);
  add_var(ncid, "l2_z", "km", "altitude", NC_DOUBLE, dimid, &z_id, 3);
  add_var(ncid, "l2_lon", "deg", "longitude", NC_DOUBLE, dimid, &lon_id, 2);
  add_var(ncid, "l2_lat", "deg", "latitude", NC_DOUBLE, dimid, &lat_id, 2);
  add_var(ncid, "l2_press", "hPa", "pressure",
	  NC_DOUBLE, &dimid[2], &p_id, 1);
  add_var(ncid, "l2_temp", "K", "temperature", NC_DOUBLE, dimid, &t_id, 3);

  /* Leave define mode... */
  NC(nc_enddef(ncid));

  /* Write data... */
  NC(nc_put_var_double(ncid, time_id, l2->time[0]));
  NC(nc_put_var_double(ncid, z_id, l2->z[0][0]));
  NC(nc_put_var_double(ncid, lon_id, l2->lon[0]));
  NC(nc_put_var_double(ncid, lat_id, l2->lat[0]));
  NC(nc_put_var_double(ncid, p_id, l2->p));
  NC(nc_put_var_double(ncid, t_id, l2->t[0][0]));

  /* Close file... */
  NC(nc_close(ncid));
}
