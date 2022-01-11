#include <netcdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_spline.h>
#include "coda.h"
#include "jurassic.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Number of IASI radiance channels (don't change). */
#define L1_NCHAN 33

/*! Maximum along-track size of IASI radiance granule (don't change). */
#define L1_NTRACK 1800

/*! Across-track size of IASI radiance granule (don't change). */
#define L1_NXTRACK 60

/*! Number of AIRS pressure layers (don't change). */
#define L2_NLAY 27

/*! Maximum along-track size of AIRS retrieval granule (don't change). */
#define L2_NTRACK 1800

/*! Across-track size of AIRS retrieval granule (don't change). */
#define L2_NXTRACK 60

/*! Number of channels of IASI radiance granule. */
#define IASI_L1_NCHAN 8700

/*! Raw data across-track size of IASI radiance granule. */
#define IASI_NXTRACK 30

/*! Raw data size of measurement matrix (2x2). */
#define IASI_PM 4

/*! Expected value for the computation of the first wavenumber. */
#define IASI_IDefNsfirst1b 2581

/*! Expected value for the computation of the last wavenumber.  */
#define IASI_IDefNslast1b 11041

/*! Expected value for the interval of the IASI wavenumbers [m^-1]. */
#define IASI_IDefSpectDWn1b 25

/*! Along-track size of perturbation data. */
#define PERT_NTRACK 132000

/*! Across-track size of perturbation data. */
#define PERT_NXTRACK 360

/*! Across-track size of wave analysis data. */
#define WX 300

/*! Along-track size of wave analysis data. */
#define WY 33000

/* ------------------------------------------------------------
   Macros...
   ------------------------------------------------------------ */

/*! Execute CODA library command and check result. */
#define CODA(cmd) {					\
    int coda_result=(cmd);				\
    if(coda_result!=0)					\
      ERRMSG("%s", coda_errno_to_string(coda_errno));	\
  }

/*! Execute netCDF library command and check result. */
#define NC(cmd) {				     \
  int nc_result=(cmd);				     \
  if(nc_result!=NC_NOERR)			     \
    ERRMSG("%s", nc_strerror(nc_result));	     \
}

/* ------------------------------------------------------------
   Structs...
   ------------------------------------------------------------ */

/*! IASI Level-1 data. */
typedef struct {

  /*! Number of along-track values. */
  size_t ntrack;

  /*! Time (seconds since 2000-01-01T00:00Z). */
  double time[L1_NTRACK][L1_NXTRACK];

  /*! Footprint longitude [deg]. */
  double lon[L1_NTRACK][L1_NXTRACK];

  /*! Footprint latitude [deg]. */
  double lat[L1_NTRACK][L1_NXTRACK];

  /*! Satellite altitude [km]. */
  double sat_z[L1_NTRACK];

  /*! Satellite longitude [deg]. */
  double sat_lon[L1_NTRACK];

  /*! Satellite latitude [deg]. */
  double sat_lat[L1_NTRACK];

  /*! Channel frequencies [cm^-1]. */
  double nu[L1_NCHAN];

  /*! Radiance [W/(m^2 sr cm^-1)]. */
  float rad[L1_NTRACK][L1_NXTRACK][L1_NCHAN];

} iasi_l1_t;

/*! IASI Level-2 data. */
typedef struct {

  /*! Number of along-track values. */
  size_t ntrack;

  /*! Time (seconds since 2000-01-01T00:00Z). */
  double time[L2_NTRACK][L2_NXTRACK];

  /*! Geopotential height [km]. */
  double z[L2_NTRACK][L2_NXTRACK][L2_NLAY];

  /*! Longitude [deg]. */
  double lon[L2_NTRACK][L2_NXTRACK];

  /*! Latitude [deg]. */
  double lat[L2_NTRACK][L2_NXTRACK];

  /*! Pressure [hPa]. */
  double p[L2_NLAY];

  /*! Temperature [K]. */
  double t[L2_NTRACK][L2_NXTRACK][L2_NLAY];

} iasi_l2_t;

/*! Perturbation data. */
typedef struct {

  /*! Number of along-track values. */
  int ntrack;

  /*! Number of across-track values. */
  int nxtrack;

  /*! Time (seconds since 2000-01-01T00:00Z). */
  double time[PERT_NTRACK][PERT_NXTRACK];

  /*! Longitude [deg]. */
  double lon[PERT_NTRACK][PERT_NXTRACK];

  /*! Latitude [deg]. */
  double lat[PERT_NTRACK][PERT_NXTRACK];

  /*! Brightness temperature (8 micron) [K]. */
  double dc[PERT_NTRACK][PERT_NXTRACK];

  /*! Brightness temperature (4 or 15 micron) [K]. */
  double bt[PERT_NTRACK][PERT_NXTRACK];

  /*! Brightness temperature perturbation (4 or 15 micron) [K]. */
  double pt[PERT_NTRACK][PERT_NXTRACK];

  /*! Brightness temperature variance (4 or 15 micron) [K]. */
  double var[PERT_NTRACK][PERT_NXTRACK];

} pert_t;

/*! IASI raw Level-1 data. */
typedef struct {

  /*! Number of along-track samples. */
  long ntrack;

  /*! Constants for radiation spectrum (must be equal to the expected constant).  */
  float IDefSpectDWn1b[L1_NTRACK];

  /*! Constants for radiation spectrum (must be equal to the expected constant).  */
  int32_t IDefNsfirst1b[L1_NTRACK];

  /*! Constants for radiation spectrum (must be equal to the expected constant).  */
  int32_t IDefNslast1b[L1_NTRACK];

  /*! Time (seconds since 2000-01-01T00:00Z). */
  double Time[L1_NTRACK][IASI_NXTRACK];

  /*! Location of the sounder pixel (long,lat). */
  double Loc[L1_NTRACK][IASI_NXTRACK][IASI_PM][2];

  /*! Wavenumbers are computed with the expected values. */
  float Wavenumber[IASI_L1_NCHAN];

  /*! Radiance [W/(m^2 sr m^-1)]. */
  short int Radiation[L1_NTRACK][IASI_NXTRACK][IASI_PM][IASI_L1_NCHAN];

  /*! Satellite altitude [m]. */
  unsigned int Sat_z[L1_NTRACK];

} iasi_raw_t;

/*! IASI converted Level-1 radiation data. */
typedef struct {

  /*! Number of along-track samples. */
  int ntrack;

  /*! channel wavenumber [cm^-1] */
  double freq[IASI_L1_NCHAN];

  /*! Seconds since 2000-01-01 for each sounder pixel. */
  double Time[L1_NTRACK][L1_NXTRACK];

  /*! Longitude of the sounder pixel. */
  double Longitude[L1_NTRACK][L1_NXTRACK];

  /*! Latitude of the sounder pixel. */
  double Latitude[L1_NTRACK][L1_NXTRACK];

  /*! Radiance [W/(m^2 sr cm^-1)]. */
  float Rad[L1_NTRACK][L1_NXTRACK][IASI_L1_NCHAN];

  /*! Altitude of the satellite. */
  double Sat_z[L1_NTRACK];

  /*! Estimated longitude of the satellite. */
  double Sat_lon[L1_NTRACK];

  /*! Estimated latitude of the satellite. */
  double Sat_lat[L1_NTRACK];

} iasi_rad_t;

/*! Wave analysis data. */
typedef struct {

  /*! Number of across-track values. */
  int nx;

  /*! Number of along-track values. */
  int ny;

  /*! Time (seconds since 2000-01-01T00:00Z). */
  double time;

  /*! Altitude [km]. */
  double z;

  /*! Longitude [deg]. */
  double lon[WX][WY];

  /*! Latitude [deg]. */
  double lat[WX][WY];

  /*! Across-track distance [km]. */
  double x[WX];

  /*! Along-track distance [km]. */
  double y[WY];

  /*! Temperature [K]. */
  double temp[WX][WY];

  /*! Background [K]. */
  double bg[WX][WY];

  /*! Perturbation [K]. */
  double pt[WX][WY];

  /*! Variance [K]. */
  double var[WX][WY];

} wave_t;

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Add variable to netCDF file. */
void add_var(
  int ncid,
  const char *varname,
  const char *unit,
  const char *longname,
  int type,
  int dimid[],
  int *varid,
  int ndims);

/*! Get background based on polynomial fits. */
void background_poly(
  wave_t * wave,
  int dim_x,
  int dim_y);

/*! Get background based on polynomial fits. */
void background_poly_help(
  double *xx,
  double *yy,
  int n,
  int dim);

/*! Get closest channel for a wavenumber [cm] (uses expected min wavenumber).  */
int get_chan_for_wavenumber(
  float wavenumber);

/*! Read IASI Level-1 data and convert to radiation type. */
void iasi_read(
  char *filename,
  iasi_rad_t * iasi_rad);

/*! Estimate noise. */
void noise(
  wave_t * wave,
  double *mu,
  double *sig);

/*! Convert radiance perturbation data to wave analysis struct. */
void pert2wave(
  pert_t * pert,
  wave_t * wave,
  int track0,
  int track1,
  int xtrack0,
  int xtrack1);

/*! Compute local variance. */
void variance(
  wave_t * wave,
  double dh);

/*! Calculate Earth radius according to WGS-84 reference ellipsoid. */
double wgs84(
  double lat);

/*! Write IASI Level-1 data. */
void write_l1(
  char *filename,
  iasi_l1_t * l1);

/*! Write IASI Level-2 data. */
void write_l2(
  char *filename,
  iasi_l2_t * l2);
