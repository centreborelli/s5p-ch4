// build the inverse of the localization function, i.e. the projection function

// Each S5P product is annotated by its LOCALIZATION function, described by
// three two-dimensional images with the following names:
//
//	/.../GEODATA/latitude
//	/.../GEODATA/longitude
//	/.../INSTRUMENT/nominal_wavelength
//
// The data itself is stored on the image ".../radiance".
//
// The data is a 3D array whose (x,y,z) dimensions correspond rhoughly to
// longitude, latitude and wavelength.  This is just an indication of the
// overall orientation and is not a linear mapping by any means.  The correct
// mapping is defined by the localization function.
//
// The localization function answers the following question.  Given a sample at
// position (x,y,z) in the 3D array, to what (longitude,latitude,wavelength)
// does it correspond?  The answer is:
//
//	(p,q,r) = ( latitude[x,y] , longitude[x,y] , nominal_wavelength[x,z] )
//
// Then, evaluating the radiance 3D array at the position [p,q,r] we obtain the
// value of the radiance.
//
// While the localization function gives a complete description for the
// geometric meaning of the data, sometimes you wannt to access this same
// information in the opposite direction.  This is the goal of the projection
// function, which is the inverse of the localization function.
//
// The PROJECTION function answers the following question.  Given a geographic
// coordinate (longitude, latitude) and a target wavelength, what is the
// observed radiance?  The answer is obtained by the following procedure:
//
//	x = Px [ longitude, latitude ]
//	y = Py [ longitude, latitude ]
//	f = F  [ x, wavelength ]
//
// Notice that the three 2D images "Px, Py, F" are not directly provided with
// the product metadata.  However, they can be easily computed by re-sampling
// the localization function.  The present file provides a simple C
// implenentation of this resampling.
//
struct projection_data {
	// spatio-frequential grid
	float lon_min, lon_max;    // longitude bounds
	float lat_min, lat_max;    // latitude bounds
	float freq_min, freq_max;  // wavelength bounds
	int w;                     // number of longitude samples
	int h;                     // number of latitude samples
	int nf;                    // number of wavelenght samples

	// further geometrical context
	int nx;                    // pixel width of the spectrometer
	int ny;                    // number of lines in the 3D image
	int nz;                    // number of spectral samples in the sensor
	float *Lon;                // localization (longitudes) [nx*ny]
	float *Lat;                // localization (latitudes) [nx*ny]
	float *Wav;                // localization (wavelengths) [nz*nx]

	// geometric projection
	float *Px;                 // image of size w*h with x-projections
	float *Py;                 // image of size w*h with y-projections
	float *P;                  // vector field of projections

	// frequential projection
	float *F;                  // image of size nf*nx with f-projections

	// radiance data (not always available)
	float *rad;
};


// convert geodesic coordinates (longitude, latitude) into a raster position
static void raster_from_geo(float pq[2], struct projection_data *P,
		float lon, float lat)
{
	pq[0] = (P->w - 1) * (lon - P->lon_min) / (P->lon_max - P->lon_min);
	pq[1] = (P->h - 1) * (lat - P->lat_min) / (P->lat_max - P->lat_min);
	pq[1] = P->h - 1 - pq[1];
}

// index of the given wavelength in the table of wavelengths
static float wlpos_from_wl(struct projection_data *P, float w)
{
	return (P->nf - 1) * (w - P->freq_min) / (P->freq_max - P->freq_min);
}


// auxiliary function to evaluate an image anywhere (constant extrapolation)
static float getsamplec(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return x[i+j*w];
}

// evaluate the localization function:
// from image int pixels to orthoimage float pixels
static void Lpixpix(float out[2], struct projection_data *P, int i, int j)
{
	// extract localization values in angular units
	float lon = getsamplec(P->Lon, P->nx, P->ny, i, j);
	float lat = getsamplec(P->Lat, P->nx, P->ny, i, j);

	// scale the localization values to the current grid
	raster_from_geo(out, P, lon, lat);
}

static void Wpixpix(float out[2], struct projection_data *P, int i, int j)
{
	float w = getsamplec(P->Wav, P->nz, P->nx, i, j);
	float wlpos = wlpos_from_wl(P, w);
	out[0] = wlpos;
	out[1] = j;
}


// auxiliary function with the formula for bilinear interpolation
static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	float r = 0;
	r += a * (1-x) * (1-y);
	r += b * ( x ) * (1-y);
	r += c * (1-x) * ( y );
	r += d * ( x ) * ( y );
	return r;
}

// auxiliary function to evaluate an image using bilinear interpolation
static float bilinear_interpolation_at(float *x, int w, int h, float p, float q)
{
	int ip = p;// - 0.5;
	int iq = q;// - 0.5;
	float a = getsamplec(x, w, h, ip  , iq  );
	float b = getsamplec(x, w, h, ip+1, iq  );
	float c = getsamplec(x, w, h, ip  , iq+1);
	float d = getsamplec(x, w, h, ip+1, iq+1);
	float r = evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
	return r;
}

//static float trilinear_interpolation_at(float *f, int w, int h, int d,
//		float x, float y, float z)
//{
//	int i = floor(x);
//	int j = floor(y);
//	int k = floor(z);
//	double X = x - i;
//	double Y = y - j;
//	double Z = z - k;
//	double 
//}

// evaluate the localization function using bilinear interpolation
// from image float pixels to orthoimage float pixels
static void Lpixpix_bil(float out[2], struct projection_data *P,
		float p, float q)
{
	float lon = bilinear_interpolation_at(P->Lon, P->nx, P->ny, p, q);
	float lat = bilinear_interpolation_at(P->Lat, P->nx, P->ny, p, q);
	raster_from_geo(out, P, lon, lat);
	// note: the bilinear interpolation is linear along each coordinate,
	// thus it commutes with the rescaling.  Yet, for other interpolators
	// we really should perform this computation in the opposite order
}

static void Wpixpix_bil(float out[2], struct projection_data *P, float p, float q)
{
	float w = bilinear_interpolation_at(P->Wav, P->nz, P->nx, p, q);
	float wlpos = wlpos_from_wl(P, w);
	out[0] = wlpos;
	out[1] = q;
}




#include <math.h>    // INFINITY, hypot
#include <stdlib.h>  // malloc
#include <stdio.h>   // fprintf












//static void find_centralest_pixel(int c[2], struct projection_data *P)
//{
//	// central raster pixel
//	float C[2] = { (P->w - 1) / 2.0 , (P->h - 1) / 2.0 };
//
//	// find c that projects closest to C
//	float min_distance_so_far = INFINITY;
//	for (int j = 0; j < P->ny; j++)
//	for (int i = 0; i < P->nx; i++)
//	{
//		float p[2];
//		Lpixpix(p, P, i, j);
//		if (p[0] < 0 || p[1] < 0 || p[0] >= P->w || p[1] >= P->h)
//			continue;
//		if (hypot(p[0] - C[0], p[1] - C[1]) < min_distance_so_far)
//		{
//			min_distance_so_far = hypot(p[0] - C[0], p[1] - C[1]);
//			c[0] = i;
//			c[1] = j;
//		}
//	}
//	// TODO: find and store a pixel-domain bounding box to ease further
//	// extractions
//}

#include "iio.h"
#include "fill_bill.c"

static int insideP(int w, int h, int i, int j)
{
	return i >= 0 && j >= 0 && i < w && j < h;
}


static void splat_pixels_into_raster(struct projection_data *P)
{
	// initialize background
	for (int i = 0; i < 2 * P->w * P->h; i++)
		P->P[i] = NAN;

	// splat each sample (has a certain bias when samples overlap)
	for (int j = 0; j < P->ny; j++)
	for (int i = 0; i < P->nx; i++)
	{
		float p[2];
		Lpixpix(p, P, i, j);
		int ip = p[0];
		int iq = p[1];
		if (!insideP(P->w, P->h, ip, iq))
			continue;
		P->P[2*(ip + iq * P->w) + 0] = i;
		P->P[2*(ip + iq * P->w) + 1] = j;
	}

	iio_write_image_float_vec("/tmp/dbg_llsplat.npy", P->P, P->w, P->h, 2);
}

static void interpolate_splatted_raster(struct projection_data *P)
{
	float *tmp = malloc(2 * P->w * P->h * sizeof*tmp);
	fill_bill_vec(tmp, P->P, P->w, P->h, 2);
	free(P->P);
	P->P = tmp;

	iio_write_image_float_vec("/tmp/dbg_illsplat.npy", P->P, P->w, P->h, 2);
}


// compute the differential matrix of the pixel->raster localization at c
static void compute_deriv_L(float M[2][2], struct projection_data *P, int c[2])
{
	// scheme: central differences
	float pq_10[2]; // right
	float pq_m0[2]; // left
	float pq_01[2]; // up
	float pq_0m[2]; // down
	int d = 1;      // difference step
	Lpixpix(pq_10, P, c[0] + d, c[1]    );
	Lpixpix(pq_m0, P, c[0] - d, c[1]    );
	Lpixpix(pq_01, P, c[0]    , c[1] + d);
	Lpixpix(pq_0m, P, c[0]    , c[1] - d);
	int D = 2.0/d;  // scaling factor
	M[0][0] = D * (pq_10[0] - pq_m0[0]);
	M[0][1] = D * (pq_10[1] - pq_m0[1]);
	M[1][0] = D * (pq_01[0] - pq_0m[0]);
	M[1][1] = D * (pq_01[1] - pq_0m[1]);
}

// compute the differential matrix of the frequency localization at c
static void compute_deriv_W(float M[2][2], struct projection_data *P, int c[2])
{
	// scheme: central differences
	float pq_10[2]; // right
	float pq_m0[2]; // left
	float pq_01[2]; // up
	float pq_0m[2]; // down
	int d = 1;      // difference step
	Wpixpix(pq_10, P, c[0] + d, c[1]    );
	Wpixpix(pq_m0, P, c[0] - d, c[1]    );
	Wpixpix(pq_01, P, c[0]    , c[1] + d);
	Wpixpix(pq_0m, P, c[0]    , c[1] - d);
	int D = 2.0/d;  // scaling factor
	M[0][0] = D * (pq_10[0] - pq_m0[0]);
	M[0][1] = D * (pq_10[1] - pq_m0[1]);
	M[1][0] = D * (pq_01[0] - pq_0m[0]);
	M[1][1] = D * (pq_01[1] - pq_0m[1]);
}

static void invert_matrix_2x2(float invM[2][2], float M[2][2])
{
	// a b     \_   d -b
	// c d     /   -c  a
	float a = M[0][0];
	float b = M[0][1];
	float c = M[1][0];
	float d = M[1][1];
	float D = a*d - b*c;
	invM[0][0] =  d / D;
	invM[0][1] = -b / D;
	invM[1][0] = -c / D;
	invM[1][1] =  a / D;
}

static void inv_deriv_L(float A[2][2], struct projection_data *P, float c[2])
{
	float M[2][2];
	int ic[2] = { round(c[0]), round(c[1]) };
	compute_deriv_L(M, P, ic);
	invert_matrix_2x2(A, M);
}

static void inv_deriv_W(float A[2][2], struct projection_data *P, float c[2])
{
	float M[2][2];
	int ic[2] = { round(c[0]), round(c[1]) };
	compute_deriv_W(M, P, ic);
	invert_matrix_2x2(A, M);
}

static void matrix_2x2_times_vector(float y[2], float A[2][2], float x[2])
{
	y[0] = A[0][0] * x[0] + A[0][1] * x[1];
	y[1] = A[1][0] * x[0] + A[1][1] * x[1];
}


static void build_L_from_P(struct projection_data *P)
{
	// Field inversion algorithm:
	// 1. splat the raster positions of each image pixel
	// 2. interpolate the splats using naive bilinear multiscale;
	//    these interpolated splats are a first approximation of P
	// 3. run a few Newton iterations to refine P towards the inverse of L:
	//    Let L:pixels->lonlats be the localization function
	//    given a position x in the surface of the earth, we want to
	//    find a pixel x such that L(x)=y.  Then we will define x=P(y).
	//    To find y, we solve the equation L(x)=y by semi-newton method,
	//    which is given by the fixed point iterations
	//    y -> y - A * (L(y) - x)
	//    where A is the inverse of the differential of L at a base point
	//    (quasi-Newton) or at y (Newton).  Notice that if we reach a fixed
	//    point it is necessary a solution of the equation.

	// 1. initialize the projection table with splatted data
	// 2. interpolate the splats using bilinear multiscale
	splat_pixels_into_raster(P);
	interpolate_splatted_raster(P);

	// 3. run a few newton or quasi-newton iterations
	fprintf(stderr, "building spatial projection...\n");
	for (int j = 0; j < P->h; j++)
	for (int i = 0; i < P->w; i++)
	{
		// compute the derivative at y0 for a quasi-newton method
		//float y0[2] = { P->P[2*(i+j*P->w)+0], P->P[2*(i+j*P->w)+1] };
		for (int iter = 0; iter < 150; iter++)
		{
			//    y -> y - A * (L(y) - x)
			float x[2] = { i, j };
			float *y = P->P + 2 * (i + j*P->w);
			float Ly[2];
			Lpixpix_bil(Ly, P, y[0], y[1]);
			float dx[2] = { Ly[0] - x[0], Ly[1] - x[1] };
			float A[2][2];
			inv_deriv_L(A, P, y);
			float Adx[2];
			matrix_2x2_times_vector(Adx, A, dx);
			y[0] -= Adx[0];
			y[1] -= Adx[1];
			if (!isfinite(y[0]) || !isfinite(y[1])) break;
		}
		//char f[FILENAME_MAX];
		//snprintf(f, FILENAME_MAX, "/tmp/dbg_P%d.npy", iter);
		//iio_write_image_float_vec(f, P->P, P->w, P->h, 2);
		//to debug, you must invert the loops and remove the breaks
	}
	fprintf(stderr, "\t...done\n");
	iio_write_image_float_vec("/tmp/dbg_P.npy", P->P, P->w, P->h, 2);

	// visual sanity check that the inverse is correctly computed
	float *Q = malloc(2 * P->w * P->h * sizeof*Q);
	for (int j = 0; j < P->h; j++)
	for (int i = 0; i < P->w; i++)
	{
		float p = P->P[2*(j*P->w+i)+0];
		float q = P->P[2*(j*P->w+i)+1];
		if (!isfinite(p) || !isfinite(q)) {
			Q[2*(j*P->w+i)+0] = NAN;
			Q[2*(j*P->w+i)+1] = NAN;
		}
		else
			Lpixpix_bil(Q+2*(j*P->w+i),P,p,q);
	}
	iio_write_image_float_vec("/tmp/dbg_Q.npy", Q, P->w, P->h, 2);
	free(Q);

}

//#include "smapa.h"
//SMART_PARAMETER(NEWTON_W,40)

static void build_F_from_L(struct projection_data *P)
{
	// check that we have a wavelength localization function
	if (!P->Wav) return;

	// we build a P->F is supposed to be the (scaled) inverse of P->Wav

	// init to NAN
	for (int i = 0; i < P->nx * P->nf; i++)
		P->F[i] = NAN;

	// splat
	for (int j = 0; j < P->nx; j++) // pixel index
	for (int i = 0; i < P->nz; i++) // wavelength index in the sensor
	{
		float w = P->Wav[j*P->nz + i];  // wavelength
		float f = wlpos_from_wl(P, w); // frequency position
		int f_idx = round(f);
		if (f_idx < 0 || f_idx >= P->nf) continue;
		P->F[j*P->nf + f_idx] = i;
	}
	iio_write_image_float("/tmp/dbg_F.npy", P->F, P->nf, P->nx);

	// fill-in interior NANS
	int nfirst[P->nx], nlast[P->nx];
	for (int j = 0; j < P->nx; j++)
	{
		float *f = P->F + j*P->nf;
		nfirst[j] = nlast[j] = 0;
		for (int i = 0; i < P->nf; i++)
			if (isfinite(f[i]))
				break;
			else
				nfirst[j] += 1;
		for (int i = P->nf - 1; i >= 0; i--)
			if (isfinite(f[i]))
				break;
			else
				nlast[j] += 1;
	}
	float *tmp = malloc(P->nx * P->nf * sizeof*tmp);
	fill_bill_vec(tmp, P->F, P->nf, P->nx, 1);
	free(P->F);
	P->F = tmp;
	//for (int j = 0; j < P->nx; j++)
	//{
	//	float *f = P->F + j*P->nf;
	//	for (int i = 0; i < nfirst[j]; i++)
	//		f[i] = NAN;
	//	for (int i = P->nf - 1; i >= P->nf - nlast[j]; i--)
	//		f[i] =NAN;
	//}
	iio_write_image_float("/tmp/dbg_Fi.npy", P->F, P->nf, P->nx);

	fprintf(stderr, "building frequential projection...\n");
	// refine this rhough (i.e. pixellic) data using a few Newton iterations
	for (int j = 0; j < P->nx; j++)
	for (int i = 0; i < P->nf; i++)
	{
		//for (int iter = 0; iter < NEWTON_W(); iter++)
		for (int iter = 0; iter < 40; iter++)
		{
			float x[2] = { i, j };
			float y[2] = { P->F[j*P->nf + i], j };
			float Wy[2];
			Wpixpix_bil(Wy, P, y[0], y[1]);
			float dx[2] = { Wy[0] - x[0], Wy[1] - x[1] };
			float A[2][2];
			inv_deriv_W(A, P, y);
			float Adx[2];
			matrix_2x2_times_vector(Adx, A, dx);
			y[0] -= Adx[0];
			y[1] -= Adx[1];
			P->F[j*P->nf+i] = y[0];
			// y[1] is not used, presumed to be "j" or very close
			if (!isfinite(y[0]) || !isfinite(y[1])) break;
		}
	}
	fprintf(stderr, "\t...done\n");
	iio_write_image_float("/tmp/dbg_FiN.npy", P->F, P->nf, P->nx);

	// sanity check
	float *FQ = malloc(P->nx * P->nf * sizeof*FQ);
	float *FQ2 = malloc(P->nx * P->nf * sizeof*FQ);
	for (int j = 0; j < P->nx; j++)
	for (int i = 0; i < P->nf; i++)
	{
		FQ[j*P->nf + i] = FQ2[j*P->nf + i] = NAN;
		float f = P->F[j*P->nf + i];
		if (!isfinite(f)) continue;
		int w_idx = round(f);
		if (w_idx < 0 || w_idx >= P->nz) continue;
		float tmp[2];
		Wpixpix_bil(tmp, P, f, j);
		FQ[j*P->nf+i] = tmp[0];
		FQ2[j*P->nf+i] = tmp[1];
		//float w = P->Wav[j*P->nz + w_idx];
		//FQ[j*P->nf + i] = wlpos_from_wl(P, w);
	}
	iio_write_image_float("/tmp/dbg_FQ.npy", FQ, P->nf, P->nx);
	iio_write_image_float("/tmp/dbg_FQ2.npy", FQ2, P->nf, P->nx);
}



#include <assert.h>


static float *read_subdataset(char *fname, char *dset, int *w, int *h, int *pd)
{
	char f[FILENAME_MAX];
	snprintf(f, FILENAME_MAX, "%s,%s", fname, dset);
	float *r = iio_read_image_float_vec(f, w, h, pd);
	if (r)
	fprintf(stderr, "subdataset %s w=%d h=%d pd=%d\n", dset, *w, *h, *pd);
	return r;
}

static void projection_data_init(struct projection_data *P, char *fname)
{
	// read the three images
	int w[3], h[3], pd[3];
	P->Lon = read_subdataset(fname, "/longitude",          w+0, h+0, pd+0);
	P->Lat = read_subdataset(fname, "/latitude",           w+1, h+1, pd+1);
	P->Wav = read_subdataset(fname, "/nominal_wavelength", w+2, h+2, pd+2);
	fprintf(stderr, "P->Lon = %p\n", (void*)P->Lon);
	fprintf(stderr, "P->Lat = %p\n", (void*)P->Lat);
	fprintf(stderr, "P->Wav = %p\n", (void*)P->Wav);

	if (!P->Lon || !P->Lat)
		exit(fprintf(stderr, "cannot find localization data\n"));

	iio_write_image_float_vec("/tmp/dbg_lon.npy", P->Lon, w[0],h[0],pd[0]);
	iio_write_image_float_vec("/tmp/dbg_lat.npy", P->Lat, w[1],h[1],pd[1]);
	if (P->Wav)
	iio_write_image_float_vec("/tmp/dbg_wav.npy", P->Wav, w[2],h[2],pd[2]);

	// copy sizes
	P->nx = w[0]; // == w[1] == h[2]
	P->ny = h[0]; // == h[1]
	P->nz = w[2];

	// check sizing consistency
	if (pd[0] != 1 || pd[1] != 1 || (P->Wav && pd[2] != 1))
		exit(fprintf(stderr, "bad pds = %d %d %d\n", *pd,pd[1],pd[2]));
	if (w[0] != w[1] || h[0] != h[1])
		exit(fprintf(stderr,"bad lls %d %d %d %d\n",*w,w[1],*h,h[1]));
	if (P->Wav && w[0] != h[2])
		exit(fprintf(stderr,"inconsistent spec %d %d\n",*w,h[2]));

	// assign frequency range (if missing)
	if (P->Wav && (!isfinite(P->freq_min) || !isfinite(P->freq_max)) )
	{
		P->freq_min = INFINITY;
		P->freq_max = -INFINITY;
		for (int i = 0; i < w[2] * h[2]; i++)
			P->freq_min = fmin(P->freq_min, P->Wav[i]);
		for (int i = 0; i < w[2] * h[2]; i++)
			P->freq_max = fmax(P->freq_max, P->Wav[i]);
		//fprintf(stderr, "set freq_min = %g\n", P->freq_min);
		//fprintf(stderr, "set freq_max = %g\n", P->freq_max);
	}
	fprintf(stderr, "freq_min = %g\n", P->freq_min);
	fprintf(stderr, "freq_max = %g\n", P->freq_max);
	fprintf(stderr, "nf       = %d\n", P->nf);

	// alloc the three images
	P->Px = 0;//malloc(P->w * P->h * sizeof(float));
	P->Py = 0;//malloc(P->w * P->h * sizeof(float));
	P->P  = malloc(2 * P->w * P->h * sizeof(float));
	P->F  = malloc(P->nx * P->nf * sizeof(float));

	// fill-in the tables
	build_L_from_P(P); // geometrical projection
	build_F_from_L(P); // spectral projection

	P->rad = NULL;
}

static void projection_data_free(struct projection_data *P)
{
	void *t[6] = {P->Lon, P->Lat, P->Wav, P->Px, P->Py, P->F};
	for (int i = 0; i < 6; i++)
		free(t[i]);
}

static bool global_debug_getradiance = false;
static float get_radiance_at_xy(struct projection_data *P,
		int x, int y, float w)
{
	if (global_debug_getradiance) {
		fprintf(stderr, "get_radiance_at %d %d %g\n", x, y, w);
		fprintf(stderr, "\tP->rad = %p\n", (void*)P->rad);
		fprintf(stderr, "\tP->Wav = %p\n", (void*)P->Wav);
	}
	if (!P->rad) return NAN;
	if (!P->Wav) return NAN;
	if (!insideP(P->nx, P->ny, x, y)) return NAN;

	float f = wlpos_from_wl(P, w);
	int f_idx = floor(f);
	if (global_debug_getradiance)
		fprintf(stderr, "\tf=%g f_idx=%d\n", f, f_idx);
	if (f_idx < 0 || f_idx >= P->nf) return NAN;

	float eff_0 =                            P->F[x*P->nf + f_idx+0];
	float eff_1 = f_idx == P->nf-1 ? eff_0 : P->F[x*P->nf + f_idx+1];
	if (global_debug_getradiance)
		fprintf(stderr, "\teff_0=%g eff_1=%g\n", eff_0, eff_1);

	float beta = f - f_idx;
	assert(beta >= 0);
	assert(beta < 1);
	float z = (1-beta) * eff_0 + beta * eff_1;
	if (global_debug_getradiance)
		fprintf(stderr, "\tbeta=%g z=%g\n", beta, z);

	if (!isfinite(z)) return NAN;
	int iz = floor(z);
	if (iz < 0 || iz >= P->nz) return NAN;

	float rad_z0 =                          P->rad[(y*P->nx+x)*P->nz+iz+0];
	float rad_z1 = iz == P->nz-1 ? rad_z0 : P->rad[(y*P->nx+x)*P->nz+iz+1];
	float alpha = z - iz;
	assert(alpha >= 0);
	assert(alpha < 1);
	if (global_debug_getradiance)
		fprintf(stderr, "\tiz=%d rad_z0=%g rad_z1=%g alpha=%g\n",
				iz, rad_z0, rad_z1, alpha);

	return (1-alpha) * rad_z0 + alpha * rad_z1;
}

static float get_radiance_at_lonlat(struct projection_data *P,
		float lon, float lat, float w)
{
	float pq[2];
	raster_from_geo(pq, P, lon, lat);
	if (!isfinite(pq[0]) || !isfinite(pq[1])) return NAN; // never reached

	int ip = round(pq[0]);
	int iq = round(pq[1]);
	if (!insideP(P->w, P->h, ip, iq)) return NAN;

	// TODO: better  interpolate P->P bilinearly instead of
	// through quantized (ip,iq)
	//
	// NOTE: this is irrespective of wether we interpolate
	// the radiances by nearest neighbor or some other
	// interpolater
	float x = P->P[2*(P->w*iq+ip)+0];
	float y = P->P[2*(P->w*iq+ip)+1];
	if (!isfinite(x) || !isfinite(y)) return NAN;

	int ix = round(x);
	int iy = round(y);
	if (!insideP(P->nx, P->ny, ix, iy)) return NAN;

	return get_radiance_at_xy(P, ix, iy, w);
}

static void load_radiance_data(struct projection_data *P, char *filename)
{
	fprintf(stderr, "loading radiance datacube...\n");
	int ww, hh, pd;
	float *img = read_subdataset(filename, "/radiance", &ww, &hh, &pd);
	fprintf(stderr, "\t...done ww=%d hh=%d pd=%d\n", ww, hh, pd);
	assert(ww == P->nx);
	assert(hh == P->ny);
	assert(pd == P->nz);
	assert(P->Wav);

	P->rad = img;
}

//
// Note about naming: the localization function is used to map a pixel into a
// geographic position, and can be used to project geo-referenced information
// into the image domain.  The projection function maps geographic positions
// into the image domain, thus it can be used to localize image-based data
// (like the radiances or the methane mixing ratios) into an orthoreferenced
// grid.  Thus, the projection function localizes and the localization function
// projects.  This is not an unfortunate convention, it is supposed to be like
// that.
//
//            a
//       A -------> R
//       |      -/
//       |      /             a : A -> R    // data on A
//       |     /              b : B -> R    // data on B
//     f |    / b             f : A -> B    // bijective mapping from A to B
//       |   /                a = b * f     // a = pull-back of b by f
//       |  /                 b = a * f^-1  // b = push-forward of a by f
//       V
//       B
//
// The pull-back is composition with a function and the push-forward is
// defined by composition with its inverse.
//
// Thus the program below, that composes an image with the projection function,
// is called "localize" and not "project".
//
#ifndef DISABLE_LOCALIZE_MAIN


#include "pickopt.c"
#include "parsenumbers.c"
int main(int c, char *v[])
{
	// extract named arguments
	char *dataset_id     = pick_option(&c, &v, "d",    "/radiance");
	char *filename_splat = pick_option(&c, &v, "s",    "");
	float freq_min  = atof(pick_option(&c, &v, "fmin", "nan"));
	float freq_max  = atof(pick_option(&c, &v, "fmax", "nan"));
	int nf          = atoi(pick_option(&c, &v, "nf",   "1000"));
	char *sfrequencies   = pick_option(&c, &v, "sf",   "");
	char *llpixels       = pick_option(&c, &v, "p",    "");

	// extract positional arguments
	if (c != 8)
		return fprintf(stderr, "usage:\n\t"
			"%s lon0 lonf lat0 latf w h in.nc "
			//0      1    2    3    4    5 6 7
			"[-d dset] [-s splatavg.npy] [-sf \"freqs\"]\n", *v);
	float lon_min = atof(v[1]);
	float lon_max = atof(v[2]);
	float lat_min = atof(v[3]);
	float lat_max = atof(v[4]);
	int w = atoi(v[5]);
	int h = atoi(v[6]);
	char *filename_in = v[7];

	// default settings for the wavelength
	//float freq_min = 200;
	//float freq_max = 2400;
	//int nf = 4000;

	// do stuff
	struct projection_data P[1] = {{
		lon_min,lon_max, lat_min,lat_max, freq_min,freq_max,
		w, h, nf, 0
	}};
	projection_data_init(P, filename_in);

	int ww, hh, pd;
	float *img = NULL;
	if (*filename_splat || *llpixels)
	{
		img = read_subdataset(filename_in, dataset_id, &ww, &hh, &pd);
		assert(ww == P->nx);
		assert(hh == P->ny);
		if (P->Wav)
			assert(pd == P->nz);
		P->rad = img;
	}

	if (*filename_splat && !*sfrequencies)
	{
		fprintf(stderr, "splatting the average datacube\n");

		float *splat = malloc(w*h*sizeof*splat);
		for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		{
			splat[j*w+i] = NAN;
			float x = P->P[2*(w*j+i)+0];
			float y = P->P[2*(w*j+i)+1];
			if (!isfinite(x) || !isfinite(y)) continue;
			int ix = round(x);
			int iy = round(y);
			if (!insideP(ww, hh, ix, iy)) continue;
			long double a = 0;
			for (int l = 0; l < pd; l++)
				a += img[(iy*ww+ix)*pd+l];
			splat[j*w+i] = a / pd;
		}

		iio_write_image_float(filename_splat, splat, w, h);
		free(splat);
	}


	if (*filename_splat && *sfrequencies)
	{
		fprintf(stderr, "splatting individual frequencies \"%s\"\n",
				sfrequencies);

		if (!P->Wav)
			exit(fprintf(stderr,"no wavelength localization\n"));

		// read frequencies
		int nsfreqs;
		double *frequencies = alloc_read_some_doubles_from_string(
				sfrequencies, &nsfreqs);
		if (nsfreqs==2 && isnan(*frequencies))
		{ // build linspace frequency range for this band
			nsfreqs = frequencies[1];
			free(frequencies);
			frequencies = malloc(nsfreqs * sizeof*frequencies);
			for (int i = 0; i < nsfreqs; i++)
				frequencies[i] = P->freq_min
				+ i*(P->freq_max - P->freq_min)/(nsfreqs-1.0);
		}
		fprintf(stderr, "nsfreqs = %d\n", nsfreqs);
		for (int i = 0; i < nsfreqs; i++)
			fprintf(stderr, "sfreq[%d] = %g (%g)\n",
					i, frequencies[i],
					wlpos_from_wl(P, frequencies[i]));

		// build un-smiled datacube
		float *out = malloc(w * h * nsfreqs * sizeof*out);
		for (int j = 0; j < h; j++)
		for (int i = 0; i < w; i++)
		for (int l = 0; l < nsfreqs; l++)
		{
			if (i==w/2 && j==h/2)
				fprintf(stderr, "ijl = %d %d %d\n", i, j, l);

			out[nsfreqs*(j*w+i)+l] = NAN;
			// TODO: write the following computation as a function
			float x = P->P[2*(w*j+i)+0];
			float y = P->P[2*(w*j+i)+1];
			if (i==w/2 && j==h/2)
				fprintf(stderr, "\txy = %g %g\n", x, y);
			if (!isfinite(x) || !isfinite(y)) continue;


			int ix = round(x);
			int iy = round(y);
			if (i==w/2 && j==h/2)
				fprintf(stderr, "\tixy = %d %d\n", ix, iy);
			if (!insideP(ww, hh, ix, iy)) continue;


			if (i==w/2 && j==h/2)
				global_debug_getradiance = true;
			float r = get_radiance_at_xy(P, ix, iy, frequencies[l]);
			global_debug_getradiance = false;
			if (i==w/2 && j==h/2)
				fprintf(stderr, "\tr = %g\n", r);
			out[nsfreqs*(j*w+i)+l] = r;
		}
		iio_write_image_float_vec(filename_splat, out, w, h, nsfreqs);
		free(out);
	}

	if (*llpixels)
	{
		fprintf(stderr, "values at ll-pixels: \"%s\"\n", llpixels);

		int npixels;
		double *llpix = alloc_read_some_doubles_from_string(
				llpixels, &npixels);
		npixels /= 2;
		fprintf(stderr, "got %d ll-pixels:\n", npixels);
		for (int i = 0; i < npixels; i++)
		{
			fprintf(stderr, "lonlat = %g %g\n",
					llpix[2*i], llpix[2*i+1]);

			float pq[2];
			raster_from_geo(pq, P, llpix[2*i], llpix[2*i+1]);
			fprintf(stderr, "\tgridpos {%g %g}\n", pq[0], pq[1]);
			if (!isfinite(pq[0]) || !isfinite(pq[1])) {
				fprintf(stderr, "ERROR: strange point!\n");
				continue;
			}

			int ip = round(pq[0]);
			int iq = round(pq[1]);
			fprintf(stderr, "\tqgridpos [%d %d]\n", ip, iq);
			if (!insideP(P->w, P->h, ip, iq)) {
				fprintf(stderr, "WARN: outside of ROI\n");
				continue;
			}

			// TODO: better  interpolate P->P bilinearly instead of
			// through quantized (ip,iq)
			//
			// NOTE: this is irrespective of wether we interpolate
			// the radiances by nearest neighbor or some other
			// interpolater
			float x = P->P[2*(w*iq+ip)+0];
			float y = P->P[2*(w*iq+ip)+1];
			fprintf(stderr, "\timgpos « %g %g »\n", x, y);
			if (!isfinite(x) || !isfinite(y)) {
				fprintf(stderr, "WARN: non-covered site\n");
				continue;
			}

			int ix = round(x);
			int iy = round(y);
			fprintf(stderr, "\tqimgpos ( %d %d )\n", ix, iy);
			if (!insideP(P->nx, P->ny, ix, iy)) {
				fprintf(stderr, "WARN: outside of img\n");
				continue;
			}

			float out[P->nz][2];
			for (int i = 0; i < P->nz; i++)
			{
				int idx_img = (iy * P->nx + ix)*P->nz + i;
				int idx_wav = ix * P->nz + i;
				out[i][0] = P->Wav[idx_wav];
				out[i][1] = img[idx_img];
			}
			for (int i = 0; i < P->nz; i++)
				fprintf(stdout, "%g %g\n", out[i][0],out[i][1]);
		}
	}

	if (img) free(img);
	projection_data_free(P);

	return 0;
}



//int main_old(int c, char *v[])
//{
//	char *dataset_id     = pick_option(&c, &v, "d",    "/radiance");
//	char *filename_splat = pick_option(&c, &v, "s",    "");
//	float freq_min  = atof(pick_option(&c, &v, "fmin", "nan"));
//	float freq_max  = atof(pick_option(&c, &v, "fmax", "nan"));
//	int nf          = atoi(pick_option(&c, &v, "nf",   "1000"));
//	char *sfrequencies   = pick_option(&c, &v, "sf",   "");
//	char *llpixels       = pick_option(&c, &v, "p",    "");
//	// extract input arguments
//	if (c != 8)
//		return fprintf(stderr, "usage:\n\t"
//			"%s lon0 lonf lat0 latf w h in.nc "
//			//0      1    2    3    4    5 6 7
//			"[-d dset] [-s splatavg.npy] [-sf \"freqs\"]\n", *v);
//	float lon_min = atof(v[1]);
//	float lon_max = atof(v[2]);
//	float lat_min = atof(v[3]);
//	float lat_max = atof(v[4]);
//	int w = atoi(v[5]);
//	int h = atoi(v[6]);
//	char *filename_in = v[7];
//
//	// default settings for the wavelength
//	//float freq_min = 200;
//	//float freq_max = 2400;
//	//int nf = 4000;
//
//	// do stuff
//	struct projection_data P[1] = {{
//		lon_min,lon_max, lat_min,lat_max, freq_min,freq_max,
//		w, h, nf, 0
//	}};
//	projection_data_init(P, filename_in);
//
//	int ww, hh, pd;
//	float *img = NULL;
//	if (*filename_splat || *llpixels)
//	{
//		char fname_radiance[FILENAME_MAX];
//		snprintf(fname_radiance, FILENAME_MAX, "%s,%s",
//						filename_in, dataset_id);
//		fprintf(stderr, "loading radiance datacube...\n");
//		img = iio_read_image_float_vec(fname_radiance, &ww, &hh, &pd);
//		fprintf(stderr, "\t...done ww=%d hh=%d pd=%d\n", ww, hh, pd);
//		assert(ww == P->nx);
//		assert(hh == P->ny);
//		if (P->Wav)
//			assert(pd == P->nz);
//	}
//
//	if (*filename_splat && !*sfrequencies)
//	{
//		fprintf(stderr, "splatting the average datacube\n");
//
//		float *splat = malloc(w*h*sizeof*splat);
//		for (int j = 0; j < h; j++)
//		for (int i = 0; i < w; i++)
//		{
//			splat[j*w+i] = NAN;
//			float x = P->P[2*(w*j+i)+0];
//			float y = P->P[2*(w*j+i)+1];
//			if (!isfinite(x) || !isfinite(y)) continue;
//			int ix = round(x);
//			int iy = round(y);
//			if (!insideP(ww, hh, ix, iy)) continue;
//			long double a = 0;
//			for (int l = 0; l < pd; l++)
//				a += img[(iy*ww+ix)*pd+l];
//			a /= pd;
//			if (pd == 1 && a > 1e30) a = NAN;
//			splat[j*w+i] = a;
//		}
//
//		iio_write_image_float(filename_splat, splat, w, h);
//		free(splat);
//	}
//
//
//	if (*filename_splat && *sfrequencies)
//	{
//		fprintf(stderr, "splatting individual frequencies \"%s\"\n",
//				sfrequencies);
//
//		if (!P->Wav)
//			exit(fprintf(stderr,"no wavelength localization\n"));
//
//
//		// read frequencies
//		int nsfreqs;
//		double *frequencies = alloc_read_some_doubles_from_string(
//				sfrequencies, &nsfreqs);
//		fprintf(stderr, "nsfreqs = %d\n", nsfreqs);
//		for (int i = 0; i < nsfreqs; i++)
//			fprintf(stderr, "sfreq[%d] = %g (%g)\n",
//					i, frequencies[i],
//					wlpos_from_wl(P, frequencies[i]));
//
//		// build un-smiled datacube
//		float *out = malloc(w * h * nsfreqs * sizeof*out);
//		for (int j = 0; j < h; j++)
//		for (int i = 0; i < w; i++)
//		for (int l = 0; l < nsfreqs; l++)
//		{
//			out[nsfreqs*(j*w+i)+l] = NAN;
//			// TODO: write the following computation as a function
//			float x = P->P[2*(w*j+i)+0];
//			float y = P->P[2*(w*j+i)+1];
//			if (!isfinite(x) || !isfinite(y)) continue;
//			int ix = round(x);
//			int iy = round(y);
//			if (!insideP(ww, hh, ix, iy)) continue;
//			float f = wlpos_from_wl(P, frequencies[l]);
//			int f_idx = floor(f);
//			if (f_idx < 0 || f_idx >= P->nf) continue;
//
//			//float z = P->F[ix*P->nf + f_idx];
//			float eff_0 = P->F[ix*P->nf + f_idx+0];
//			float eff_1 = f_idx == P->nf-1 ? eff_0 : P->F[ix*P->nf + f_idx+1];
//			float beta = f - f_idx;
//			assert(beta >= 0);
//			assert(beta < 1);
//			float z = (1-beta) * eff_0 + beta * eff_1;
//
//			if (!isfinite(z)) continue;
//			int iz = floor(z);
//			if (iz < 0 || iz >= P->nz) continue;
//			float rad_z0 =                          img[(iy*ww+ix)*pd+iz+0];
//			float rad_z1 = iz == P->nz-1 ? rad_z0 : img[(iy*ww+ix)*pd+iz+1];
//			float alpha = z - iz;
//			assert(alpha >= 0);
//			assert(alpha < 1);
//			out[nsfreqs*(j*w+i)+l] = (1-alpha) * rad_z0 + alpha * rad_z1;
//		}
//		iio_write_image_float_vec(filename_splat, out, w, h, nsfreqs);
//		free(out);
//	}
//
//	if (*llpixels)
//	{
//		fprintf(stderr, "values at ll-pixels: \"%s\"\n", llpixels);
//
//		int npixels;
//		double *llpix = alloc_read_some_doubles_from_string(
//				llpixels, &npixels);
//		npixels /= 2;
//		fprintf(stderr, "got %d ll-pixels:\n", npixels);
//		for (int i = 0; i < npixels; i++)
//		{
//			fprintf(stderr, "lonlat = %g %g\n",
//					llpix[2*i], llpix[2*i+1]);
//
//			float pq[2];
//			raster_from_geo(pq, P, llpix[2*i], llpix[2*i+1]);
//			fprintf(stderr, "\tgridpos {%g %g}\n", pq[0], pq[1]);
//			if (!isfinite(pq[0]) || !isfinite(pq[1])) {
//				fprintf(stderr, "ERROR: strange point!\n");
//				continue;
//			}
//
//			int ip = round(pq[0]);
//			int iq = round(pq[1]);
//			fprintf(stderr, "\tqgridpos [%d %d]\n", ip, iq);
//			if (!insideP(P->w, P->h, ip, iq)) {
//				fprintf(stderr, "WARN: outside of ROI\n");
//				continue;
//			}
//
//			// TODO: better  interpolate P->P bilinearly instead of
//			// through quantized (ip,iq)
//			//
//			// NOTE: this is irrespective of wether we interpolate
//			// the radiances by nearest neighbor or some other
//			// interpolater
//			float x = P->P[2*(w*iq+ip)+0];
//			float y = P->P[2*(w*iq+ip)+1];
//			fprintf(stderr, "\timgpos « %g %g »\n", x, y);
//			if (!isfinite(x) || !isfinite(y)) {
//				fprintf(stderr, "WARN: non-covered site\n");
//				continue;
//			}
//
//			int ix = round(x);
//			int iy = round(y);
//			fprintf(stderr, "\tqimgpos ( %d %d )\n", ix, iy);
//			if (!insideP(P->nx, P->ny, ix, iy)) {
//				fprintf(stderr, "WARN: outside of img\n");
//				continue;
//			}
//
//			float out[P->nz][2];
//			for (int i = 0; i < P->nz; i++)
//			{
//				int idx_img = (iy * P->nx + ix)*P->nz + i;
//				int idx_wav = ix * P->nz + i;
//				out[i][0] = P->Wav[idx_wav];
//				out[i][1] = img[idx_img];
//			}
//			for (int i = 0; i < P->nz; i++)
//				fprintf(stdout, "%g %g\n", out[i][0],out[i][1]);
//		}
//	}
//
//	if (img) free(img);
//	projection_data_free(P);
//
//	return 0;
//}
#endif//DISABLE_LOCALIZE_MAIN
