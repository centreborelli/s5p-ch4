// this is a kind of "plyflatten" for S5P images

#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include "iio.h"

#include "pickopt.c"

static int insideP(int w, int h, int i, int j)
{
	return i >= 0 && j >= 0 && i < w && j < h;
}

int main(int c, char *v[])
{
	// extract input arguments
	bool compute_all = pick_option(&c, &v, "a", 0);
	if (c != 9)
		return fprintf(stderr, "usage:\n\t"
			"%s [-a] lon0 lonf lat0 latf w h in.nc out.npy\n", *v);
			//0      1    2    3    4    5 6 7     8
	float lon_min = atof(v[1]);
	float lon_max = atof(v[2]);
	float lat_min = atof(v[3]);
	float lat_max = atof(v[4]);
	int w = atoi(v[5]);
	int h = atoi(v[6]);
	char *filename_in = v[7];
	char *filename_out = v[8];

	// read input data
	char filename_lon[FILENAME_MAX];
	char filename_lat[FILENAME_MAX];
	snprintf(filename_lon, FILENAME_MAX, "%s,A/longitude", filename_in);
	snprintf(filename_lat, FILENAME_MAX, "%s,A/latitude", filename_in);
	char fname_radiance[FILENAME_MAX];
	snprintf(fname_radiance, FILENAME_MAX, "%s,radiance",filename_in);
	int ww[3], hh[3], pd;
	float *lon = iio_read_image_float(filename_lon, ww+0, hh+0);
	float *lat = iio_read_image_float(filename_lat, ww+1, hh+1);
	if (ww[0] != ww[1] || hh[0] != hh[1])
		return fprintf(stderr, "bad lonlat %d,%d != %d,%d\n",
				ww[0], hh[0], ww[1], hh[1]);

	// initialize output images
	float *y_cnt = malloc(4*w*h*sizeof*y_cnt);
	float *y_avg = y_cnt + 1*w*h;//malloc(w*h*sizeof*y);
	float *y_min = y_cnt + 2*w*h;//malloc(w*h*sizeof*y);
	float *y_max = y_cnt + 3*w*h;//malloc(w*h*sizeof*y);
	for (int i = 0; i < w*h; i++)
	{
		y_cnt[i] = 0;
		y_avg[i] = NAN;
		y_min[i] = INFINITY;
		y_max[i] = -INFINITY;
	}

	float *x = 0;
	if (compute_all)
	{
		char fname_radiance[FILENAME_MAX];
		snprintf(fname_radiance, FILENAME_MAX, "%s,/radiance",
								filename_in);
		x = iio_read_image_float_vec(fname_radiance, ww+2, hh+2, &pd);
		if (ww[0] != ww[2] || hh[0] != hh[2])
			return fprintf(stderr, "bad radiance %d,%d != %d,%d\n",
					ww[0], hh[0], ww[2], hh[2]);

	}

	// splat lon,lat points into output image
	for (int j = 0; j < *hh; j++)
	for (int i = 0; i < *ww; i++)
	{
		float p = lon[j**ww+i];
		float q = lat[j**ww+i];
		int P = (w-1) * (p - lon_min) / (lon_max - lon_min);
		int Q = (h-1) * (q - lat_min) / (lat_max - lat_min);
		Q = h-1-Q; // because fuck you
		//fprintf(stderr,"i,j=%d,%d p,q=%g,%g P,Q=%d,%d\n",i,j,p,q,P,Q);
		if (insideP(w, h, P, Q))
		{
			y_cnt[Q*w+P] += 1;
			if (y_cnt[Q*w+P]==1 && compute_all)
				y_avg[Q*w+P] = 0;
			if (compute_all)
			{
				long double a = 0;
				for (int l = 0; l < pd; l++)
					a += x[(j**ww+i)*pd+l];
				a /= pd;
				y_avg[Q*w+P] += a;
				y_min[Q*w+P] = fmin(y_min[Q*w+P], a);
				y_max[Q*w+P] = fmax(y_max[Q*w+P], a);
			}
		}
	}
	for (int i = 0; i < w*h; i++)
		if (y_cnt[i])
			y_avg[i] /= y_cnt[i];

	iio_write_image_float_split(filename_out, y_cnt, w, h, 1+3*compute_all);

	return 0;
}
