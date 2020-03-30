// match l2 data and band 8 spectra
//
// given an image "L2-CH4" and its corresponding "L1-B8" image, produce a
// rectangular matrix of size npoints * (1 + nfreqs),
// where "npoints" is the number of good data points found on the L2 product
// and "nfreqs" is the number of user-given wavelengths specified on the
// command line.  The first two columns of the matrix contain the L2 methane
// mixing ratio reading~~s (both raw and "bias-corrected")~~, and the rest of
// the columns the selected part of the spectrum (interpolated linearly between
// contiguous wavelength samples

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define DISABLE_LOCALIZE_MAIN
#include "localize.c"

#include "pickopt.c"
#include "parsenumbers.c"
#include "iio.h"
int main(int c, char *v[])
{
	// extract input arguments
	if (c != 5)
		return fprintf(stderr, "usage:\n\t"
				"%s L2 L1 \"freqs\" out.npy\n", *v);
				//0 1  2    3       4
	char *filename_L2  = v[1];
	char *filename_L1  = v[2];
	char *sfreqs       = v[3];
	char *filename_out = v[4];

	// read L1 product and its localization data
	struct projection_data P[1] = {{
		-180, 180, -90, 90, NAN, NAN, 100, 50, 1000, 0
	}};
	projection_data_init(P, filename_L1);
	load_radiance_data(P, filename_L1);

	// read L2 product (mixing ratio and quality)
	int ww, hh, pd;
	float *mmr = read_subdataset(filename_L2, "mixing_ratio", &ww,&hh,&pd);
	assert(ww == P->nx); assert(hh == P->ny); assert(pd == 1);
	float *qav = read_subdataset(filename_L2, "qa_value", &ww,&hh,&pd);
	assert(ww == P->nx); assert(hh == P->ny); assert(pd == 1);

	// parse selected frequencies
	int nfreqs;
	double *freqs = alloc_read_some_doubles_from_string(sfreqs, &nfreqs);
	fprintf(stderr, "nfreqs = %d\n", nfreqs);
	for (int i = 0; i < nfreqs; i++)
		fprintf(stderr, "freq[%d] = %g (%g)\n",
				i, freqs[i], wlpos_from_wl(P, freqs[i]));

	// count number of valid points
	int npoints = 0;
	for (int i = 0; i < ww*hh; i++)
		npoints += qav[i] > 60;
	fprintf(stderr, "got %d valid methane points\n", npoints);

	// alloc output array
	int ow = nfreqs + 1;
	float *out = malloc(npoints * ow * sizeof*out);

	// fill-in output array
	int idx = 0;
	for (int j = 0; j < hh; j++)
	for (int i = 0; i < ww; i++)
	{
		int ij = j*ww + i;
		if (qav[ij] <= 60) continue;
		out[ow*idx+0] = mmr[ij];
		for (int l = 0; l < nfreqs; l++)
			out[ow*idx+1+l] = get_radiance_at_xy(P, i, j, freqs[l]);
		idx += 1;
	}
	assert(idx == npoints);
	fprintf(stderr, "got here without much ado\n");

	// save output array
	iio_write_image_float(filename_out, out, ow, npoints);

	return 0;
}
