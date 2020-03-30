// scalar product of an image by a vector

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
int main(int c, char *v[])
{
	// extract command line arguments
	if (c != 3)
		return fprintf(stderr, "usage:\n\t%s img vector > out\n", *v);
		//                                 0 1   2
	char *filename_img = v[1];
	char *filename_vec = v[2];

	// read input data
	int w, h, pd, n, m;
	float *x = iio_read_image_float_vec(filename_img, &w, &h, &pd);
	float *y = iio_read_image_float(filename_vec, &n, &m);
	fprintf(stderr, "w=%d h=%d pd=%d n=%d m=%d\n", w, h, pd, n, m);
	assert(n == pd);
	assert(m == 1);

	// compute output image
	float *z = malloc(w * h * sizeof*z);
	for (int j = 0; j < w*h; j++)
	{
		z[j] = 0;
		for (int i = 0; i < n; i++)
			z[j] += y[i] * x[j*pd+i];
	}

	// save and exit
	iio_write_image_float("-", z, w, h);
	return 0;
}
