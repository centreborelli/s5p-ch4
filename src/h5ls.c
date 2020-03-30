#include <hdf5.h>

// callback function for the iterator
static herr_t print_shit(hid_t o, const char *n, const H5O_info_t *i, void *d)
{
	if (i->type == H5O_TYPE_DATASET)
		printf("/%s\n", n);
	return 0;
}

// print the names of all datasets
static void h5ls(char *filename)
{
	hid_t f = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

	H5O_iterate_t u = print_shit;
	herr_t e = H5Ovisit(f, H5_INDEX_NAME, H5_ITER_NATIVE, u, NULL);
}

int main(int c, char *v[])
{
	if (c != 2)
		return fprintf(stderr, "usage:\n\t%s file.nc\n", *v);
	char *filename_in = v[1];
	h5ls(filename_in);
	return 0;
}
