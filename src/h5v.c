// a fancier version of h5ls

#include <stdbool.h>
#include <string.h>
#include <hdf5/serial/hdf5.h>

//struct twostrings { char *a, *b; };
//
//static bool string_suffix(const char *s, const char *suf)
//{
//	int len_s = strlen(s);
//	int len_suf = strlen(suf);
//	if (len_s < len_suf)
//		return false;
//	return 0 == strcmp(suf, s + (len_s - len_suf));
//}
//
//static herr_t find_suffix(hid_t o, const char *n, const H5O_info_t *i, void *d)
//{
//	struct twostrings *p = d;
//	if (i->type == H5O_TYPE_DATASET && string_suffix(n, p->a))
//		return strncpy(p->b, n, FILENAME_MAX),1;
//	return 0;
//}
//
//static hid_t my_hd5open(hid_t f, char *suffix)
//{
//	char dset[FILENAME_MAX] = {0};
//
//	H5O_iterate_t u = find_suffix;
//	struct twostrings p = {suffix, dset};
//	p.a = suffix;
//	herr_t e = H5Ovisit(f, H5_INDEX_NAME, H5_ITER_NATIVE, u, &p);
//	if (*dset) fprintf(stderr, "HDF5_DSET = /%s\n", dset);
//	return *dset ? H5Dopen2(f, dset, H5P_DEFAULT) : -1;
//}

// callback function for the iterator
static herr_t print_shit(hid_t o, const char *n, const H5O_info_t *i, void *no)
{
	char dset[FILENAME_MAX];
	snprintf(dset, FILENAME_MAX, "/%s", n);

	if (i->type != H5O_TYPE_DATASET)
		return 0;
	printf("%ld %ld %s\n", (long)i->fileno, (long)i->addr, dset);

	hid_t fff = 0;//*(hid_t*)no;
	printf("\tfff = %ld\n", fff);
	hid_t d = H5Dopen2(fff, dset, H5P_DEFAULT);
	printf("\td = %d\n", (int)d);
	if (d < 0) return 0;

	hid_t        t = H5Dget_type(d);
	H5T_class_t  c = H5Tget_class(t);
	H5T_sign_t sgn = H5Tget_sign(t); // fml
	size_t Bps     = H5Tget_size(t);
	printf("\tt = %d\n", (int)t);
	printf("\tc = %d\n", (int)c);
	printf("\tsgn = %d\n", (int)sgn);
	printf("\tBps = %zu\n", Bps);
	return 0;
}

// print the names of all datasets
static void h5ls(char *filename)
{
	hid_t f = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	printf("f = %d\n", (int)f);

	H5O_iterate_t u = print_shit;
	herr_t e = H5Ovisit(f, H5_INDEX_NAME, H5_ITER_NATIVE, u, &f);
}

int main(int c, char *v[])
{
	if (c != 2)
		return fprintf(stderr, "usage:\n\t%s file.nc\n", *v);
	char *filename_in = v[1];
	h5ls(filename_in);
	return 0;
}
