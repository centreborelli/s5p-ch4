#include <string.h>
#include <stdbool.h>
#include <hdf5.h>

static bool strsuffix(const char *s, const char *suf)
{
	int len_s = strlen(s);
	int len_suf = strlen(suf);
	if (len_s < len_suf)
		return false;
	return 0 == strcmp(suf, s + (len_s - len_suf));
}

struct twostrings { char *a, b[FILENAME_MAX]; };

static herr_t find_suffix(hid_t o, const char *n, const H5O_info_t *i, void *d)
{
	struct twostrings *p = d;
	if (i->type == H5O_TYPE_DATASET && strsuffix(n, p->a))
		return strncpy(p->b, n, FILENAME_MAX),1;
	return 0;
}

static void h5get(char *filename, char *suffix)
{
	hid_t f = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

	H5O_iterate_t u = find_suffix;
	struct twostrings p;
	p.a = suffix;
	p.b[0] = 0;
	herr_t e = H5Ovisit(f, H5_INDEX_NAME, H5_ITER_NATIVE, u, &p);
	if (*p.b)
		printf("FOUND(e=%d) \"%s\"\n", (int)e, p.b);
	else
		printf("NOT FOUND(e=%d) %s\n", (int)e, suffix);
}

int main(int c, char *v[])
{
	if (c != 3)
		return fprintf(stderr, "usage:\n\t%s file.nc suffix\n", *v);
	char *filename_in = v[1];
	char *suffix = v[2];
	h5get(filename_in, suffix);
	return 0;
}
