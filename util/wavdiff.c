#include <sndfile.h>
#include <stdlib.h>
#include <stdio.h>

int cmp(const void *a , const void *b)
{
	if (*(long long *)a > *(long long *)b)
		return 1;
	else if (*(long long *)a < *(long long *)b)
		return -1;
	else
		return 0;

}

void dump(long long *list, sf_count_t length)
{
	sf_count_t prev = 0;
	for (sf_count_t i = 1; i < length; i++) {
		if (list[i] > list[i - 1]) {
			printf("%lld\t%d\n", list[i - 1], i - prev);
			prev = i;
		}
	}
	printf("%lld\t%d\n\n", list[length - 1], length - prev);
}

int main(int argc, char *argv[])
{
	SF_INFO info;
	SNDFILE* wave1 = sf_open(argv[1], SFM_READ, &info);
	SNDFILE* wave2 = sf_open(argv[2], SFM_READ, &info);

	sf_count_t length = info.frames;
	short *data1 = malloc(sizeof(short) * 2 * length);
	short *data2 = malloc(sizeof(short) * 2 * length);

	sf_readf_short(wave1, data1, length);
	sf_readf_short(wave2, data2, length);

	sf_close(wave1);
	sf_close(wave2);

	long long *diffl = malloc(sizeof(long long) * length);
	long long *diffr = malloc(sizeof(long long) * length);
	long long sigmal = 0;
	long long sigmar = 0;
	for (sf_count_t i = 0; i < length; i++) {
		diffl[i] = data1[2 * i] - data2[2 * i];
		diffr[i] = data1[2 * i + 1] - data2[2 * i + 1];
		sigmal += diffl[i];
		sigmar += diffr[i];
	}

	free(data1);
	free(data2);

	for (sf_count_t i = 0; i < length; i++) {
		diffl[i] *= diffl[i];
		diffr[i] = diffr[i] * diffr[i];
	}

	qsort(diffl, length, sizeof(long long), cmp);
	qsort(diffr, length, sizeof(long long), cmp);

	long long suml = 0;
	long long sumr = 0;
	for (sf_count_t i = 0; i < length; i++) {
		suml += diffl[i];
		sumr += diffr[i];
	}

	printf("diff.left = %#.10g - %#.10g\n",
	       (double)suml / length,
	       (double)sigmal / length * sigmal / length);
	dump(diffl, length);
	printf("diff.right = %#.10g - %#.10g\n",
	       (double)sumr / length,
	       (double)sigmar / length * sigmar / length);
	dump(diffr, length);

	free(diffl);
	free(diffr);
}
