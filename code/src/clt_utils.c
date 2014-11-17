#include "clt_utils.h"
#include "utils.h"

#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void
clt_mul_mats(mpz_t *result, const mpz_t *left, const mpz_t *right,
             const mpz_t q, long m, long n, long p)
{
    mpz_t *tmparray;
    double start, end;

    start = current_time();
    tmparray = (mpz_t *) malloc(sizeof(mpz_t) * m * p);
    for (int i = 0; i < m * p; ++i) {
        mpz_init(tmparray[i]);
    }
#pragma omp parallel for
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < p; ++j) {
            mpz_t tmp, sum;
            mpz_inits(tmp, sum, NULL);
            for (int k = 0; k < n; ++k) {
                mpz_mul(tmp,
                        left[k * m + (i * m + j) % m],
                        right[k + n * ((i * m + j) / m)]);
                mpz_add(sum, sum, tmp);
                mpz_mod(sum, sum, q);
            }
            mpz_set(tmparray[i * n + j], sum);
            mpz_clears(tmp, sum, NULL);
        }
    }
    for (int i = 0; i < m * p; ++i) {
        mpz_swap(result[i], tmparray[i]);
        mpz_clear(tmparray[i]);
    }
    free(tmparray);
    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, " Multiplying took: %f\n", end - start);
}

void
clt_mul_vect_by_mat(mpz_t *v, const mpz_t *m, const mpz_t q, int size,
                    mpz_t *tmparray)
{
#pragma omp parallel for
    for (int i = 0; i < size; ++i) {
        mpz_t tmp, sum;
        mpz_inits(tmp, sum, NULL);
        for (int j = 0; j < size; ++j) {
            mpz_mul(tmp, v[j], m[i * size + j]);
            mpz_add(sum, sum, tmp);
            mpz_mod(sum, sum, q);
        }
        mpz_set(tmparray[i], sum);
        mpz_clears(tmp, sum, NULL);
    }
    for (int i = 0; i < size; ++i) {
        mpz_swap(v[i], tmparray[i]);
    }
}

void
clt_mul_vect_by_vect(mpz_t out, const mpz_t *v, const mpz_t *u, const mpz_t q,
                     int size)
{
    mpz_set_ui(out, 0);
#pragma omp parallel for
    for (int i = 0; i < size; ++i) {
        mpz_t tmp;
        mpz_init(tmp);
        mpz_mul(tmp, v[i], u[i]);
#pragma omp critical
        {
            mpz_add(out, out, tmp);
            mpz_mod(out, out, q);
        }
        mpz_clears(tmp, NULL);
    }
}
