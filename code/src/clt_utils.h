#ifndef __CLT_UTILS_H__
#define __CLT_UTILS_H__

#include <gmp.h>

void
clt_mul_mats(mpz_t *result, const mpz_t *left, const mpz_t *right,
             const mpz_t q, long m, long n, long p);

void
clt_mul_vect_by_mat(mpz_t *v, const mpz_t *m, const mpz_t q, int size,
                    mpz_t *tmparray);

void
clt_mul_vect_by_vect(mpz_t out, const mpz_t *v, const mpz_t *u, const mpz_t q,
                     int size);

#endif
