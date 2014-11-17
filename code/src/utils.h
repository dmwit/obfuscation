#ifndef __UTILS_H__
#define __UTILS_H__

// stdio.h needed for gmp to define prototypes for mpz_inp_raw and mpz_out_raw
#include <stdio.h>
#include <gmp.h>

extern int g_verbose;

double
current_time(void);

int
seed_rng(gmp_randstate_t *rng);

int
load_mpz_scalar(const char *fname, mpz_t x);

int
save_mpz_scalar(const char *fname, const mpz_t x);

int
load_mpz_vector(const char *fname, mpz_t *m, const int len);

int
save_mpz_vector(const char *fname, const mpz_t *m, const int len);

void
mpz_genrandom(mpz_t rnd, gmp_randstate_t *rng, const long nbits);

#endif
