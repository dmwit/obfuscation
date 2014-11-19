#include "utils.h"

#include <fcntl.h>
#include <sys/time.h>
#include <unistd.h>

#include <oz/flint-addons.h>

int g_verbose;

double
current_time(void)
{
    struct timeval t;
    (void) gettimeofday(&t, NULL);
    return (double) (t.tv_sec + (double) (t.tv_usec / 1000000.0));
}

// XXX: The use of /dev/urandom is not secure; however, the supercomputer we run
// on doesn't appear to have enough entropy, and blocks for long periods of
// time.  Thus, we use /dev/urandom instead.
#ifndef RANDFILE
#  define RANDFILE "/dev/urandom"
#endif

static unsigned long
get_seed(void)
{
    int file;
    unsigned long seed;
    if ((file = open(RANDFILE, O_RDONLY)) == -1) {
        (void) fprintf(stderr, "Error opening %s\n", RANDFILE);
        return 1;
    } else {
        if (read(file, &seed, sizeof seed) == -1) {
            (void) fprintf(stderr, "Error reading from %s\n", RANDFILE);
            (void) close(file);
            return 0;
        } else {
            if (g_verbose)
                (void) fprintf(stderr, "  Seed: %lu\n", seed);
        }
    }
    if (file != -1)
        (void) close(file);
    return seed;
}

int
seed_gmp_rng(gmp_randstate_t *rng)
{
    unsigned long seed;
    seed = get_seed();
    if (seed == 0)
        return 1;
    gmp_randinit_default(*rng);
    gmp_randseed_ui(*rng, seed);
    return 0;
}

int
seed_flint_rng(flint_rand_t *rng)
{
    unsigned long seed;
    seed = get_seed();
    if (seed == 0)
        return 1;
    flint_randinit_seed(*rng, 0, 1); /* XXX: set seed to 0 for now */
    return 0;
}


int
load_mpz_scalar(const char *fname, mpz_t x)
{
    FILE *f;
    if ((f = fopen(fname, "r")) == NULL) {
        perror(fname);
        return 1;
    }
    (void) mpz_inp_raw(x, f);
    (void) fclose(f);
    return 0;
}

int
save_mpz_scalar(const char *fname, const mpz_t x)
{
    FILE *f;
    if ((f = fopen(fname, "w")) == NULL) {
        perror(fname);
        return 1;
    }
    if (mpz_out_raw(f, x) == 0) {
        (void) fclose(f);
        return 1;
    }
    (void) fclose(f);
    return 0;
}

int
load_mpz_vector(const char *fname, mpz_t *m, int len)
{
    FILE *f;
    if ((f = fopen(fname, "r")) == NULL) {
        perror(fname);
        return 1;
    }
    for (int i = 0; i < len; ++i) {
        (void) mpz_inp_raw(m[i], f);
    }
    (void) fclose(f);
    return 0;
}

int
save_mpz_vector(const char *fname, const mpz_t *v, int len)
{
    FILE *f;
    if ((f = fopen(fname, "w")) == NULL) {
        perror(fname);
        return 1;
    }
    for (int i = 0; i < len; ++i) {
        if (mpz_out_raw(f, v[i]) == 0) {
            (void) fclose(f);
            return 1;
        }
    }
    (void) fclose(f);
    return 0;
}

int
save_fmpz_mod_poly_vector(const char *fname, const fmpz_mod_poly_t *v, int len)
{
    FILE *f;
    if ((f = fopen(fname, "w")) == NULL) {
        perror(fname);
        return 1;
    }
    /* for (int i = 0; i < len; ++i) { */
        fprintf(stderr, "NOT SAVING VECTOR YET\n");
    /* } */
    (void) fclose(f);
    return 0;
}

void
mpz_genrandom(mpz_t rnd, gmp_randstate_t *rng, long nbits)
{
    mpz_t one;
    mpz_init_set_ui(one, 1 << (nbits - 1));
    mpz_urandomb(rnd, *rng, nbits);
    mpz_clear(one);
}
