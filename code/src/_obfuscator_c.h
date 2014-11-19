#ifndef ___OBFUSCATOR_C_H__
#define ___OBFUSCATOR_C_H__

#include <gmp.h>
#include <gghlite/gghlite.h>

int
write_clt_vector(const char *dir, const mpz_t *vector, long size,
                 const char *name);

int
write_ggh_vector(const char *dir, const gghlite_enc_t *vector, long size,
                 const char *name);

int
write_clt_layer(const char *dir, int inp, long idx, const mpz_t *zero,
                const mpz_t *one, long nrows, long ncols);

int
evaluate_agis(const char *dir, const char *input, long bplen);

int
evaluate_sz(const char *dir, const char *input, long bplen);

#endif
