#include "_obfuscator_c.h"

#include "clt_mlm.h"
#include "clt_utils.h"
#include "utils.h"

#include <string.h>

int
write_clt_vector(const char *dir, const mpz_t *vector, long size,
                 const char *name)
{
    char *fname;
    int fnamelen;

    fnamelen = strlen(dir) + strlen(name) + 2;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return 1;
    (void) snprintf(fname, fnamelen, "%s/%s", dir, name);
    (void) save_mpz_vector(fname, vector, size);
    free(fname);
    return 0;
}

int
write_ggh_vector(const char *dir, const gghlite_enc_t *vector, long size,
                 const char *name)
{
    char *fname;
    int fnamelen;

    fnamelen = strlen(dir) + strlen(name) + 2;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return 1;
    (void) snprintf(fname, fnamelen, "%s/%s", dir, name);
    (void) save_fmpz_mod_poly_vector(fname, vector, size);
    free(fname);
    return 0;
}

int
write_clt_layer(const char *dir, int inp, long idx, const mpz_t *zero,
                const mpz_t *one, long nrows, long ncols)
{
    mpz_t tmp;
    char *fname;
    int fnamelen;

    if (idx < 0)
        return 1;
    fnamelen = strlen(dir) + sizeof idx + 7;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return 1;
    mpz_init_set_ui(tmp, inp);
    (void) snprintf(fname, fnamelen, "%s/%ld.input", dir, idx);
    (void) save_mpz_scalar(fname, tmp);
    (void) snprintf(fname, fnamelen, "%s/%ld.zero", dir, idx);
    (void) save_mpz_vector(fname, zero, nrows * ncols);
    (void) snprintf(fname, fnamelen, "%s/%ld.one", dir, idx);
    (void) save_mpz_vector(fname, one, nrows * ncols);
    mpz_set_ui(tmp, nrows);
    (void) snprintf(fname, fnamelen, "%s/%ld.nrows", dir, idx);
    (void) save_mpz_scalar(fname, tmp);
    mpz_set_ui(tmp, ncols);
    (void) snprintf(fname, fnamelen, "%s/%ld.ncols", dir, idx);
    (void) save_mpz_scalar(fname, tmp);
    free(fname);
    mpz_clear(tmp);
    return 0;
}

int
evaluate_agis(const char *dir, const char *input, long bplen)
{
    char *fname;
    int fnamelen;
    int iszero = -1;
    mpz_t *comp, *s, *t;
    mpz_t tmp, q;
    long size;
    int err = 0;
    double start, end;

    fnamelen = strlen(dir) + sizeof bplen + 7;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return 1;

    mpz_inits(tmp, q, NULL);

    // Get the size of the matrices
    (void) snprintf(fname, fnamelen, "%s/size", dir);
    (void) load_mpz_scalar(fname, tmp);
    size = mpz_get_ui(tmp);

    // Load q
    (void) snprintf(fname, fnamelen, "%s/q", dir);
    (void) load_mpz_scalar(fname, q);

    comp = (mpz_t *) malloc(sizeof(mpz_t) * size * size);
    s = (mpz_t *) malloc(sizeof(mpz_t) * size);
    t = (mpz_t *) malloc(sizeof(mpz_t) * size);
    if (!comp || !s || !t) {
        err = 1;
        goto cleanup;
    }
    for (int i = 0; i < size; ++i) {
        mpz_inits(s[i], t[i], NULL);
    }
    for (int i = 0; i < size * size; ++i) {
        mpz_init(comp[i]);
    }
    for (int layer = 0; layer < bplen; ++layer) {
        unsigned int input_idx;

        start = current_time();
        // find out the input bit for the given layer
        (void) snprintf(fname, fnamelen, "%s/%d.input", dir, layer);
        (void) load_mpz_scalar(fname, tmp);
        input_idx = mpz_get_ui(tmp);
        if (input_idx >= strlen(input)) {
            err = 1;
            break;
        }
        if (input[input_idx] != '0' && input[input_idx] != '1') {
            err = 1;
            break;
        }

        // load in appropriate matrix for the given input value
        if (input[input_idx] == '0') {
            (void) snprintf(fname, fnamelen, "%s/%d.zero", dir, layer);
        } else {
            (void) snprintf(fname, fnamelen, "%s/%d.one", dir, layer);
        }
        (void) load_mpz_vector(fname, comp, size * size);

        // for the first matrix, multiply 'comp' by 's' to get a vector
        if (layer == 0) {
            (void) snprintf(fname, fnamelen, "%s/s_enc", dir);
            (void) load_mpz_vector(fname, s, size);
        }
        clt_mul_vect_by_mat(s, (const mpz_t *) comp, q, size, t);
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, " Multiplying matrices: %f\n",
                           end - start);
    }

    if (!err) {
        start = current_time();
        (void) snprintf(fname, fnamelen, "%s/t_enc", dir);
        (void) load_mpz_vector(fname, t, size);
        clt_mul_vect_by_vect(tmp, (const mpz_t *) s, (const mpz_t *) t, q, size);
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, " Multiplying vectors: %f\n",
                           end - start);

        start = current_time();
        {
            mpz_t pzt, nu;
            mpz_inits(pzt, nu, NULL);
            (void) snprintf(fname, fnamelen, "%s/pzt", dir);
            (void) load_mpz_scalar(fname, pzt);
            (void) snprintf(fname, fnamelen, "%s/nu", dir);
            (void) load_mpz_scalar(fname, nu);
            iszero = clt_mlm_is_zero(tmp, pzt, q, mpz_get_ui(nu));
            mpz_clears(pzt, nu, NULL);
        }
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, " Zero test: %f\n", end - start);
    }
    for (int i = 0; i < size; ++i) {
        mpz_clears(s[i], t[i], NULL);
    }
    for (int i = 0; i < size * size; ++i) {
        mpz_clear(comp[i]);
    }
cleanup:
    mpz_clears(tmp, q, NULL);
    if (comp)
        free(comp);
    if (s)
        free(s);
    if (t)
        free(t);
    if (fname)
        free(fname);

    return err ? -1 : iszero;
}

int
evaluate_sz(const char *dir, const char *input, long bplen)
{
    char *fname;
    int fnamelen;
    int iszero = -1;
    mpz_t tmp, q;
    mpz_t *result = NULL;
    long nrows, ncols = -1, nrows_prev = -1;
    int err = 0;
    double start, end;

    fnamelen = strlen(dir) + sizeof bplen + 7;
    fname = (char *) malloc(sizeof(char) * fnamelen);
    if (fname == NULL)
        return 1;

    mpz_inits(tmp, q, NULL);

    // Load q
    (void) snprintf(fname, fnamelen, "%s/q", dir);
    (void) load_mpz_scalar(fname, q);

    for (int layer = 0; layer < bplen; ++layer) {
        unsigned int input_idx;
        mpz_t *left, *right;

        start = current_time();

        // determine the size of the matrix
        (void) snprintf(fname, fnamelen, "%s/%d.nrows", dir, layer);
        (void) load_mpz_scalar(fname, tmp);
        nrows = mpz_get_ui(tmp);
        (void) snprintf(fname, fnamelen, "%s/%d.ncols", dir, layer);
        (void) load_mpz_scalar(fname, tmp);
        ncols = mpz_get_ui(tmp);

        // find out the input bit for the given layer
        (void) snprintf(fname, fnamelen, "%s/%d.input", dir, layer);
        (void) load_mpz_scalar(fname, tmp);
        input_idx = mpz_get_ui(tmp);
        if (input_idx >= strlen(input)) {
            err = 1;
            break;
        }
        if (input[input_idx] != '0' && input[input_idx] != '1') {
            err = 1;
            break;
        }
        // load in appropriate matrix for the given input value
        if (input[input_idx] == '0') {
            (void) snprintf(fname, fnamelen, "%s/%d.zero", dir, layer);
        } else {
            (void) snprintf(fname, fnamelen, "%s/%d.one", dir, layer);
        }

        if (layer == 0) {
            result = (mpz_t *) malloc(sizeof(mpz_t) * nrows * ncols);
            for (int i = 0; i < nrows * ncols; ++i) {
                mpz_init(result[i]);
            }
            (void) load_mpz_vector(fname, result, nrows * ncols);
            mpz_set(tmp, result[0]);
            nrows_prev = nrows;
        } else {
            left = result;
            right = (mpz_t *) malloc(sizeof(mpz_t) * nrows * ncols);
            for (int i = 0; i < nrows * ncols; ++i) {
                mpz_init(right[i]);
            }
            (void) load_mpz_vector(fname, right, nrows * ncols);
            result = (mpz_t *) malloc(sizeof(mpz_t) * nrows_prev * ncols);
            for (int i = 0; i < nrows_prev * ncols; ++i) {
                mpz_init(result[i]);
            }
            clt_mul_mats(result, (const mpz_t *) left, (const mpz_t *) right,
                         q, nrows_prev, nrows, ncols);
            for (int i = 0; i < nrows_prev * nrows; ++i) {
                mpz_clear(left[i]);
            }
            for (int i = 0; i < nrows * ncols; ++i) {
                mpz_clear(right[i]);
            }
            free(left);
            free(right);
        }
        end = current_time();

        if (g_verbose)
            (void) fprintf(stderr, "  Multiplying matrices: %f\n",
                           end - start);
    }

    if (!err) {
        mpz_t pzt, nu;

        start = current_time();
        mpz_inits(pzt, nu, NULL);
        (void) snprintf(fname, fnamelen, "%s/pzt", dir);
        (void) load_mpz_scalar(fname, pzt);
        (void) snprintf(fname, fnamelen, "%s/nu", dir);
        (void) load_mpz_scalar(fname, nu);
        iszero = clt_mlm_is_zero(result[1], pzt, q, mpz_get_ui(nu));
        mpz_clears(pzt, nu, NULL);
        end = current_time();
        if (g_verbose)
            (void) fprintf(stderr, "  Zero test: %f\n", end - start);
    }

    for (int i = 0; i < nrows_prev * ncols; ++i) {
        mpz_clear(result[i]);
    }
    free(result);

    mpz_clears(tmp, q, NULL);

    if (fname)
        free(fname);

    return err ? -1 : iszero;
}
