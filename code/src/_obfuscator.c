#include "_obfuscator_c.h"

#include "clt_utils.h"
#include "state.h"
#include "utils.h"
#include "pyutils.h"

#include <omp.h>

static int
extract_indices(PyObject *py_list, int *idx1, int *idx2)
{
    *idx1 = -1;
    *idx2 = -1;
    switch (PyList_GET_SIZE(py_list)) {
    case 2:
        *idx2 = PyLong_AsLong(PyList_GET_ITEM(py_list, 1));
        /* fallthrough */
    case 1:
        *idx1 = PyLong_AsLong(PyList_GET_ITEM(py_list, 0));
        break;
    default:
        return 1;
    }
    return 0;
}

static void
encode_vectors_clt(struct state *s, const PyObject *py_vectors,
                   const int indices[2], const int pows[2], const char *name)
{
    mpz_t *vector;
    ssize_t length;

    // We assume that all vectors have the same length, and thus just grab the
    // length of the first vector
    length = PyList_GET_SIZE(PyList_GET_ITEM(py_vectors, 0));
    vector = (mpz_t *) malloc(sizeof(mpz_t) * length);

#pragma omp parallel for
    for (ssize_t i = 0; i < length; ++i) {
        mpz_t *elems;
        mpz_init(vector[i]);
        elems = (mpz_t *) malloc(sizeof(mpz_t) * clt(s).secparam);
        for (unsigned long j = 0; j < clt(s).secparam; ++j) {
            mpz_init(elems[j]);
            py_to_mpz(elems[j],
                      PyList_GET_ITEM(PyList_GET_ITEM(py_vectors, j), i));
        }
        clt_mlm_encode(&clt(s), vector[i], clt(s).secparam,
                       (const mpz_t *) elems, 2, indices, pows);
        for (unsigned long j = 0; j < clt(s).secparam; ++j) {
            mpz_clear(elems[j]);
        }
        free(elems);
    }
    (void) write_clt_vector(s->dir, (const mpz_t *) vector, length, name);
    for (ssize_t i = 0; i < length; ++i) {
        mpz_clear(vector[i]);
    }
    free(vector);
}

static void
encode_vectors_ggh(struct state *s, const PyObject *py_vectors,
                   const int indices[2], const int pows[2], const char *name)
{
    gghlite_enc_t *vector;
    ssize_t length;

    length = PyList_GET_SIZE(PyList_GET_ITEM(py_vectors, 0));
    vector = (gghlite_enc_t *) malloc(sizeof(gghlite_enc_t) * length);

#pragma omp parallel for
    for (ssize_t i = 0; i < length; ++i) {
        gghlite_enc_t elem;
        gghlite_enc_init(elem, ggh(s).params->pk);
        gghlite_enc_init(vector[i], ggh(s).params->pk);
        py_to_fmpz_mod_poly(elem, PyList_GET_ITEM(py_vectors, i));
        ggh_mlm_encode(&ggh(s), vector[i], elem);
        gghlite_enc_clear(elem);
    }
    (void) write_ggh_vector(s->dir, (const gghlite_enc_t *) vector, length,
                            name);
    for (ssize_t i = 0; i < length; ++i) {
        gghlite_enc_clear(vector[i]);
    }
    free(vector);
}

static void
encode_layers_clt(struct state *s, long inp, long idx, long nrows, long ncols,
                  PyObject *py_zero_ms, PyObject *py_one_ms,
                  int zero_indices[2], int one_indices[2],
                  const int pows[2])
{
    mpz_t *zero, *one;

    zero = (mpz_t *) malloc(sizeof(mpz_t) * nrows * ncols);
    one = (mpz_t *) malloc(sizeof(mpz_t) * nrows * ncols);

#pragma omp parallel for
    for (Py_ssize_t ctr = 0; ctr < 2 * nrows * ncols; ++ctr) {
        mpz_t *elems;
        PyObject *py_array;
        int *indices;
        mpz_t *val;
        size_t i;

        if (ctr < nrows * ncols) {
            i = ctr;
            val = &zero[i];
            py_array = py_zero_ms;
            indices = zero_indices;
        } else {
            i = ctr - nrows * ncols;
            val = &one[i];
            py_array = py_one_ms;
            indices = one_indices;
        }

        mpz_init(*val);
        elems = (mpz_t *) malloc(sizeof(mpz_t) * clt(s).secparam);
        for (unsigned long j = 0; j < clt(s).secparam; ++j) {
            mpz_init(elems[j]);
            py_to_mpz(elems[j],
                      PyList_GET_ITEM(PyList_GET_ITEM(py_array, j), i));
        }
        clt_mlm_encode(&clt(s), *val, clt(s).secparam, (const mpz_t *) elems, 2,
                       indices, pows);
        for (unsigned long j = 0; j < clt(s).secparam; ++j) {
            mpz_clear(elems[j]);
        }
        free(elems);
    }

    (void) write_clt_layer(s->dir, inp, idx, (const mpz_t *) zero,
                           (const mpz_t *) one, nrows, ncols);

    for (int i = 0; i < nrows * ncols; ++i) {
        mpz_clears(zero[i], one[i], NULL);
    }
    free(zero);
    free(one);
}

//
//
// Python functions
//
//

static PyObject *
obf_setup_clt(PyObject *self, PyObject *args)
{
    long kappa, size, nthreads;
    struct state *s;

    s = (struct state *) malloc(sizeof(struct state));
    if (s == NULL)
        return NULL;
    s->mlm.choice = CLT;

    if (!PyArg_ParseTuple(args, "llllsl", &clt(s).secparam, &kappa, &size,
                          &clt(s).nzs, &s->dir, &nthreads)) {
        free(s);
        return NULL;
    }

    (void) omp_set_num_threads(nthreads);

    {
        long *pows;

        pows = (long *) malloc(sizeof(long) * clt(s).nzs);
        for (unsigned long i = 0; i < clt(s).nzs; ++i) {
            pows[i] = 1L;
        }

        (void) clt_mlm_setup(&clt(s), s->dir, pows, kappa, size, g_verbose);

        free(pows);
    }

    /* Convert g_i values to python objects */
    {
        PyObject *py_gs, *py_state;

        py_gs = PyList_New(clt(s).secparam);

        //
        // Only convert the first secparam g_i values since we only need to fill in
        // the first secparam slots of the plaintext space.
        //
        for (unsigned long i = 0; i < clt(s).secparam; ++i) {
            PyList_SetItem(py_gs, i, mpz_to_py(clt(s).gs[i]));
        }

        /* Encapsulate state as python object */
        py_state = PyCapsule_New((void *) s, NULL, NULL);

        return PyTuple_Pack(2, py_state, py_gs);
    }
}

static PyObject *
obf_setup_ggh(PyObject *self, PyObject *args)
{
    long secparam, kappa, nzs, size, nthreads;
    struct state *s;

    s = (struct state *) malloc(sizeof(struct state));
    if (s == NULL)
        return NULL;
    s->mlm.choice = GGH;

    if (!PyArg_ParseTuple(args, "lllsl", &secparam, &kappa, &size, &nzs,
                          &s->dir, &nthreads)) {
        free(s);
        return NULL;
    }

    (void) omp_set_num_threads(nthreads);

    ggh_mlm_setup(&ggh(s), secparam, kappa);

    {
        PyObject *py_state, *py_prime;

        fprintf(stderr, "g = ");
        (void) fmpz_poly_fprint_pretty(stderr, ggh(s).params->g, "x");
        fprintf(stderr, "\n");

        py_prime = PyCapsule_New((void *) ggh(s).params->g, NULL, NULL);
        py_state = PyCapsule_New((void *) s, NULL, NULL);

        return PyTuple_Pack(2, py_state, py_prime);
    }
}

//
// Encode N vectors across all slots of the MLM
//
static PyObject *
obf_encode_vectors(PyObject *self, PyObject *args)
{
    PyObject *py_state, *py_vectors, *py_list;
    char *name;
    int indices[2];
    int pows[] = {1, 1};
    double start, end;
    struct state *s;

    if (!PyArg_ParseTuple(args, "OOOs", &py_state, &py_vectors, &py_list, &name))
        return NULL;

    s = (struct state *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;

    (void) extract_indices(py_list, &indices[0], &indices[1]);

    start = current_time();
    switch (s->mlm.choice) {
    case CLT:
        encode_vectors_clt(s, py_vectors, indices, pows, name);
        break;
    case GGH:
        encode_vectors_ggh(s, py_vectors, indices, pows, name);
        break;
    default:
        return NULL;
    }

    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Encoding %ld elements: %f\n",
                       PyList_GET_SIZE(PyList_GET_ITEM(py_vectors, 0)),
                       end - start);

    Py_RETURN_NONE;
}

//
// Encode N layers across all slots of the MLM
//
static PyObject *
obf_encode_layers(PyObject *self, PyObject *args)
{
    PyObject *py_zero_ms, *py_one_ms;
    PyObject *py_zero_set, *py_one_set;
    PyObject *py_state;
    int zero_indices[2], one_indices[2];
    int pows[] = {1, 1};
    long inp, idx, nrows, ncols;
    double start, end;
    struct state *s;

    if (!PyArg_ParseTuple(args, "OllllOOOO", &py_state, &idx, &nrows, &ncols,
                          &inp, &py_zero_ms, &py_one_ms, &py_zero_set,
                          &py_one_set))
        return NULL;

    s = (struct state *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;

    start = current_time();

    (void) extract_indices(py_zero_set, &zero_indices[0], &zero_indices[1]);
    (void) extract_indices(py_one_set, &one_indices[0], &one_indices[1]);

    switch (s->mlm.choice) {
    case CLT:
        encode_layers_clt(s, inp, idx, nrows, ncols, py_zero_ms, py_one_ms,
                          zero_indices, one_indices, pows);
        break;
    case GGH:
        return NULL;
    default:
        return NULL;
    }

    end = current_time();
    if (g_verbose)
        (void) fprintf(stderr, "  Encoding %ld elements: %f\n",
                       2 * nrows * ncols, end - start);

    Py_RETURN_NONE;
}

static PyObject *
obf_evaluate_agis(PyObject *self, PyObject *args)
{
    char *dir, *input;
    long bplen;
    int iszero;

    if (!PyArg_ParseTuple(args, "ssl", &dir, &input, &bplen))
        return NULL;

    iszero = evaluate_agis(dir, input, bplen);

    if (iszero == -1)
        return NULL;
    else
        return Py_BuildValue("i", iszero ? 0 : 1);
}

static PyObject *
obf_evaluate_sz(PyObject *self, PyObject *args)
{
    char *dir, *input;
    long bplen;
    int iszero;

    if (!PyArg_ParseTuple(args, "ssl", &dir, &input, &bplen))
        return NULL;

    iszero = evaluate_sz(dir, input, bplen);

    if (iszero == -1)
        return NULL;
    else
        return Py_BuildValue("i", iszero ? 0 : 1);
}

static PyObject *
obf_cleanup(PyObject *self, PyObject *args)
{
    PyObject *py_state;
    struct state *s;

    if (!PyArg_ParseTuple(args, "O", &py_state))
        return NULL;

    s = (struct state *) PyCapsule_GetPointer(py_state, NULL);
    if (s == NULL)
        return NULL;

    switch (s->mlm.choice) {
    case CLT:
        clt_mlm_cleanup(&clt(s));
        break;
    case GGH:
        ggh_mlm_cleanup(&ggh(s));
        break;
    }

    free(s);

    Py_RETURN_NONE;
}

static PyMethodDef
ObfMethods[] = {
    {"verbose", obf_verbose, METH_VARARGS,
     "Set verbosity."},
    {"setup_clt", obf_setup_clt, METH_VARARGS,
     "Set up obfuscator using CLT multilinear map."},
    {"setup_ggh", obf_setup_ggh, METH_VARARGS,
     "Set up obfuscator using GGH multilinear map."},
    {"encode_vectors", obf_encode_vectors, METH_VARARGS,
     "Encode a vector in each slot."},
    {"encode_layers", obf_encode_layers, METH_VARARGS,
     "Encode a branching program layer in each slot."},
    {"max_mem_usage", obf_max_mem_usage, METH_VARARGS,
     "Print out the maximum memory usage."},
    {"cleanup", obf_cleanup, METH_VARARGS,
     "Clean up objects created during setup."},
    {"evaluate_agis", obf_evaluate_agis, METH_VARARGS,
     "evaluate the obfuscation using AGIS."},
    {"evaluate_sz", obf_evaluate_sz, METH_VARARGS,
     "evaluate the obfuscation using SZ."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
init_obfuscator(void)
{
    (void) Py_InitModule("_obfuscator", ObfMethods);
}
