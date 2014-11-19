#include "ggh_mlm.h"

#include "utils.h"

void
ggh_mlm_setup(struct ggh_mlm_state *s, long secparam, long kappa)
{
    seed_flint_rng(&s->rng);
    gghlite_init(s->params, secparam, kappa, 0,
                 GGHLITE_FLAGS_VERBOSE, s->rng);
    gghlite_print_params(s->params->pk);
}

void
ggh_mlm_cleanup(struct ggh_mlm_state *s)
{
    gghlite_clear(s->params, 1);
    flint_randclear(s->rng);
}

void
ggh_mlm_encode(struct ggh_mlm_state *s, gghlite_enc_t out, gghlite_enc_t in)
{
    gghlite_enc(out, s->params->pk, in, 0,
                GGHLITE_FLAGS_VERBOSE, s->rng);
}

void
ggh_mlm_add(const struct ggh_mlm_state *s, gghlite_enc_t out,
            const gghlite_enc_t a, const gghlite_enc_t b)
{
    gghlite_add(out, s->params->pk, a, b);
}

void
ggh_mlm_mul(const struct ggh_mlm_state *s, gghlite_enc_t out,
            const gghlite_enc_t a, const gghlite_enc_t b)
{
    gghlite_mul(out, s->params->pk, a, b);
}

int
ggh_mlm_is_zero(struct ggh_mlm_state *s, const gghlite_enc_t c)
{
    return gghlite_is_zero(s->params->pk, c);
}
