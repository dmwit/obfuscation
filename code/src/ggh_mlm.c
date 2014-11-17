#include "ggh_mlm.h"

void
ggh_mlm_setup(struct ggh_mlm_state *s, long secparam, long kappa)
{
    flint_randinit(s->rng);

    gghlite_init(s->ggh, secparam, kappa, 0, GGHLITE_FLAGS_DEFAULT, s->rng);
    gghlite_print_params(s->ggh->pk);
}

void
ggh_mlm_cleanup(struct ggh_mlm_state *s)
{
    flint_randclear(s->rng);
    gghlite_clear(s->ggh, 1);
}

void
ggh_mlm_encode(struct ggh_mlm_state *s, gghlite_enc_t out, gghlite_enc_t in)
{
    gghlite_enc(out, s->ggh->pk, in, 0, GGHLITE_FLAGS_DEFAULT, s->rng);
}

void
ggh_mlm_add(const struct ggh_mlm_state *s, gghlite_enc_t out,
            const gghlite_enc_t a, const gghlite_enc_t b)
{
    gghlite_add(out, s->ggh->pk, a, b);
}

void
ggh_mlm_mul(const struct ggh_mlm_state *s, gghlite_enc_t out,
            const gghlite_enc_t a, const gghlite_enc_t b)
{
    gghlite_mul(out, s->ggh->pk, a, b);
}

int
ggh_mlm_is_zero(struct ggh_mlm_state *s, const gghlite_enc_t c)
{
    return gghlite_is_zero(s->ggh->pk, c);
}
