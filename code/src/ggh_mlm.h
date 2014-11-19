#ifndef __GGH_MLM_H__
#define __GGH_MLM_H__

#include "gghlite/gghlite.h"

struct ggh_mlm_state {
    flint_rand_t rng;
    gghlite_t params;
};

void
ggh_mlm_setup(struct ggh_mlm_state *s, long secparam, long kappa);

void
ggh_mlm_cleanup(struct ggh_mlm_state *s);

void
ggh_mlm_encode(struct ggh_mlm_state *s, gghlite_enc_t out, gghlite_enc_t in);

void
ggh_mlm_add(const struct ggh_mlm_state *s, gghlite_enc_t out,
            const gghlite_enc_t a, const gghlite_enc_t b);

void
ggh_mlm_mul(const struct ggh_mlm_state *s, gghlite_enc_t out,
            const gghlite_enc_t a, const gghlite_enc_t b);

int
ggh_mlm_is_zero(struct ggh_mlm_state *s, const gghlite_enc_t c);

#endif
