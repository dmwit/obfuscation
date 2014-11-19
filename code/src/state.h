#ifndef __STATE_H__
#define __STATE_H__

#include "clt_mlm.h"
#include "ggh_mlm.h"

enum mlm_type { CLT, GGH };

struct mlm_state {
    enum mlm_type choice;
    union {
        struct clt_mlm_state clt;
        struct ggh_mlm_state ggh;
    };
};

struct state {
    struct mlm_state mlm;
    char *dir;
};

#define clt(s) (s)->mlm.clt
#define ggh(s) (s)->mlm.ggh

#endif
