#ifndef BWTINDEX_H_
#define BWTINDEX_H_

#include "bwt.h"

int64_t bwa_seq_len(const char *fn_pac);
bwt_t *bwt_pac2bwt(const char *fn_pac, int use_is);

#endif
