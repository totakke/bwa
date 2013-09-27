#include "bwa_java.h"

#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "bntseq.h"
#include "bwt.h"
#include "bwtindex.h"
#include "utils.h"

#ifdef _DIVBWT
#include "divsufsort.h"
#endif

#ifdef USE_MALLOC_WRAPPERS
#include "malloc_wrap.h"
#endif

int bwa_java_index(char *in_fa, char *out_prefix, char *algo, int is_64)
{
    int algo_type = -1;
    char *str, *str2, *str3;
	int64_t l_pac;

    /* Validate arguments */
    if (strcmp(algo, "auto") == 0) algo_type = 0;
    else if (strcmp(algo, "div") == 0) algo_type = 1;
    else if (strcmp(algo, "bwtsw") == 0) algo_type = 2;
    else if (strcmp(algo, "is") == 0) algo_type = 3;
    else return 1;

    if (is_64 < 0 || 1 < is_64) return 1;

    /* Preparation */
    if (is_64) strcat(out_prefix, ".64");
    str  = (char *)calloc(strlen(out_prefix) + 10, 1);
    str2 = (char *)calloc(strlen(out_prefix) + 10, 1);
    str3 = (char *)calloc(strlen(out_prefix) + 10, 1);

    /* Indexing nucleotide */
    {
        gzFile fp = xzopen(in_fa, "r");
        l_pac = bns_fasta2bntseq(fp, out_prefix, 0);
        err_gzclose(fp);
    }

    if (algo_type == 0) algo_type = l_pac > 50000000 ? 2 : 3; // set the algorithm for generating BWT

    {
        strcpy(str, out_prefix); strcat(str, ".pac");
        strcpy(str2, out_prefix); strcat(str2, ".bwt");
        if (algo_type == 2) bwt_bwtgen(str, str2);
        else if (algo_type == 1 || algo_type == 3) {
            bwt_t *bwt;
            bwt = bwt_pac2bwt(str, algo_type == 3);
            bwt_dump_bwt(str2, bwt);
            bwt_destroy(bwt);
        }
    }

    {
        bwt_t *bwt;
        strcpy(str, out_prefix); strcat(str, ".bwt");
        bwt = bwt_restore_bwt(str);
        bwt_bwtupdate_core(bwt);
        bwt_dump_bwt(str, bwt);
        bwt_destroy(bwt);
    }

    {
        gzFile fp = xzopen(in_fa, "r");
        l_pac = bns_fasta2bntseq(fp, out_prefix, 1);
        err_gzclose(fp);
    }

    {
        bwt_t *bwt;
        strcpy(str, out_prefix); strcat(str, ".bwt");
        strcpy(str3, out_prefix); strcat(str3, ".sa");
        bwt = bwt_restore_bwt(str);
        bwt_cal_sa(bwt, 32);
        bwt_dump_sa(str3, bwt);
        bwt_destroy(bwt);
    }

    free(str3);
    free(str2);
    free(str);

    return 0;
}

int bwa_java_mem()
{
    return 0;
}

int bwa_java_aln()
{
    return 0;
}

int bwa_java_samse()
{
    return 0;
}

int bwa_java_sampe()
{
    return 0;
}

int bwa_java_bwasw()
{
    return 0;
}
