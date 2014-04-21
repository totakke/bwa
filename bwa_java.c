#include "bwa_java.h"

#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "bwa.h"
#include "bntseq.h"
#include "bwt.h"
#include "bwtindex.h"
#include "bwamem.h"
#include "fastmap.h"
#include "kseq.h"
#include "utils.h"

#ifdef _DIVBWT
#include "divsufsort.h"
#endif

#ifdef USE_MALLOC_WRAPPERS
#include "malloc_wrap.h"
#endif

KSEQ_DECLARE(gzFile)

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

/* From fastmap.c/main_mem */
int bwa_java_mem(char *ref_fa, char *in_fq, char *out_sam,
                 int n_threads, int min_seed_len, int w, int zdrop,
                 float split_factor, int max_occ, int skip_mate_rescue,
                 int skip_pairing, int a, int b,
                 int o_del, int o_ins, int e_del, int e_ins,
                 int pen_clip5, int pen_unpaired, int p,
                 int T, int append_comment, int mark_shorter_split)
{
    mem_opt_t *opt;
    int fd, i, n, copy_comment = 0;
    gzFile fp, fp2 = 0;
    kseq_t *ks, *ks2 = 0;
    bseq1_t *seqs;
    bwaidx_t *idx;
    char *rg_line = 0;
    void *ko = 0, *ko2 = 0;
    FILE *fpo;
    int64_t n_processed = 0;

    /* Validate arguments */
    if (skip_mate_rescue < 0 || 1 < skip_mate_rescue) return 1;
    if (skip_pairing < 0 || 1 < skip_pairing) return 1;
    if (p < 0 || 1 < p) return 1;
    if (append_comment < 0 || 1 < append_comment) return 1;
    if (mark_shorter_split < 0 || 1 < mark_shorter_split) return 1;

    opt = mem_opt_init();
    opt->min_seed_len = min_seed_len;
    opt->w = w;
    opt->a = a;
    opt->b = b;
    opt->o_del = o_del;
    opt->o_ins = o_ins;
    opt->e_del = e_del;
    opt->e_ins = e_ins;
    opt->T = T;
    opt->pen_unpaired = pen_unpaired;
    opt->n_threads = n_threads > 1 ? n_threads : 1;
    if (skip_pairing) opt->flag |= MEM_F_NOPAIRING;
    if (p) opt->flag |= MEM_F_PE;
    if (mark_shorter_split) opt->flag |= MEM_F_NO_MULTI;
    if (skip_mate_rescue) opt->flag |= MEM_F_NO_RESCUE;
    opt->max_occ = max_occ;
    opt->zdrop = zdrop;
    opt->split_factor = split_factor;
    if (append_comment) copy_comment = 1;
    opt->pen_clip5 = opt->pen_clip3 = pen_clip5;

    bwa_fill_scmat(opt->a, opt->b, opt->mat);
    if ((idx = bwa_idx_load(ref_fa, BWA_IDX_ALL)) == 0) return 1; // FIXME: memory leak

    ko = kopen(in_fq, &fd);
    if (ko == 0) {
        if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, in_fq);
        return 1;
    }
    fp = gzdopen(fd, "r");
    ks = kseq_init(fp);

    fpo = xopen(out_sam, "w");

    bwa_fprint_sam_hdr(fpo, idx->bns, rg_line);
    while ((seqs = bseq_read(opt->chunk_size * opt->n_threads, &n, ks, ks2)) != 0) {
        int64_t size = 0;
        if ((opt->flag & MEM_F_PE) && (n&1) == 1) {
            if (bwa_verbose >= 2)
                fprintf(stderr, "[W::%s] odd number of reads in the PE mode; last read dropped\n", __func__);
            n = n>>1<<1;
        }
        if (!copy_comment)
            for (i = 0; i < n; ++i) {
                free(seqs[i].comment); seqs[i].comment = 0;
            }
        for (i = 0; i < n; ++i) size += seqs[i].l_seq;
        if (bwa_verbose >= 3)
            fprintf(stderr, "[M::%s] read %d sequences (%ld bp)...\n", __func__, n, (long)size);
        mem_process_seqs(opt, idx->bwt, idx->bns, idx->pac, n_processed, n, seqs, 0);
        n_processed += n;
        for (i = 0; i < n; ++i) {
            err_fputs(seqs[i].sam, fpo);
            free(seqs[i].name); free(seqs[i].comment); free(seqs[i].seq); free(seqs[i].qual); free(seqs[i].sam);
        }
        err_fflush(fpo);
        free(seqs);
    }

    free(opt);
    bwa_idx_destroy(idx);
    kseq_destroy(ks);
    err_fclose(fpo);
    err_gzclose(fp); kclose(ko);
    if (ks2) {
        kseq_destroy(ks2);
        err_gzclose(fp2); kclose(ko2);
    }

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
