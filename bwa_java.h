#ifndef BWA_JAVA_H_
#define BWA_JAVA_H_

int bwa_java_index(char *in_fa, char *out_prefix, char *algo, int is_64);

int bwa_java_mem(char *ref_fa, char *in_fq, char *out_sam,
                 int n_threads, int min_seed_len, int w, int zdrop,
                 float fplit_factor, int max_occ, int skip_mate_rescue,
                 int skip_pairing, int a, int b, int q, int r,
                 int pen_clip5, int pen_unpaired, int p,
                 int T, int append_comment, int mark_shorter_split);

int bwa_java_aln();
int bwa_java_samse();
int bwa_java_sampe();
int bwa_java_bwasw();

#endif
