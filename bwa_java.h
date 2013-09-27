#ifndef BWA_JAVA_H_
#define BWA_JAVA_H_

int bwa_java_index(char *in_fa, char *out_prefix, char *algo, int is_64);
int bwa_java_mem();
int bwa_java_aln();
int bwa_java_samse();
int bwa_java_sampe();
int bwa_java_bwasw();

#endif
