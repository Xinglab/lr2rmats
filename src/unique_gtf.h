#ifndef _UNIQUE_GTF_H
#define _UNIQUE_GTF_H
#include "htslib/sam.h"

typedef struct {
    int input_mode; samFile *gtf_bam;
    int min_exon, min_intron, deletion_max, ss_dis, end_dis;
    float single_exon_ovlp_frac;
    FILE *out_gtf_fp;
    char source[1024];
} unique_gtf_para;


int unique_gtf(int argc, char *argv[]);


#endif
