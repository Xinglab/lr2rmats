#ifndef _UPDATE_GTF_H
#define _UPDATE_GTF_H
#include "htslib/sam.h"


typedef struct {
    uint8_t input_mode; samFile *gtf_bam;
    FILE *sj_fp;
    int min_exon, min_intron, ss_dis, full_level, split_trans;
    float single_exon_ovlp_frac;
    FILE *out_gtf_fp, *bam_gtf_fp, *bam_detail_fp, *known_gtf_fp, *novel_gtf_fp, *unrecog_gtf_fp, *summary_fp;
    char source[1024];
} update_gtf_para;

int update_gtf(int argc, char *argv[]);

#endif
