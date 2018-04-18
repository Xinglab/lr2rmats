#ifndef _UPDATE_GTF_H
#define _UPDATE_GTF_H
#include "htslib/sam.h"
#include "gtf.h"

#define MIN_SJ_CNT 1

typedef struct {
    uint8_t input_mode; samFile *gtf_bam;
    FILE *sj_fp; uint8_t use_multi; int min_sj_cnt;
    int min_exon, min_intron, ss_dis, end_dis, full_level, split_trans;
    float single_exon_ovlp_frac;
    FILE *out_gtf_fp, *exon_bed_fp, *bam_gtf_fp, *bam_detail_fp, *known_gtf_fp, *novel_gtf_fp, *unrecog_gtf_fp, *summary_fp;
    char source[1024];
} update_gtf_para;

int merge_trans(trans_t *t, read_trans_t *T, int ss_dis, int end_dis, float single_exon_ovlp_frac);

int update_gtf(int argc, char *argv[]);


#endif
