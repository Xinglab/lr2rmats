#ifndef _BAM2GTF_H
#define _BAM2GTF_H
#include "htslib/sam.h"
#include "gtf.h"

#define bam_unmap(b) ((b)->core.flag & BAM_FUNMAP)

typedef struct {
    uint8_t input_mode, uncla, full_len_level, only_bam;
    char in_bam[1024], source[1024];
    FILE *intron_fp, *out_gtf_fp;
    int min_exon, min_intron, ss_dis;
} update_gtf_para;

int read_bam_trans(samFile *in, bam_hdr_t *h, bam1_t *b, update_gtf_para *ugp, read_trans_t *T);
int read_intron_group(intron_group_t *I, FILE *fp);
int read_anno_trans1(read_trans_t *T, FILE *fp);

int bam2gtf(int argc, char *argv[]);

#endif
