#ifndef _BAM2GTF_H
#define _BAM2GTF_H
#include "htslib/sam.h"
#include "gtf.h"

#define bam_unmap(b) ((b)->core.flag & BAM_FUNMAP)
int read_bam_trans(samFile *in, bam_hdr_t *h, bam1_t *b, int min_exon, int min_intron, read_trans_t *T);
int bam2gtf(int argc, char *argv[]);

#endif
