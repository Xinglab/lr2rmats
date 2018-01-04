#ifndef _BAM_FUSION_H
#define _BAM_FUSION_H

int bam_fusion(int argc, char *argv[]);

typedef struct {
    float ovlp_frac, each_cov, all_cov;
    int dis;
} bam_fusion_para;

#define OVLP_FRAC 0.1
#define EACH_COV 0.1
#define ALL_COV 0.99
#define FUSION_DIS 100000
#define FUSION_DIS_STR "100k"

#endif
