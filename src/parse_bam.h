#ifndef _BAM_SJ_H
#define _BAM_SJ_H
#include <stdlib.h>
#include "htslib/sam.h"
#include "gtf.h"
#include "kseq.h"
#include "utils.h"

KSEQ_INIT(gzFile, gzread)

typedef struct {
    uint32_t tid:30, is_uniq:1, is_splice:1;
    int intv_n, intv_m;
    int32_t start, end, rlen;
    int32_t *exon_end, *intr_end; // exon_end[intv_n]: exonic, intr_end[intv_n-1]: intronic
} ad_t;   // alignment details: start, end, intv_n, intv[]

typedef struct {
    int tid, is_rev, score, ed;
    int read_start, read_end;
    int ref_start, ref_end; char gname[1024];
    bam1_t *b;
} bam_seg_t;


#define PAIR "paried"
#define SING "single"
#define PAIR_T 1
#define SING_T 0


#define bam_is_prop(b) (((b)->core.flag&BAM_FPROPER_PAIR) != 0)
#define bam_unmap(b) ((b)->core.flag & BAM_FUNMAP)


typedef struct {
    char fn[1024];
    hts_idx_t *idx;
    samFile *in;
    bam_hdr_t *h;

    bam1_t *b;
    hts_itr_t *itr;
} bam_aux_t;

int bam_ref_len(bam1_t *b);
int bam_query_len(bam1_t *b);
int ad_sim_comp(ad_t *ad1, ad_t *ad2);
int ad_comp(ad_t *ad1, ad_t *ad2);
ad_t *ad_init(int n);
void ad_copy(ad_t *dest, ad_t *src);
int push_exon_coor(exon_t **e, int *e_n, int *e_m, ad_t *ad);
int push_sj(int **don, int *don_n, int *don_m, ad_t *ad);
exon_t *infer_exon_coor(int *infer_e_n, exon_t *e, int e_n, int *don, int don_n);

kseq_t *kseq_load_genome(gzFile genome_fp, int *_seq_n, int *_seq_m);
int bam2sj(int argc, char *argv[]);
void free_sj_group(sj_t *sj_g, int sj_n);
void free_ad_group(ad_t *ad_g, int ad_n);
uint8_t bam_is_uniq_NH(bam1_t *b);

bam_aux_t *bam_aux_init();
void bam_aux_destroy(bam_aux_t *aux);


bam_seg_t *bam_seg_init(int n);
void bam_seg_free(bam_seg_t *s, int n);
bam_seg_t *push_bam_seg(bam_seg_t *seg, int *seg_n, int *seg_m, bam_seg_t *s);
void move_bam_seg(bam_seg_t *seg, int dest, int src);
int bam2seg(bam1_t *b, bam_seg_t *s);

#define err_sam_open(in, fn) { if ((in = sam_open(fn, "rb")) == NULL) err_fatal(__func__, "fail to open \"%s\"\n", fn); }
#define err_sam_hdr_read(h, in, fn) { if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "fail to read header for \"%s\"\n", fn); }
#define err_sam_idx_load(idx, in, fn) { if ((idx = sam_index_load(in, fn)) == NULL) err_fatal(__func__, "fail to load the BAM index for \"%s\"\n", fn); }

#endif
