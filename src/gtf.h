#ifndef _GTF_H
#define _GTF_H

#include <stdint.h>
#include <stdio.h>
#include "htslib/sam.h"

#define MAX_SITE 2147483647
#define DON_SITE_F 0
#define ACC_SITE_F 1

typedef struct {
    int start, end; //1-base
} coor_t;

typedef struct {
    int32_t tid; uint8_t is_rev;
    int32_t start, end; //1-based, ref
                        //0: start of init exon
                        //MAX: end of term exon
    int sg_node_id;
} exon_t;

typedef struct {
    int32_t tid; int32_t don, acc;
    uint8_t strand:2, is_anno:2, motif:4;// strand: 0:undefined, 1:+, 2:-
    int uniq_c, multi_c, max_over;
} sj_t;


#define set_l_iden(map) (map |= 0x4)
#define set_r_iden(map) (map |= 0x2)
#define set_b_iden(map) (map |= 0x1)

#define check_l_iden(map) (map & 0x4)
#define check_r_iden(map) (map & 0x2)
#define check_b_iden(map) (map & 0x1)

typedef struct {
    exon_t *exon; int exon_n, exon_m;
    uint8_t *novel_exon_map, *novel_sj_map; // 3-bit map: l-iden | r-iden | both-iden
    int tid; uint8_t is_rev;
    int start, end;
    char tname[100], trans_id[100];
    char gname[100], gid[100];
    int novel_gene_flag, cov;
    uint8_t lfull:2, lnoth:2, rfull:2, rnoth:2;
    uint8_t full:2, novel:2, all_novel:2, all_iden:2;
} trans_t;

typedef struct {
    int tid;
    int start, end; //1-based, ref
    uint8_t is_rev:2,is_canon:2, is_anno:4;
    int uniq_c, multi_c;
} intron_t;

typedef struct {
    intron_t *intron; int intron_n, intron_m;
} intron_group_t;

typedef struct {
    trans_t *t; int trans_n, trans_m;
} read_trans_t;

typedef struct {
    trans_t *trans; int trans_n, anno_tran_n, trans_m;
    int tid; uint8_t is_rev;
    int start, end;
    char gname[1024], gid[1024];
} gene_t;

typedef struct {
    gene_t *g; int gene_n, gene_m;
} gene_group_t;

typedef struct {
    char **chr_name;
    int chr_n, chr_m;
} chr_name_t;

exon_t *exon_init(int n);
void exon_free(exon_t *e);


chr_name_t *chr_name_init(void);
void chr_name_free(chr_name_t *cname);
int read_sj_group(FILE *sj_fp, chr_name_t *cname, sj_t **sj_group, int sj_m);
int bam_set_cname(bam_hdr_t *h, chr_name_t *cname);

trans_t *trans_init(int n);
int add_exon(trans_t *t, int tid, int start, int end, uint8_t is_rev);
void sort_exon(trans_t *t);
int set_trans_name(trans_t *t, char *gid, char *gname, char *trans_id, char *tname);
trans_t *exon_realloc(trans_t *t);
void trans_free(trans_t *t);

read_trans_t *read_trans_init(void);
void add_read_trans(read_trans_t *r, trans_t t);
read_trans_t *read_trans_realloc(read_trans_t *r);
void novel_read_trans_free(read_trans_t *r);
void read_trans_free(read_trans_t *r);
//int set_read_trans(read_trans_t *r);

intron_t *intron_init(int n);
intron_group_t *intron_group_init(void);
void add_intron(intron_group_t *i, intron_t i1);
int read_intron_group(intron_group_t *I, FILE *fp);

void intron_group_free(intron_group_t *i);

gene_t *gene_init(void);
gene_t *copy_gene(gene_t *g);
void add_trans(gene_t *g, trans_t t, int novel_gene_flag);
gene_t *trans_realloc(gene_t *g);
void gene_free(gene_t *g);
void gtf_add_info(char add_info[], char tag[], char *info);

gene_group_t *gene_group_init(void);
gene_group_t *gene_group_realloc(gene_group_t *gg);
void add_gene(gene_group_t *gg, gene_t g, int novel_gene_flag);
void set_gene_group(gene_group_t *gg);
void gene_group_free(gene_group_t *gg);
int read_gene_group(char *fn, chr_name_t *cname, gene_group_t *gg);

int print_exon(exon_t e, FILE *out);
int print_trans(trans_t t, bam_hdr_t *h, char *src, FILE *out);
int print_read_trans(read_trans_t *anno_T, read_trans_t *novel_T, bam_hdr_t *h, char *src, FILE *out);
void print_gene(FILE* out, char *src, gene_t *g, char **cname);
void print_gene_group(gene_group_t gg, bam_hdr_t *h, char *src, FILE *out, char **group_line, int *group_line_n);
void print_gtf_trans(gene_t g, bam_hdr_t *h, char *src, FILE *out);

// for filtering splice-junction
#define INTRON_MIN_LEN 3
#define INTER_EXON_MIN_LEN 3
#define SPLICE_DISTANCE 0
#define MIN_INTRON_NUM 0
#define ANCHOR_MIN_LEN 1 // for annotated sj
#define UNIQ_MIN 0 // for annotated sj
#define ALL_MIN 0 // for annotated sj
// for novel sj
#define NON_ANCHOR 30 // non-canonical anchor-len
#define ANCHOR1 12    // GT/AG, CT/AC anchor-len
#define ANCHOR2 12    // GC/AG, CT/GC anchor-len
#define ANCHOR3 12    // AT/AC, GT/AT anchor-len
#define NON_UNIQ_MIN 3 // non-canonical uniq-map
#define UNIQ_MIN1 1    // GT/AG, CT/AC uniq-map
#define UNIQ_MIN2 1    // GC/AG, CT/GC uniq-map
#define UNIQ_MIN3 1    // AT/AC, GT/AT uniq-map
#define NON_ALL_MIN 3 // non-canonical all-map
#define ALL_MIN1 1    // GT/AG, CT/AC all-map 
#define ALL_MIN2 1    // GC/AG, CT/GC all-map
#define ALL_MIN3 1    // AT/AC, GT/AT all-map


int check_sub_iden(trans_t *t1, trans_t *t2, int dis);
int check_iden(trans_t *t1, trans_t *t2, int dis);

#endif
