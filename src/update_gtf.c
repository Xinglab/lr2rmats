/* update_gtf.c
 *   generate junction information based on sam/bam file
 *   then, update existing GTF file
 *   currently, only work for single-end long read data
 * 
 * Author:  Yan Gao
 * Contact: yangao07@hit.edu.cn                             */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "htslib/sam.h"
#include "utils.h"
#include "gtf.h"
#include "bam2gtf.h"
#include "update_gtf.h"
#include "parse_bam.h"
#include "utils.h"

extern const char PROG[20];

update_gtf_para *update_gtf_init_para(void) {
    update_gtf_para *ugp = (update_gtf_para*)_err_malloc(sizeof(update_gtf_para));
    ugp->input_mode = 0/*bam*/; ugp->gtf_bam = NULL;
    ugp->sj_fp = NULL; ugp->use_multi = 0; ugp->min_sj_cnt = MIN_SJ_CNT;
    ugp->min_exon = INTER_EXON_MIN_LEN, ugp->min_intron = INTRON_MIN_LEN, ugp->max_delet = DELETION_MAX_LEN, ugp->ss_dis = SPLICE_DISTANCE; ugp->full_level = 5/*most relax*/; ugp->split_trans = 0; ugp->end_dis = END_DISTANCE;
    ugp->single_exon_ovlp_frac = SING_OVLP_FRAC;
    ugp->keep_min_set = 0;
    ugp->out_gtf_fp = stdout; ugp->exon_bed_fp = NULL; ugp->bam_gtf_fp = NULL; ugp->bam_detail_fp = NULL; ugp->known_gtf_fp = NULL; ugp->novel_gtf_fp = NULL; ugp->unrecog_gtf_fp = NULL; ugp->summary_fp = NULL;
    strcpy(ugp->source, PROG);

    return ugp;
}

int update_gtf_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s update-gtf [option] <in.bam/in.gtf> <old.gtf> > new.gtf\n\n", PROG);
    err_printf("Notice:  the BAM and GTF files should be sorted in advance.\n\n");
    err_printf("Input options:\n\n");
    err_printf("         -m --input-mode   [STR]    format of input file <in.bam/in.gtf>, BAM file(b) or GTF file(g). [b]\n");
    err_printf("         -b --bam          [STR]    for GTF input <in.gtf>, BAM file is needed to obtain BAM header information. [NULL]\n");
    err_printf("         -j --sj           [STR]    junction information file output by STAR(*.out.tab). [NULL]\n");
    err_printf("\n");

    err_printf("Function options:\n\n");
    err_printf("         -e --min-exon     [INT]    minimum length of internal exon. [%d]\n", INTER_EXON_MIN_LEN);
    err_printf("         -i --min-intron   [INT]    minimum length of intron. [%d]\n", INTRON_MIN_LEN);
    err_printf("         -t --max-delet    [INT]    maximum length of deletion, longer deletion will be considered as intron. [%d]\n", DELETION_MAX_LEN);
    err_printf("         -d --distance     [INT]    consider same if distance between two splice site is not bigger than d. [%d]\n", SPLICE_DISTANCE);
    err_printf("         -D --DISTANCE     [INT]    consider same if distance between two start/end site is not bigger than D. [%d]\n", END_DISTANCE);
    err_printf("         -f --frac         [INT]    consider same if overlapping between two single-exon transcript is bigger than f. [%.2f]\n", SING_OVLP_FRAC);
    err_printf("         -s --split-trans           split read on unreliable junctions. [False]\n");
    //err_printf("         -s --match-strand         only transcript of matched strand will be 
    err_printf("         -M --use-multi             use junction information of multi-mapped read. [False]\n");
    err_printf("         -J --min-junc-cnt [INT]    minimum short-read junction count of novel junction. [%d]\n", MIN_SJ_CNT);
    err_printf("         -l --full-length  [INT]    level of strict criterion for considering full-length transcript. \n");
    err_printf("                                    (1->5, most strict->most relaxed) [%d]\n", 5);
    err_printf("\n");

    err_printf("Output options:\n\n");
    err_printf("         -o --output       [STR]    updated GTF file. [stdout]\n");
    err_printf("         -n --min-output            only keep the minimal set of novel transcripts in the updated GTF file. [False]\n");
    err_printf("         -E --exon-bed     [STR]    updated novel exon file in bed format. [NULL]\n");
    err_printf("         -a --bam-gtf      [STR]    bam-derived transcript GTF file. [NULL]\n");
    err_printf("         -A --bam-detial   [STR]    detailed information of each bam-derived transcript. [NULL]\n");
    err_printf("         -k --known-gtf    [STR]    bam-derived known transcript GTF file. [NULL]\n");
    err_printf("         -v --novel-gtf    [STR]    bam-derived novel transcript GTF file. [NULL]\n");
    err_printf("         -u --unrecog      [STR]    bam-derived unrecognized transcript GTF file. [NULL]\n");
    err_printf("         -y --summary      [STR]    Staticstic summary of bam-derived transcript. [NULL]\n");
    err_printf("         -S --source       [STR]    \'source\' field in GTF: program, database or project name. [%s]\n", PROG);
    err_printf("\n");
    return 1;
}

// overlap_len / min(len1, len2)
float exon_overlap_frac(exon_t e1, exon_t e2)
{
    if (e1.start > e2.end || e2.start > e1.end) return 0.0;
    int start1 = e1.start, end1 = e1.end;
    int start2 = e2.start, end2 = e2.end;
    //int overlap_len = end1 - start2 + 1 > 0 ? end1 - start2 + 1 : end2 - start1 + 1;
    int overlap_len = MIN_OF_TWO(end1, end2) - MAX_OF_TWO(start1, start2) + 1; 
    int min_len = MIN_OF_TWO(end1 - start1 + 1, end2 - start2 + 1);
    return (overlap_len / (min_len + 0.0));
}

int exon_overlap(exon_t e1, exon_t e2)
{
    if (e1.start > e2.end || e2.start > e1.end) return 0;
    return 1;
}

// merge contained and identical multi-exon trans
int merge_trans1(trans_t *t, trans_t *T, int ss_dis, int end_dis)
{
    int i = t->exon_n-1, j = T->exon_n-1;
    int ret = check_iden(t, T, ss_dis, end_dis);
    if (ret == 0) { // fully identical
        T->cov++;
        if (t->exon[0].start < T->exon[0].start)  {
            T->exon[0].start = t->exon[0].start;
            T->start = t->exon[0].start;
        }
        if (t->exon[i].end > T->exon[j].end) {
            T->exon[j].end = t->exon[i].end;
            T->end = t->exon[i].end;
        }
        return 1;
    } else if (ret == 2) return 1;
    else if (ret == 1) { // t fully contains T: change T to t
        modify_read_trans(T, *t);
        return 1;
    }
    else return 0;
}

// merge contained and identical single-exon trans
int merge_trans2(trans_t *t, trans_t *T, int end_dis, float single_exon_ovlp_frac)
{
    if (abs(t->exon[0].start - T->exon[0].start) > end_dis) return 0;
    if (abs(t->exon[0].end - T->exon[0].end) > end_dis) return 0;

    if (exon_overlap_frac(t->exon[0], T->exon[0]) >= single_exon_ovlp_frac) {
        T->cov++;

        if (t->exon[0].start < T->exon[0].start)  {
            T->exon[0].start = t->exon[0].start;
            T->start = t->exon[0].start;
        }
        if (t->exon[0].end > T->exon[0].end) {
            T->exon[0].end = t->exon[0].end;
            T->end = t->exon[0].end;
        }
        return 1;
    } else return 0;
}

// TODO : remove all the redundant
// what if t fully contains some transcripts in T ???
int merge_trans(trans_t *t, read_trans_t *T, int ss_dis, int end_dis, float single_exon_ovlp_frac)
{
    int i; 
    for (i = T->trans_n-1; i >= 0; --i) {
        if (t->tid > T->t[i].tid || t->start > T->t[i].end) return 0;
        if (t->exon_n == 1 && T->t[i].exon_n == 1) {
            if (merge_trans2(t, T->t+i, end_dis, single_exon_ovlp_frac)) return 1;
        } else if (t->exon_n > 1 && T->t[i].exon_n > 1) {
            if (merge_trans1(t, T->t+i, ss_dis, end_dis)) return 1;
        }
    }
    return 0;
}

typedef struct {
    int tid; uint8_t is_rev;
    int start, end;
    char gene_id[1024];
} simp_gene_t;

typedef struct {
    int tid; uint8_t is_rev; int site;
} simp_site_t;

// merge contained and identical trans
int merge_gene1(simp_gene_t *g1, simp_gene_t *g2)
{
    return (strcmp(g1->gene_id, g2->gene_id) == 0);
}

int merge_gene(simp_gene_t *G, int G_n, simp_gene_t *g)
{
    int i; 
    for (i = G_n-1; i >= 0; --i) {
        if (merge_gene1(g, G+i)) return 1;
        if (g->tid > G[i].tid) return 0;
    }
    return 0;
}

simp_gene_t *add_simp_gene(simp_gene_t *G, simp_gene_t *g, int *G_n, int *G_m) {
    if (! merge_gene(G, *G_n, g)) {
        if (*G_n == *G_m) {
            *G_m <<= 1;
            G = (simp_gene_t*)_err_realloc(G, *G_m * sizeof(simp_gene_t));
        }
        G[*G_n].tid = g->tid; G[*G_n].is_rev = G[*G_n].is_rev;
        G[*G_n].start = g->start; G[*G_n].end = g->end;
        strcpy(G[*G_n].gene_id, g->gene_id);
        (*G_n)++;
    }
    return G;
}

int merge_exon1(exon_t *e1, exon_t *e2)
{
    if (e1->tid != e2->tid || e1->start != e2->start || e1->end != e2->end) return 0;
    else return 1;
}

int merge_exon(exon_t *E, int E_n, exon_t *e, int cov)
{
    int i; 
    for (i = E_n-1; i >= 0; --i) {
        if (merge_exon1(e, E+i)) {
            E[i].score += cov;
            return 1;
        }
        if (e->tid > E[i].tid) return 0;
    }
    return 0;
}

exon_t *add_simp_exon(exon_t *E, exon_t *e, int cov, int *E_n, int *E_m) {
    if (! merge_exon(E, *E_n, e, cov)) {
        if (*E_n == *E_m) {
            *E_m <<= 1;
            E = (exon_t*)_err_realloc(E, *E_m * sizeof(exon_t));
        }
        E[*E_n].tid = e->tid; E[*E_n].is_rev = e->is_rev; E[*E_n].exon_type = e->exon_type;
        E[*E_n].start = e->start; E[*E_n].end = e->end; E[*E_n].score = cov;
        (*E_n)++;
    }
    return E;
}

int merge_site1(simp_site_t *s1, simp_site_t *s2)
{
    if (s1->tid != s2->tid || s1->site != s2->site) return 0;
    else return 1;
}

int merge_site(simp_site_t *S, int S_n, simp_site_t *s)
{
    int i; 
    for (i = S_n-1; i >= 0; --i) {
        if (merge_site1(s, S+i)) return 1;
        if (s->tid > S[i].tid) return 0;
    }
    return 0;
}

simp_site_t *add_simp_site(simp_site_t *S, simp_site_t *s, int *S_n, int *S_m) {
    if (! merge_site(S, *S_n, s)) {
        if (*S_n == *S_m) {
            (*S_m) <<= 1;
            S = (simp_site_t*)_err_realloc(S, *S_m * sizeof(simp_site_t));
        }
        S[*S_n].tid = s->tid; S[*S_n].is_rev = s->is_rev;
        S[*S_n].site = s->site;
        (*S_n)++;
    }
    return S;
}

int merge_sj1(sj_t *s1, sj_t *s2)
{
    if (s1->tid != s2->tid || s1->don != s2->don || s1->acc != s2->acc) return 0;
    else return 1;
}

int merge_sj(sj_t *S, int S_n, sj_t *s)
{
    int i; 
    for (i = S_n-1; i >= 0; --i) {
        if (merge_sj1(s, S+i)) {
            S[i].score++;
            return 1;
        } if (s->tid > S[i].tid) return 0;
    }
    return 0;
}

sj_t *add_simp_sj(sj_t *S, sj_t *s, int *S_n, int *S_m) {
    if (! merge_sj(S, *S_n, s)) {
        if (*S_n == *S_m) {
            *S_m <<= 1;
            S = (sj_t*)_err_realloc(S, *S_m * sizeof(sj_t));
        }
        S[*S_n].tid = s->tid; S[*S_n].is_rev = s->is_rev;
        S[*S_n].don = s->don; S[*S_n].acc = s->acc; S[*S_n].score = 1;
        (*S_n)++;
    }
    return S;
}

int print_bam_detail_trans(read_trans_t *bam_T, chr_name_t *cname, FILE *fp)
{
    int i, j; char na[10] = "NA";
    int tot_known_n = 0, tot_novel_n=0, tot_unrecog_n=0, tot_novel_exon_n=0, tot_novel_site_n=0, tot_novel_junction_n=0, tot_unreliable_junction_n=0;
    int first = 1;
    // header
    // Novel: 0:known, 1:novel, 2:unrecognized
    // GeneID/GeneName: NA if unrecognized
    // NovelExonIndex: NA if NovelExonCount is 0
    //           0         1    2       3      4       5         6          7          8        9               10              11              12              13                  14                  15                       16
    fprintf(fp, "ReadName\tchr\tstrand\tNovel\tGeneID\tGeneName\tExonCount\tExonStart\tExonEnd\tNovelExonCount\tNovelExonIndex\tNovelSiteCount\tNovelSiteIndex\tNovelJunctionCount\tNovelJunctionIndex\tUnreliableJunctionCount\tUnreliableJunctionIndex\n");

    for (i = 0; i < bam_T->trans_n; ++i) {
        trans_t *bam_t = bam_T->t+i;
        int novel = 0;
        if (bam_t->known) {
            tot_known_n += 1;
            novel = 0;
        } else if (bam_t->has_known_site) {
            tot_novel_n += 1;
            novel = 1;
        } else {
            tot_unrecog_n += 1;
            novel = 2;
        }
        //           0   1   2   3   4   5   6
        fprintf(fp, "%s\t%s\t%c\t%d\t%s\t%s\t%d\t", bam_t->trans_name, cname->chr_name[bam_t->tid], "+-"[bam_t->is_rev], novel, bam_t->gene_id, bam_t->gene_name, bam_t->exon_n);

        first = 1;
        for (j = 0; j < bam_t->exon_n; ++j) {
            if (first == 1) first = 0; else fprintf(fp, ",");
            //           7
            fprintf(fp, "%d", bam_t->exon[j].start);
        } fprintf(fp, "\t");
        first = 1;
        for (j = 0; j < bam_t->exon_n; ++j) {
            if (first == 1) first = 0; else fprintf(fp, ",");
            //           8
            fprintf(fp, "%d", bam_t->exon[j].end);
        } fprintf(fp, "\t");

        int novel_exon_n = 0;
        for (j = 0; j < bam_t->exon_n; ++j) {
            novel_exon_n += bam_t->novel_exon_flag[j];
        }
        //           9
        fprintf(fp, "%d\t", novel_exon_n);
        tot_novel_exon_n += novel_exon_n;
        //                                 10
        if(novel_exon_n == 0) fprintf(fp, "%s\t", na);
        else {
            first = 1;
            for (j = 0; j < bam_t->exon_n; ++j) {
                if (bam_t->novel_exon_flag[j]) {
                    if (first == 1) first = 0; else fprintf(fp, ",");
                    //           10
                    fprintf(fp, "%d", j);
                }
            } fprintf(fp, "\t");
        }

        int novel_site_n = 0;
        for (j = 0; j < (bam_t->exon_n-1)*2; ++j) {
            novel_site_n += bam_t->novel_site_flag[j];
        }
        //           11
        fprintf(fp, "%d\t", novel_site_n);
        tot_novel_site_n += novel_site_n;
        //                                 12
        if(novel_site_n == 0) fprintf(fp, "%s\t", na);
        else {
            first = 1;
            for (j = 0; j < (bam_t->exon_n-1)*2; ++j) {
                if (bam_t->novel_site_flag[j]) {
                    if (first == 1) first = 0; else fprintf(fp, ",");
                    //           12
                    fprintf(fp, "%d", j);
                }
            } fprintf(fp, "\t");
        }
        int novel_junction_n = 0;
        for (j = 0; j < bam_t->exon_n-1; ++j) {
            novel_junction_n += bam_t->novel_junction_flag[j];
        }
        //           13
        fprintf(fp, "%d\t", novel_junction_n);
        tot_novel_junction_n += novel_junction_n;
        //                                     14
        if(novel_junction_n == 0) fprintf(fp, "%s\t", na);
        else {
            first = 1;
            for (j = 0; j < bam_t->exon_n-1; ++j) {
                if (bam_t->novel_junction_flag[j]) {
                    if (first == 1) first = 0; else fprintf(fp, ",");
                    //           14 
                    fprintf(fp, "%d", j);
                }
            } fprintf(fp, "\t");
        }
        int unreliable_junction_n = 0;
        for (j = 0; j < bam_t->exon_n-1; ++j) {
            unreliable_junction_n += bam_t->unreliable_junction_flag[j];
        }
        //           15
        fprintf(fp, "%d\t", unreliable_junction_n);
        tot_unreliable_junction_n += unreliable_junction_n;
        //                                          16
        if(unreliable_junction_n == 0) fprintf(fp, "%s\t", na);
        else {
            first = 1;
            for (j = 0; j < bam_t->exon_n-1; ++j) {
                if (bam_t->unreliable_junction_flag[j]) {
                    if (first == 1) first = 0; else fprintf(fp, ",");
                    //           16
                    fprintf(fp, "%d", j);
                }
            } // fprintf(fp, "\t");
        }
        fprintf(fp, "\n");
    }

    return 0;
}

int print_trans_summary(bam_hdr_t *h, read_trans_t *anno_T, read_trans_t *updated_T, read_trans_t *bam_T, update_gtf_para *ugp, FILE *summary_fp, FILE *novel_exon_fp)
{
    int i, j;
    // annotation
    int anno_trans_n = anno_T->trans_n, anno_gene_n = anno_T->gene_n;
    // 0. updated_T
    // 0.1. total count, total_gene
    // 0.2. total novel exon/site/junction
    // 0.3. total updated gene
    int updated_trans_n=0, updated_gene_n=0, updated_full_trans_n=0, updated_partial_trans_n=0;
    int updated_novel_exon_n=0, updated_novel_don_site_n=0, updated_novel_acc_site_n=0, updated_novel_junction_n=0;

    simp_gene_t *G = (simp_gene_t*)_err_malloc(1 * sizeof(simp_gene_t)); int gene_m = 1;
    simp_gene_t *g = (simp_gene_t*)_err_malloc(sizeof(simp_gene_t));
    exon_t *novel_exon = (exon_t*)_err_malloc(1 * sizeof(exon_t)); int novel_exon_m = 1;
    simp_site_t *novel_don_site = (simp_site_t*)_err_malloc(1 * sizeof(simp_site_t)); int novel_don_site_m = 1;
    simp_site_t *novel_acc_site = (simp_site_t*)_err_malloc(1 * sizeof(simp_site_t)); int novel_acc_site_m = 1;
    simp_site_t *s = (simp_site_t*)_err_malloc(1 * sizeof(simp_site_t));
    sj_t *novel_junction = (sj_t*)_err_malloc(1 * sizeof(sj_t)); int novel_junction_m = 1;
    sj_t *sj = (sj_t*)_err_malloc(1 * sizeof(sj_t));
    updated_trans_n = updated_T->trans_n;
    for (i = 0; i < updated_T->trans_n; ++i) {
        trans_t *t = updated_T->t + i;
        g->tid = t->tid, g->is_rev = t->is_rev, g->start = t->start, g->end = t->end;
        strcpy(g->gene_id, t->gene_id);
        // gene
        G = add_simp_gene(G, g, &updated_gene_n, &gene_m);
        // partial/full-read
        updated_partial_trans_n += t->partial_read;
        // novel_exon
        for (j = 0; j < t->exon_n; ++j) {
            if (t->novel_exon_flag[j]) {
                if (t->exon_n > 1) {
                    if (j == 0 || j == t->exon_n-1)
                        t->exon[j].exon_type = 0;
                    else
                        t->exon[j].exon_type = 1;
                } else t->exon[j].exon_type = 2;
                novel_exon = add_simp_exon(novel_exon, t->exon+j, t->cov, &updated_novel_exon_n, &novel_exon_m);
            }
        }
        // novel don site
        for (j = 0; j < t->exon_n-1; ++j) {
            if (t->novel_site_flag[j*2]) {
                s->tid = t->tid, s->is_rev = t->is_rev, s->site = t->exon[j].end;
                novel_don_site = add_simp_site(novel_don_site, s, &updated_novel_don_site_n, &novel_don_site_m);
            }
        }
        // novel acc site
        for (j = 0; j < t->exon_n-1; ++j) {
            if (t->novel_site_flag[j*2+1]) {
                s->tid = t->tid, s->is_rev = t->is_rev, s->site = t->exon[j+1].start;
                novel_acc_site = add_simp_site(novel_acc_site, s, &updated_novel_acc_site_n, &novel_acc_site_m);
            }
        }
        // novel junction
        for (j = 0; j < t->exon_n-1; ++j) {
            if (t->novel_junction_flag[j]) {
                sj->tid = t->tid; sj->is_rev = t->is_rev;
                sj->don = t->exon[j].end, sj->acc = t->exon[j+1].start;
                novel_junction = add_simp_sj(novel_junction, sj, &updated_novel_junction_n, &novel_junction_m);
            }
        }
    }
    updated_full_trans_n = updated_trans_n - updated_partial_trans_n;

    // 1. known_T
    // 1.1 total_count, total_gene
    // 1.2 uniq_count
    // 2. novel_T(reliable/unreliable junction)
    // 2.1 total_count, total_gene(updated_gene_n)
    // 2.2 uniq_count
    // 3. unrecognized_T
    // 3.1 total_count
    // 3.2 uniq_count
    int known_trans_n=0, known_gene_n=0, uniq_known_trans_n=0;
    int reliable_novel_trans_n=0, uniq_reliable_novel_trans_n=0, unreliable_novel_trans_n=0, uniq_unreliable_novel_trans_n=0;
    int unrecog_trans_n=0, uniq_unrecog_trans_n=0;
    read_trans_t *uniq_known_T=read_trans_init(1), *uniq_unreliable_T=read_trans_init(1), *uniq_reliable_T=read_trans_init(1), *uniq_unrecog_T=read_trans_init(1);
    //int tot_novel_exon_n=0, tot_novel_site_n=0, tot_novel_junction_n=0, tot_unreliable_junction_n=0;
    for (i = 0; i < bam_T->trans_n; ++i) {
        trans_t *t = bam_T->t + i;
        if (t->known) { //known
            known_trans_n += 1;
            g->tid = t->tid, g->is_rev = t->is_rev, g->start = t->start, g->end = t->end; strcpy(g->gene_id, t->gene_id);
            G = add_simp_gene(G, g, &known_gene_n, &gene_m);

            if (merge_trans(t, uniq_known_T, ugp->ss_dis, ugp->end_dis, ugp->single_exon_ovlp_frac) == 0)
                add_read_trans(uniq_known_T, *t);
        } else if (t->has_known_site) { // reliable novel and unreliable novel
            if (t->has_unreliable_junction) {
                unreliable_novel_trans_n++;

                if (merge_trans(t, uniq_unreliable_T, ugp->ss_dis, ugp->end_dis, ugp->single_exon_ovlp_frac) == 0)
                    add_read_trans(uniq_unreliable_T, *t);
            } else {
                reliable_novel_trans_n++;
                if (merge_trans(t, uniq_reliable_T, ugp->ss_dis, ugp->end_dis, ugp->single_exon_ovlp_frac) == 0)
                    add_read_trans(uniq_reliable_T, *t);
            }
        } else { // unrecognized
            unrecog_trans_n++;
            if (merge_trans(t, uniq_unrecog_T, ugp->ss_dis, ugp->end_dis, ugp->single_exon_ovlp_frac) == 0)
                add_read_trans(uniq_unrecog_T, *t);
        }
    }
    uniq_known_trans_n = uniq_known_T->trans_n; uniq_reliable_novel_trans_n = uniq_reliable_T->trans_n;
    uniq_unreliable_novel_trans_n = uniq_unreliable_T->trans_n; uniq_unrecog_trans_n = uniq_unrecog_T->trans_n;

    // 2. novel_site/exon/junction, unreliable_junction
    // 2.1 total count of novel_site/exon/junction, unreliable_junction
    // 2.2 unique count of novel_site/exon/junction, unreliable_junction
    // 2.3 novel_site/exon/junction count, unreliable_junction -> trans_count

    if (summary_fp) { // print summary information
        // annotation
        fprintf(summary_fp, "==== Annotaion ====\n");
        fprintf(summary_fp, "%s\t%d\n", "Genes_of_annotation_GTF", anno_gene_n);
        fprintf(summary_fp, "%s\t%d\n", "Transcripts_of_annotation_GTF", anno_trans_n);
        fprintf(summary_fp, "\n===================\n");
        // updated information
        fprintf(summary_fp, "\n==== Updated information ====\n");
        fprintf(summary_fp, "%s\t%d\n", "Updated_Genes", updated_gene_n);
        fprintf(summary_fp, "%s\t%d\n", "Added_Novel_Transcripts", updated_full_trans_n+updated_partial_trans_n);
        fprintf(summary_fp, "%s\t%d\n", "Added_Novel_Full-read_Transcripts", updated_full_trans_n);
        fprintf(summary_fp, "%s\t%d\n", "Added_Novel_Partial-read_Transcripts", updated_partial_trans_n);
        fprintf(summary_fp, "%s\t%d\n", "Added_Novel_Exons", updated_novel_exon_n);
        fprintf(summary_fp, "%s\t%d\n", "Added_Novel_Sites", updated_novel_don_site_n+updated_novel_acc_site_n);
        fprintf(summary_fp, "%s\t%d\n", "Added_Novel_Splice_Junctions", updated_novel_junction_n);
        fprintf(summary_fp, "\n=============================\n");
        // known information
        fprintf(summary_fp, "\n==== Known information ====\n");
        fprintf(summary_fp, "%s\t%d\n", "Known_Transcripts_from_BAM", known_trans_n);
        fprintf(summary_fp, "%s\t%d\n", "Genes_of_Known_Transcripts_from_BAM", known_gene_n);
        fprintf(summary_fp, "%s\t%d\n", "Uniq_Known_Transcripts_from_BAM", uniq_known_trans_n);
        fprintf(summary_fp, "\n===========================\n");
        // novel information
        fprintf(summary_fp, "\n==== Novel information ====\n");
        fprintf(summary_fp, "%s\t%d\n", "Novel_Transcript_from_BAM", reliable_novel_trans_n+unreliable_novel_trans_n);
        fprintf(summary_fp, "%s\t%d\n", "Novel_Transcript_from_BAM_with_All_Reliable_Junction", reliable_novel_trans_n);
        fprintf(summary_fp, "%s\t%d\n", "Uniq_Novel_Transcript_from_BAM_with_All_Reliable_Junction", uniq_reliable_novel_trans_n);
        fprintf(summary_fp, "%s\t%d\n", "Novel_Transcript_from_BAM_with_Unreliable_Junction", unreliable_novel_trans_n);
        fprintf(summary_fp, "%s\t%d\n", "Uniq_Novel_Transcript_from_BAM_with_Unreliable_Junction", uniq_unreliable_novel_trans_n);
        fprintf(summary_fp, "\n===========================\n");
        // unrecognized information
        fprintf(summary_fp, "\n==== Unrecognized information ====\n");
        fprintf(summary_fp, "%s\t%d\n", "Unrecognized_Transcript_from_BAM", unrecog_trans_n);
        fprintf(summary_fp, "%s\t%d\n", "Uniq_Unrecognized_Transcript_from_BAM", uniq_unrecog_trans_n);
        fprintf(summary_fp, "\n==================================\n");
    }
    if (novel_exon_fp) { // print novel exon information
        // chrom    start0base  end1base    name    count   strand
        for (i = 0; i < updated_novel_exon_n; ++i) {
            fprintf(novel_exon_fp, "%s\t%d\t%d\t%c_exon\t%d\t%c\n", h->target_name[novel_exon[i].tid], novel_exon[i].start-1, novel_exon[i].end, "TIS"[novel_exon[i].exon_type], novel_exon[i].score, "+-"[novel_exon[i].is_rev]);
        }
    }
    // TODO
    /*if (novel_sj_fp) { // print novel splice-junction information
    }
    if (novel_site_fp) { // print novel splice-site information
    }*/

    free(G); free(novel_exon); free(novel_don_site); free(novel_acc_site); free(novel_junction);
    free(g); free(s); free(sj);
    read_trans_free(uniq_known_T), read_trans_free(uniq_unreliable_T), read_trans_free(uniq_reliable_T), read_trans_free(uniq_unrecog_T);
    return 0;
}

int check_short_sj1(int tid, int start, int end, sj_t *sj_group, int sj_n, int i_start, update_gtf_para *ugp)
{
    int i = i_start, sj_cnt;
    int dis = ugp->ss_dis, min_cnt = ugp->min_sj_cnt;
    while (i < sj_n) {
        if (sj_group[i].tid > tid || (sj_group[i].tid == tid && sj_group[i].don >= end)) return 0;
        if (abs(sj_group[i].don-start)<=dis && abs(sj_group[i].acc-end)<=dis) {
            if (ugp->use_multi) sj_cnt = sj_group[i].uniq_c + sj_group[i].multi_c;
            else sj_cnt = sj_group[i].uniq_c;
            if (sj_cnt >= min_cnt) return 1;
        }
        i++;
    }
    return 0;
}


// return:
// 1: support by short sj
// 0: not support by short sj
int check_short_sj(trans_t *bam_t, int *sj_map, sj_t *sj_group, int sj_n, int *sj_i, update_gtf_para *ugp)
{
    int i = *sj_i, j, ret= 1;
    while (i < sj_n) {
        if (sj_group[i].tid < bam_t->tid || (sj_group[i].tid == bam_t->tid && sj_group[i].acc <= bam_t->start)) {
            i++; *sj_i = i;
        } else if (sj_group[i].tid > bam_t->tid || (sj_group[i].tid == bam_t->tid && sj_group[i].don >= bam_t->end)) return 0;
        else {
            for (j = 0; j < bam_t->exon_n-1; ++j) {
                if (sj_map[j] == 0 && check_short_sj1(bam_t->tid, bam_t->exon[j].end+1, bam_t->exon[j+1].start-1, sj_group, sj_n, i, ugp) == 0) {
                    bam_t->unreliable_junction_flag[j] = 1;
                    ret = 0;
                }
            }
            return ret;
        }
    }
    return 0;
}

int check_full(trans_t *t, trans_t anno_t, int level)
{
    if (t->lfull && t->rfull) return 0;
    int i = t->exon_n-1, j = anno_t.exon_n-1;
    if (level == 1) { // identical first and last splice-site
        if (t->lfull == 0) {
            if (t->exon[0].end == anno_t.exon[0].end) t->lfull = 1;
        }
        if (t->rfull == 0) {
            if (t->exon[i].start == anno_t.exon[j].start) t->rfull = 1;
        }
    } else if (level == 2) { // overlapping first and last exon
        if (t->lfull == 0) {
            if (exon_overlap(t->exon[0], anno_t.exon[0])) t->lfull = 1;
        }
        if (t->rfull == 0) {
            if (exon_overlap(t->exon[i], anno_t.exon[j])) t->rfull = 1;
        }
    } else if (level == 3) { // overlapping first and last exon, or overlapping nothing
        int ii;
        if (t->lfull == 0) {
            i = 0; j = 0;
            int ii;
            if (exon_overlap(t->exon[i], anno_t.exon[j])) t->lfull = 1;
            else {
                for (ii = 0; ii < anno_t.exon_n; ++ii) {
                    if (exon_overlap(t->exon[i], anno_t.exon[ii])) { t->lnoth = 0; break; }
                }
            }
        }
        if (t->rfull == 0) {
            i = t->exon_n-1, j = anno_t.exon_n-1;
            if (exon_overlap(t->exon[i], anno_t.exon[j])) t->rfull = 1;
            else { // overlap nothing
                for (ii = 0; ii < anno_t.exon_n; ++ii) {
                    if (exon_overlap(t->exon[i], anno_t.exon[ii])) { t->rnoth=0; break; }
                }
            }
        }
    } else if (level == 4) { // XXX most 5' exon meets #3, most 3' exon has a polyA+ tail of 15bp or longer
        if (t->lfull == 0) {
            i = 0; j = 0;
            int ii;
            if (exon_overlap(t->exon[i], anno_t.exon[j])) t->lfull = 1;
            else {
                for (ii = 0; ii < anno_t.exon_n; ++ii) {
                    if (exon_overlap(t->exon[i], anno_t.exon[ii])) { t->lnoth = 0; break; }
                }
            }
        }
    }
    return 0;
}

void set_full(trans_t *t, int l)
{
    if (l == 5) t->full = 1;
    else if (l == 4) {
        if (t->lfull == 1 || t->lnoth == 1) t->full = 1;
        else t->full = 0;
    } else if (l == 3){
        if ((t->lfull == 1 || t->lnoth == 1) && (t->rfull == 1 || t->rnoth == 1)) t->full = 1;
        else t->full = 0;
    } else {
        if (t->lfull == 1 && t->rfull == 1) t->full = 1;
        else t->full = 0;
    }
}

int check_with_short_sj(trans_t *bam_t, sj_t *sj_group, int sj_n, int *last_sj_i, update_gtf_para *ugp)
{
    int *sj_map = (int*)_err_calloc((bam_t->exon_n-1), sizeof(int));
    int i;
    for (i = 0; i < bam_t->exon_n-1; ++i) {
        sj_map[i] = 1 - bam_t->novel_junction_flag[i];
    }
    int ret = check_short_sj(bam_t, sj_map, sj_group, sj_n, last_sj_i, ugp);
    free(sj_map);
    bam_t->has_unreliable_junction = 1-ret;
    return ret;
}

// @function:
// compare bam_t with anno_t
// @return:
// 0: novel, has no known site
// 1: known
// 2: novel, has known site
int check_splice_site(trans_t *bam_t, trans_t anno_t, int dis)
{
    int bam_ovlp_site_n=0, anno_ovlp_site_n=0, identical_site_n = 0;
    int bam_all_site_n = (bam_t->exon_n - 1) * 2;// anno_all_site_n = (anno_t.exon_n - 1) * 2;
    int ovlp_start = MAX_OF_TWO(bam_t->start, anno_t.start), ovlp_end = MIN_OF_TWO(bam_t->end, anno_t.end);
    int i, j;
    for (i = 0; i < bam_t->exon_n-1; ++i) {
        if (bam_t->exon[i].end >= ovlp_start && bam_t->exon[i].end  <= ovlp_end) {
            bam_ovlp_site_n++;
        }
        if (bam_t->exon[i+1].start >= ovlp_start && bam_t->exon[i+1].start <= ovlp_end) {
            bam_ovlp_site_n++;
        }
    }
    for (i = 0; i < anno_t.exon_n-1; ++i) {
        if (anno_t.exon[i].end >= ovlp_start && anno_t.exon[i].end  <= ovlp_end) {
            anno_ovlp_site_n++;
            // search in bam_don_site => remove novel-stie
            for (j = 0; j < bam_t->exon_n-1; ++j) {
                if (abs(anno_t.exon[i].end - bam_t->exon[j].end) <= dis) {
                    identical_site_n++;
                    bam_t->novel_site_flag[j*2] = 0;
                }
            }
        }
        if (anno_t.exon[i+1].start >= ovlp_start && anno_t.exon[i+1].start <= ovlp_end) {
            anno_ovlp_site_n++;
            // search in bam_acc_site => remove novel-site
            for (j = 0; j < bam_t->exon_n-1; ++j) {
                if (abs(anno_t.exon[i+1].start - bam_t->exon[j].start) <= dis) {
                    identical_site_n++;
                    bam_t->novel_site_flag[j*2+1] = 0;
                }
            }
        }
    }
    for (i = 0; i < anno_t.exon_n; ++i) {
        // search in bam_exon  => remove novel-exon
        for (j = 0; j < bam_t->exon_n; ++j) {
            if (abs(anno_t.exon[i].start - bam_t->exon[j].start) <= dis && abs(anno_t.exon[i].end - bam_t->exon[j].end) <= dis) {
                bam_t->novel_exon_flag[j] = 0;
            }
        }
    }
    for (i = 0; i < anno_t.exon_n-1; ++i) {
        // search in bam_junction => remove novel-junction
        for (j = 0; j < bam_t->exon_n-1; ++j) {
            if (abs(anno_t.exon[i].end - bam_t->exon[j].end) <= dis && abs(anno_t.exon[i+1].start - bam_t->exon[j+1].start) <= dis) {
                bam_t->novel_junction_flag[j] = 0;
            }
        }
    }
    int ret = 0;
    if (bam_all_site_n == bam_ovlp_site_n && bam_ovlp_site_n == identical_site_n) {
        bam_t->known = 1;
        ret = 1;
    } else if (identical_site_n > 0) {
        bam_t->has_known_site = 1;
        ret = 2;
    }

    return ret;
}


//@ return:
//0: overlap
//-1: t1 ... t2
//1:  t2 ... t1
int comp_trans(trans_t t1, trans_t t2) {
    if (t1.tid < t2.tid || (t1.tid == t2.tid && t1.end <= t2.start)) return -1;
    else if (t2.tid < t1.tid || (t2.tid == t1.tid && t2.end <= t1.start)) return 1;
    else return 0;
}

void check_with_anno_trans(trans_t *bam_t, read_trans_t *anno_T, int *last_anno_i, update_gtf_para *ugp) {
    int i = *last_anno_i, single_exon = 0, ret = 0;
    int ref_anno_i = -1; char na[10] = "NA";
    if (bam_t->exon_n == 1) single_exon = 1;
    for (; i < anno_T->trans_n; ++i) {
        trans_t anno_t = anno_T->t[i];
        ret = comp_trans(*bam_t, anno_t);
        if (ret < 0) {
            break;
        } else if (ret > 0) {
            if (*last_anno_i == i) ++(*last_anno_i);
        } else {
            // full-length transcript
            check_full(bam_t, anno_t, ugp->full_level);
            if (single_exon && anno_t.exon_n == 1) {
                if (exon_overlap_frac(bam_t->exon[0], anno_t.exon[0]) >= ugp->single_exon_ovlp_frac) {
                    ref_anno_i = i;
                    bam_t->known = 1;
                    break;
                }
            } else if (single_exon == 0 && anno_t.exon_n > 1) { // check splice site/junction
                ret = check_splice_site(bam_t, anno_t, ugp->ss_dis);
                if (ret == 1) { // known
                    ref_anno_i = i;
                    break;
                } else if (ret == 2) { // has_known_site
                    ref_anno_i = i;
                }
            }
        }
    }
    if (ref_anno_i != -1) { // set gene_id, gene_name, gene_is_rev
        trans_t *anno_t = anno_T->t + ref_anno_i;
        uint8_t anno_is_rev = anno_t->is_rev;
        if (anno_is_rev != bam_t->is_rev) {
            for (i = 0; i < bam_t->exon_n; ++i) {
                bam_t->exon[i].is_rev = anno_is_rev;
            }
            bam_t->is_rev = anno_is_rev;
        }
        set_trans_name(bam_t, anno_t->gene_id, anno_t->gene_name, NULL, NULL);
    } else set_trans_name(bam_t, na, na, NULL, NULL);
    set_full(bam_t, ugp->full_level);
}

read_trans_t *split_trans(trans_t *bam_t)
{
    int i, j;
    int split_trans_n = 0, last_exon_i = 0, trans_i=0, has_novel = 0, has_known = 0;
    for (i = 0; i < bam_t->exon_n-1; ++i)
        split_trans_n += bam_t->unreliable_junction_flag[i];
    read_trans_t *split_trans = read_trans_init(split_trans_n+1);
    trans_t *t;
    for (i = 0; i < bam_t->exon_n-1; ++i) {
        if (bam_t->novel_junction_flag[i]) has_novel = 1;
        else has_known = 1;
        if (bam_t->unreliable_junction_flag[i]) {
            if (has_novel && has_known && i - last_exon_i >= 1) {// discard single-exon trans
                // copy [last_exon_i, i] to split_trans
                t = split_trans->t + trans_i;
                t->exon_n = 0, t->cov = 1;
                for (j = last_exon_i; j <= i; ++j)
                    add_exon(t, bam_t->exon[j].tid, bam_t->exon[j].start, bam_t->exon[j].end, bam_t->exon[j].is_rev);
                t->full = 0, t->lfull = 0, t->lnoth = 1, t->rfull = 0, t->rnoth = 1;
                t->known = 0; t->has_known_site = 0; t->has_unreliable_junction = 0; t->partial_read = 1;
                // set novel flag
                t->novel_exon_flag = (uint8_t*)_err_malloc(t->exon_n * sizeof(uint8_t)); memset(t->novel_exon_flag, 1, t->exon_n);
                for (j = last_exon_i; j <= i; ++j) 
                    t->novel_exon_flag[j-last_exon_i] =bam_t->novel_exon_flag[j];
                t->novel_site_flag = (uint8_t*)_err_malloc((t->exon_n-1)*2 * sizeof(uint8_t)); memset(t->novel_site_flag, 1, (t->exon_n-1)*2);
                for (j = last_exon_i; j < i; ++j) {
                    t->novel_site_flag[(j-last_exon_i)*2] = bam_t->novel_site_flag[j*2];
                    t->novel_site_flag[(j-last_exon_i)*2+1] = bam_t->novel_site_flag[j*2+1];
                }
                t->novel_junction_flag = (uint8_t*)_err_malloc((t->exon_n-1) * sizeof(uint8_t)); memset(t->novel_junction_flag, 1, t->exon_n-1);
                for (j = last_exon_i; j < i; ++j) 
                    t->novel_junction_flag[j-last_exon_i] = bam_t->novel_junction_flag[j];
                t->unreliable_junction_flag = (uint8_t*)_err_malloc((t->exon_n-1) * sizeof(uint8_t)); memset(t->unreliable_junction_flag, 0, t->exon_n-1);
                
                sprintf(t->trans_id, "%s.split.%d", bam_t->trans_id, trans_i);
                sprintf(t->trans_name, "%s.split.%d", bam_t->trans_name, trans_i);
                strcpy(t->gene_id, bam_t->gene_id);
                strcpy(t->gene_name, bam_t->gene_name);

                split_trans->trans_n = ++trans_i;
            }
            last_exon_i = i+1;
            has_novel = 0, has_known = 0;
        }
    }
    if (has_novel && has_known && i - last_exon_i >= 1) {// discard single-exon trans
        // copy [last_exon_i, i] to split_trans
        t = split_trans->t + trans_i;
        t->exon_n = 0, t->cov = 1;
        for (j = last_exon_i; j <= i; ++j)
            add_exon(t, bam_t->exon[j].tid, bam_t->exon[j].start, bam_t->exon[j].end, bam_t->exon[j].is_rev);
        t->known = 0; t->has_known_site = 0; t->has_unreliable_junction = 0; t->partial_read = 1;
        // set novel flag
        t->full = 0, t->lfull = 0, t->lnoth = 1, t->rfull = 0, t->rnoth = 1;
        t->novel_exon_flag = (uint8_t*)_err_malloc(t->exon_n * sizeof(uint8_t)); memset(t->novel_exon_flag, 1, t->exon_n);
        for (j = last_exon_i; j <= i; ++j) 
            t->novel_exon_flag[j-last_exon_i] =bam_t->novel_exon_flag[j];
        t->novel_site_flag = (uint8_t*)_err_malloc((t->exon_n-1)*2 * sizeof(uint8_t)); memset(t->novel_site_flag, 1, (t->exon_n-1)*2);
        for (j = last_exon_i; j < i; ++j) {
            t->novel_site_flag[(j-last_exon_i)*2] = bam_t->novel_site_flag[j*2];
            t->novel_site_flag[(j-last_exon_i)*2+1] = bam_t->novel_site_flag[j*2+1];
        }
        t->novel_junction_flag = (uint8_t*)_err_malloc((t->exon_n-1) * sizeof(uint8_t)); memset(t->novel_junction_flag, 1, t->exon_n-1);
        for (j = last_exon_i; j < i; ++j) 
            t->novel_junction_flag[j-last_exon_i] = bam_t->novel_junction_flag[j];
        t->unreliable_junction_flag = (uint8_t*)_err_malloc((t->exon_n-1) * sizeof(uint8_t)); memset(t->unreliable_junction_flag, 0, t->exon_n-1);

        sprintf(t->trans_id, "%s.split.%d", bam_t->trans_id, trans_i);
        sprintf(t->trans_name, "%s.split.%d", bam_t->trans_name, trans_i);
        strcpy(t->gene_id, bam_t->gene_id);
        strcpy(t->gene_name, bam_t->gene_name);

        split_trans->trans_n = ++trans_i;
    }

    return split_trans;
}

//@ input:
//bam_T: input bam/gtf transcript
//anno_T: annotation transcript
//bam_T and anno_T are sorted by chr and start
//sj_group: short-read sj information
//@ output:
//updated_T: reliable novel transcript (may be cut into small ones)
//           full-transcript level
//           poly A contained
//           full/partial read
//known_read: known bam/gtf transcript
//           full-transcript level
//           poly A contained
//novel_read: novel bam/gtf transcript
//           full-transcript level
//           poly A contained
//novel_read_additional:
//           novel_site/exon/junction
//           unreliable site/exon/junction
//           poly A contained
//unrecog_read: unrecognized bam/gtf transcript
void check_trans(read_trans_t *bam_T, read_trans_t *anno_T, sj_t *sj_group, int sj_n, read_trans_t *updated_T, read_trans_t *known_T, read_trans_t *novel_T, read_trans_t *unrecog_T, update_gtf_para *ugp) 
{
    int i, j, last_anno_i=0, last_sj_i=0;
    for (i = 0; i < bam_T->trans_n; ++i) {
        trans_t *bam_t = bam_T->t+i;

        check_with_anno_trans(bam_t, anno_T, &last_anno_i, ugp);
        if (bam_t->full == 0) continue;
        if (bam_t->known) {
            add_read_trans(known_T, *bam_t);
        } else if (bam_t->has_known_site) {
            if (sj_n == 0 || check_with_short_sj(bam_t, sj_group, sj_n, &last_sj_i, ugp)) { // novel junction are supported by short read
                add_read_trans(novel_T, *bam_t);
                if (merge_trans(bam_t, updated_T, ugp->ss_dis, ugp->end_dis, ugp->single_exon_ovlp_frac) == 0)
                    add_read_trans(updated_T, *bam_t);
            } else if (ugp->split_trans) { // has unreliable novel splice junction
                // split into short transcripts
                read_trans_t *split_read_trans = split_trans(bam_t);
                for (j = 0; j < split_read_trans->trans_n; ++j) {
                    add_read_trans(novel_T, split_read_trans->t[j]);
                    if (merge_trans(split_read_trans->t+j, updated_T, ugp->ss_dis, ugp->end_dis, ugp->single_exon_ovlp_frac) == 0) 
                        add_read_trans(updated_T, split_read_trans->t[j]);
                }
                read_trans_free(split_read_trans);
            }
        } else { // novel and no_known_site
            add_read_trans(unrecog_T, *bam_t);
        }
    }
}

const struct option update_long_opt [] = {
    { "input-mode", 1, NULL, 'm' },
    { "bam", 1, NULL, 'b' },

    { "sj", 1, NULL, 'j' },
    { "min-exon", 1, NULL, 'e' },
    { "min-intron", 1, NULL, 'i' },
    { "distance", 1, NULL, 'd' },
    { "DISTANCE", 1, NULL, 'D' },

    { "frac", 1, NULL, 'f' },

    { "full-gtf", 1, NULL, 'l' },
    { "use-multi", 0, NULL, 'M' },
    { "min_sj_cnt", 1, NULL, 'J' },

    { "output", 1, NULL, 'o' }, 
    { "bam-gtf", 1, NULL, 'a' },
    { "known-gtf", 1, NULL, 'k' },
    { "novel-gtf", 1, NULL, 'v' },
    { "unrecog", 1, NULL, 'u' },
    { "source", 1, NULL, 's' },


    { 0, 0, 0, 0}
};

int update_gtf(int argc, char *argv[])
{
    int c; 
    update_gtf_para *ugp = update_gtf_init_para();
    while ((c = getopt_long(argc, argv, "m:b:j:J:M:e:i:t:sd:D:f:l:o:nE:a:A:k:v:u:y:S:", update_long_opt, NULL)) >= 0) {
        switch(c)
        {
            case 'm': if (optarg[0] == 'b') ugp->input_mode=0; 
                          else if (optarg[0] == 'g') ugp->input_mode=1; 
                          else return update_gtf_usage(); 
                          break;
            case 'b': if ((ugp->gtf_bam = sam_open(optarg, "rb")) == NULL) {
                          err_fatal(__func__, "Cannot open \"%s\"\n", optarg); 
                          return update_gtf_usage();
                      }
                      break;

            case 'j': if ((ugp->sj_fp = fopen(optarg, "r")) == NULL) {
                          err_fatal(__func__, "Can not open splice-junction file \"%s\"\n", optarg);
                          return update_gtf_usage();
                      } 
                      break;
            case 'e': ugp->min_exon = atoi(optarg); break;
            case 'i': ugp->min_intron = atoi(optarg); break;
            case 't': ugp->max_delet = atoi(optarg); break;
            case 'd': ugp->ss_dis = atoi(optarg); break;
            case 'D': ugp->end_dis = atoi(optarg); break;
            case 'f': ugp->single_exon_ovlp_frac = atof(optarg); break;

            case 's': ugp->split_trans = 1; break;
            case 'l': ugp->full_level = atoi(optarg); break;
            case 'M': ugp->use_multi = 1; break;
            case 'J': ugp->min_sj_cnt = atoi(optarg); break;

            case 'o': ugp->out_gtf_fp = fopen(optarg, "w"); break;
            case 'n': ugp->keep_min_set = 1; break;
            case 'E': ugp->exon_bed_fp = fopen(optarg, "w"); break;
            case 'a': ugp->bam_gtf_fp = fopen(optarg, "w"); break;
            case 'A': ugp->bam_detail_fp = fopen(optarg, "w"); break;
            case 'k': ugp->known_gtf_fp = fopen(optarg, "w"); break;
            case 'v': ugp->novel_gtf_fp = fopen(optarg, "w"); break;
            case 'u': ugp->unrecog_gtf_fp = fopen(optarg, "w"); break;
            case 'y': ugp->summary_fp = fopen(optarg, "w"); break;
            case 'S': strcpy(ugp->source, optarg); break;
            default:
                      err_printf("Error: unknown option: %s.\n", optarg);
                      return update_gtf_usage();
                      break;
        }
    }
    if (argc - optind != 2) return update_gtf_usage();

    chr_name_t *cname = chr_name_init();
    read_trans_t *anno_T, *bam_T, *updated_T, *known_T, *unique_known_T, *novel_T, *unique_novel_T, *unrecog_T, *unique_unrecog_T;
    gene_group_t *gg = gene_group_init();
    anno_T = read_trans_init(1); bam_T = read_trans_init(1); 
    updated_T = read_trans_init(1); 
    known_T = read_trans_init(1); 
    novel_T = read_trans_init(1); 
    unrecog_T = read_trans_init(1);
    unique_known_T = read_trans_init(1); 
    unique_novel_T = read_trans_init(1); 
    unique_unrecog_T = read_trans_init(1);
    sj_t *sj_group = (sj_t*)_err_malloc(10000 * sizeof(sj_t)); int sj_m = 10000;

    // read all input-transcript
    samFile *in = NULL; bam_hdr_t *h; 
    if (ugp->input_mode == 0) { // bam input
        bam1_t *b;
        if ((in = sam_open(argv[optind], "rb")) == NULL) err_fatal(__func__, "Can not open \"%s\"\n", argv[optind]);
        if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", argv[optind]);
        bam_set_cname(h, cname);
        b = bam_init1(); 
        read_bam_trans(in, h, b, ugp->min_exon, ugp->min_intron, ugp->max_delet, bam_T);
        bam_destroy1(b);
    } else { // gtf input
        if ((h = sam_hdr_read(ugp->gtf_bam)) == NULL) err_fatal(__func__, "Couldn't read header of provided BAM file.\n");
        bam_set_cname(h, cname);
        FILE *fp = xopen(argv[optind], "r");
        read_gtf_trans(fp, h, bam_T);
    }

    FILE *gtf_fp = xopen(argv[optind+1], "r");
    // read all anno-transcript
    read_anno_trans(gtf_fp, h, anno_T);
    // read intron file
    int sj_n = read_sj_group(ugp->sj_fp, cname, &sj_group, sj_m);

    // identify novel transcript
    check_trans(bam_T, anno_T, sj_group, sj_n, updated_T, known_T, novel_T, unrecog_T, ugp);
    // if (ugp->keep_min_set) updated_T = min_trans_set(updated_T);

    // print novel transcript
    print_read_trans(updated_T, cname, ugp->source, ugp->out_gtf_fp);
    if (ugp->bam_gtf_fp) print_read_trans(bam_T, cname, ugp->source, ugp->bam_gtf_fp); 
    // novel info of each read
    if (ugp->bam_detail_fp) print_bam_detail_trans(bam_T, cname, ugp->bam_detail_fp);
    if (ugp->known_gtf_fp) print_read_trans(known_T, cname, ugp->source, ugp->known_gtf_fp);
    if (ugp->novel_gtf_fp) print_read_trans(novel_T, cname, ugp->source, ugp->novel_gtf_fp);
    if (ugp->unrecog_gtf_fp) print_read_trans(unrecog_T, cname, ugp->source, ugp->unrecog_gtf_fp);
    // summary of novel statistics
    if (ugp->summary_fp || ugp->exon_bed_fp) print_trans_summary(h, anno_T, updated_T, bam_T, ugp, ugp->summary_fp, ugp->exon_bed_fp);

    chr_name_free(cname);
    read_trans_free(bam_T); read_trans_free(updated_T); 
    trans_free(anno_T); err_fclose(gtf_fp);
    read_trans_free(novel_T); read_trans_free(known_T); read_trans_free(unrecog_T); 
    read_trans_free(unique_novel_T); read_trans_free(unique_known_T); read_trans_free(unique_unrecog_T); 
    free(sj_group); gene_group_free(gg);
    err_fclose(ugp->out_gtf_fp); 
    bam_hdr_destroy(h); 
    if (in) sam_close(in); 
    if (ugp->gtf_bam) sam_close(ugp->gtf_bam);
    if (ugp->sj_fp) err_fclose(ugp->sj_fp);
    if (ugp->exon_bed_fp) err_fclose(ugp->exon_bed_fp);
    if (ugp->bam_gtf_fp) err_fclose(ugp->bam_gtf_fp);
    if (ugp->bam_detail_fp) err_fclose(ugp->bam_detail_fp);
    if (ugp->known_gtf_fp) err_fclose(ugp->known_gtf_fp);
    if (ugp->novel_gtf_fp) err_fclose(ugp->novel_gtf_fp);
    if (ugp->unrecog_gtf_fp) err_fclose(ugp->unrecog_gtf_fp);
    if (ugp->summary_fp) err_fclose(ugp->summary_fp);
    free(ugp);
    return 0;
}
