#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gtf.h"
#include "utils.h"
#include "htslib/sam.h"

extern int gen_trans(bam1_t *b, trans_t *t, int exon_min);

// exon
exon_t *exon_init(int n) {
    exon_t *e = (exon_t*)_err_malloc(n * sizeof(exon_t));
    return e;
}
void exon_free(exon_t *e) { free(e); }

//transcript
trans_t *trans_init(int n) {
    trans_t *t = (trans_t*)_err_calloc(n, sizeof(trans_t));
    strcpy(t->trans_name, "");
    t->exon_n = 0; t->exon_m = 0; t->cov = 1;
    //t->exon = exon_init(2);
    return t;
}

int add_exon(trans_t *t, int tid, int start, int end, uint8_t is_rev)
{
    if (t->exon_n == t->exon_m) t = exon_realloc(t);
    t->exon[t->exon_n].tid = tid;
    t->exon[t->exon_n].start = start;
    t->exon[t->exon_n].end = end;
    t->exon[t->exon_n].is_rev = is_rev;
    t->exon_n++;
    return 0;
}

int trans_exon_comp(const void *_a, const void *_b)
{
    exon_t *a = (exon_t*)_a, *b = (exon_t*)_b;
    if (a->is_rev != b->is_rev) err_fatal_simple("Strands of exons do NOT match.\n");
    if (a->start != b->start) return a->start - b->start;
    else return a->end - b->end;
}
// '+': s->e => S->E
// '-': S->E => s->e
void sort_exon(trans_t *t)
{
    qsort(t->exon, t->exon_n, sizeof(exon_t), trans_exon_comp);
}

// check if t1 is identical to t2 or fully contains t2
int check_iden(trans_t *t1, trans_t *t2, int ss_dis, int end_dis) {
    //if (t1->is_rev != t2->is_rev) return 0;
    trans_t *l, *s; int full_match = -1, partial_match = -1;
    if (t1->exon_n > t2->exon_n) {
        l = t1; s = t2;
        partial_match = 1;
    } else if (t1->exon_n < t2->exon_n) {
        l = t2; s = t1;
        partial_match = 2;
    } else {
        l = t1, s = t2;
        full_match = 0;
    }
    int i, j, un_iden=-1;
    if (full_match == 0) {
        if (abs(l->exon[0].start - s->exon[0].start) > end_dis) return un_iden;
        for (i = 0; i < l->exon_n-1; ++i) {
            if (abs(l->exon[i].end - s->exon[i].end) > ss_dis) return un_iden;
            if (abs(l->exon[i+1].start - s->exon[i+1].start) > ss_dis) return un_iden;
        }
        if (abs(l->exon[l->exon_n-1].end - s->exon[s->exon_n-1].end) > end_dis) return un_iden;
        return full_match;
    } else {
        partial_match = un_iden;
        if (abs(l->exon[0].start - s->exon[0].start) > end_dis) return un_iden;
        for (i = 0; i < l->exon_n-1; ++i) {
            if ((abs(l->exon[i].end - s->exon[0].end) <= ss_dis) && (abs(l->exon[i+1].start - s->exon[1].start) <= ss_dis)) {
                partial_match = 2;
                for (i = i+1, j = 1; i < l->exon_n-1 && j < s->exon_n-1; ++i, ++j) {
                    if (abs(l->exon[i].end - s->exon[j].end) > ss_dis) return un_iden;
                    if (abs(l->exon[i+1].start - s->exon[j+1].start) > ss_dis) return un_iden;
                }
                break;
            }
        }
        if (abs(l->exon[l->exon_n-1].end - s->exon[s->exon_n-1].end) > end_dis) return un_iden;
        return partial_match;
    }
}

int set_trans_name(trans_t *t, char *gene_id, char *gene_name, char *trans_id, char *trans_name)
{
    sort_exon(t);
    t->tid = t->exon[0].tid;
    t->is_rev = t->exon[0].is_rev;
    t->start = t->exon[0].start;
    t->end = t->exon[t->exon_n-1].end;
    if (gene_id) strcpy(t->gene_id, gene_id);
    if (gene_name) strcpy(t->gene_name, gene_name);
    if (trans_id) strcpy(t->trans_id, trans_id);
    if (trans_name) strcpy(t->trans_name, trans_name);
    return 0;
}

int set_gene(gene_t *g, char *gene_name)
{
    int i;
    g->tid = g->trans[0].tid; g->is_rev = g->trans[0].is_rev;
    g->start = g->trans[0].start; g->end = g->trans[0].end;
    for (i = 1; i < g->trans_n; ++i) {
        if (g->start > g->trans[i].start) g->start = g->trans[i].start;
        if (g->end < g->trans[i].end) g->end= g->trans[i].end;
    }
    if (gene_name) strcpy(g->gene_name, gene_name);
    return 0;
}

trans_t *exon_realloc(trans_t *t) {
    if (t->exon_m == 0) {
        t->exon_m = 2;
        t->exon = (exon_t*)_err_malloc(t->exon_m * sizeof(exon_t));
    } else {
        t->exon_m <<= 1;
        t->exon = (exon_t*)_err_realloc(t->exon, t->exon_m * sizeof(exon_t));
    }
    return t;
}


//for one read: multi-alignments => multi-transcripts
read_trans_t *read_trans_init(int trans_m)
{
    read_trans_t *r = (read_trans_t*)_err_malloc(sizeof(read_trans_t));
    r->trans_n = 0, r->trans_m = trans_m; r->gene_n = 0;
    r->t = trans_init(trans_m);
    return r;
}

void add_read_trans(read_trans_t *r, trans_t t)
{
    if (r->trans_n == r->trans_m) r = read_trans_realloc(r);
    int i;
    trans_t *new_t = r->t + r->trans_n;
    new_t->exon_n = 0;
    new_t->cov = t.cov;
    for (i = 0; i < t.exon_n; ++i)
        add_exon(new_t, t.exon[i].tid, t.exon[i].start, t.exon[i].end, t.exon[i].is_rev);
    new_t->tid = t.tid; new_t->is_rev = t.is_rev;
    new_t->start = t.start; new_t->end = t.end;

    strcpy(new_t->gene_id, t.gene_id); strcpy(new_t->gene_name, t.gene_name);
    strcpy(new_t->trans_id, t.trans_id); strcpy(new_t->trans_name, t.trans_name);

    new_t->full = t.full, new_t->lfull = t.lfull, new_t->lnoth = t.lnoth, new_t->rfull = t.rfull, new_t->rnoth = t.rnoth; 
    new_t->known = t.known, new_t->has_known_site = t.has_known_site, new_t->has_unreliable_junction = t.has_unreliable_junction, new_t->partial_read = t.partial_read;
    new_t->novel_exon_flag = (uint8_t*)_err_malloc(new_t->exon_n * sizeof(uint8_t)); memcpy(new_t->novel_exon_flag, t.novel_exon_flag, new_t->exon_n);
    new_t->novel_site_flag = (uint8_t*)_err_malloc((new_t->exon_n-1) * 2 * sizeof(uint8_t)); memcpy(new_t->novel_site_flag, t.novel_site_flag, (new_t->exon_n-1)*2);
    new_t->novel_junction_flag = (uint8_t*)_err_malloc((new_t->exon_n-1) * sizeof(uint8_t)); memcpy(new_t->novel_junction_flag, t.novel_junction_flag, new_t->exon_n-1);
    new_t->unreliable_junction_flag = (uint8_t*)_err_malloc((new_t->exon_n-1) * sizeof(uint8_t)); memcpy(new_t->unreliable_junction_flag, t.unreliable_junction_flag, new_t->exon_n-1);
    r->trans_n++;
}

void modify_read_trans(trans_t *T, trans_t t)
{
    int i;
    T->exon_n = 0;
    T->cov = t.cov;
    for (i = 0; i < t.exon_n; ++i)
        add_exon(T, t.exon[i].tid, t.exon[i].start, t.exon[i].end, t.exon[i].is_rev);
    T->tid = t.tid; T->is_rev = t.is_rev;
    T->start = t.start; T->end = t.end;

    strcpy(T->gene_id, t.gene_id); strcpy(T->gene_name, t.gene_name);
    strcpy(T->trans_id, t.trans_id); strcpy(T->trans_name, t.trans_name);


    T->full = t.full, T->lfull = t.lfull, T->lnoth = t.lnoth, T->rfull = t.rfull, T->rnoth = t.rnoth; 
    T->known = t.known, T->has_known_site = t.has_known_site, T->has_unreliable_junction = t.has_unreliable_junction, T->partial_read = t.partial_read;
    T->novel_exon_flag = (uint8_t*)_err_realloc(T->novel_exon_flag, T->exon_n * sizeof(uint8_t)); memcpy(T->novel_exon_flag, t.novel_exon_flag, T->exon_n);
    T->novel_site_flag = (uint8_t*)_err_realloc(T->novel_site_flag, (T->exon_n-1) * 2 * sizeof(uint8_t)); memcpy(T->novel_site_flag, t.novel_site_flag, (T->exon_n-1)*2);
    T->novel_junction_flag = (uint8_t*)_err_realloc(T->novel_junction_flag, (T->exon_n-1) * sizeof(uint8_t)); memcpy(T->novel_junction_flag, t.novel_junction_flag, T->exon_n-1);
    T->unreliable_junction_flag = (uint8_t*)_err_realloc(T->unreliable_junction_flag, (T->exon_n-1) * sizeof(uint8_t)); memcpy(T->unreliable_junction_flag, t.unreliable_junction_flag, T->exon_n-1);
}

void add_anno_trans(read_trans_t *r, trans_t t)
{
    if (r->trans_n == r->trans_m) r = read_trans_realloc(r);
    int i;
    trans_t *new_t = r->t + r->trans_n;
    new_t->exon_n = 0;
    new_t->cov = 1;
    for (i = 0; i < t.exon_n; ++i)
        add_exon(new_t, t.exon[i].tid, t.exon[i].start, t.exon[i].end, t.exon[i].is_rev);
    new_t->tid = t.tid; new_t->is_rev = t.is_rev;
    new_t->start = t.start; new_t->end = t.end;

    strcpy(new_t->gene_id, t.gene_id); strcpy(new_t->gene_name, t.gene_name);
    strcpy(new_t->trans_id, t.trans_id); strcpy(new_t->trans_name, t.trans_name);

    r->trans_n++;
}

read_trans_t *read_trans_realloc(read_trans_t *r)
{
    r->trans_m <<= 1;
    r->t = (trans_t*)_err_realloc(r->t, r->trans_m * sizeof(trans_t));

    int i;
    for (i = (r->trans_m >> 1); i < r->trans_m; ++i) {
        r->t[i].exon_n = 0, r->t[i].exon_m = 0;
    //    r->t[i].exon = exon_init(2);
    }
    return r;
}

void read_trans_free1(trans_t *t) { 
    free(t->exon); 
    free(t->novel_exon_flag);
    free(t->novel_site_flag);
    free(t->novel_junction_flag);
    free(t->unreliable_junction_flag);
    free(t); 
}

void read_trans_free(read_trans_t *r)
{
    int i;
    for (i = 0; i < r->trans_n; ++i) {
        free(r->t[i].exon);
        free(r->t[i].novel_exon_flag);
        free(r->t[i].novel_site_flag);
        free(r->t[i].novel_junction_flag);
        free(r->t[i].unreliable_junction_flag);
    }
    free(r->t); free(r);
}

void trans_free(read_trans_t *r)
{
    int i;
    for (i = 0; i < r->trans_n; ++i) free(r->t[i].exon);
    free(r->t); free(r);
}

// sj_group
sj_t *sj_init(int n)
{
    sj_t *i = (sj_t *)_err_calloc(n, sizeof(sj_t));
    return i;
}


//gene
gene_t *gene_init(void) {
    gene_t *g = (gene_t*)_err_malloc(sizeof(gene_t));
    g->trans_n = 0, g->trans_m = 1;
    g->trans = trans_init(1);
    return g;
}

gene_t *copy_gene(gene_t *g) {
    gene_t *r_g = gene_init();
    r_g->tid = g->tid; r_g->is_rev = g->is_rev;
    r_g->start = g->start, r_g->end = g->end;
    strcpy(r_g->gene_name, g->gene_name), strcpy(r_g->gene_id, g->gene_id);
    int i;
    for (i = 0; i < g->trans_n; ++i) {
        add_trans(r_g, g->trans[i]);
    }

    return r_g;
}

void add_trans(gene_t *g, trans_t t)
{
    if (g->trans_n == g->trans_m) g = trans_realloc(g);
    int i;
    g->trans[g->trans_n].tid = t.tid;
    g->trans[g->trans_n].is_rev = t.is_rev;
    g->trans[g->trans_n].start = t.start;
    g->trans[g->trans_n].end = t.end;
    g->trans[g->trans_n].cov = 1;
    strcpy(g->trans[g->trans_n].trans_name, t.trans_name);
    strcpy(g->trans[g->trans_n].trans_id, t.trans_id);

    g->trans[g->trans_n].exon_n = 0;
    for (i = 0; i < t.exon_n; ++i) {
        add_exon(g->trans+g->trans_n, t.exon[i].tid, t.exon[i].start, t.exon[i].end, t.exon[i].is_rev);
    }
    // update gene start and end
    g->trans_n++;
}

gene_t *trans_realloc(gene_t *g) {
    g->trans_m <<= 1;
    g->trans = (trans_t*)_err_realloc(g->trans, g->trans_m * sizeof(trans_t));
    int i;
    for (i = (g->trans_m >> 1); i < g->trans_m; ++i) {
        g->trans[i].exon_n = 0, g->trans[i].exon_m = 2;
        g->trans[i].exon = exon_init(2);
    }
    return g;
}

void gene_free(gene_t *g) {
    int i;
    for (i = 0; i < g->trans_m; ++i) {
        free(g->trans[i].exon);
    }
    free(g->trans); free(g);
}

// gtf additional information
void gtf_add_info(char add_info[], char tag[], char *info)
{
    int i;
    for (i=0; add_info[i] != '\0'; ++i) {
        if (strncmp(add_info+i, tag, strlen(tag)) == 0) {
            sscanf(add_info+i+strlen(tag)+2, "%[^\"]", info);
            return;
        }
    }
}

// gene_group
gene_group_t *gene_group_init(void)
{
    gene_group_t *gg = (gene_group_t*)_err_malloc(sizeof(gene_group_t));
    gg->g = gene_init(); gg->gene_n = 0; gg->gene_m = 1;
    return gg;
}

chr_name_t *chr_name_init(void)
{
    chr_name_t *cname = (chr_name_t*)_err_malloc(sizeof(chr_name_t));
    cname->chr_name = (char**)_err_malloc(30 * sizeof(char*));
    int i;
    for (i = 0; i < 30; ++i) cname->chr_name[i] = (char*)_err_malloc(100 * sizeof(char));
    cname->chr_n = 0; cname->chr_m = 30;
    return cname;
}

void chr_name_free(chr_name_t *cname)
{
    int i; for (i = 0; i < cname->chr_m; ++i) free(cname->chr_name[i]);
    free(cname->chr_name); free(cname);
}

gene_group_t *gene_group_realloc(gene_group_t *gg)
{
    int i;
    gg->gene_m <<= 1;
    gg->g = (gene_t*)_err_realloc(gg->g, gg->gene_m * sizeof(gene_t));
    for (i=gg->gene_m>>1; i < gg->gene_m; ++i) {
        gg->g[i].trans_n = 0; gg->g[i].trans_m = 1;
        gg->g[i].trans = trans_init(1);
    }
    return gg;
}

void add_gene(gene_group_t *gg, gene_t g)
{
    if (gg->gene_n == gg->gene_m) gg = gene_group_realloc(gg);
    int i;
    for (i = 0; i < g.trans_n; ++i)
        add_trans(gg->g+gg->gene_n, g.trans[i]);
    strcpy(gg->g[gg->gene_n].gene_name, g.gene_name);
    gg->g[gg->gene_n].tid = g.tid;
    gg->g[gg->gene_n].start = g.start;
    gg->g[gg->gene_n].end = g.end;
    gg->g[gg->gene_n].is_rev = g.is_rev;
    gg->gene_n++;
}

void gene_group_free(gene_group_t *gg)
{
    int i, j;
    for (i = 0; i < gg->gene_m; ++i) {
        for (j = 0; j < gg->g[i].trans_m; ++j)
            free(gg->g[i].trans[j].exon);
        free(gg->g[i].trans);
    }
    free(gg->g); free(gg);
}

int get_chr_id(chr_name_t *cname, char *chr)
{
    int i;
    for (i = 0; i < cname->chr_n; ++i) {
        if (strcmp(cname->chr_name[i], chr) == 0) return i;
    }
    if (cname->chr_n == cname->chr_m) {
        cname->chr_m <<= 1;
        cname->chr_name = (char**)_err_realloc(cname->chr_name, cname->chr_m * sizeof(char*));
        for (i = cname->chr_m>>1; i < cname->chr_m; ++i)
            cname->chr_name[i] = (char*)_err_malloc(100 * sizeof(char));
    }
    strcpy(cname->chr_name[cname->chr_n], chr);
    return cname->chr_n++;
}

int bam_set_cname(bam_hdr_t *h, chr_name_t *cname)
{
    int i;
    for (i = 0; i < h->n_targets; ++i) {
        get_chr_id(cname, h->target_name[i]);
    }
    return 0;
}

int sj_group_comp(const void *_a, const void *_b)
{
    sj_t *a = (sj_t*)_a, *b = (sj_t*)_b;
    if (a->tid != b->tid) return a->tid - b->tid;
    else if (a->don != b->don) return a->don - b->don;
    else return a->acc - b->acc;
}

int gene_group_comp(const void *_a, const void *_b)
{
    gene_t *a = (gene_t*)_a, *b = (gene_t*)_b;
    if (a->tid != b->tid) return a->tid - b->tid;
    else if (a->start != b->start) return a->start - b->start;
    else return a->end - b->end;
}

// read splice-junction
int read_sj_group(FILE *sj_fp, chr_name_t *cname, sj_t **sj_group, int sj_m)
{
    if (sj_fp == NULL) return 0;
    char line[1024], ref[100]="\0";
    int sj_n = 0;
    int _strand, _motif, _is_anno;
    while (fgets(line, 1024, sj_fp) != NULL) {
        if (sj_n == sj_m) _realloc(*sj_group, sj_m, sj_t)
        sj_t *sj = (*sj_group)+sj_n;

        sscanf(line, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d", ref, &(sj->don), &(sj->acc), &(_strand), &(_motif), &(_is_anno), &(sj->uniq_c), &(sj->multi_c), &(sj->max_over));
        sj->strand = _strand, sj->is_rev = (_strand == 1 ? 0 : 1), sj->motif = _motif, sj->is_anno = _is_anno;
        int tid = get_chr_id(cname, ref);
        (*sj_group)[sj_n++].tid = tid;
    }
    // sort with cname
    qsort(*sj_group, sj_n, sizeof(sj_t), sj_group_comp);
    return sj_n;
}

void reverse_exon_order(gene_group_t *gg) {
    int i, j, k; exon_t tmp;
    for (i = 0; i < gg->gene_n; ++i) {
        for (j = 0; j < gg->g[i].trans_n; ++j) {
            if (gg->g[i].trans[j].is_rev == 0) continue;
            int exon_n = gg->g[i].trans[j].exon_n;
            for (k = 0; k < exon_n >> 1; ++k) {
                tmp = gg->g[i].trans[j].exon[k];
                gg->g[i].trans[j].exon[k] = gg->g[i].trans[j].exon[exon_n-1-k];
                gg->g[i].trans[j].exon[exon_n-1-k] = tmp;
            }
        }
    }
}

// from annotation gtf file extract transcript-exon structure
int read_anno_trans(FILE *fp, bam_hdr_t *h, read_trans_t *T)
{
    char line[1024], ref[100]="\0", type[20]="\0"; int start, end; char strand, add_info[1024], name[100], last_gene[100]="";
    trans_t *t = trans_init(1); int new_gene = 0;
    while (fgets(line, 1024, fp) != NULL) {
        sscanf(line, "%s\t%*s\t%s\t%d\t%d\t%*s\t%c\t%*s\t%[^\n]", ref, type, &start, &end, &strand, add_info);
        uint8_t is_rev = (strand == '-' ? 1 : 0);
        if (strcmp(type, "transcript") == 0) {
            if (t->exon_n >= 1) {
                set_trans_name(t, NULL, NULL, NULL, NULL);
                add_anno_trans(T, *t);
                if (new_gene) T->gene_n++;
                new_gene = 0;
                read_trans_free1(t);
                t = trans_init(1);
            }
        } else if (strcmp(type, "exon") == 0) { // exon
            add_exon(t, bam_name2id(h, ref), start, end, is_rev);
            char tag[20]="gene_id";
            gtf_add_info(add_info, tag, name); strcpy(t->gene_id, name);
            strcpy(tag, "gene_name");
            gtf_add_info(add_info, tag, name); strcpy(t->gene_name, name);
            if (strcmp(name, last_gene) != 0) {
                new_gene = 1;
                strcpy(last_gene, name);
            }
            strcpy(tag, "transcript_name");
            gtf_add_info(add_info, tag, name); strcpy(t->trans_name, name);
            strcpy(tag, "transcript_id");
            gtf_add_info(add_info, tag, name); strcpy(t->trans_id, name);
        }
    }
    if (t->exon_n != 0) {
        set_trans_name(t, NULL, NULL, NULL, NULL);
        add_anno_trans(T, *t);
        if (new_gene) T->gene_n++;
    }
    read_trans_free1(t);
    return T->trans_n;
}

// from gtf file extract transcript-exon structure
int read_gtf_trans(FILE *fp, bam_hdr_t *h, read_trans_t *T)
{
    char line[1024], ref[100]="\0", type[20]="\0"; int start, end; char strand, add_info[1024], name[100];
    trans_t *t = trans_init(1);
    while (fgets(line, 1024, fp) != NULL) {
        sscanf(line, "%s\t%*s\t%s\t%d\t%d\t%*s\t%c\t%*s\t%[^\n]", ref, type, &start, &end, &strand, add_info);
        uint8_t is_rev = (strand == '-' ? 1 : 0);
        if (strcmp(type, "transcript") == 0) {
            if (t->exon_n >= 1) {
                // for bam_trans
                t->full = 0, t->lfull = 0, t->lnoth = 1, t->rfull = 0, t->rnoth = 1;
                t->known = 0; t->has_known_site = 0; t->has_unreliable_junction = 0; t->partial_read = 0; //t->polyA = 0;
                t->novel_exon_flag = (uint8_t*)_err_malloc(t->exon_n * sizeof(uint8_t)); memset(t->novel_exon_flag, 1, t->exon_n);
                t->novel_site_flag = (uint8_t*)_err_malloc((t->exon_n-1)*2 * sizeof(uint8_t)); memset(t->novel_site_flag, 1, (t->exon_n-1)*2);
                t->novel_junction_flag = (uint8_t*)_err_malloc((t->exon_n-1) * sizeof(uint8_t)); memset(t->novel_junction_flag, 1, t->exon_n-1);
                t->unreliable_junction_flag = (uint8_t*)_err_malloc((t->exon_n-1) * sizeof(uint8_t)); memset(t->unreliable_junction_flag, 0, t->exon_n-1);
                set_trans_name(t, NULL, NULL, NULL, NULL);
                add_read_trans(T, *t);
                read_trans_free1(t);
                t = trans_init(1);
            }
        } else if (strcmp(type, "exon") == 0) { // exon
            add_exon(t, bam_name2id(h, ref), start, end, is_rev);
            char tag[20]="gene_id";
            gtf_add_info(add_info, tag, name); strcpy(t->gene_id, name);
            strcpy(tag, "gene_name");
            gtf_add_info(add_info, tag, name); strcpy(t->gene_name, name);
            strcpy(tag, "transcript_name");
            gtf_add_info(add_info, tag, name); strcpy(t->trans_name, name);
            strcpy(tag, "transcript_id");
            gtf_add_info(add_info, tag, name); strcpy(t->trans_id, name);
        }
    }
    if (t->exon_n != 0) {
        // for bam_trans
        t->full = 0, t->lfull = 0, t->lnoth = 1, t->rfull = 0, t->rnoth = 1;
        t->known = 0; t->has_known_site = 0; t->has_unreliable_junction = 0; t->partial_read = 0; //t->polyA = 0;
        t->novel_exon_flag = (uint8_t*)_err_malloc(t->exon_n * sizeof(uint8_t)); memset(t->novel_exon_flag, 1, t->exon_n);
        t->novel_site_flag = (uint8_t*)_err_malloc((t->exon_n-1)*2 * sizeof(uint8_t)); memset(t->novel_site_flag, 1, (t->exon_n-1)*2);
        t->novel_junction_flag = (uint8_t*)_err_malloc((t->exon_n-1) * sizeof(uint8_t)); memset(t->novel_junction_flag, 1, t->exon_n-1);
        t->unreliable_junction_flag = (uint8_t*)_err_malloc((t->exon_n-1) * sizeof(uint8_t)); memset(t->unreliable_junction_flag, 0, t->exon_n-1);

        set_trans_name(t, NULL, NULL, NULL, NULL);
        add_read_trans(T, *t);
    }
    read_trans_free1(t);
    return T->trans_n;
}

int print_trans(trans_t t, chr_name_t *cname, char *src, FILE *out)
{
    int i;
    fprintf(out, "%s\t%s\t%s\t%d\t%d\t.\t%c\t.\tgene_id \"%s\"; transcript_id \"%s\";\n", cname->chr_name[t.tid], src, "transcript", t.start, t.end, "+-"[t.is_rev], t.gene_id, t.trans_id);
    for (i = 0; i < t.exon_n; ++i)
        fprintf(out, "%s\t%s\t%s\t%d\t%d\t.\t%c\t.\tgene_id \"%s\"; transcript_id \"%s\";\n", cname->chr_name[t.tid], src, "exon", t.exon[i].start, t.exon[i].end, "+-"[t.exon[i].is_rev], t.gene_id, t.trans_id);
    return 0;
}

// tid source feature start end score(.) strand phase(.) additional
int print_read_trans(read_trans_t *read_trans, chr_name_t *cname, char *src, FILE *out)
{
    int i, j; char tmp[1024], exon_attr[1024], trans_attr[1024];
    //int score_min = 450, score_step=50;

    for (i = 0; i < read_trans->trans_n; ++i) {
        memset(tmp, 0, 1024); memset(exon_attr, 0, 1024);
        if (strlen(read_trans->t[i].gene_id) > 0) sprintf(tmp, " gene_id \"%s\";", read_trans->t[i].gene_id), strcat(exon_attr, tmp);
        if (strlen(read_trans->t[i].trans_id) > 0) sprintf(tmp, " transcript_id \"%s\";", read_trans->t[i].trans_id), strcat(exon_attr, tmp);
        if (strlen(read_trans->t[i].gene_name) > 0) sprintf(tmp, " gene_name \"%s\";", read_trans->t[i].gene_name), strcat(exon_attr, tmp);
        if (strlen(read_trans->t[i].trans_name) > 0) sprintf(tmp, " transcript_name \"%s\";", read_trans->t[i].trans_name), strcat(exon_attr, tmp);
        strcpy(trans_attr, exon_attr); sprintf(tmp, " transcript_cov \"%d\";", read_trans->t[i].cov), strcat(trans_attr, tmp);

        fprintf(out, "%s\t%s\t%s\t%d\t%d\t.\t%c\t.\t%s\n", cname->chr_name[read_trans->t[i].tid], src, "transcript", read_trans->t[i].start, read_trans->t[i].end, "+-"[read_trans->t[i].is_rev], trans_attr+1);

        if (read_trans->t[i].is_rev) { // '-' strand
            for (j = read_trans->t[i].exon_n-1; j >= 0; --j)
                fprintf(out, "%s\t%s\t%s\t%d\t%d\t.\t%c\t.\t%s\n", cname->chr_name[read_trans->t[i].exon[j].tid], src, "exon", read_trans->t[i].exon[j].start, read_trans->t[i].exon[j].end, "+-"[read_trans->t[i].exon[j].is_rev], exon_attr+1);
        } else { // '+' strand
            for (j = 0; j < read_trans->t[i].exon_n; ++j)
                fprintf(out, "%s\t%s\t%s\t%d\t%d\t.\t%c\t.\t%s\n", cname->chr_name[read_trans->t[i].exon[j].tid], src, "exon", read_trans->t[i].exon[j].start, read_trans->t[i].exon[j].end, "+-"[read_trans->t[i].exon[j].is_rev], exon_attr+1);
        }
    }
    //err_func_format_printf(__func__, "Total transcript: %d\n", read_trans->trans_n);
    return 0;
}

void print_gene(FILE *out, char *src, gene_t *gene, char **cname) {
    int i, j; char tmp[1024], name[1024];
    // print gene line
    memset(tmp, 0, 1024); memset(name, 0, 1024);
    if (strlen(gene->gene_id) > 0) sprintf(tmp, " gene_id \"%s\";", gene->gene_id), strcat(name, tmp);
    if (strlen(gene->gene_name) > 0) sprintf(tmp, " gene_name \"%s\";", gene->gene_name), strcat(name, tmp);
    fprintf(out, "%s\t%s\t%s\t%d\t%d\t.\t%c\t.\t%s\n", cname[gene->tid], src, "gene", gene->start, gene->end, "+-"[gene->is_rev], name+1);

    for (i = 0; i < gene->trans_n; ++i) {
        trans_t *t = gene->trans + i;
        memset(tmp, 0, 1024); memset(name, 0, 1024);
        if (strlen(gene->gene_id) > 0) sprintf(tmp, " gene_id \"%s\";", gene->gene_id), strcat(name, tmp);
        if (strlen(t->trans_id) > 0) sprintf(tmp, " transcript_id \"%s\";", t->trans_id), strcat(name, tmp);
        if (strlen(gene->gene_name) > 0) sprintf(tmp, " gene_name \"%s\";", gene->gene_name), strcat(name, tmp);
        if (strlen(t->trans_name) > 0) sprintf(tmp, " transcript_name \"%s\";", t->trans_name), strcat(name, tmp);

        // print transcript line
        fprintf(out, "%s\t%s\t%s\t%d\t%d\t.\t%c\t.\t%s\n", cname[t->tid], src, "transcript", t->start, t->end, "+-"[t->is_rev], name+1);

        // print exon line
        if (t->is_rev) { // '-' strand
            for (j = t->exon_n-1; j >= 0; --j)
                fprintf(out, "%s\t%s\t%s\t%d\t%d\t0\t%c\t.\t%s\n", cname[t->exon[j].tid], src, "exon", t->exon[j].start, t->exon[j].end, "+-"[t->exon[j].is_rev], name+1);
        } else { // '+' strand
            for (j = 0; j < t->exon_n; ++j)
                fprintf(out, "%s\t%s\t%s\t%d\t%d\t0\t%c\t.\t%s\n", cname[t->exon[j].tid], src, "exon", t->exon[j].start, t->exon[j].end, "+-"[t->exon[j].is_rev], name+1);
        }
    }
}

read_trans_t *min_trans_set(read_trans_t *trans_T) {
    read_trans_t *min_trans_T;
    // for each trans cluster, build one splice-graph
    // generate minimal trans set from splice-graph with Heaviest-Bundling
    return min_trans_T;
}

/*void print_gtf_trans(gene_t g, bam_hdr_t *h, char *src, FILE *out)
{
    if (g.trans_n <= g.anno_tran_n) return;
    int i, j; char gene_name[100]; int score_min=450, score_step=50;
    for (i = g.anno_tran_n; i < g.trans_n; ++i) {
        if (g.trans[i].novel_gene_flag) strcpy(gene_name, "UNCLASSIFIED");
        else strcpy(gene_name, g.gene_name);
        trans_t t = g.trans[i];
        fprintf(out, "%s\t%s\t%s\t%d\t%d\t.\t%c\t.\tgene_id \"%s\"; transcript_id \"%s\";\n", h->target_name[t.tid], src, "transcript", t.start, t.end, "+-"[t.is_rev], gene_name, t.trans_name);
        for (j = 0; j < t.exon_n; ++j)
            fprintf(out, "%s\t%s\t%s\t%d\t%d\t%d\t%c\t.\tgene_id \"%s\"; transcript_id \"%s\";\n", h->target_name[t.tid], src, "exon", t.exon[j].start, t.exon[j].end, score_min+t.cov*score_step, "+-"[t.exon[j].is_rev], gene_name, t.trans_name);
    }
}*/

/*void print_gene_group(gene_group_t gg, bam_hdr_t *h, char *src, FILE *out, char **group_line, int *group_line_n)
{
    int l_i = 0, i, j;
    for (i = 0; i < gg.gene_n; ++i) {
        // print anno
        for (j = 0; j < group_line_n[i]; ++j)
            fprintf(out, "%s", group_line[l_i++]);
        // print novel trans
        print_gtf_trans(gg.g[i], h, src, out);
    }
}*/
