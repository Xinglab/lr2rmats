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
    trans_t *t = (trans_t*)_err_malloc(n * sizeof(trans_t));
    strcpy(t->tname, "");
    t->exon_n = 0; t->exon_m = 2;
    t->exon = exon_init(2);
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
    if (a->is_rev && !b->is_rev) return -1;
    else if (a->is_rev == b->is_rev) return a->start - b->start;
    else return a->end - b->end;
}
// '+': s->e => S->E
// '-': S->E => s->e
void sort_exon(trans_t *t)
{
    qsort(t->exon, t->exon_n, sizeof(exon_t), trans_exon_comp);
}

int check_sub_iden(trans_t *t1, trans_t *t2, int dis) {
    if (t1->is_rev != t2->is_rev) return 0;
    trans_t *l, *s;
    if (t1->exon_n > t2->exon_n) {
        l = t1; s = t2;
    } else {
        l = t2; s = t1;
    }
    int i, j, match=0;
    for (i = 0; i < l->exon_n-1; ++i) {
        if ((abs(l->exon[i].end - s->exon[0].end) <= dis) && (abs(l->exon[i+1].start - s->exon[1].start) <= dis)) {
            match = 1;
            for (i = i+1, j = 1; i < l->exon_n-1 && j < s->exon_n-1; ++i, ++j) {
                if (abs(l->exon[i].end - s->exon[j].end) > dis) return 0;
                if (abs(l->exon[i+1].start - s->exon[j+1].start) > dis) return 0;
            }
            break;
        }
    }
    return match;
}

int check_iden(trans_t *t1, trans_t *t2, int dis)
{
    if (t1->is_rev != t2->is_rev || t1->exon_n != t2->exon_n) return 0;
    int i;
    for (i = 0; i < t1->exon_n-1; ++i) {
        if (abs(t1->exon[i].end - t2->exon[i].end) > dis) return 0;
        if (abs(t1->exon[i+1].start - t2->exon[i+1].start) > dis) return 0;
    }
    return 1;
}

int set_trans_name(trans_t *t, char *gid, char *gname, char *tname, char *trans_id)
{
    sort_exon(t);
    t->tid = t->exon[0].tid;
    t->is_rev = t->exon[0].is_rev;
    t->start = t->exon[0].start;
    t->end = t->exon[t->exon_n-1].end;
    if (gid) strcpy(t->gid, gid);
    if (gname) strcpy(t->gname, gname);
    if (trans_id) strcpy(t->trans_id, trans_id);
    if (tname) strcpy(t->tname, tname);
    return 0;
}

int set_gene(gene_t *g, char *gname)
{
    int i;
    g->tid = g->trans[0].tid; g->is_rev = g->trans[0].is_rev;
    g->start = g->trans[0].start; g->end = g->trans[0].end;
    for (i = 1; i < g->trans_n; ++i) {
        if (g->start > g->trans[i].start) g->start = g->trans[i].start;
        if (g->end < g->trans[i].end) g->end= g->trans[i].end;
    }
    if (gname) strcpy(g->gname, gname);
    return 0;
}

trans_t *exon_realloc(trans_t *t) {
    t->exon_m <<= 1;
    t->exon = (exon_t*)_err_realloc(t->exon, t->exon_m * sizeof(exon_t));
    return t;
}

void trans_free(trans_t *t) { free(t->exon); free(t); }

//for one read: multi-alignments => multi-transcripts
read_trans_t *read_trans_init(void)
{
    read_trans_t *r = (read_trans_t*)_err_malloc(sizeof(read_trans_t));
    r->trans_n = 0, r->trans_m = 1;
    r->t = trans_init(1);
    return r;
}

void add_read_trans(read_trans_t *r, trans_t t)
{
    if (r->trans_n == r->trans_m) r = read_trans_realloc(r);
    int i;
    r->t[r->trans_n].exon_n = 0;
    r->t[r->trans_n].cov = 1;
    for (i = 0; i < t.exon_n; ++i)
        add_exon(r->t+r->trans_n, t.exon[i].tid, t.exon[i].start, t.exon[i].end, t.exon[i].is_rev);
    strcpy(r->t[r->trans_n].trans_id, t.trans_id);
    strcpy(r->t[r->trans_n].tname, t.tname);
    strcpy(r->t[r->trans_n].gid, t.gid);
    strcpy(r->t[r->trans_n].gname, t.gname);
    r->trans_n++;
}

read_trans_t *read_trans_realloc(read_trans_t *r)
{
    r->trans_m <<= 1;
    r->t = (trans_t*)_err_realloc(r->t, r->trans_m * sizeof(trans_t));
    int i;
    for (i = (r->trans_m >> 1); i < r->trans_m; ++i) {
        r->t[i].exon_n = 0, r->t[i].exon_m = 2;
        r->t[i].exon = exon_init(2);
    }
    return r;
}

void novel_read_trans_free(read_trans_t *r)
{
    int i;
    for (i = 0; i < r->trans_m; ++i) { 
        free(r->t[i].exon); 
        free(r->t[i].novel_exon_map);
        free(r->t[i].novel_sj_map);
    }
    free(r->t); free(r);
}

void read_trans_free(read_trans_t *r)
{
    int i;
    for (i = 0; i < r->trans_m; ++i) free(r->t[i].exon);
    free(r->t); free(r);
}

// intron_group
intron_t *intron_init(int n)
{
    intron_t *i = (intron_t *)_err_malloc(n * sizeof(intron_t));
    return i;
}

intron_group_t *intron_group_init(void)
{
    intron_group_t *i = (intron_group_t*)_err_malloc(sizeof(intron_group_t));
    i->intron = intron_init(2);
    i->intron_n = 0, i->intron_m = 2;
    return i;
}

intron_group_t *intron_group_realloc(intron_group_t *i){
    i->intron_m <<= 1;
    i->intron = (intron_t*)_err_realloc(i->intron, i->intron_m * sizeof(intron_t));
    return i;
}

void add_intron(intron_group_t *i, intron_t i1)
{
    if (i->intron_n == i->intron_m) {
        i = intron_group_realloc(i);
    }
    i->intron[i->intron_n].tid = i1.tid;
    i->intron[i->intron_n].is_rev = i1.is_rev;
    i->intron[i->intron_n].start = i1.start;
    i->intron[i->intron_n].end = i1.end;
    i->intron[i->intron_n].is_canon = i1.is_canon;
    i->intron[i->intron_n].is_anno = i1.is_anno;
    i->intron[i->intron_n].uniq_c = i1.uniq_c;
    i->intron[i->intron_n].multi_c = i1.multi_c;
    i->intron_n++;
}

//XXX
int name2id(char ref[])
{
    if (ref[3] == 'X') return 22;
    else if (ref[3] == 'Y') return 23;
    else if (ref[3] == 'M') return 24;
    else return atoi(ref+3)-1;
}

int read_intron_group(intron_group_t *I, FILE *fp)
{
    if (fp == NULL) return 0;
    char line[1024], ref[100]; int start, end, nstrand, canon, anno, uniq_map, multi_map, overlang;
    intron_t *i = intron_init(1);
    while (fgets(line, 1024, fp) != NULL) {
        sscanf(line, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d", ref, &start, &end, &nstrand, &canon, &anno, &uniq_map, &multi_map, &overlang);
        i->tid = name2id(ref); i->start = start, i->end = end;
        i->is_rev = (nstrand == 1 ? 0 : (nstrand == 2 ? 1 : -1)); i->is_canon = canon;
        i->uniq_c = uniq_map; i->multi_c = multi_map;
        add_intron(I, *i);
    }
    free(i);
    return I->intron_n;
}

void intron_group_free(intron_group_t *i) { free(i->intron); free(i); }

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
    strcpy(r_g->gname, g->gname), strcpy(r_g->gid, g->gid);
    int i;
    for (i = 0; i < g->trans_n; ++i) {
        add_trans(r_g, g->trans[i], 0);
    }

    return r_g;
}

void add_trans(gene_t *g, trans_t t, int novel_gene_flag)
{
    if (g->trans_n == g->trans_m) g = trans_realloc(g);
    int i;
    g->trans[g->trans_n].tid = t.tid;
    g->trans[g->trans_n].is_rev = t.is_rev;
    g->trans[g->trans_n].start = t.start;
    g->trans[g->trans_n].end = t.end;
    g->trans[g->trans_n].novel_gene_flag = novel_gene_flag;
    g->trans[g->trans_n].cov = 1;
    strcpy(g->trans[g->trans_n].tname, t.tname);
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

void add_gene(gene_group_t *gg, gene_t g, int novel_gene_flag)
{
    if (gg->gene_n == gg->gene_m) gg = gene_group_realloc(gg);
    int i;
    for (i = 0; i < g.trans_n; ++i)
        add_trans(gg->g+gg->gene_n, g.trans[i], novel_gene_flag);
    strcpy(gg->g[gg->gene_n].gname, g.gname);
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
    char line[1024], ref[100]="\0";
    int sj_n = 0;
    int _strand, _motif, _is_anno;
    while (fgets(line, 1024, sj_fp) != NULL) {
        if (sj_n == sj_m) _realloc(*sj_group, sj_m, sj_t)
        sj_t *sj = (*sj_group)+sj_n;

        sscanf(line, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d", ref, &(sj->don), &(sj->acc), &(_strand), &(_motif), &(_is_anno), &(sj->uniq_c), &(sj->multi_c), &(sj->max_over));
        sj->strand = _strand, sj->motif = _motif, sj->is_anno = _is_anno;
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

// merge overlapping gene into one complex
// sorted with chr and start
// exon sorted by start
int read_gene_group(char *fn, chr_name_t *cname, gene_group_t *gg)
{
    err_func_format_printf(__func__, "read gene annotation from GTF file ...\n");
    FILE *gtf = fopen(fn, "r");
    char line[1024], ref[100]="\0", type[20]="\0"; int start, end; char strand, add_info[1024], gname[1024], gid[1024], trans_name[1024], trans_id[1024], tag[20];
    gene_t *cur_g=0; trans_t *cur_t=0; exon_t *cur_e=0;
    gg->gene_n = 0;
    int last_tid=-1, last_start=-1, last_end=-1;

    while (fgets(line, 1024, gtf) != NULL) {
        if (line[0] == '#') continue;
        sscanf(line, "%s\t%*s\t%s\t%d\t%d\t%*s\t%c\t%*s\t%[^\n]", ref, type, &start, &end, &strand, add_info);
        int tid = get_chr_id(cname, ref);
        uint8_t is_rev = (strand == '-' ? 1 : 0);
        strcpy(tag, "gene_id"); gtf_add_info(add_info, tag, gid);
        strcpy(tag, "gene_name"); gtf_add_info(add_info, tag, gname);
        strcpy(tag, "transcript_id"); gtf_add_info(add_info, tag, trans_id);
        strcpy(tag, "transcript_name"); gtf_add_info(add_info, tag, trans_name);
 
        if (strcmp(type, "gene") == 0) { // new gene starts old gene ends
            if (tid == last_tid &&  start < last_end) {
                if (start < last_start) {
                    last_start = start;
                    cur_g->start = start;
                }
                if (end > last_end) {
                    last_end = end;
                    cur_g->end = end;
                }
                continue;
            }
            last_tid = tid, last_start = start, last_end = end;
            if (++gg->gene_n == gg->gene_m) gg = gene_group_realloc(gg);
            cur_g = gg->g + gg->gene_n-1;
            cur_g->tid = tid; cur_g->is_rev = is_rev;
            cur_g->start = start; cur_g->end = end;
            strcpy(cur_g->gname, gname); strcpy(cur_g->gid, gid);
            cur_g->trans_n = 0;
        } else if (strcmp(type, "transcript") == 0) { // new trans starts, old trans ends
            if (cur_g == 0) err_fatal_core(__func__, "GTF format error in %s.\n", fn);
            if (++cur_g->trans_n == cur_g->trans_m) cur_g = trans_realloc(cur_g);
            cur_t = cur_g->trans + cur_g->trans_n-1;
            cur_t->tid = tid; cur_t->is_rev = is_rev;
            cur_t->start = start; cur_t->end = end;
            strcpy(cur_t->tname, trans_name); strcpy(cur_t->trans_id, trans_id);
            cur_t->exon_n = 0;
        } else if (strcmp(type, "exon") == 0) { // new exon starts, old exon ends
            if (cur_t == 0) err_fatal_core(__func__, "GTF format error in %s.\n", fn);
            // add exon to gg
            if (++cur_t->exon_n == cur_t->exon_m) cur_t = exon_realloc(cur_t);
            cur_e = cur_t->exon + cur_t->exon_n-1;
            cur_e->tid = tid; cur_e->is_rev = is_rev;
            cur_e->start = start; cur_e->end = end;
        }
    }
    // reverse '-' transcript
    reverse_exon_order(gg);
    // sort with cname
    qsort(gg->g, gg->gene_n, sizeof(gene_t), gene_group_comp);
    err_fclose(gtf);
    err_func_format_printf(__func__, "read gene annotation from GTF file done!\n");
    return gg->gene_n;
}

//print
int print_exon(exon_t e, FILE *out)
{
    return 0;
}

int print_trans(trans_t t, bam_hdr_t *h, char *src, FILE *out)
{
    int i;
    fprintf(out, "%s\t%s\t%s\t%d\t%d\t.\t%c\t.\tgene_id \"%s\"; transcript_id \"%s\";\n", h->target_name[t.tid], src, "transcript", t.start, t.end, "+-"[t.is_rev], "UNCLASSIFIED", t.tname);
    for (i = 0; i < t.exon_n; ++i)
        fprintf(out, "%s\t%s\t%s\t%d\t%d\t.\t%c\t.\tgene_id \"%s\"; transcript_id \"%s\";\n", h->target_name[t.tid], src, "exon", t.exon[i].start, t.exon[i].end, "+-"[t.exon[i].is_rev], "UNCLASSIFIED",t.tname);
    return 0;
}

// tid source feature start end score(.) strand phase(.) additional
int print_read_trans(read_trans_t *anno_T, read_trans_t *novel_T, bam_hdr_t *h, char *src, FILE *out)
{
    int i, j; char tmp[1024], name[1024];
    int score_min = 450, score_step=50;

    for (i = 0; i < novel_T->trans_n; ++i) {
        memset(tmp, 0, 1024); memset(name, 0, 1024);
        if (strlen(novel_T->t[i].gid) > 0) sprintf(tmp, " gene_id \"%s\";", novel_T->t[i].gid), strcat(name, tmp);
        if (strlen(novel_T->t[i].trans_id) > 0) sprintf(tmp, " transcript_id \"%s\";", novel_T->t[i].trans_id), strcat(name, tmp);
        if (strlen(novel_T->t[i].gname) > 0) sprintf(tmp, " gene_name \"%s\";", novel_T->t[i].gname), strcat(name, tmp);
        if (strlen(novel_T->t[i].tname) > 0) sprintf(tmp, " transcript_name \"%s\";", novel_T->t[i].tname), strcat(name, tmp);

        fprintf(out, "%s\t%s\t%s\t%d\t%d\t.\t%c\t.\t%s\n", h->target_name[novel_T->t[i].tid], src, "transcript", novel_T->t[i].start, novel_T->t[i].end, "+-"[novel_T->t[i].is_rev], name+1);

        if (novel_T->t[i].is_rev) { // '-' strand
            for (j = novel_T->t[i].exon_n-1; j >= 0; --j)
                fprintf(out, "%s\t%s\t%s\t%d\t%d\t%d\t%c\t.\t%s\n", h->target_name[novel_T->t[i].exon[j].tid], src, "exon", novel_T->t[i].exon[j].start, novel_T->t[i].exon[j].end, score_min+score_step*novel_T->t[i].cov, "+-"[novel_T->t[i].exon[j].is_rev], name+1);
        } else { // '+' strand
            for (j = 0; j < novel_T->t[i].exon_n; ++j)
                fprintf(out, "%s\t%s\t%s\t%d\t%d\t%d\t%c\t.\t%s\n", h->target_name[novel_T->t[i].exon[j].tid], src, "exon", novel_T->t[i].exon[j].start, novel_T->t[i].exon[j].end, score_min+score_step*novel_T->t[i].cov, "+-"[novel_T->t[i].exon[j].is_rev], name+1);
        }
    }
    err_printf("Total novel transcript: %d\n", novel_T->trans_n);
    return 0;
}

void print_gene(FILE *out, char *src, gene_t *gene, char **cname) {
    int i, j; char tmp[1024], name[1024];
    // print gene line
    memset(tmp, 0, 1024); memset(name, 0, 1024);
    if (strlen(gene->gid) > 0) sprintf(tmp, " gene_id \"%s\";", gene->gid), strcat(name, tmp);
    if (strlen(gene->gname) > 0) sprintf(tmp, " gene_name \"%s\";", gene->gname), strcat(name, tmp);
    fprintf(out, "%s\t%s\t%s\t%d\t%d\t.\t%c\t.\t%s\n", cname[gene->tid], src, "gene", gene->start, gene->end, "+-"[gene->is_rev], name+1);

    for (i = 0; i < gene->trans_n; ++i) {
        trans_t *t = gene->trans + i;
        memset(tmp, 0, 1024); memset(name, 0, 1024);
        if (strlen(gene->gid) > 0) sprintf(tmp, " gene_id \"%s\";", gene->gid), strcat(name, tmp);
        if (strlen(t->trans_id) > 0) sprintf(tmp, " transcript_id \"%s\";", t->trans_id), strcat(name, tmp);
        if (strlen(gene->gname) > 0) sprintf(tmp, " gene_name \"%s\";", gene->gname), strcat(name, tmp);
        if (strlen(t->tname) > 0) sprintf(tmp, " transcript_name \"%s\";", t->tname), strcat(name, tmp);

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

/*void print_gtf_trans(gene_t g, bam_hdr_t *h, char *src, FILE *out)
{
    if (g.trans_n <= g.anno_tran_n) return;
    int i, j; char gene_name[100]; int score_min=450, score_step=50;
    for (i = g.anno_tran_n; i < g.trans_n; ++i) {
        if (g.trans[i].novel_gene_flag) strcpy(gene_name, "UNCLASSIFIED");
        else strcpy(gene_name, g.gname);
        trans_t t = g.trans[i];
        fprintf(out, "%s\t%s\t%s\t%d\t%d\t.\t%c\t.\tgene_id \"%s\"; transcript_id \"%s\";\n", h->target_name[t.tid], src, "transcript", t.start, t.end, "+-"[t.is_rev], gene_name, t.tname);
        for (j = 0; j < t.exon_n; ++j)
            fprintf(out, "%s\t%s\t%s\t%d\t%d\t%d\t%c\t.\tgene_id \"%s\"; transcript_id \"%s\";\n", h->target_name[t.tid], src, "exon", t.exon[j].start, t.exon[j].end, score_min+t.cov*score_step, "+-"[t.exon[j].is_rev], gene_name, t.tname);
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

