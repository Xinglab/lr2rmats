#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "utils.h"
#include "htslib/sam.h"
#include "gtf.h"

#define bam_unmap(b) ((b)->core.flag & BAM_FUNMAP)
#define COV_RATIO 0.67
#define MAP_QUAL  0.75
#define SEC_RATIO 0.98

extern const char PROG[20];
extern int read_anno_trans(FILE *fp, bam_hdr_t *h, read_trans_t *T);
int filter_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s filter [option] <in.bam/sam> <rRNA.gtf> | samtools sort > out.sort.bam\n\n", PROG);
    err_printf("Options:\n");
    err_printf("         -v --coverage   [FLOAT]    minimum fraction of aligned bases. [%.2f]\n", COV_RATIO);
    err_printf("         -q --map-qual   [FLOAT]    minimum fraction of identically aligned bases. [%.2f]\n", MAP_QUAL);
    err_printf("         -s --sec-rat    [FLOAT]    maximum ratio of second best and best score to retain the best\n");
    err_printf("         -i --intron     [INT]      minimum number of intron indicated by the alignment. [%d]\n", MIN_INTRON_NUM);
    err_printf("                                    alignment, or no alignments will be retained. [%.2f]\n", SEC_RATIO);

    err_printf("\n");
    return 1;
}

void add_pathid(bam1_t *b, int id)
{
    char path[20]; int l;
    sprintf(path, ".path%d", id);
    l = strlen(path);

    if (b->l_data + l > b->m_data) { //.path1
        b->m_data =  b->l_data;
        kroundup32(b->m_data);
        b->data = (uint8_t*)_err_realloc(b->data, b->m_data*sizeof(uint8_t));
    }
    memmove(b->data+b->core.l_qname+l-1, b->data+b->core.l_qname-1, b->l_data-b->core.l_qname+1);
    memcpy(b->data+b->core.l_qname-1, path, l);
    b->l_data += l;
    b->core.l_qname += l;
}

int rRNA_overlap(bam1_t *b, read_trans_t *r)
{
    int pos = b->core.pos, tid = b->core.tid;
    int rlen = bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
    int i;
    for (i = 0; i < r->trans_n; ++i) {
        if (tid == r->t[i].tid && !(pos > r->t[i].end || r->t[i].start > pos+rlen-1)) return 1;
        if (tid < r->t[i].tid) return 0;
    }
    return 0;
}

int gtf_filter(bam1_t *b, int *score, int *intron_n, float cov_rate, float map_qual, read_trans_t *r)
{
    if (bam_unmap(b)) return 1;
    uint32_t *c = bam_get_cigar(b); int n_c = b->core.n_cigar;
    // intron number
    int i, del_len=0;
    *intron_n=0;
    for (i = 0; i < n_c; ++i) {
        if (bam_cigar_op(c[i]) == BAM_CREF_SKIP) (*intron_n)++;
        else if (bam_cigar_op(c[i]) == BAM_CDEL) del_len+= bam_cigar_oplen(c[i]);
    }
    // cover len/rate
    int cigar_qlen = b->core.l_qseq;
    int op0 = bam_cigar_op(c[0]), op1 = bam_cigar_op(c[n_c-1]);
    if (op0 == BAM_CSOFT_CLIP || op0 == BAM_CHARD_CLIP) cigar_qlen -= bam_cigar_oplen(c[0]);
    if (n_c > 1 && (op1 == BAM_CSOFT_CLIP || op1 == BAM_CHARD_CLIP)) cigar_qlen -= bam_cigar_oplen(c[n_c-1]);
    if ((cigar_qlen+0.0) / b->core.l_qseq < cov_rate) return 1;
    // NM 
    uint8_t *p = bam_aux_get(b, "NM"); // Edit Distance
    int ed; 
    ed = bam_aux2i(p);
    if ((cigar_qlen - ed + del_len) < map_qual * cigar_qlen) return 1;
    if (rRNA_overlap(b, r)) return 1;
    *score = (cigar_qlen - ed + del_len);
    return 0;
}

const struct option filter_long_opt [] = {
    { "coverage", 1, NULL, 'v' },
    { "map-quality", 1, NULL, 'q' },
    { "sec-rat", 1, NULL, 's' },

    { 0, 0, 0, 0}
};

int bam_filter(int argc, char *argv[])
{
    int c; float cov_rat=COV_RATIO, map_qual = MAP_QUAL, sec_rat=SEC_RATIO; int min_intron_n = MIN_INTRON_NUM;
    int cnt=0;
    while ((c = getopt_long(argc, argv, "v:q:s:i:", filter_long_opt, NULL)) >= 0) {
        switch (c) {
            case 'v': cov_rat = atof(optarg); break;
            case 'q': map_qual = atof(optarg); break;
            case 's': sec_rat = atof(optarg); break;
            case 'i': min_intron_n = atoi(optarg); break;
            default : return filter_usage();
        }
    }

    if (argc - optind != 2) return filter_usage();

    samFile *in, *out; bam_hdr_t *h; bam1_t *b;
    bam1_t *best_b; int b_score=0, s_score=0, score, b_intron_n=0, intron_n;
    if ((in = sam_open(argv[optind], "rb")) == NULL) err_fatal(__func__, "Cannot open \"%s\"\n", argv[optind]);
    if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", argv[optind]);
    b = bam_init1();  best_b = bam_init1();

    // read rRNA gtf
    FILE *fp = xopen(argv[optind+1], "r"); read_trans_t *r = read_trans_init();
    read_anno_trans(fp, h, r);
    err_fclose(fp); 

    if ((out = sam_open_format("-", "wb", NULL)) == NULL) err_fatal_simple("Cannot open \"-\"\n");
    if (sam_hdr_write(out, h) != 0) err_fatal_simple("Error in writing SAM header\n"); //sam header
    char lqname[100]="\0"; int id=1, best_id=1;
    while (sam_read1(in, h, b) >= 0) {
        if (gtf_filter(b, &score, &intron_n, cov_rat, map_qual, r)) continue;

        if (strcmp(bam_get_qname(b), lqname) == 0) {
            id++;
            if (score > b_score) {
                bam_copy1(best_b, b);
                best_id = id;
                s_score = b_score;
                b_score = score;
                b_intron_n = intron_n;
            } else if (score > s_score)
                s_score = score;
        } else { 
            if (strcmp(lqname, "\0") != 0 && s_score < sec_rat * b_score && b_intron_n >= min_intron_n) {
                add_pathid(best_b, best_id);
                if (sam_write1(out, h, best_b) < 0) err_fatal_simple("Error in writing SAM record\n");
                cnt++;
            }
            bam_copy1(best_b, b);
            b_score = score; s_score = 0; b_intron_n=intron_n;
            best_id = id = 1;
            strcpy(lqname, bam_get_qname(b)); 
        }
    }
    if (strcmp(lqname, "\0") != 0 && s_score < sec_rat * b_score && b_intron_n >= min_intron_n) {
        add_pathid(best_b, best_id);
        if (sam_write1(out, h, best_b) < 0) err_fatal_simple("Error in writing SAM record\n");
        cnt++;
    }
    err_func_format_printf(__func__, "Filtered alignments: %d\n", cnt);
    bam_destroy1(b); bam_destroy1(best_b); bam_hdr_destroy(h); sam_close(in); sam_close(out);
    read_trans_free(r);    
    return 0;
}
