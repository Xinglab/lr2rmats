#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <getopt.h>
#include <string.h>
#include "parse_bam.h"
#include "utils.h"
#include "bam_fusion.h"
#include "../htslib/htslib/sam.h"

extern char PROG[25];

bam_fusion_para *bam_fusion_init_para(void) {
    bam_fusion_para *bfp = (bam_fusion_para*)_err_malloc(sizeof(bam_fusion_para));
    bfp->ovlp_frac = OVLP_FRAC;
    bfp->each_cov = EACH_COV;
    bfp->all_cov = ALL_COV;
    bfp->dis = FUSION_DIS;

    bfp->fs_fp = NULL;

    return bfp;
}

int fusion_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s fusion [option] <in.bam/sam> > fusion.sam\n", PROG);
    err_printf("     or: %s fusion [option] <in.bam/sam> | bedtools bamtobed -i stdin -bed12 > fusion.bed\n\n", PROG);
    err_printf("Options:\n");
    err_printf("         -o --ovlp-frac   [FLOAT]    maximum overlap fraction of each fusion part. [%.2f]\n", OVLP_FRAC);
    err_printf("         -v --each-cov    [FLOAT]    minimum fraction of each fusion part. [%.2f]\n", EACH_COV);
    err_printf("         -V --all-cov     [FLOAT]    minimum fraction of all mapped parts. [%.2f]\n", ALL_COV);
    err_printf("         -d --dis         [INT]      minimum distance of two fusion parts. [%s]\n", FUSION_DIS_STR);
    err_printf("         -f --fusion-site [STR]      output fusion site file. [NULL]\n");
    // TODO extract gene name/id from gtf
    err_printf("         -g --gtf         [STR]      gene annotation in GTF format. [None]\n");
    err_printf("\n");
    return 1;
}

const struct option fusion_long_opt [] = {
    { "ovlp-frac", 1, NULL, 'o' },
    { "each-cov", 1, NULL, 'v' },
    { "all-cov", 1, NULL, 'V' },
    { "dis", 1, NULL, 'd' },

    { 0, 0, 0, 0}
};

int str2int(char *str)
{
    double d; char *p;
    d = strtod(str, &p);
    if (*p == 'G' || *p == 'g') d *= 1e9;
    else if (*p == 'M' || *p == 'm') d *= 1e6;
    else if (*p == 'K' || *p == 'k') d *= 1e3;
    return (int)(d + 0.499); 
}

int seg_cmpfunc (const void * a, const void * b) {
    if (((bam_seg_t*)a)->score != ((bam_seg_t*)b)->score) {
        return ((bam_seg_t*)b)->score - ((bam_seg_t*)a)->score;
    } else return ((bam_seg_t*)a)->ed - ((bam_seg_t*)b)->ed;
}

float ovlp_rat(int start1, int end1, int start2, int end2) {
    if (start1 > end2 || start2 > end1) return 0.0;
    int overlap_len = end1 - start2 + 1 > 0 ? end1 - start2 + 1 : end2 - start1 + 1;
    int min_len = MIN_OF_TWO(end1 - start1 + 1, end2 - start2 + 1);
    return (overlap_len / (min_len + 0.0));
}

int check_with_exist1(bam_seg_t* s1, bam_seg_t* s2, bam_fusion_para *bfp) {
    // at most ovlp_frac overlap
    if (ovlp_rat(s1->read_start, s1->read_end, s2->read_start, s2->read_end) > bfp->ovlp_frac) return 0;
    // mapping distance >= dis
    if (s1->tid == s2->tid) { // && s1->is_rev == s2->is_rev) {
        if (ovlp_rat(s1->ref_start, s1->ref_end, s2->ref_start, s2->ref_end) > 0.0) {
            return 0;
        } else {
            if (s1->ref_start - s2->ref_end > 0 && s1->ref_start - s2->ref_end < bfp->dis) return 0;
            if (s2->ref_start - s1->ref_end > 0 && s2->ref_start - s1->ref_end < bfp->dis) return 0;
        }
    }
    return 1;
}

int check_with_exist(bam_seg_t *seg, int fusion_seg_n, int seg_i, bam_fusion_para *bfp) {
    int i;
    for (i = 0; i < fusion_seg_n; ++i) {
        if (check_with_exist1(seg+i, seg+seg_i, bfp) == 0) return 0;
    }
    return 1;
}


float bam_seg_cov(bam_seg_t *seg, int n, int rlen) {
    uint8_t *map = (uint8_t*)_err_calloc(rlen, sizeof(uint8_t));
    int i, j;
    for (i = 0; i < n; ++i) {
        for (j = seg[i].read_start; j <= seg[i].read_end; ++j) {
            map[j-1] = 1;
        }
    }
    int cov_n = 0;
    for (i = 0; i < rlen; ++i) {
        cov_n += map[i];
    }
    free(map);
    return (cov_n + 0.0) / rlen;
}

int check_fusion(bam_seg_t *seg, int rlen, int seg_n, bam_fusion_para *bfp) {
    qsort(seg, seg_n, sizeof(bam_seg_t), seg_cmpfunc);
    int fusion_seg_n = 1, i;
    for (i = 1; i < seg_n; ++i) {
        // at least each_cov of read
        if ((seg[i].read_end - seg[i].read_start + 1) / (rlen+0.0) < bfp->each_cov) continue;
        if (check_with_exist(seg, fusion_seg_n, i, bfp)) {
            move_bam_seg(seg, fusion_seg_n, i);
            fusion_seg_n++;
            if (bam_seg_cov(seg, fusion_seg_n, rlen) >= bfp->all_cov) {
                return fusion_seg_n;
            }
        }
    }
    return -1;
}

// seg_n == 2
void fusion_write(FILE *out, bam_seg_t *seg, int seg_n, char **target_name) {
    if (seg_n != 2) return;
    int left_i, right_i;
    if (seg[0].read_start < seg[1].read_start) {
        left_i = 0, right_i = 1;
    } else {
        left_i = 1, right_i = 0;
    }

    fprintf(out, "%s\t%s\t%c\t%d\t%d\t%s\t%c\t%d\t%d\n", bam_get_qname(seg[left_i].b), target_name[seg[left_i].tid], "+-"[seg[left_i].is_rev], seg[left_i].ref_start, seg[left_i].ref_end, target_name[seg[right_i].tid], "+-"[seg[right_i].is_rev], seg[right_i].ref_start, seg[right_i].ref_end);
}

int bam_fusion(int argc, char *argv[])
{
    int c;
    bam_fusion_para *bfp = bam_fusion_init_para();
    while ((c = getopt_long(argc, argv, "o:v:V:d:f:", fusion_long_opt, NULL)) >= 0) {
        switch (c) {
            case 'o': bfp->ovlp_frac = atof(optarg); break;
            case 'v': bfp->each_cov = atof(optarg); break;
            case 'V': bfp->all_cov = atof(optarg); break;
            case 's': bfp->dis = str2int(optarg); break;
            case 'f': bfp->fs_fp = xopen(optarg, "w"); break;
            default : return fusion_usage();
        }
    }

    if (argc - optind != 1) return fusion_usage();

    samFile *in, *out; bam_hdr_t *h; bam1_t *b;
    if ((in = sam_open(argv[optind], "rb")) == NULL) err_fatal(__func__, "Cannot open \"%s\"\n", argv[optind]);
    if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", argv[optind]);
    if ((out = sam_open_format("-", "wb", NULL)) == NULL) err_fatal_simple("Cannot open \"-\"\n");
    if (sam_hdr_write(out, h) != 0) err_fatal_simple("Error in writing SAM header\n");

    b = bam_init1(); 

    char lqname[1024] = "";
    bam_seg_t *s = bam_seg_init(1), *seg = bam_seg_init(2);
    int seg_m = 2, seg_n = 0, i, rlen, cnt = 0;

    if (bfp->fs_fp) fprintf(bfp->fs_fp, "#fusion_id\t1st_chr\t1st_strand\tst_start_site\t1st_end_site\t2nd_chr\t2nd_strand\t2nd_start_site\t2nd_end_site\n");

    while (sam_read1(in, h, b) >= 0) {
        if (bam2seg(b, s) < 0) continue;

        if (strcmp(bam_get_qname(b), lqname) == 0) {
            seg = push_bam_seg(seg, &seg_n, &seg_m, s);
        } else {
            if (strcmp(lqname, "\0") != 0 && seg_n >= 2) {
                seg_n = check_fusion(seg, rlen, seg_n, bfp);
                if (seg_n == 2) {
                    for (i = 0; i < seg_n; ++i) {
                        if (sam_write1(out, h, seg[i].b) < 0) err_fatal_simple("Error in writing SAM record\n");
                    }
                    if (bfp->fs_fp) fusion_write(bfp->fs_fp, seg, 2, h->target_name);
                    cnt++;
                }
            }
            strcpy(lqname, bam_get_qname(b));
            seg_n = 0; rlen = bam_query_len(b);
            seg = push_bam_seg(seg, &seg_n, &seg_m, s);
        }
    }
    if (strcmp(lqname, "\0") != 0 && seg_n >= 2) {
        seg_n = check_fusion(seg, rlen, seg_n, bfp);
        if (seg_n == 2) {
            for (i = 0; i < seg_n; ++i) {
                if (sam_write1(out, h, seg[i].b) < 0) err_fatal_simple("Error in writing SAM record\n");
            }
            cnt++;
        }
    }
    err_func_format_printf(__func__, "Candidate gene-fusion transcripts: %d\n", cnt);
    bam_destroy1(b); bam_hdr_destroy(h); sam_close(in); sam_close(out);
    bam_seg_free(s, 1); bam_seg_free(seg, seg_m);
    if (bfp->fs_fp) err_fclose(bfp->fs_fp);
    free(bfp);

    return 0;
}
