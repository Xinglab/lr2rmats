#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "utils.h"
#include "../htslib/htslib/sam.h"

extern char PROG[25];

int fusion_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s fusion [option] <in.bam/sam> | samtools sort > out.sort.bam\n\n", PROG);
    err_printf("Options:\n");
    err_printf("         -v --coverage   [FLOAT]    minimum fraction of aligned bases. [%.2f]\n", COV_RATIO);
    err_printf("         -q --map-qual   [FLOAT]    minimum fraction of identically aligned bases. [%.2f]\n", MAP_QUAL);
    err_printf("         -s --sec-rat    [FLOAT]    maximum ratio of second best and best score to retain the best\n");
    err_printf("                                    alignment, or no alignments will be retained. [%.2f]\n", SEC_RATIO);
    err_printf("         -i --intron     [INT]      minimum number of intron indicated by the alignment. [%d]\n", MIN_INTRON_NUM);
    err_printf("         -r --remove-gtf [STR]      remove all the alignment record that overlap with transcript in this GTF file. [NONE]\n");

    err_printf("\n");
    return 1;
}
const struct option fusion_long_opt [] = {
    { "coverage", 1, NULL, 'v' },
    { "map-quality", 1, NULL, 'q' },
    { "sec-rat", 1, NULL, 's' },
    { "intron", 1, NULL, 'i' },
    { "remove-gtf", 1, NULL, 'r' },

    { 0, 0, 0, 0}
};
int bam_fusion(int argc, char *argv[])
{
    int c; 
    int cnt=0;
    while ((c = getopt_long(argc, argv, "v:q:s:i:r:", fusion_long_opt, NULL)) >= 0) {
        switch (c) {
            case 'v': cov_rat = atof(optarg); break;
            case 'q': map_qual = atof(optarg); break;
            case 's': sec_rat = atof(optarg); break;
            case 'i': min_intron_n = atoi(optarg); break;
            case 'r': strcpy(remove_fn, optarg); break;
            default : return fusion_usage();
        }
    }

    if (argc - optind != 1) return fusion_usage();

    samFile *in, *out; bam_hdr_t *h; bam1_t *b;
    bam1_t *best_b; int b_score=0, s_score=0, score, b_intron_n=0, intron_n;
    if ((in = sam_open(argv[optind], "rb")) == NULL) err_fatal(__func__, "Cannot open \"%s\"\n", argv[optind]);
    if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", argv[optind]);
    b = bam_init1();  best_b = bam_init1();

    // read remove gtf
    read_trans_t *r = read_trans_init(1); 
    if (strcmp(remove_fn, "") != 0) {
        FILE *fp = xopen(remove_fn, "r"); 
        read_anno_trans(fp, h, r);
        err_fclose(fp);
    }

    if ((out = sam_open_format("-", "wb", NULL)) == NULL) err_fatal_simple("Cannot open \"-\"\n");
    if (sam_hdr_write(out, h) != 0) err_fatal_simple("Error in writing SAM header\n"); //sam header
    char lqname[100]="\0"; int id=1, best_id=1;
    while (sam_read1(in, h, b) >= 0) {
        if (gtf_fusion(b, &score, &intron_n, cov_rat, map_qual, r)) continue;

        if (strcmp(bam_get_qname(b), lqname) == 0) {
            id++;
            if (score > b_score) {
                bam_copy1(best_b, b);
                //best_id = id;
                s_score = b_score;
                b_score = score;
                b_intron_n = intron_n;
            } else if (score > s_score)
                s_score = score;
        } else {
            if (strcmp(lqname, "\0") != 0 && s_score < sec_rat * b_score && b_intron_n >= min_intron_n) {
                //add_pathid(best_b, best_id);
                if (sam_write1(out, h, best_b) < 0) err_fatal_simple("Error in writing SAM record\n");
                cnt++;
            }
            bam_copy1(best_b, b);
            b_score = score; s_score = 0; b_intron_n=intron_n;
            //best_id = id = 1;
            strcpy(lqname, bam_get_qname(b));
        }
    }
    if (strcmp(lqname, "\0") != 0 && s_score < sec_rat * b_score && b_intron_n >= min_intron_n) {
        //add_pathid(best_b, best_id);
        if (sam_write1(out, h, best_b) < 0) err_fatal_simple("Error in writing SAM record\n");
        cnt++;
    }
    err_func_format_printf(__func__, "Filtered alignments: %d\n", cnt);
    bam_destroy1(b); bam_destroy1(best_b); bam_hdr_destroy(h); sam_close(in); sam_close(out);
    read_trans_free(r);
    return 0;
}

