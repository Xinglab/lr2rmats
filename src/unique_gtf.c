#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "unique_gtf.h"
#include "update_gtf.h"
#include "utils.h"
#include "bam2gtf.h"
#include "gtf.h"
#include "htslib/sam.h"

extern const char PROG[20];

unique_gtf_para *unique_gtf_init_para(void) {
    unique_gtf_para *ugp = (unique_gtf_para*)_err_malloc(sizeof(unique_gtf_para));
    ugp->input_mode = 0/*bam*/; ugp->gtf_bam = NULL;
    ugp->force_strand = 0;
    ugp->min_exon = INTER_EXON_MIN_LEN, ugp->min_intron = INTRON_MIN_LEN, ugp->ss_dis = SPLICE_DISTANCE; ugp->end_dis = END_DISTANCE; ugp->deletion_max = DELETION_MAX_LEN;
    ugp->single_exon_ovlp_frac = SING_OVLP_FRAC;
    ugp->out_gtf_fp = stdout;
    strcpy(ugp->source, PROG);
    return ugp;
}

int unique_gtf_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s unique-gtf [option] <in.sorted.bam/in.sorted.gtf> > unique.gtf\n\n", PROG);
    err_printf("Notice:  the BAM and GTF files should be sorted in advance.\n\n");
    err_printf("Input options:\n\n");
    err_printf("         -m --input-mode  [STR]    format of input file <in.bam/in.gtf>, BAM file(b) or GTF file(g). [b]\n");
    err_printf("         -b --bam         [STR]    for GTF input <in.gtf>, BAM file is needed to obtain BAM header information. [NULL]\n");
    err_printf("\n");

    err_printf("Function options:\n\n");
    err_printf("         -s --force-strand         force to match strand when merging transcripts. [False]\n");
    err_printf("         -e --min-exon    [INT]    minimum length of internal exon. [%d]\n", INTER_EXON_MIN_LEN);
    err_printf("         -i --min-intron  [INT]    minimum length of intron. [%d]\n", INTRON_MIN_LEN);
    err_printf("         -t --max-delet   [INT]    maximum length of deletion, longer deletion will be considered as intron. [%d]\n", DELETION_MAX_LEN);
    err_printf("         -d --distance    [INT]    consider same if distance between two splice site is not bigger than d. [%d]\n", SPLICE_DISTANCE);
    err_printf("         -D --DISTANCE    [INT]    consider same if distance between two start/end site is not bigger than D. [%d]\n", END_DISTANCE);
    err_printf("         -f --frac        [INT]    consider same if overlapping between two single-exon transcript is bigger than f. [%.2f]\n", SING_OVLP_FRAC);
    err_printf("\n");

    err_printf("Output options:\n\n");
    err_printf("         -o --output      [STR]    unique GTF file. [stdout]\n");
    err_printf("         -S --source      [STR]    \'source\' field in GTF: program, database or project name. [%s]\n", PROG);
    err_printf("\n");
    return 1;
}

const struct option unique_gtf_long_opt [] = {
    { "input-mode", 1, NULL, 'm' },
    { "bam", 1, NULL, 'b' },

    { "force-strand", 0, NULL, 's' },
    { "min-exon", 1, NULL, 'e' },
    { "min-intron", 1, NULL, 'i' },
    { "distance", 1, NULL, 'd' },
    { "DISTANCE", 1, NULL, 'D' },

    { "frac", 1, NULL, 'f' },

    { "output", 1, NULL, 'o' }, 
    { "source", 1, NULL, 's' },


    { 0, 0, 0, 0}
};

int uniq_trans(read_trans_t *bam_T, read_trans_t *uniq_T, unique_gtf_para *ugp) {
    int i;
    for (i = 0; i < bam_T->trans_n; ++i) {
        trans_t *t = bam_T->t + i;
        if (merge_trans(t, uniq_T, ugp->force_strand, ugp->ss_dis, ugp->end_dis, ugp->single_exon_ovlp_frac) == 0)
            add_read_trans(uniq_T, *t);
    }

    return uniq_T->trans_n;
}

int unique_gtf(int argc, char *argv[])
{
    int c; 
    unique_gtf_para *ugp = unique_gtf_init_para();
    while ((c = getopt_long(argc, argv, "m:b:se:i:d:D:f:o:S:", unique_gtf_long_opt, NULL)) >= 0) {
        switch(c)
        {
            case 'm': if (optarg[0] == 'b') ugp->input_mode=0; 
                          else if (optarg[0] == 'g') ugp->input_mode=1; 
                          else return unique_gtf_usage(); 
                          break;
            case 'b': if ((ugp->gtf_bam = sam_open(optarg, "rb")) == NULL) {
                          err_fatal(__func__, "Cannot open \"%s\"\n", optarg); 
                          return unique_gtf_usage();
                      }
                      break;
            case 's': ugp->force_strand = 1; break;
            case 'e': ugp->min_exon = atoi(optarg); break;
            case 'i': ugp->min_intron = atoi(optarg); break;
            case 't': ugp->deletion_max = atoi(optarg); break;
            case 'd': ugp->ss_dis = atoi(optarg); break;
            case 'D': ugp->end_dis = atoi(optarg); break;
            case 'f': ugp->single_exon_ovlp_frac = atof(optarg); break;

            case 'o': ugp->out_gtf_fp = fopen(optarg, "w"); break;
            case 'S': strcpy(ugp->source, optarg); break;
            default:
                      err_printf("Error: unknown option: %s.\n", optarg);
                      return unique_gtf_usage();
                      break;
        }
    }
    if (argc - optind != 1) return unique_gtf_usage();

    chr_name_t *cname = chr_name_init();
    read_trans_t *bam_T, *unique_T;
    bam_T = read_trans_init(1);
    unique_T = read_trans_init(1); 

    // read all input-transcript
    samFile *in = NULL; bam_hdr_t *h; 
    if (ugp->input_mode == 0) { // bam input
        bam1_t *b;
        if ((in = sam_open(argv[optind], "rb")) == NULL) err_fatal(__func__, "Can not open \"%s\"\n", argv[optind]);
        if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", argv[optind]);
        bam_set_cname(h, cname);
        b = bam_init1(); 
        read_bam_trans(in, h, b, ugp->min_exon, ugp->min_intron, ugp->deletion_max, bam_T);
        bam_destroy1(b);
    } else { // gtf input
        if ((h = sam_hdr_read(ugp->gtf_bam)) == NULL) err_fatal(__func__, "Couldn't read header of provided BAM file.\n");
        bam_set_cname(h, cname);
        FILE *fp = xopen(argv[optind], "r");
        read_gtf_trans(fp, h, bam_T);
    }

    // identify novel transcript
    uniq_trans(bam_T, unique_T, ugp);

    // print novel transcript
    print_read_trans(unique_T, cname, ugp->source, ugp->out_gtf_fp);

    chr_name_free(cname);
    read_trans_free(bam_T); read_trans_free(unique_T);
    err_fclose(ugp->out_gtf_fp); 
    bam_hdr_destroy(h); 
    if (in) sam_close(in); 
    if (ugp->gtf_bam) sam_close(ugp->gtf_bam);
    free(ugp);
    return 0;
}
