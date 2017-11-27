/* bam2gtf.c
 *   generate junction information based on sam/bam file
 *   currently, only for single-end long read data
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "bam2gtf.h"
#include "htslib/sam.h"
#include "utils.h"
#include "gtf.h"

extern const char PROG[20];
int bam2gtf_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s bam2gtf [option] <in.bam> > out.gtf\n\n", PROG);
    err_printf("Options:\n\n");
    err_printf("         -e --exon-min    [INT]    minimum length of internal exon. [%d]\n", INTER_EXON_MIN_LEN);
    err_printf("         -i --intron-len  [INT]    minimum length of intron. [%d]\n", INTRON_MIN_LEN);
    err_printf("         -s --source      [STR]    source field in GTF, program, database or project name. [NONE]\n");
	err_printf("\n");
	return 1;
}

int gen_exon(trans_t *t, bam1_t *b, uint32_t *c, int n_cigar, int exon_min, int intron_len)
{
    t->exon_n = 0;
    int tid = b->core.tid; int start = b->core.pos+1, end = start-1;/*1-base*/ uint8_t is_rev, *p;
    p = bam_aux_get(b, "XS"); // strand orientation for a splice
    if (p == 0) is_rev = bam_is_rev(b);
    else is_rev = ((bam_aux2A(p) == '+' )? 0 : 1);

    int i;
    for (i = 0; i < n_cigar; ++i) {
        int l = bam_cigar_oplen(c[i]);
        switch (bam_cigar_op(c[i])) {
            case BAM_CREF_SKIP:  // N(0 1)
                if (l >= intron_len) {
                    if (t->exon_n == 0 || (end-start+1) >= exon_min)
                        add_exon(t, tid, start, end, is_rev);
                    start = end + l + 1;
                }
                end += l;
                break;
            case BAM_CDEL : // D(0 1)
                end += l;
                break;
            case BAM_CMATCH: // 1 1
            case BAM_CEQUAL:
            case BAM_CDIFF:
                end += l;
                break;
            case BAM_CINS: // 1 0
            case BAM_CSOFT_CLIP:
            case BAM_CHARD_CLIP:
                break;
            case BAM_CPAD: // 0 0
            case BAM_CBACK:
                break;
            default:
                err_printf("Error: unknown cigar type: %d.\n", bam_cigar_op(c[i]));
                break;
        }
    }
    add_exon(t, tid, start, end, is_rev);
    return 0;
}

int gen_trans(bam1_t *b, trans_t *t, int exon_min, int intron_len)
{
    if (bam_unmap(b)) return 0;

    uint32_t *c = bam_get_cigar(b); int n_cigar = b->core.n_cigar;
    gen_exon(t, b, c, n_cigar, exon_min, intron_len);
    return 1;
}

int read_bam_trans(samFile *in, bam_hdr_t *h, bam1_t *b, update_gtf_para *ugp, read_trans_t *T)
{
    trans_t *t = trans_init(1);
    int sam_ret = sam_read1(in, h, b) ;
    while (sam_ret >= 0) {
        gen_trans(b, t, ugp->min_exon, ugp->min_intron); set_trans_name(t, NULL, NULL, NULL, bam_get_qname(b));
        add_read_trans(T, *t); set_trans_name(T->t+T->trans_n-1, NULL, NULL, NULL, bam_get_qname(b));
        // for bam_trans
        T->t[T->trans_n-1].novel_exon_map = (uint8_t*)calloc(t->exon_n, sizeof(uint8_t));
        T->t[T->trans_n-1].novel_sj_map = (uint8_t*)calloc(t->exon_n-1, sizeof(uint8_t));
        //strcpy(T->t[T->trans_n-1].gname, "UNCLASSIFIED");
        T->t[T->trans_n-1].lfull = 0, T->t[T->trans_n-1].lnoth = 1, T->t[T->trans_n-1].rfull = 0, T->t[T->trans_n-1].rnoth = 1;
        T->t[T->trans_n-1].novel = 0, T->t[T->trans_n-1].all_novel=0, T->t[T->trans_n-1].all_iden=0;
        sam_ret = sam_read1(in, h, b) ;
    }
    trans_free(t);
    return T->trans_n;
}

const struct option bam2gtf_long_opt [] = {
    { "exon-min", 1, NULL, 'e' },
    { "intron-len", 1, NULL, 'i' },
    { "source", 1, NULL, 's' },

    { 0, 0, 0, 0}
};

int bam2gtf(int argc, char *argv[])
{
    int c, exon_min=INTER_EXON_MIN_LEN, intron_len=INTRON_MIN_LEN;
    char src[100]="NONE";
	while ((c = getopt_long(argc, argv, "s:e:i:", bam2gtf_long_opt, NULL)) >= 0)
    {
        switch(c)
        {
            case 'e': exon_min = atoi(optarg); break;
            case 's': strcpy(src, optarg); break;
            case 'i': intron_len = atoi(optarg); break;
            default: err_printf("Error: unknown option: %s.\n", optarg);
                     return bam2gtf_usage();
        }
    }
    if (argc - optind != 1) return bam2gtf_usage();

    samFile *in; bam_hdr_t *h; bam1_t *b;

    in = sam_open(argv[optind], "rb");
    if (in == NULL) err_fatal(__func__, "Cannot open \"%s\"\n", argv[optind]);
    h = sam_hdr_read(in);
    if (h == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", argv[optind]);
    b = bam_init1();

    trans_t *t = trans_init(1);

    while (sam_read1(in, h, b) >= 0) {
        if (gen_trans(b, t, exon_min, intron_len)) set_trans_name(t, NULL, NULL, NULL, bam_get_qname(b));
        print_trans(*t, h, src, stdout);
    }

    trans_free(t);
    bam_destroy1(b); bam_hdr_destroy(h); sam_close(in);
    return 0;
}
