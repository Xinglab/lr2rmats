#include <stdio.h>
#include <stdlib.h>
#include "../htslib/htslib/sam.h"
#include "utils.h"
#include "gtf.h"
#include "parse_bam.h"

extern const char PROG[20];

int update_gtf_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s sr-valid [option] <long.bam/long.gtf> <short.sj/bam> > out.bam/gtf\n\n", PROG);
    err_printf("Notice:  the BAM/GTF files should be sorted in advance.\n\n");
    err_printf("Input options:\n\n");
    err_printf("         -m --input-mode   [STR]    format of input file <long.bam/gtf>, BAM file(b) or GTF file(g). [b]\n");
    err_printf("         -M --output-mode  [STR]    format of output file <out.bam/gtf>, BAM file(b) or GTF file(g). [b]\n");
    err_printf("         -j --input-sj     [STR]    format of input splice-junction file <short.sj/bam>, BAM file(b) or STAR.SJ.out.tab(t). [t]\n");
    err_printf("         -b --bam          [STR]    for GTF input <in.gtf>, BAM file is needed to obtain BAM header information. [NULL]\n");
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
    //err_printf("         -J --min-junc-cnt [INT]    minimum short-read junction count of novel junction. [%d]\n", MIN_SJ_CNT);
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


