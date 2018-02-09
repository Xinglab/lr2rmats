#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "bam_filter.h"
#include "bam_fusion.h"
#include "update_gtf.h"
#include "bam2gtf.h"
#include "parse_bam.h"

const char PROG[20] = "lr2rmats";
const char VERSION[20] = "1.0.0";
const char DATE[20] = "2017-11-22";
const char CONTACT[30] = "yangaoucla@gmail.com";


static int usage(void)
{
    err_printf("\n");
	err_printf("Program: %s (long read to GTF)\n", PROG);
    err_printf("Version: %s, Date: %s\n", VERSION, DATE);
    err_printf("Contact: %s\n", CONTACT);
    err_printf("Usage:   %s <command> [options]\n\n", PROG);
	err_printf("Commands: \n");
    err_printf("         filter       filter out alignment records with low confidence\n");
    err_printf("         fusion       generate candidate gene-fusion transcripts\n");
	err_printf("         update-gtf   generate new GTF file based on BAM/SAM and existing GTF file\n");
	err_printf("         bam2gtf      generate transcript and exon information based on BAM/SAM file\n");
	err_printf("         bam2sj       generate splice-junction information based on BAM/SAM file\n");
	err_printf("\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage();

    if (strcmp(argv[1], "filter") == 0) return bam_filter(argc-1, argv+1);
    else if (strcmp(argv[1], "fusion") == 0) return bam_fusion(argc-1, argv+1);
	else if (strcmp(argv[1], "update-gtf") == 0) return update_gtf(argc-1, argv+1);
	else if (strcmp(argv[1], "bam2gtf") == 0) return bam2gtf(argc-1, argv+1);
    else if (strcmp(argv[1], "bam2sj") == 0) return bam2sj(argc-1, argv+1);
	else { fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]); return 1; }
    return 0;
}
