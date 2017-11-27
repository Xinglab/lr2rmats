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

#define bam_unmap(b) ((b)->core.flag & BAM_FUNMAP)

extern const char PROG[20];

update_gtf_para *update_gtf_init_para(void) {
    update_gtf_para *ugp = (update_gtf_para*)_err_malloc(sizeof(update_gtf_para));
    ugp->input_mode = 0/*bam*/, ugp->full_len_level = 5/*most relax*/, ugp->uncla = 0, ugp->only_bam = 0;
    ugp->intron_fp = NULL, ugp->out_gtf_fp = stdout; strcpy(ugp->source, PROG);
    ugp->min_exon = INTER_EXON_MIN_LEN, ugp->min_intron = INTRON_MIN_LEN, ugp->ss_dis = SPLICE_DISTANCE;

    return ugp;
}

int update_gtf_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s update-gtf [option] <in.bam/in.gtf> <old.gtf> > new.gtf\n\n", PROG);
    err_printf("Notice:  the BAM and GTF files should be sorted in advance.\n\n");
    err_printf("Options:\n\n");
    err_printf("         -m --input-mode  [STR]    format of input file <in.bam/in.gtf>, BAM file(b) or GTF file(g). [b]\n");
    err_printf("         -b --bam         [STR]    for GTF input <in.gtf>, BAM file is needed to obtain BAM header information. [NULL]\n");
    err_printf("         -I --intron      [STR]    intron information file output by STAR(*.out.tab). [NONE]\n");
    err_printf("         -e --min-exon    [INT]    minimum length of internal exon. [%d]\n", INTER_EXON_MIN_LEN);
    err_printf("         -i --intron-len  [INT]    minimum length of intron. [%d]\n", INTRON_MIN_LEN);
    err_printf("         -d --distance    [INT]    consider same if distance between two splice site is not bigger than d. [%d]\n", SPLICE_DISTANCE);
    err_printf("         -l --full-length [INT]    level of strict criterion for considering full-length transcript. \n");
    err_printf("                                   (1->5, most strict->most relaxed) [%d]\n", 5);

    err_printf("         -u --unclassified         output UNCLASSIFIED novel transcript. [False]\n");
    err_printf("         -s --source      [STR]    source field in GTF, program, database or project name. [%s]\n", PROG);
    err_printf("         -n --only-bam             only output bam-derived transcript. [False]\n");
    err_printf("         -o --output               output GTF file. [stdout]\n");
    err_printf("\n");
    return 1;
}

void check_exon_junction(trans_t *bam_t, trans_t anno_t, int dis, int *iden_n, int *iden_intron_n, int *not_iden_iden, int *intron_map)
{
    int i,j, left=0,right=0, last_j=-1;
    *iden_n=0, *iden_intron_n = 0, *not_iden_iden=0;
    for (i = 0; i < bam_t->exon_n; ++i) {
        for (j = 0; j < anno_t.exon_n; ++j) {
            // for exon
            if (abs(bam_t->exon[i].start-anno_t.exon[j].start) <= dis || i==0) left = 1;
            if (abs(bam_t->exon[i].end - anno_t.exon[j].end) <= dis || i==bam_t->exon_n-1) right = 1;
            if (left) set_l_iden(bam_t->novel_exon_map[i]);
            if (right) set_r_iden(bam_t->novel_exon_map[i]);
            if (left && right) set_b_iden(bam_t->novel_exon_map[i]);
            left = right = 0;
            // for junction
            if (i != bam_t->exon_n-1 && j != anno_t.exon_n-1) {
                if (abs(bam_t->exon[i].end-anno_t.exon[j].end) <= dis) left = 1;
                if (abs(bam_t->exon[i+1].start-anno_t.exon[j+1].start) <= dis) right= 1;
                if (left) set_l_iden(bam_t->novel_sj_map[i]);
                if (right) set_r_iden(bam_t->novel_sj_map[i]);
                if (left && right) {
                    set_b_iden(bam_t->novel_sj_map[i]);

                    if (last_j != -1 && j != last_j+1) *not_iden_iden = 1;
                    last_j = j; *iden_intron_n += 1;
                    if (intron_map) intron_map[i] = 1;
                }
                *iden_n += (left+right);

                left = right = 0;
                if (anno_t.exon[j+1].start > bam_t->exon[i+1].start) break;
            }
        }
    }
}

int check_short_sj1(int tid, int start, int end, uint8_t is_rev, intron_group_t *I, int i_start, int dis)
{
    int i = i_start;
    while (i < I->intron_n) {
        if (I->intron[i].tid > tid || I->intron[i].start > start) return 0;
        if (I->intron[i].is_rev != is_rev) { i++; continue; }

        if (abs(I->intron[i].start-start)<=dis && abs(I->intron[i].end-end)<=dis) return 1;
        else i++;
    }
    return 0;
}

int check_short_sj(trans_t *bam_t, int *intron_map, intron_group_t *I, int *intron_i, int dis)
{
    int i = *intron_i, j;
    while (i < I->intron_n) {
        if (I->intron[i].tid < bam_t->tid || (I->intron[i].tid == bam_t->tid && I->intron[i].end <= bam_t->start)) {
            i++; *intron_i = i;
        } else if (I->intron[i].tid > bam_t->tid) return 0;
        else {
            for (j = 0; j < bam_t->exon_n-1; ++j) {
                if (intron_map[j] == 0 && check_short_sj1(bam_t->tid, bam_t->exon[j].end+1, bam_t->exon[j+1].start-1, bam_t->is_rev, I, i, dis) == 0)
                    return 0;
            }
            return 1;
        }
    }
    return 0;
}

int exon_overlap(exon_t e1, exon_t e2)
{
    if (e1.start > e2.end || e2.start > e1.end) return 0;
    return 1;
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

int check_novel_intron(trans_t *bam_t, trans_t anno_t, intron_group_t *I, int *intron_i, int dis, int l)
{
    if (bam_t->is_rev != anno_t.is_rev || bam_t->exon_n < 2) return 3;
    // check full-length
    check_full(bam_t, anno_t, l);

    // check novel
    int iden_n=0, iden_intron_n = 0, not_iden_iden=0;
    int *intron_map = (int*)_err_calloc((bam_t->exon_n-1), sizeof(int));
    check_exon_junction(bam_t, anno_t, dis, &iden_n, &iden_intron_n, &not_iden_iden, intron_map);

    // analyse check result
    if (iden_intron_n == bam_t->exon_n-1 && not_iden_iden == 0) bam_t->all_iden=1;
    else {
        if (iden_n > 0) {
            if (check_short_sj(bam_t, intron_map, I, intron_i, dis)) {
                strcpy(bam_t->gname, anno_t.gname), strcpy(bam_t->gid, anno_t.gid);
                bam_t->novel = 1;
            }
        } else {
            bam_t->all_novel=1;
        }
    }
    free(intron_map);
    return 0;
}

// check if t is novel and has identical splice site
// if t has all identical splice sites with other novel t, merge two ends
// @return value
//    0: novel, NOT share any identical splice site (unclassified)
//    1: novel, and share identical splice site (gene_id)
//    2: totally identical, can NOT be added to any anno
//    3: other cases that cannot be added to this anno(not full-length to any anno-trans)
int check_novel1(trans_t *bam_t, trans_t anno_t, int dis, int l)
{
    if (bam_t->is_rev != anno_t.is_rev || bam_t->exon_n < 2) return 0;
    // check full-length
    check_full(bam_t, anno_t, l);

    // check novel
    int iden_n=0, iden_intron_n=0, not_iden_iden=0;
    check_exon_junction(bam_t, anno_t, dis, &iden_n, &iden_intron_n, &not_iden_iden, NULL);

    // analyse check result
    if (iden_intron_n == bam_t->exon_n-1 && not_iden_iden == 0) bam_t->all_iden=1;
    else if (iden_n > 0) {
        strcpy(bam_t->gname, anno_t.gname), strcpy(bam_t->gid, anno_t.gid);
        bam_t->novel = 1;
    } else {
        bam_t->all_novel = 1;
    }
    return 0;
}

// merge contained and identical trans
int merge_trans1(trans_t *t1, trans_t *t2, int dis)
{
    int i = t1->exon_n-1, j = t2->exon_n-1;
    if (check_iden(t1, t2, dis)) {
        t2->cov++;

        if (t1->exon[0].start < t2->exon[0].start)  {
            t2->exon[0].start = t1->exon[0].start;
            t2->start = t1->exon[0].start;
        }
        if (t1->exon[i].end > t2->exon[j].end) {
            t2->exon[j].end = t1->exon[i].end;
            t2->end = t1->exon[i].end;
        }
        return 1;
    } else if (check_sub_iden(t1, t2, dis)) return 1;
    else return 0;
}

int merge_trans(trans_t *t, read_trans_t *T, int dis)
{
    int i; 
    for (i = T->trans_n-1; i >= 0; --i) {
        if (merge_trans1(t, T->t+i, dis)) return 1;
        if (t->tid > T->t[i].tid || t->start > T->t[i].end) return 0;
    }
    return 0;
}

int check_novel_trans(read_trans_t *bam_T, read_trans_t *anno_T, intron_group_t *I, read_trans_t *novel_T, update_gtf_para *ugp)
{
    int i=0, j=0, last_j=0, k=0;
    int unclassify = 0; char uncla[1024];
    while (i < bam_T->trans_n && j < anno_T->trans_n) {
        // check if redundant
        if (merge_trans(bam_T->t+i, novel_T, ugp->ss_dis)) { 
            //err_printf("merge: %s\n", bam_T->t[i].tname);
            i++; continue; 
        } else if (bam_T->t[i].tid > anno_T->t[j].tid || (bam_T->t[i].tid == anno_T->t[j].tid && bam_T->t[i].start > anno_T->t[j].end)) {
            if (j == last_j) last_j++;
            j++;
        } // compare with all anno trans done
        else if (anno_T->t[j].tid > bam_T->t[i].tid || (anno_T->t[j].tid == bam_T->t[i].tid && anno_T->t[j].start > bam_T->t[i].end)) {
            //if (bam_T->t[i].novel != 1) err_printf("unrecover: %s\n", bam_T->t[i].tname);
            set_full(bam_T->t+i, ugp->full_len_level);
            if (bam_T->t[i].full && bam_T->t[i].novel) {
                if (bam_T->t[i].all_iden) {
                    err_printf("Error: all iden %s.\n", bam_T->t[i].tname);
                    exit(0);
                } else if (bam_T->t[i].all_novel && ugp->uncla) {
                    add_read_trans(novel_T, bam_T->t[i]);
                    sprintf(uncla, "UNCLA_%d", unclassify++);
                    set_trans_name(novel_T->t+novel_T->trans_n-1, NULL, uncla, NULL, NULL);
                } else if (bam_T->t[i].all_novel == 0) {
                    add_read_trans(novel_T, bam_T->t[i]);
                    set_trans_name(novel_T->t+novel_T->trans_n-1, NULL, NULL, NULL, NULL);
                }
            }

            i++; j = last_j;
        } else {
            if (I->intron_n > 0) check_novel_intron(bam_T->t+i, anno_T->t[j], I, &k, ugp->ss_dis, ugp->full_len_level);
            else check_novel1(bam_T->t+i, anno_T->t[j], ugp->ss_dis, ugp->full_len_level);

            if (bam_T->t[i].all_iden == 1) {
                //err_printf("all-iden: %s\n", bam_T->t[i].tname);
                i++; j = last_j;
            } else {
                j++;
            }
        }
    }
    if (i < bam_T->trans_n) {
        set_full(bam_T->t+i, ugp->full_len_level);
        if (bam_T->t[i].full && bam_T->t[i].novel) {
            if (bam_T->t[i].all_iden) {
                err_printf("Error: all iden %s.\n", bam_T->t[i].tname);
                exit(0);
            } else if (bam_T->t[i].all_novel && ugp->uncla) {
                add_read_trans(novel_T, bam_T->t[i]);
                sprintf(uncla, "UNCLA_%d", unclassify++);
                set_trans_name(novel_T->t+novel_T->trans_n-1, NULL, uncla, NULL, NULL);
            } else if (bam_T->t[i].all_novel == 0) {
                add_read_trans(novel_T, bam_T->t[i]);
                set_trans_name(novel_T->t+novel_T->trans_n-1, NULL, NULL, NULL, NULL);
            }
        }
    }
    return 0;
}

// from annotation gtf file extract transcript-exon structure
int read_anno_trans(FILE *fp, bam_hdr_t *h, read_trans_t *T)
{
    char line[1024], ref[100]="\0", type[20]="\0"; int start, end; char strand, add_info[1024], gname[100];
    trans_t *t = trans_init(1);
    while (fgets(line, 1024, fp) != NULL) {
        sscanf(line, "%s\t%*s\t%s\t%d\t%d\t%*s\t%c\t%*s\t%[^\n]", ref, type, &start, &end, &strand, add_info);
        uint8_t is_rev = (strand == '-' ? 1 : 0);
        if (strcmp(type, "transcript") == 0) {
            if (t->exon_n > 1) {
                add_read_trans(T, *t);
                set_trans_name(T->t+T->trans_n-1, NULL, NULL, NULL, NULL);
                // for bam_trans
                T->t[T->trans_n-1].novel_exon_map = (uint8_t*)calloc(t->exon_n, sizeof(uint8_t));
                T->t[T->trans_n-1].novel_sj_map = (uint8_t*)calloc(t->exon_n-1, sizeof(uint8_t));
                T->t[T->trans_n-1].lfull = 0, T->t[T->trans_n-1].lnoth = 1, T->t[T->trans_n-1].rfull = 0, T->t[T->trans_n-1].rnoth = 1;
                T->t[T->trans_n-1].novel = 0, T->t[T->trans_n-1].all_novel=0, T->t[T->trans_n-1].all_iden=0;
            }
            t->exon_n = 0;
        } else if (strcmp(type, "exon") == 0) { // exon
            add_exon(t, bam_name2id(h, ref), start, end, is_rev);
            char tag[20]="gene_id";
            gtf_add_info(add_info, tag, gname); strcpy(t->gid, gname);
            strcpy(tag, "gene_name");
            gtf_add_info(add_info, tag, gname); strcpy(t->gname, gname);
            strcpy(tag, "transcript_name");
            gtf_add_info(add_info, tag, gname); strcpy(t->tname, gname);
            strcpy(tag, "transcript_id");
            gtf_add_info(add_info, tag, gname); strcpy(t->trans_id, gname);
        }
    }
    if (t->exon_n != 0) {
        add_read_trans(T, *t);
        set_trans_name(T->t+T->trans_n-1, NULL, NULL, NULL, NULL);
        // for bam_trans
        T->t[T->trans_n-1].novel_exon_map = (uint8_t*)calloc(t->exon_n, sizeof(uint8_t));
        T->t[T->trans_n-1].novel_sj_map = (uint8_t*)calloc(t->exon_n-1, sizeof(uint8_t));
        T->t[T->trans_n-1].lfull = 0, T->t[T->trans_n-1].lnoth = 1, T->t[T->trans_n-1].rfull = 0, T->t[T->trans_n-1].rnoth = 1;
        T->t[T->trans_n-1].novel = 0, T->t[T->trans_n-1].all_novel=0, T->t[T->trans_n-1].all_iden=0;
    }
    trans_free(t);
    return T->trans_n;
}

const struct option update_long_opt [] = {
    { "input-mode", 1, NULL, 'm' },
    { "bam", 1, NULL, 'b' },
    { "intron", 1, NULL, 'I' },
    { "min-exon", 1, NULL, 'e' },
    { "intron-len", 1, NULL, 'i' },
    { "distance", 1, NULL, 'd' },
    { "full-gtf", 1, NULL, 'l' },

    { "unclassified", 0, NULL, 'u' },
    { "source", 1, NULL, 's' },
    { "only-bam", 0, NULL, 'n' },
    { "full-bam", 0, NULL, 'f' },

    { 0, 0, 0, 0}
};

int update_gtf(int argc, char *argv[])
{
    int c; 
    update_gtf_para *ugp = update_gtf_init_para();
    while ((c = getopt_long(argc, argv, "m:b:i:I:e:d:l:us:no:", update_long_opt, NULL)) >= 0) {
        switch(c)
        {
            case 'm': if (optarg[0] == 'b') ugp->input_mode=0; else if (optarg[0] == 'g') ugp->input_mode=1; else return update_gtf_usage();
            case 'b': strcpy(ugp->in_bam, optarg); break;
            case 'I': if ((ugp->intron_fp = fopen(optarg, "r")) == NULL) {
                          err_fatal(__func__, "Can not open intron file \"%s\"\n", optarg);
                          return update_gtf_usage();
                      } 
                      break;
            case 'e': ugp->min_exon = atoi(optarg); break;
            case 'i': ugp->min_intron = atoi(optarg); break;
            case 'd': ugp->ss_dis = atoi(optarg); break;
            case 'l': ugp->full_len_level = atoi(optarg); break;
            case 'u': ugp->uncla = 1; break;
            case 's': strcpy(ugp->source, optarg); break;
            case 'n': ugp->only_bam = 1; break;
            case 'o': ugp->out_gtf_fp = fopen(optarg, "w"); break;
            default:
                      err_printf("Error: unknown option: %s.\n", optarg);
                      return update_gtf_usage();
                      break;
        }
    }
    if (argc - optind != 2) return update_gtf_usage();

    chr_name_t *cname = chr_name_init();
    read_trans_t *anno_T, *bam_T, *novel_T; gene_group_t *gg = gene_group_init();
    anno_T = read_trans_init(); bam_T = read_trans_init(); novel_T = read_trans_init(); 
    intron_group_t *I; I = intron_group_init();

    // read all input-transcript
    samFile *in; bam_hdr_t *h; 
    if (ugp->input_mode == 0) { // bam input
        bam1_t *b;     
        if ((in = sam_open(argv[optind], "rb")) == NULL) err_fatal(__func__, "Cannot open \"%s\"\n", argv[optind]);
        if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", argv[optind]);
        b = bam_init1(); 
        read_bam_trans(in, h, b, ugp, bam_T);
        bam_destroy1(b);
    } else { // gtf input
        if ((in = sam_open(ugp->in_bam, "rb")) == NULL) err_fatal(__func__, "Cannot open \"%s\"\n", ugp->in_bam);
        if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", ugp->in_bam);
        FILE *fp = xopen(argv[optind], "r");
        read_anno_trans(fp, h, bam_T);
    }

    FILE *gfp = xopen(argv[optind+1], "r");
    // read all anno-transcript
    read_anno_trans(gfp, h, anno_T);
    // read intron file
    read_intron_group(I, ugp->intron_fp);

    // identify novel transcript
    check_novel_trans(bam_T, anno_T, I, novel_T, ugp);

    // print novel transcript
    print_read_trans(anno_T, novel_T, h, ugp->source, ugp->out_gtf_fp);

    chr_name_free(cname);
    novel_read_trans_free(bam_T); novel_read_trans_free(anno_T); 
    read_trans_free(novel_T); intron_group_free(I); gene_group_free(gg);
    bam_hdr_destroy(h); sam_close(in); err_fclose(ugp->out_gtf_fp); if (ugp->intron_fp) err_fclose(ugp->intron_fp);
    return 0;
}
