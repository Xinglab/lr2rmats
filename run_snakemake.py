import os, sys, re
import yaml
import argparse


def parse_argv(args, cfg_temp):
    cfg = cfg_temp
    # genome
    ## genome fasta and GTF
    cfg['genome']['fasta'] = args.genome
    cfg['genome']['fasta_dir'] = os.path.dirname(args.genome)
    cfg['genome']['gtf'] = args.gtf
    cfg['genome']['rm_gtf'] = '-r ' + args.rm_gtf if args.rm_gtf else ''

    ## genome index
    cfg['genome']['minimap_idx'] = args.minimap2_idx if args.minimap2_idx else args.genome + '.smmi'
    cfg['genome']['star_idx'] = args.STAR_idx if args.STAR_idx else args.genome + '.STAR'

    ## long reads
    cfg['sample']['long_read'] = {}
    with open(args.long_read_list) as in_list:
        n_samp_long = int(in_list.readline().split()[0])
        samp_i = 1
        for i in range(0, n_samp_long):
            n_rep = int(in_list.readline().split()[0])
            for j in range(0, n_rep):
                cfg['sample']['long_read']['samp' + str(samp_i)] = os.path.abspath(in_list.readline().split()[0])
                samp_i += 1

    ## short reads
    cfg['sample']['short_read'] = {}
    with open(args.short_read_list) as in_list:
        n_samp_short = int(in_list.readline().split()[0])
        samp_i = 1
        for i in range(0, n_samp_short):
            n_rep = int(in_list.readline().split()[0])
            for j in range(0, n_rep):
                line = in_list.readline()
                if '#' in line: line = line[:line.index('#')]

                cfg['sample']['short_read']['samp' + str(samp_i)] = {}
                cfg['sample']['short_read']['samp' + str(samp_i)]['first'] = os.path.abspath(line.split()[0])
                cfg['sample']['short_read']['samp' + str(samp_i)]['second'] = os.path.abspath(line.split()[1]) if len(
                    line.split()) >= 2 else []
                samp_i += 1

    # running parameter
    cfg['lr2rmats'] = {}
    cfg['lr2rmats']['aln_cov'] = args.aln_cov
    cfg['lr2rmats']['iden_frac'] = args.iden_frac
    cfg['lr2rmats']['sec_rat'] = args.sec_rat
    cfg['lr2rmats']['sup_cnt'] = args.sup_cnt
    cfg['lr2rmats']['split_trans'] = '-s' if args.split_trans == True else ''

    # default running resources
    cfg['__default__']['threads'] = args.default_threads
    cfg['__default__']['h_data'] = args.default_mem
    cfg['__default__']['h_rt'] = args.default_time

    # running resources on local machine
    cfg['star_idx']['threads'] = args.star_idx_threads
    cfg['star_map']['threads'] = args.star_map_threads
    cfg['minimap_idx']['threads'] = args.minimap_idx_threads
    cfg['minimap_map']['threads'] = args.minimap_map_threads
    cfg['novel_gtf']['threads'] = args.novel_gtf_threads
    return cfg


def add_parse_argv(parser):
    genome_par = parser.add_argument_group('Genome sequence and gene annotation')
    genome_par.add_argument('--genome', metavar='FILE', type=str, required=True,
                            help='Genome sequence in FASTA format.')
    genome_par.add_argument('--STAR-idx', type=str,
                            help='STAR index directory for genome sequence. Index will be built if not set.')
    genome_par.add_argument('--minimap2-idx', type=str,
                            help='Minimap2 index directory for genome sequence. Index will be built if not set.')

    genome_par.add_argument('--gtf', metavar='FILE', type=str, required=True, help='Gene annotation in GTF format.')
    genome_par.add_argument('--rm-gtf', type=str, default='',
                            help='Transcript file in GTF format that need to be removed from long-read data, e.g., rRNA.')

    sample_par = parser.add_argument_group('Input data (long read and short read)')
    sample_par.add_argument('--long-read-list', metavar='LIST', type=str, required=True,
                            help='Long-read input file list. Refer to \'long_read_fa_input.example.list\'')
    sample_par.add_argument('--short-read-list', metavar='LIST', type=str, required=True,
                            help='Corresponding short-read input file list. Refer to \'short_read_fa_input.example.list\'')

    long_argv = parser.add_argument_group('Running parameters')
    long_argv.add_argument('--aln-cov', type=float, default=0.67,
                           help='Minimum fraction of long-read\'s aligned bases.')
    long_argv.add_argument('--iden-frac', type=float, default=0.75,
                           help='Minimum fraction of long-read\'s identically aligned bases.')
    long_argv.add_argument('--sec-rat', type=float, default=0.98,
                           help='Maximum ratio of second best and best score to retain the best alignment.')
    long_argv.add_argument('--sup-cnt', type=int, default=1,
                           help='Minimum supporting count of novel-junction spanning short-read to make novel-junction reliable.')
    long_argv.add_argument('--split-trans', default=False, action='store_true',
                           help='Split long-read derived transcript on unreliable junctions.')

    def_running_res = parser.add_argument_group(
        'Default running resources for all the jobs\n (overwritten if specified for specific job)')
    def_running_res.add_argument('--out-dir', metavar='DIR', type=str, help='output folder', default='.')
    def_running_res.add_argument('--default-threads', metavar='INT', type=int, help='default threads number', default=4)
    def_running_res.add_argument('--default-mem', metavar='MEM', type=str, help='default memory allocation',
                                 default='8G')
    def_running_res.add_argument('--default-time', metavar='TIME', type=str, help='default time allocation',
                                 default='02:00:00')

    running_res = parser.add_argument_group("Running resources")

    running_res.add_argument('--cores', metavar='INT', type=int, help='number of available cores', default=8)
    running_res.add_argument('--star-idx-threads', metavar='INT', type=int, help='threads for STAR indexing', default=8)
    running_res.add_argument('--star-map-threads', metavar='INT', type=int, help='threads for STAR mapping', default=8)
    running_res.add_argument('--minimap-idx-threads', metavar='INT', type=int, help='threads for minimap2 indexing',
                             default=8)
    running_res.add_argument('--minimap-map-threads', metavar='INT', type=int, help='threads for minimap2 mapping',
                             default=8)
    running_res.add_argument('--novel-gtf-threads', metavar='INT', type=int,
                             help='threads for generating novel GTF', default=4)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    add_parse_argv(parser)
    args = parser.parse_args()

    dir = os.path.dirname(os.path.abspath(sys.argv[0]))
    with open(dir + '/config.template.yaml') as cfg_temp_fp:
        cfg_temp = yaml.load(cfg_temp_fp)
    cfg_temp_fp.close()
    cfg = parse_argv(args, cfg_temp)
    cfg_fn = args.out_dir + '/config.yaml'
    with open(cfg_fn, 'w') as cfg_fp:
        yaml.dump(cfg, cfg_fp, default_flow_style=False)
    cfg_fp.close()

    snake = dir + '/Snakefile'
    cmd = "snakemake -p --snakefile " + snake + " --configfile " + cfg_fn + " --jobs " + str(
        args.cores) + " --directory " + args.out_dir
    print cmd
    os.system(cmd)
