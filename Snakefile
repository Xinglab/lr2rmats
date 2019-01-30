rule all:
    input:
        config["output"]["updated_gtf"]

# STAR build index file
rule star_idx:
    input:
        config["genome"]["fasta"]
    output:
        directory(config["genome"]["star_idx_dir"])
    threads:
        config["star_idx"]["threads"]
    log:
        "logs/star_idx"
    benchmark:
        "benchmark/star_idx.benchmark.txt"
    params:
        star=config["exe_files"]["star"]
    shell:
        "mkdir {output}; "
        "{params.star} --runMode genomeGenerate --runThreadN {threads} --genomeFastaFiles {input} --genomeDir {output} --outFileNamePrefix {log} >> {log}"

# minimap build index file
rule minimap_idx:
    input:
        config["genome"]["fasta"]
    output:
        config["genome"]["minimap_idx_file"]
    threads:
        config["minimap_idx"]["threads"]
    log:
        "logs/minimap_idx.log"
    benchmark:
        "benchmark/minimap_idx.benchmark.txt"
    params:
        minimap=config["exe_files"]["minimap2"]
    shell:
        "{params.minimap} -x splice {input} -d {output} -t {threads} 2> {log}"


# minimap mapping for long reads
rule minimap_map:
    input:
        genome=config["genome"]["minimap_idx_file"],
        reads=lambda wildcards: config["sample"]["long_read"][wildcards.sample]
    output:
        sam="alignment/{sample}.minimap.sam",
        bam=temp("alignment/{sample}.minimap.bam"),
        bed="alignment/{sample}.minimap.bed"
    threads:
        config["minimap_map"]["threads"]
    log:
        "logs/minimap_map/{sample}.log"
    benchmark:
        "benchmark/{sample}.minimap.benchmark.txt"
    params:
        minimap=config["exe_files"]["minimap2"],
        samtools=config["exe_files"]["samtools"],
        bedtools=config["exe_files"]["bedtools"]
    shell:
        "{params.minimap} -ax splice -ub -t {threads} {input.genome} {input.reads} > {output.sam} 2> {log}; "
        "{params.samtools} view -b {output.sam} > {output.bam}; "
        "{params.bedtools} bamtobed -bed12 -i {output.bam} > {output.bed} "


# filter long read alignment and generate novel GTF file
rule sam_novel_gtf:
    input:
        sam="alignment/{sample}.minimap.sam",
        gtf=config["genome"]["gtf"]
    output:
        filtered_bam="alignment/{sample}.filtered.bam",
        sam_gtf="gtf/{sample}_sam_novel.gtf"
    threads:
        config["novel_gtf"]["threads"]
    log:
        "logs/sam_novel_gtf/{sample}.log"
    benchmark:
        "benchmark/{sample}.novel_gtf.benchmark.txt"
    params:
        rm_gtf=config["lr2rmats"]["rm_gtf"],
        aln_cov=config["lr2rmats"]["aln_cov"],
        iden_frac=config["lr2rmats"]["iden_frac"],
        sec_rat=config["lr2rmats"]["sec_rat"],
        full_level=config["lr2rmats"]["full_level"],
        lr2rmats=config["exe_files"]["lr2rmats"],
        samtools=config["exe_files"]["samtools"]
    run:
        if {params.rm_gtf} == '':
            shell("{params.lr2rmats} filter {input.sam} -v {params.aln_cov} -q {params.iden_frac} -s {params.sec_rat} 2> {log} | {params.samtools} sort -@ {threads} > {output.filtered_bam} 2>> {log}; ")
        else:
            shell("{params.lr2rmats} filter {input.sam} -r {params.rm_gtf} -v {params.aln_cov} -q {params.iden_frac} -s {params.sec_rat} 2> {log} | {params.samtools} sort -@ {threads} > {output.filtered_bam} 2>> {log}; ")
        shell("{params.lr2rmats} update-gtf {output.filtered_bam} {input.gtf} -l {params.full_level} 2>> {log} > {output.sam_gtf}")

# merge and sort gtf
rule new_gtf:
    input:
        gtf=config["genome"]["gtf"],
        novel_gtf="gtf/{sample}_sam_novel.gtf"
    output:
        tmp=temp("gtf/{sample}_tmp.gtf"),
        gtf="gtf/{sample}_new.gtf"
    log:
        "logs/new_gtf/{sample}.log"
    benchmark:
        "benchmark/{sample}_new_gtf.benchmark.txt"
    params:
        sort_gtf=config["exe_files"]["sort_gtf"]
    shell:
        "cp {input.gtf} {output.tmp} 2> {log}; "
        "cat {input.novel_gtf} >> {output.tmp} 2> {log}; "
        "{params.sort_gtf} {output.tmp} {output.gtf} 2>> {log}; "


# STAR mapping for short reads
rule star_map:
    input:
        read1=lambda wildcards: config["sample"]["short_read"][wildcards.sample]["first"],
        read2=lambda wildcards: config["sample"]["short_read"][wildcards.sample]["second"],
        genome=config["genome"]["star_idx_dir"],
        gtf="gtf/{sample}_new.gtf"
    output:
        bam="alignment/{sample}.STARAligned.out.bam",
        SJ="alignment/{sample}.STARSJ.out.tab"
    threads:
        config["star_map"]["threads"]
    log:
        "logs/star_map/{sample}.log"
    benchmark:
        "benchmark/{sample}.star.benchmark.txt"
    params:
        star=config["exe_files"]["star"],
        samtools=config["exe_files"]["samtools"],
        prefix="alignment/{sample}"
    shell:
        "{params.star} --runThreadN {threads}  --genomeDir {input.genome} --readFilesIn {input.read1} {input.read2} "
        "--outFileNamePrefix {params.prefix}.STAR --outSAMtype BAM Unsorted "
        "--outFilterType BySJout   --outFilterMultimapNmax 20 "
        "--outFilterMismatchNmax 999   --alignIntronMin 25   --alignIntronMax 1000000   --alignMatesGapMax 1000000 "
        "--alignSJoverhangMin 8   --alignSJDBoverhangMin 5   --sjdbGTFfile {input.gtf}  --sjdbOverhang 100 > {log} 2 >& 1; "

rule gtf_novel_gtf:
    input:
        gtf=config["genome"]["gtf"],
        #novel_gtf="gtf/{sample}_sam_novel.gtf",
        filtered_bam="alignment/{sample}.filtered.bam",
        #bam="alignment/{sample}.STARAligned.out.bam",
        SJ="alignment/{sample}.STARSJ.out.tab"
    output:
        update_gtf="gtf/{sample}_gtf_novel.gtf",
        known_gtf="output/{sample}.known.gtf",
        novel_gtf="output/{sample}.novel.gtf",
        unrecog_gtf="output/{sample}.unrecog.gtf",
        bam_gtf="output/{sample}.bam.gtf",
        detail="output/{sample}.detail.txt",
        summary="output/{sample}.summary.txt",
        exon_bed="output/{sample}.novel_exon.bed"
    log:
        "logs/gtf_novel_gtf/{sample}.log"
    benchmark:
        "benchmark/{sample}_gtf_novel_gtf.benchmark.txt"
    params:
        sup_cnt=config["lr2rmats"]["sup_cnt"],
        split_trans=config["lr2rmats"]["split_trans"],
        full_level=config["lr2rmats"]["full_level"],
        lr2rmats=config["exe_files"]["lr2rmats"],
        sort_gtf=config["exe_files"]["sort_gtf"],
        samtools=config["exe_files"]["samtools"]
    shell:
        "{params.lr2rmats} update-gtf {params.split_trans} -l {params.full_level} -J {params.sup_cnt} -j {input.SJ} {input.filtered_bam} {input.gtf} -y {output.summary} -a {output.bam_gtf} -A  {output.detail} -k {output.known_gtf} -v {output.novel_gtf} -u {output.unrecog_gtf} -E {output.exon_bed}  > {output.update_gtf} 2> {log}"

rule update_gtf:
    input:
        gtf=config["genome"]["gtf"],
        novel_gtf=expand("gtf/{sample}_gtf_novel.gtf", sample=config["sample"]["long_read"]),
        sam=expand("alignment/{sample}.minimap.sam", sample=config["sample"]["long_read"]),
    output:
        uniq_gtf=temp("gtf/uniq.gtf"),
        updated_gtf=config["output"]["updated_gtf"]
    log:
        "logs/update_gtf.log"
    benchmark:
        "benchmark/update_gtf.benchmark.txt"
    params:
        lr2rmats=config["exe_files"]["lr2rmats"],
        sort_gtf=config["exe_files"]["sort_gtf"],
        samtools=config["exe_files"]["samtools"]
    shell:
        "cat {input.novel_gtf} >> gtf/tmp.gtf; "
        "{params.lr2rmats} unique-gtf -mg -b {input.sam[0]} gtf/tmp.gtf > {output.uniq_gtf} 2> {log}; "
        "cat {input.gtf} {output.uniq_gtf} > gtf/tmp.gtf 2>> {log}; "
        "{params.sort_gtf} gtf/tmp.gtf {output.updated_gtf} 2>> {log}; rm gtf/tmp.gtf"
