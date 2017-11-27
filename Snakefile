rule all:
    input:
        config["output"]["gtf"]

# STAR build index file
rule star_idx:
    input:
        config["genome"]["fasta"]
    output:
        config["genome"]["star_idx"]
    threads:
        config["star_idx"]["threads"]
    benchmark:
        "benchmark/star_idx.benchmark.txt"
    params:
        star=config["exe_files"]["star"]
    shell:
        "mkdir {output}; "
        "{params.star} --runMode genomeGenerate --runThreadN {threads} --genomeFastaFiles {input} --genomeDir {output}"

# minimap build index file
rule minimap_idx:
    input:
        config["genome"]["fasta"]
    output:
        config["genome"]["minimap_idx"]
    threads:
        config["minimap_idx"]["threads"]
    benchmark:
        "benchmark/minimap_idx.benchmark.txt"
    params:
        minimap=config["exe_files"]["minimap2"]
    shell:
        "{params.minimap} -x splice {input} -d {output} -t {threads}"


# minimap mapping for long reads
rule minimap_map:
    input:
        genome=config["genome"]["minimap_idx"],
        reads=lambda wildcards: config["sample"]["long_read"][wildcards.sample]
    output:
        "alignment/{sample}.minimap.sam"
    threads:
        config["minimap_map"]["threads"]
    benchmark:
        "benchmark/{sample}.minimap.benchmark.txt"
    params:
        minimap=config["exe_files"]["minimap2"]
    shell:
        "{params.minimap} -ax splice -ub -t {threads} {input.genome} {input.reads} > {output}"

# filter long read alignment and generate novel GTF file
rule sam_novel_gtf:
    input:
        sam="alignment/{sample}.minimap.sam",
        rRNA=config["genome"]["rRNA"],
        gtf=config["genome"]["gtf"]
    output:
        "gtf/{sample}_sam_novel.gtf"
    threads:
        config["novel_gtf"]["threads"]
    benchmark:
        "benchmark/{sample}.novel_gtf.benchmark.txt"
    params:
        lr2gtf=config["exe_files"]["lr2gtf"],
        samtools=config["exe_files"]["samtools"]
    shell:
        "{params.lr2gtf} filter {input.sam} {input.rRNA} | {params.samtools} sort -@ {threads} | {params.lr2gtf} update-gtf - {input.gtf} > {output}"

# merge and sort gtf
rule new_gtf:
    input:
        gtf=config["genome"]["gtf"],
        novel_gtf="gtf/{sample}_sam_novel.gtf"
    output:
        tmp=temp("gtf/{sample}_tmp.gtf"),
        gtf="gtf/{sample}_new.gtf"
    benchmark:
        "benchmark/{sample}_new_gtf.benchmark.txt"
    params:
        sort_gtf=config["exe_files"]["sort_gtf"]
    shell:
        "cp {input.gtf} {output.tmp}; "
        "cat {input.novel_gtf} >> {output.tmp}; "
        "{params.sort_gtf} {output.tmp} {output.gtf}; "


# STAR mapping for short reads
rule star_new_map:
    input:
        genome=config["genome"]["star_idx"],
        gtf="gtf/{sample}_new.gtf",
        read1=lambda wildcards: config["sample"]["short_read"][wildcards.sample]["first"],
        read2=lambda wildcards: config["sample"]["short_read"][wildcards.sample]["second"]
    output:
        bam="alignment/{sample}.new.STAR.sort.bam",
        SJ="alignment/{sample}.new.STARSJ.out.tab"
    threads:
        config["star_new_map"]["threads"]
    benchmark:
        "benchmark/{sample}.star.new.benchmark.txt"
    params:
        star=config["exe_files"]["star"],
        samtools=config["exe_files"]["samtools"],
        prefix="alignment/{sample}",
        b="alignment/{sample}.new.STAR.sort.bam.txt"
    shell:
        "{params.star} --runThreadN {threads}  --genomeDir {input.genome} --readFilesIn {input.read1} {input.read2} "
        "--outFileNamePrefix {params.prefix}.new.STAR --outSAMtype BAM Unsorted "
        "--outFilterType BySJout   --outFilterMultimapNmax 20 "
        "--outFilterMismatchNmax 999   --alignIntronMin 25   --alignIntronMax 1000000   --alignMatesGapMax 1000000 "
        "--alignSJoverhangMin 8   --alignSJDBoverhangMin 5   --sjdbGTFfile {input.gtf}  --sjdbOverhang 100 ; "
        "{params.samtools} sort {params.prefix}.new.STARAligned.out.bam -@ {threads} > {output.bam} ;"
        "ls -1 {output.bam} > {params.b}"

rule gtf_novel_gtf:
    input:
        gtf=config["genome"]["gtf"],
        novel_gtf="gtf/{sample}_sam_novel.gtf",
        bam="alignment/{sample}.new.STAR.sort.bam",
        SJ="alignment/{sample}.new.STARSJ.out.tab"
    output:
        "gtf/{sample}_gtf_novel.gtf"
    benchmark:
        "benchmark/{sample}_gtf_novel_gtf.benchmark.txt"
    params:
        lr2gtf=config["exe_files"]["lr2gtf"],
        sort_gtf=config["exe_files"]["sort_gtf"],
        samtools=config["exe_files"]["samtools"]
    shell:
        "{params.lr2gtf} update-gtf -mg -b {input.bam} -I {input.SJ} {input.novel_gtf} {input.gtf} > {output}"

rule update_gtf:
    input:
        gtf=config["genome"]["gtf"],
        novel_gtf=expand("gtf/{sample}_gtf_novel.gtf", sample=config["sample"]["long_read"])
    output:
        config["output"]["gtf"]
    benchmark:
        "benchmark/update_gtf.benchmark.txt"
    params:
        lr2gtf=config["exe_files"]["lr2gtf"],
        sort_gtf=config["exe_files"]["sort_gtf"],
        samtools=config["exe_files"]["samtools"]
    shell:
        "cp {input.gtf} gtf/tmp.gtf; "
        "cat {input.novel_gtf} >> gtf/tmp.gtf; "
        "{params.sort_gtf} gtf/tmp.gtf {output}; rm gtf/tmp.gtf"

