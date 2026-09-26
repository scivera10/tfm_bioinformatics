rule bwa_index:
    input:
        fasta="results/08_purge_dups/{hap}/purged.fa"
    output:
        idx=multiext("results/11_arima_mapping/{hap}/purged", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/bwa_index_{hap}.log"
    conda:
        "../envs/arima_mapping.yaml"
    shell:
        "mkdir -p results/11_arima_mapping/{wildcards.hap}; "
        "bwa index -p results/11_arima_mapping/{wildcards.hap}/purged {input.fasta} &> {log}"


rule bwa_mem:
    input:
        idx=multiext("results/11_arima_mapping/{hap}/purged", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        reads=lambda wildcards: f"results/02_fastp/{HIC_R1_BASE if wildcards.read == '1' else HIC_R2_BASE}.trimmed.fastq.gz"
    output:
        bam="results/11_arima_mapping/{hap}/RAW_DIR/hic_vs_contigs_{read}.bam"
    log:
        "logs/bwa_mem_{hap}_R{read}.log"
    threads: config["threads"]
    conda:
        "../envs/arima_mapping.yaml"
    shell:
        "mkdir -p results/11_arima_mapping/{wildcards.hap}/RAW_DIR; "
        "bwa mem -5SPM -t {threads} results/11_arima_mapping/{wildcards.hap}/purged {input.reads} 2> {log} "
        "| samtools view -@ {threads} -Sb - > {output.bam}"

rule filter_five_end:
    input:
        bam="results/11_arima_mapping/{hap}/RAW_DIR/hic_vs_contigs_{read}.bam",
        script="scripts/arima_mapping_pipeline/filter_five_end.pl"
    output:
        bam="results/11_arima_mapping/{hap}/FILT_DIR/hic_vs_contigs_filt_{read}.bam"
    log:
        "logs/filter_five_end_{hap}_R{read}.log"
    conda:
        "../envs/arima_mapping.yaml"
    shell:
        "mkdir -p results/11_arima_mapping/{wildcards.hap}/FILT_DIR; "
        "samtools view -h {input.bam} | perl {input.script} 2> {log} | samtools view -Sb - > {output.bam}"

rule faidx:
    input:
        fasta="results/08_purge_dups/{hap}/purged.fa"
    output:
        fai="results/08_purge_dups/{hap}/purged.fa.fai"
    log:
        "logs/faidx_{hap}.log"
    conda:
        "../envs/arima_mapping.yaml"
    shell:
        "samtools faidx {input.fasta} &> {log}"

rule combine_pairs:
    input:
        r1="results/11_arima_mapping/{hap}/FILT_DIR/hic_vs_contigs_filt_1.bam",
        r2="results/11_arima_mapping/{hap}/FILT_DIR/hic_vs_contigs_filt_2.bam",
        fai="results/08_purge_dups/{hap}/purged.fa.fai",
        script="scripts/arima_mapping_pipeline/two_read_bam_combiner.pl"
    output:
        bam="results/11_arima_mapping/{hap}/TMP_DIR/hic_vs_contigs_filt_paired.bam"
    log:
        "logs/two_read_bam_combiner_{hap}.log"
    threads: config["threads"]
    params:
        mapq=10
    conda:
        "../envs/arima_mapping.yaml"
    shell:
        "mkdir -p results/11_arima_mapping/{wildcards.hap}/TMP_DIR; "
        "perl {input.script} {input.r1} {input.r2} samtools {params.mapq} 2> {log} "
        "| samtools view -bS -t {input.fai} - "
        "| samtools sort -@ {threads} -o {output.bam} - 2>> {log}"

rule add_read_group:
    input:
        bam="results/11_arima_mapping/{hap}/TMP_DIR/hic_vs_contigs_filt_paired.bam"
    output:
        bam="results/11_arima_mapping/{hap}/PAIR_DIR/paired_add_read_group.bam"
    log:
        "logs/picard_addreadgroup_{hap}.log"
    params:
        sample="SRR32416950",
        mem="4G"
    conda:
        "../envs/arima_mapping.yaml"
    shell:
        "mkdir -p results/11_arima_mapping/{wildcards.hap}/PAIR_DIR; "
        "picard -Xmx{params.mem} AddOrReplaceReadGroups INPUT={input.bam} OUTPUT={output.bam} "
        "ID={params.sample} LB={params.sample} SM={params.sample} PL=ILLUMINA PU=none &> {log}"


rule mark_duplicates:
    input:
        bam="results/11_arima_mapping/{hap}/PAIR_DIR/paired_add_read_group.bam"
    output:
        bam="results/11_arima_mapping/{hap}/REP_DIR/paired_mark_dups_final.bam",
        metrics="results/11_arima_mapping/{hap}/REP_DIR/metrics.final.txt"
    log:
        "logs/picard_markdup_{hap}.log"
    params:
        mem="4G"
    conda:
        "../envs/arima_mapping.yaml"
    shell:
        "mkdir -p results/11_arima_mapping/{wildcards.hap}/REP_DIR; "
        "picard -Xmx{params.mem} -XX:-UseGCOverheadLimit MarkDuplicates "
        "INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.metrics} "
        "ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE &> {log}"

rule index_bam:
    input:
        bam="results/11_arima_mapping/{hap}/REP_DIR/paired_mark_dups_final.bam"
    output:
        bai="results/11_arima_mapping/{hap}/REP_DIR/paired_mark_dups_final.bam.bai"
    log:
        "logs/index_bam_{hap}.log"
    conda:
        "../envs/arima_mapping.yaml"
    shell:
        "samtools index {input.bam} &> {log}"


rule bam_stats:
    input:
        bam="results/11_arima_mapping/{hap}/REP_DIR/paired_mark_dups_final.bam",
        script="scripts/arima_mapping_pipeline/get_stats.pl"
    output:
        stats="results/11_arima_mapping/{hap}/REP_DIR/final.bam.stats"
    log:
        "logs/get_stats_{hap}.log"
    conda:
        "../envs/arima_mapping.yaml"
    shell:
        "perl {input.script} {input.bam} > {output.stats} 2> {log}"


rule yahs:
    input:
        fasta="results/08_purge_dups/{hap}/purged.fa",
        fai="results/08_purge_dups/{hap}/purged.fa.fai",
        bam="results/11_arima_mapping/{hap}/REP_DIR/paired_mark_dups_final.bam"
    output:
        fasta="results/12_yahs/{hap}/asm_yahs_scaffolds_final.fa",
        agp="results/12_yahs/{hap}/asm_yahs_scaffolds_final.agp",
        bin="results/12_yahs/{hap}/asm_yahs.bin"
    log:
        "logs/yahs_{hap}.log"
    params:
        enzyme="AAGCTT"
    conda:
        "../envs/yahs.yaml"
    shell:
        "mkdir -p results/12_yahs/{wildcards.hap}; "
        "yahs {input.fasta} {input.bam} -o results/12_yahs/{wildcards.hap}/asm_yahs -e {params.enzyme} &> {log}"


rule busco_scaffolds:
    input:
        assembly="results/12_yahs/{hap}/asm_yahs_scaffolds_final.fa"
    output:
        dir=directory("results/13_busco_scaffolds/{hap}")
    log:
        "logs/busco_scaffolds_{hap}.log"
    threads: config["threads"]
    params:
        lineage="eudicots_odb10"
    conda:
        "../envs/busco.yaml"
    shell:
        "busco -i {input.assembly} -o {wildcards.hap} --out_path results/13_busco_scaffolds "
        "-l {params.lineage} -m genome -c {threads} -f &> {log}"


rule quast_scaffolds:
    input:
        assemblies=["results/12_yahs/hap1/asm_yahs_scaffolds_final.fa", "results/12_yahs/hap2/asm_yahs_scaffolds_final.fa"]
    output:
        dir=directory("results/14_quast_scaffolds")
    log:
        "logs/quast_scaffolds.log"
    threads: config["threads"]
    params:
        labels="hap1_scaffolds,hap2_scaffolds"
    conda:
        "../envs/quast.yaml"
    shell:
        "quast {input.assemblies} -l {params.labels} -o {output.dir} -t {threads} --fragmented --large &> {log}"
