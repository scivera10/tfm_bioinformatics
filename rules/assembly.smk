rule hifiasm:
    input:
        hifi=config["hifi"],
        hic_r1=f"results/02_fastp/{HIC_R1_BASE}.trimmed.fastq.gz",
        hic_r2=f"results/02_fastp/{HIC_R2_BASE}.trimmed.fastq.gz"
    output:
        gfa_hap1="results/05_hifiasm/asm.hic.hap1.p_ctg.gfa",
        gfa_hap2="results/05_hifiasm/asm.hic.hap2.p_ctg.gfa",
        fasta_hap1="results/05_hifiasm/asm.hic.hap1.p_ctg.fa",
        fasta_hap2="results/05_hifiasm/asm.hic.hap2.p_ctg.fa"
    log:
        "logs/hifiasm.log"
    threads: config["hifiasm_threads"]
    params:
        prefix="results/05_hifiasm/asm",
        opts="-f0"
    conda:
        "../envs/hifiasm.yaml"
    shell:
        "hifiasm {params.opts} -t {threads} -o {params.prefix} "
        "--h1 {input.hic_r1} --h2 {input.hic_r2} {input.hifi} &> {log}; "
        "awk -f scripts/hifiasm/gfa_to_fasta.awk < {output.gfa_hap1} > {output.fasta_hap1}; "
        "awk -f scripts/hifiasm/gfa_to_fasta.awk < {output.gfa_hap2} > {output.fasta_hap2}"

rule busco:
    input:
        assembly="results/05_hifiasm/asm.hic.{hap}.p_ctg.fa"
    output:
        dir=directory("results/06_busco/{hap}_busco")
    log:
        "logs/busco_{hap}.log"
    threads: config["threads"]
    params:
        lineage="eudicots_odb10"
    conda:
        "../envs/busco.yaml"
    shell:
        "busco -i {input.assembly} -o {wildcards.hap}_busco --out_path results/06_busco "
        "-l {params.lineage} -m genome -c {threads} -f &> {log}"

rule quast:
    input:
        assemblies=["results/05_hifiasm/asm.hic.hap1.p_ctg.fa", "results/05_hifiasm/asm.hic.hap2.p_ctg.fa"]
    output:
        dir=directory("results/07_quast")
    log:
        "logs/quast.log"
    threads: config["threads"]
    params:
        labels="hap1,hap2"
    conda:
        "../envs/quast.yaml"
    shell:
        "quast {input.assemblies} -l {params.labels} -o {output.dir} -t {threads} &> {log}"


rule purge_dups:
    input:
        hifi=config["hifi"],
        assembly="results/05_hifiasm/asm.hic.{hap}.p_ctg.fa"
    output:
        dir=directory("results/08_purge_dups/{hap}"),
        stat="results/08_purge_dups/{hap}/PB.stat",
        cutoffs="results/08_purge_dups/{hap}/cutoffs",
        fasta="results/08_purge_dups/{hap}/purged.fa"
    log:
        "logs/purge_dups_{hap}.log"
    threads: config["threads"]
    conda:
        "../envs/purge_dups.yaml"
    shell:
        "LOG=$(realpath {log}); "
        "mkdir -p {output.dir}; "
        "minimap2 -xasm20 {input.assembly} {input.hifi} -t {threads} 2> \"$LOG\" | gzip -c - > {output.dir}/hifi_vs_{wildcards.hap}.paf.gz; "
        "cd {output.dir} && "
        "pbcstat hifi_vs_{wildcards.hap}.paf.gz >> \"$LOG\" 2>&1 && "
        "calcuts PB.stat > cutoffs 2>> \"$LOG\" && "
        "split_fa ../../05_hifiasm/asm.hic.{wildcards.hap}.p_ctg.fa > asm.hic.{wildcards.hap}.p_ctg.split 2>> \"$LOG\" && "
        "minimap2 -xasm5 -DP asm.hic.{wildcards.hap}.p_ctg.split asm.hic.{wildcards.hap}.p_ctg.split -t {threads} 2>> \"$LOG\" | gzip -c - > asm.hic.{wildcards.hap}.p_ctg.split.self.paf.gz && "
        "purge_dups -2 -T cutoffs -c PB.base.cov asm.hic.{wildcards.hap}.p_ctg.split.self.paf.gz > dups.bed 2>> \"$LOG\" && "
        "get_seqs -e dups.bed ../../05_hifiasm/asm.hic.{wildcards.hap}.p_ctg.fa > /dev/null 2>> \"$LOG\""


rule hist_plot:
    input:
        stat="results/08_purge_dups/{hap}/PB.stat",
        cutoffs="results/08_purge_dups/{hap}/cutoffs",
        script="scripts/purge_dups/hist_plot.py"
    output:
        png="results/08_purge_dups/{hap}/hist.out.png"
    log:
        "logs/hist_plot_{hap}.log"
    conda:
        "../envs/purge_dups.yaml"
    shell:
        "S=$(realpath {input.script}); "
        "LOG=$(realpath {log}); "
        "cd results/08_purge_dups/{wildcards.hap} && "
        "MPLBACKEND=Agg python3 $S -c cutoffs PB.stat hist.out.png &> $LOG"

