
rule nanoplot:
    input:
        hifi=config["hifi"]
    output:
        report="results/01_nanoplot/NanoPlot-report.html",
        stats="results/01_nanoplot/NanoStats.txt"
    log:
        "logs/nanoplot.log"
    threads: config["threads"]
    conda:
        "../envs/nanoplot.yaml"
    shell:
        "NanoPlot -t {threads} --fastq {input.hifi} -o results/01_nanoplot &> {log}"


rule fastp:
    input:
        r1=config["hic_r1"],
        r2=config["hic_r2"]
    output:
        r1=f"results/02_fastp/{HIC_R1_BASE}.trimmed.fastq.gz",
        r2=f"results/02_fastp/{HIC_R2_BASE}.trimmed.fastq.gz",
        json="results/02_fastp/fastp.json",
        html="results/02_fastp/fastp.html"
    log:
        "logs/fastp.log"
    threads: config["threads"]
    conda:
        "../envs/fastp.yaml"
    shell:
        "fastp -p -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
        "--detect_adapter_for_pe --thread {threads} --json {output.json} --html {output.html} &> {log}"


rule kmc:
    input:
        hifi=config["hifi"]
    output:
        kmc_pre="results/03_kmc/output.kmc_pre",
        kmc_suf="results/03_kmc/output.kmc_suf",
        hist="results/03_kmc/output_hist"
    log:
        "logs/kmc.log"
    threads: config["threads"]
    params:
        k=config["kmc"]["k"],
        ci=config["kmc"]["ci"],
        cs=config["kmc"]["cs"],
        cx=config["kmc"]["cx"]
    conda:
        "../envs/kmc.yaml"
    shell:
        "mkdir -p results/03_kmc/tmp; "
        "kmc -k{params.k} -t{threads} -ci{params.ci} -cs{params.cs} -fq {input.hifi} "
        "results/03_kmc/output results/03_kmc/tmp &> {log}; "
        "kmc_tools transform results/03_kmc/output histogram {output.hist} -cx{params.cx} &>> {log}"


rule genomescope2:
    input:
        hist="results/03_kmc/output_hist"
    output:
        outdir=directory("results/04_genomescope2")
    log:
        "logs/genomescope2.log"
    params:
        k=config["genomescope2"]["k"],
        p=config["genomescope2"]["p"]
    conda:
        "../envs/genomescope2.yaml"
    shell:
        "genomescope2 -i {input.hist} -o {output.outdir} -k {params.k} -p {params.p} &> {log}"
