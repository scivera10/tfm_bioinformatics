import os

configfile: "config/config.yaml"

HIC_R1_BASE = os.path.basename(config["hic_r1"]).replace(".fastq.gz", "")
HIC_R2_BASE = os.path.basename(config["hic_r2"]).replace(".fastq.gz", "")

wildcard_constraints:
    hap="hap[12]",
    read="[12]"

include: "rules/preprocessing.smk"
include: "rules/assembly.smk"
include: "rules/arima_yahs.smk"
include: "rules/repeat_annotation.smk"

rule all:
    input:
        "results/01_nanoplot/NanoPlot-report.html",
        "results/01_nanoplot/NanoStats.txt",
        f"results/02_fastp/{HIC_R1_BASE}.trimmed.fastq.gz",
        f"results/02_fastp/{HIC_R2_BASE}.trimmed.fastq.gz",
        "results/02_fastp/fastp.json",
        "results/02_fastp/fastp.html",
        "results/03_kmc/output_hist",
        "results/04_genomescope2",
        "results/05_hifiasm/asm.hic.hap1.p_ctg.fa",
        "results/05_hifiasm/asm.hic.hap2.p_ctg.fa",
        expand("results/06_busco/{hap}_busco", hap=["hap1", "hap2"]),
        "results/07_quast",
        expand("results/08_purge_dups/{hap}", hap=["hap1", "hap2"]),
        expand("results/08_purge_dups/{hap}/hist.out.png", hap=["hap1", "hap2"]),
        expand("results/11_arima_mapping/{hap}/RAW_DIR/hic_vs_contigs_{read}.bam", hap=["hap1", "hap2"], read=["1", "2"]),
        expand("results/11_arima_mapping/{hap}/FILT_DIR/hic_vs_contigs_filt_{read}.bam", hap=["hap1", "hap2"], read=["1", "2"]),
        expand("results/08_purge_dups/{hap}/purged.fa.fai", hap=["hap1", "hap2"]),
        expand("results/11_arima_mapping/{hap}/REP_DIR/paired_mark_dups_final.bam", hap=["hap1", "hap2"]),
        expand("results/11_arima_mapping/{hap}/REP_DIR/paired_mark_dups_final.bam.bai", hap=["hap1", "hap2"]),
        expand("results/11_arima_mapping/{hap}/REP_DIR/final.bam.stats", hap=["hap1", "hap2"]),
        expand("results/12_yahs/{hap}/asm_yahs_scaffolds_final.fa", hap=["hap1", "hap2"]),
        expand("results/13_busco_scaffolds/{hap}", hap=["hap1", "hap2"]),
        "results/14_quast_scaffolds",
        expand("results/15_ltrharvest/{hap}/{hap}.ltrharvest.scn", hap=["hap1", "hap2"]),
        expand("results/16_ltrfinder/{hap}/{hap}.ltrfinder.scn", hap=["hap1", "hap2"]),
        expand("results/17_ltrretriever/{hap}/{hap}.LTRlib.fa", hap=["hap1", "hap2"]),
        expand("results/17_ltrretriever/{hap}/{hap}.pass.list", hap=["hap1", "hap2"]),
        expand("results/18_lai/{hap}/{hap}.out.LAI", hap=["hap1", "hap2"]),
        expand("results/19_repeatmodeler/{hap}/{hap}-families.fa", hap=["hap1", "hap2"]),
        expand("results/20_tesorter/{hap}/{hap}.rexdb.cls.lib", hap=["hap1", "hap2"]),
        expand("results/21_combine_lib/{hap}/{hap}.combined_lib.fa", hap=["hap1", "hap2"]),
        expand("results/22_seqkit_clean/{hap}/{hap}.combined_lib.clean.fa", hap=["hap1", "hap2"]),
        expand("results/23_repeatmasker/{hap}/{hap}.masked.fa", hap=["hap1", "hap2"]),
        expand("results/23_repeatmasker/{hap}/{hap}.out.gff", hap=["hap1", "hap2"]),
        expand("results/24_misa/{hap}/{hap}.fa.misa", hap=["hap1", "hap2"]),
        expand("results/25_trf/{hap}/{hap}.trf.dat", hap=["hap1", "hap2"])
