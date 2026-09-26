rule ltrharvest:
    input:
        fasta="results/12_yahs/{hap}/asm_yahs_scaffolds_final.fa"
    output:
        gff3="results/15_ltrharvest/{hap}/{hap}.ltrharvest.gff3",
        fasta="results/15_ltrharvest/{hap}/{hap}.ltrharvest.fa",
        scn="results/15_ltrharvest/{hap}/{hap}.ltrharvest.scn"
    log:
        "logs/ltrharvest_{hap}.log"
    params:
        indexname="results/15_ltrharvest/{hap}/index"
    conda:
        "../envs/ltr.yaml"
    shell:
        "mkdir -p results/15_ltrharvest/{wildcards.hap}; "
        "gt suffixerator -db {input.fasta} -indexname {params.indexname} "
        "-tis -suf -lcp -des -ssp -sds -dna &> {log}; "
        "gt ltrharvest -index {params.indexname} -seqids yes "
        "-similar 90 -vic 10 -seed 20 "
        "-minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 "
        "-motif TGCA -motifmis 1 "
        "-out {output.fasta} -gff3 {output.gff3} "
        "> {output.scn} 2>> {log}"

rule ltrfinder:
    input:
        fasta="results/12_yahs/{hap}/asm_yahs_scaffolds_final.fa"
    output:
        scn="results/16_ltrfinder/{hap}/{hap}.ltrfinder.scn"
    log:
        "logs/ltrfinder_{hap}.log"
    threads: config["threads"]
    conda:
        "../envs/ltr.yaml"
    shell:
        "mkdir -p results/16_ltrfinder/{wildcards.hap}; "
        "LTR_FINDER_parallel -seq {input.fasta} -threads {threads} "
        "-harvest_out -output results/16_ltrfinder/{wildcards.hap} &> {log}; "
        "mv results/16_ltrfinder/{wildcards.hap}/$(basename {input.fasta}).finder.combine.scn {output.scn}"

rule ltrretriever:
    input:
        fasta="results/12_yahs/{hap}/asm_yahs_scaffolds_final.fa",
        harvest="results/15_ltrharvest/{hap}/{hap}.ltrharvest.scn",
        finder="results/16_ltrfinder/{hap}/{hap}.ltrfinder.scn"
    output:
        lib="results/17_ltrretriever/{hap}/{hap}.LTRlib.fa",
        passlist="results/17_ltrretriever/{hap}/{hap}.pass.list",
        out="results/17_ltrretriever/{hap}/{hap}.out"
    log:
        "logs/ltrretriever_{hap}.log"
    threads: config["threads"]
    conda:
        "../envs/ltr.yaml"
    shell:
        "mkdir -p results/17_ltrretriever/{wildcards.hap}; "
        "cd results/17_ltrretriever/{wildcards.hap} && "
        "LTR_retriever -genome $(realpath ../../../{input.fasta}) "
        "-inharvest $(realpath ../../../{input.harvest}) "
        "-infinder $(realpath ../../../{input.finder}) "
        "-threads {threads} &> $(realpath ../../../{log}); "
        "cd results/17_ltrretriever/{wildcards.hap} && "
        "mv $(basename {input.fasta}).LTRlib.fa {wildcards.hap}.LTRlib.fa && "
        "mv $(basename {input.fasta}).pass.list {wildcards.hap}.pass.list && "
        "mv $(basename {input.fasta}).out {wildcards.hap}.out"

rule lai:
    input:
        fasta="results/12_yahs/{hap}/asm_yahs_scaffolds_final.fa",
        passlist="results/17_ltrretriever/{hap}/{hap}.pass.list",
        out="results/17_ltrretriever/{hap}/{hap}.out"
    output:
        lai="results/18_lai/{hap}/{hap}.out.LAI"
    log:
        "logs/lai_{hap}.log"
    threads: config["threads"]
    conda:
        "../envs/ltr.yaml"
    shell:
        "mkdir -p results/18_lai/{wildcards.hap}; "
        "cd results/18_lai/{wildcards.hap} && "
        "LAI -genome $(realpath ../../../{input.fasta}) "
        "-intact $(realpath ../../../{input.passlist}) "
        "-all $(realpath ../../../{input.out}) "
        "-t {threads} &> $(realpath ../../../{log})"

rule repeatmodeler:
    input:
        fasta="results/12_yahs/{hap}/asm_yahs_scaffolds_final.fa"
    output:
        families_fa="results/19_repeatmodeler/{hap}/{hap}-families.fa",
        families_stk="results/19_repeatmodeler/{hap}/{hap}-families.stk"
    log:
        "logs/repeatmodeler_{hap}.log"
    threads: config["threads"]
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        "mkdir -p results/19_repeatmodeler/{wildcards.hap}; "
        "cd results/19_repeatmodeler/{wildcards.hap} && "
        "BuildDatabase -name {wildcards.hap} $(realpath ../../../{input.fasta}) "
        "&> $(realpath ../../../{log}); "
        "cd results/19_repeatmodeler/{wildcards.hap} && "
        "RepeatModeler -database {wildcards.hap} -pa {threads} "
        "&>> $(realpath ../../../{log})"

rule tesorter:
    input:
        lib="results/17_ltrretriever/{hap}/{hap}.LTRlib.fa"
    output:
        cls_lib="results/20_tesorter/{hap}/{hap}.rexdb.cls.lib",
        cls_tsv="results/20_tesorter/{hap}/{hap}.rexdb.cls.tsv"
    log:
        "logs/tesorter_{hap}.log"
    threads: config["threads"]
    conda:
        "../envs/tesorter.yaml"
    shell:
        "mkdir -p results/20_tesorter/{wildcards.hap}; "
        "cd results/20_tesorter/{wildcards.hap} && "
        "TEsorter $(realpath ../../../{input.lib}) -db rexdb-plantv3 "
        "-pre {wildcards.hap} -p {threads} &> $(realpath ../../../{log})"

rule combine_lib:
    input:
        repeatmodeler="results/19_repeatmodeler/{hap}/{hap}-families.fa",
        tesorter="results/20_tesorter/{hap}/{hap}.rexdb.cls.lib"
    output:
        combined="results/21_combine_lib/{hap}/{hap}.combined_lib.fa"
    log:
        "logs/combine_lib_{hap}.log"
    shell:
        "mkdir -p results/21_combine_lib/{wildcards.hap}; "
        "cat {input.repeatmodeler} {input.tesorter} > {output.combined} 2> {log}"

rule seqkit_clean:
    input:
        combined="results/21_combine_lib/{hap}/{hap}.combined_lib.fa"
    output:
        clean="results/22_seqkit_clean/{hap}/{hap}.combined_lib.clean.fa"
    log:
        "logs/seqkit_clean_{hap}.log"
    conda:
        "../envs/seqkit.yaml"
    shell:
        "mkdir -p results/22_seqkit_clean/{wildcards.hap}; "
        "seqkit rmdup -s {input.combined} -o {output.clean} &> {log}"

rule repeatmasker:
    input:
        fasta="results/12_yahs/{hap}/asm_yahs_scaffolds_final.fa",
        lib="results/22_seqkit_clean/{hap}/{hap}.combined_lib.clean.fa"
    output:
        masked="results/23_repeatmasker/{hap}/{hap}.masked.fa",
        out="results/23_repeatmasker/{hap}/{hap}.out",
        gff="results/23_repeatmasker/{hap}/{hap}.out.gff"
    log:
        "logs/repeatmasker_{hap}.log"
    threads: config["threads"]
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        "mkdir -p results/23_repeatmasker/{wildcards.hap}; "
        "RepeatMasker -pa {threads} -lib {input.lib} -xsmall -gff "
        "-dir results/23_repeatmasker/{wildcards.hap} {input.fasta} &> {log}; "
        "cd results/23_repeatmasker/{wildcards.hap} && "
        "mv $(basename {input.fasta}).masked {wildcards.hap}.masked.fa && "
        "mv $(basename {input.fasta}).out {wildcards.hap}.out && "
        "mv $(basename {input.fasta}).out.gff {wildcards.hap}.out.gff"

rule misa:
    input:
        fasta="results/12_yahs/{hap}/asm_yahs_scaffolds_final.fa",
        script="scripts/misa/misa.pl",
        ini="scripts/misa/misa.ini"
    output:
        misa="results/24_misa/{hap}/{hap}.fa.misa",
        stats="results/24_misa/{hap}/{hap}.fa.statistics"
    log:
        "logs/misa_{hap}.log"
    shell:
        "mkdir -p results/24_misa/{wildcards.hap}; "
        "cp {input.ini} results/24_misa/{wildcards.hap}/misa.ini; "
        "cp {input.fasta} results/24_misa/{wildcards.hap}/{wildcards.hap}.fa; "
        "cd results/24_misa/{wildcards.hap} && "
        "perl $(realpath ../../../{input.script}) {wildcards.hap}.fa "
        "&> $(realpath ../../../{log})"

rule trf:
    input:
        fasta="results/12_yahs/{hap}/asm_yahs_scaffolds_final.fa"
    output:
        dat="results/25_trf/{hap}/{hap}.trf.dat"
    log:
        "logs/trf_{hap}.log"
    conda:
        "../envs/tandem.yaml"
    shell:
        "mkdir -p results/25_trf/{wildcards.hap}; "
        "cp {input.fasta} results/25_trf/{wildcards.hap}/{wildcards.hap}.fa; "
        "cd results/25_trf/{wildcards.hap} && "
        "(trf {wildcards.hap}.fa 2 7 7 80 10 50 2000 -d -h "
        "&> $(realpath ../../../{log}) || true); "
        "mv {wildcards.hap}.fa.2.7.7.80.10.50.2000.dat {wildcards.hap}.trf.dat"
