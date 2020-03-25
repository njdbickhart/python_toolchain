import os

KING = ["full", "euk"]

rule all_binned:
    input:
        "FinishedBinning"

rule temp_completion:
    input:
        expand("binning/metabat2/{assembly_group}/metabat_bin_full.tab", assembly_group=getAssemblyBaseName(config["assemblies"])),
        expand("binning/metabat2/{assembly_group}/metabat_bin_euk.tab", assembly_group=getAssemblyBaseName(config["assemblies"])),
        expand("binning/concoct/{assembly_group}/clustering_merged_{king}.csv", assembly_group=getAssemblyBaseName(config["assemblies"]), king=KING)
    output:
        temp(touch("FinishedBinning"))


rule bwa_index:
    input:
        "assembly/{assembly_group}.fa"
    output:
        "assembly/{assembly_group}.fa.amb",
        "assembly/{assembly_group}.fa.ann",
        "assembly/{assembly_group}.fa.bwt",
        "assembly/{assembly_group}.fa.pac",
        "assembly/{assembly_group}.fa.sa"
    log:
        config['logdir'] + "/{assembly_group}.log"
    conda:
        "../envs/metabat.yaml"
    shell:
        """
        bwa index {input} 2> {log}
        """

rule bwa_mem:
    input:
        amb = "assembly/{assembly_group}.fa.amb",
        ann = "assembly/{assembly_group}.fa.ann",
        bwt = "assembly/{assembly_group}.fa.bwt",
        pac = "assembly/{assembly_group}.fa.pac",
        sa = "assembly/{assembly_group}.fa.sa",
        reference = "assembly/{assembly_group}.fa",
        r1 = lambda wildcards: config["samples"][wildcards.sample][0],
        r2 = lambda wildcards: config["samples"][wildcards.sample][1],
    output:
        "mapping/{assembly_group}/{sample}.bam"
    log:
        config['logdir'] + "/bwa_mem/{assembly_group}_{sample}.log"
    params:
        extra="",
        #pipe_cmd = "samtools sort -o {output} -",
        threads = 8
    conda:
        "../envs/metabat.yaml"
    shell:
        """
        bwa mem -t {params.threads} {params.extra} {input.reference} {input.r1} {input.r2} | samtools sort -o {output} - >> {log} 2>&1
        """

# Metabat

rule metabat_abundance:
    input:
        expand("mapping/{assembly_group}/{sample}.bam", assembly_group=getAssemblyBaseName(config["assemblies"]), sample=config["samples"])
    output:
        "binning/metabat2/{assembly_group}/jgi_abund.txt"
    conda:
         "../envs/metabat.yaml"
    log:
        config['logdir'] + "/binning/{assembly_group}.abun.log"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output} {input} > {log} 2>&1
        """

rule metabat_binning:
    input:
        assembly = "assembly/{assembly_group}.fa",
        depth = "binning/metabat2/{assembly_group}/jgi_abund.txt"
    output:
        "binning/metabat2/{assembly_group}/metabat_bin_full.tab"
    conda:
         "../envs/metabat.yaml"
    params:
        other = "--saveCls --minContig 3000 --noBinOut",
        threads = 8
    log:
        config['logdir'] + "/binning/{assembly_group}.metabat.log"
    shell:
        """
        metabat2 {params.other} --numThreads {params.threads} -i {input.assembly} -a {input.depth} -o {output} > {log} 2>&1
        """

rule eukrep:
    input:
        assembly = "assembly/{assembly_group}.fa"
    output:
        "eukrep/{assembly_group}/euk.final.contigs.fa"
    conda:
        "../envs/EukRep.yaml"
    log:
        config['logdir'] + "/eukrep/{assembly_group}.eukrep.log"
    params:
        prok = "eukrep/{assembly_group}/prok.final.contigs.fa",
        min_contig = 1000
    shell:
        """
        EukRep -i {input} -o {output} --prokarya {params.prok} --min {params.min_contig} > {log} 2>&1
        """

rule metabat_binning_euk:
    input:
        assembly = "eukrep/{assembly_group}/euk.final.contigs.fa",
        depth = "binning/metabat2/{assembly_group}/jgi_abund.txt"
    output:
        "binning/metabat2/{assembly_group}/metabat_bin_euk.tab"
    conda:
         "../envs/metabat.yaml"
    params:
        other = "--saveCls --minContig 3000 --noBinOut",
        threads = 8
    log:
        config['logdir'] + "/metabat2/{assembly_group}.eukbin.log"
    shell:
        """
        metabat2 {params.other} --numThreads {params.threads} -i {input.assembly} -a {input.depth} -o {output} > {log} 2>&1
        """

# Concoct

rule concoct_ctgprep:
    input:
        assembly = "assembly/{assembly_group}.fa"
    output:
        contigs = "binning/concoct/{assembly_group}/contigs_10k_full.fa",
        bed = "binning/concoct/{assembly_group}/contigs_10k_full.bed"
    conda:
        "../envs/concoct.yaml"
    shell:
        """
        cut_up_fasta.py {input.assembly} -c 10000 -o 0 --merge_last -b {output.bed} > {output.contigs}
        """

rule concoct_ctgprep_euk:
    input:
        assembly = "eukrep/{assembly_group}/euk.final.contigs.fa"
    output:
        contigs = "binning/concoct/{assembly_group}/contigs_10k_euk.fa",
        bed = "binning/concoct/{assembly_group}/contigs_10k_euk.bed"
    conda:
        "../envs/concoct.yaml"
    shell:
        """
        cut_up_fasta.py {input.assembly} -c 10000 -o 0 --merge_last -b {output.bed} > {output.contigs}
        """

rule concoct_calc_cov:
    input:
        bed = "binning/concoct/{assembly_group}/contigs_10k_{king}.bed",
        bam = expand("mapping/{assembly_group}/{sample}.bam", assembly_group=getAssemblyBaseName(config["assemblies"]), sample=config["samples"])
    output:
        coverage = "binning/concoct/{assembly_group}/coverage_file_{king}.tab"
    conda:
        "../envs/concoct.yaml"
    shell:
        """
        concoct_coverage_table.py {input.bed} {input.bam} > {output.coverage}
        """

rule run_concoct:
    input:
        contigs = "binning/concoct/{assembly_group}/contigs_10k_{king}.fa",
        coverage = "binning/concoct/{assembly_group}/coverage_file_{king}.tab"
    output:
        bins = temp("binning/concoct/{assembly_group}/{king}/clustering_gt3000.csv"),
        fbins = "binning/concoct/{assembly_group}/clustering_merged_{king}.csv"
    params:
        basename = "binning/concoct/{assembly_group}/{king}/",
        read_length = 150,
        min_length = 3000
    conda:
        "../envs/concoct.yaml"
    shell:
        """
        concoct \
            --coverage_file {input.coverage} \
            --composition_file {input.fasta} \
            --basename {params.basename} \
            --read_length {params.read_length} \
            --length_threshold {params.min_length} \
            --converge_out

        merge_cutup_clustering.py {output.bins} > {output.fbins}
        """
