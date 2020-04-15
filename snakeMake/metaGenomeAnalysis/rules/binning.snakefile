import os

KING = ["full", "euk"]
BINS = ["metabat2", "concoct"]
if config.get("hic"):
    BINS.append("bin3c")

wildcard_constraints:
    king = "(full|euk)"

localrules: convert_concoct, metabat_convert, modify_bin3c

rule all_binned:
    input:
        "FinishedBinning"

rule temp_completion:
    input:
        expand("binning/{bins}/{assembly_group}/{bins}.{king}.clusters.tab", bins= BINS, assembly_group=getAssemblyBaseName(config["assemblies"]), king=KING),
        expand("binning/DASTool/{assembly_group}.{king}_cluster_attribution.tsv", assembly_group=getAssemblyBaseName(config["assemblies"]), king=KING),
        expand("binning/DASTool/{assembly_group}.{king}_cluster_counts.tab", assembly_group=getAssemblyBaseName(config["assemblies"]), king=KING)
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
        "assembly/{assembly_group}.fa.sa",
        "assembly/{assembly_group}.fa.fai"
    log:
        config['logdir'] + "/{assembly_group}.log"
    conda:
        "../envs/metabat.yaml"
    shell:
        """
        bwa index {input} 2> {log}
        samtools faidx {input} 2> {log}
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
rule create_bam_index:
    input:
        "mapping/{assembly_group}/{file}.bam"
    output:
        "mapping/{assembly_group}/{file}.bam.bai"
    conda:
        "../envs/metabat.yaml"
    threads:
        1
    shell:
        "samtools index {input}"
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

rule metabat_convert:
    input:
        bin = "binning/metabat2/{assembly_group}/metabat_bin_{king}.tab"
    output:
        cluster = "binning/metabat2/{assembly_group}/metabat2.{king}.clusters.tab"
    run:
        with open(input.bin, 'r') as bins, open(output.cluster, 'w') as out:
            for l in bins:
                s = l.rstrip().split()
                out.write(f'{s[0]}\tmetabat2.{s[1]}\n')

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
        bam = expand("mapping/{assembly_group}/{sample}.bam", assembly_group=getAssemblyBaseName(config["assemblies"]), sample=config["samples"]),
        bai = expand("mapping/{assembly_group}/{sample}.bam.bai", assembly_group=getAssemblyBaseName(config["assemblies"]), sample=config["samples"])
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
            --composition_file {input.contigs} \
            --basename {params.basename} \
            --read_length {params.read_length} \
            --length_threshold {params.min_length} \
            --converge_out

        merge_cutup_clustering.py {output.bins} > {output.fbins}
        """
rule convert_concoct:
    input:
        bin = "binning/concoct/{assembly_group}/clustering_merged_{king}.csv"
    output:
        cluster = "binning/concoct/{assembly_group}/concoct.{king}.clusters.tab"
    run:
        with open(input.bin, 'r') as bins, open(output.cluster, 'w') as out:
            bins.readline()
            for l in bins:
                s = l.rstrip().split(',')
                out.write(f'{s[0]}\tconcoct.{s[1]}\n')

if config.get("hic"):

    rule align_hic:
        input:
            amb = "assembly/{assembly_group}.fa.amb",
            ann = "assembly/{assembly_group}.fa.ann",
            bwt = "assembly/{assembly_group}.fa.bwt",
            pac = "assembly/{assembly_group}.fa.pac",
            sa = "assembly/{assembly_group}.fa.sa",
            reference = "assembly/{assembly_group}.fa",
            r1 = lambda wildcards: config["hic"][wildcards.enzyme][0],
            r2 = lambda wildcards: config["hic"][wildcards.enzyme][1]
        output:
            "mapping/{assembly_group}/hic_{enzyme}.bam"
        threads: 8
        conda:
            "../envs/metabat.yaml"
        shell:
            """
            bwa mem -t {threads} -5SP {input.reference} {input.r1} {input.r2} | \
            samtools view -F 0x904 -bS - | \
            samtools sort -o {output} -n -
            """

    rule bin3c_contact:
        input:
            reference = "assembly/{assembly_group}.fa",
            bam = "mapping/{assembly_group}/hic_{enzyme}.bam"
        output:
            directory("binning/bin3c/{assembly_group}/{enzyme}_full_out")
        threads: 1
        conda:
            "../envs/bin3c.yaml"
        params:
            enzyme = lambda wildcards: wildcards.enzyme,
            bin3c = config["bin3c"]
        shell:
            """
            {params.bin3c} mkmap -e {params.enzyme} -v {input.reference} {input.bam} {output}
            """

    rule bin3c_contact_euk:
        input:
            reference = "eukrep/{assembly_group}/euk.final.contigs.fa",
            bam = "mapping/{assembly_group}/hic_{enzyme}.bam"
        output:
            directory("binning/bin3c/{assembly_group}/{enzyme}_euk_out")
        threads: 1
        conda:
            "../envs/bin3c.yaml"
        params:
            enzyme = lambda wildcards: wildcards.enzyme,
            bin3c = config["bin3c"]
        shell:
            """
            {params.bin3c} mkmap -e {params.enzyme} -v {input.reference} {input.bam} {output}
            """

    rule bin3c_cluster:
        input:
            directory("binning/bin3c/{assembly_group}/{enzyme}_{king}_out")
        output:
            outclust = "binning/bin3c/{assembly_group}/{enzyme}_{king}_clust/clustering.mcl"
        threads: 1
        params:
            outfolder = "binning/bin3c/{assembly_group}/{enzyme}_{king}_temp",
            realout = "binning/bin3c/{assembly_group}/{enzyme}_{king}_clust",
            bin3c = config["bin3c"]
        conda:
            "../envs/bin3c.yaml"
        shell:
            """
            python2 {params.bin3c} cluster --no-plot -v {input}/contact_map.p.gz {params.outfolder}
            mv {params.outfolder}/clustering.mcl {params.realout}/
            rm -r {params.outfolder}
            """

    rule modify_bin3c:
        input:
            expand("binning/bin3c/{assembly_group}/{enzyme}_{king}_clust/clustering.mcl", assembly_group=getAssemblyBaseName(config["assemblies"]), enzyme=config["hic"], king=KING)
        output:
            "binning/bin3c/{assembly_group}/bin3c.{king}.clusters.tab"
        run:
            with open(output[0], 'w') as out:
                seen = set()
                bnum = 0
                for j in input:
                    with open(j, 'r') as bins:
                        for l in bins:
                            s = l.rstrip().split()
                            for i in s:
                                if not i in seen:
                                    out.write(f'{i}\tbin3c.{bnum}\n')
                                seen.add(i)
                            bnum += 1

rule das_tool:
    input:
        binfiles = expand("binning/{bins}/{assembly_group}/{bins}.full.clusters.tab", bins=BINS, assembly_group=getAssemblyBaseName(config["assemblies"])),
        reference = "assembly/{assembly_group}.fa"
    output:
        expand("binning/DASTool/{{assembly_group}}.full{postfix}",
               postfix=["_DASTool_summary.txt", "_DASTool_hqBins.pdf", "_DASTool_scores.pdf"]),
        expand("binning/DASTool/{{assembly_group}}.full_{bins}.eval",
               bins= BINS),
        cluster_attribution = "binning/DASTool/{assembly_group}.full_cluster_attribution.tsv"
    threads: 10
    log:
        config['logdir'] + "/dastool.{assembly_group}.full.log"
    conda:
        "../envs/dastool.yaml"
    params:
        binnames = ",".join(BINS),
        scaffolds2bin = lambda wildcards, input: ",".join(input.binfiles),
        output_prefix = "binning/DASTool/{assembly_group}.full"
    shell:
        """
        module load usearch/11.0.667
        echo {params.output_prefix} {params.scaffolds2bin} {params.binnames}
        DAS_Tool --outputbasename {params.output_prefix} --bins {params.scaffolds2bin} \
        --labels {params.binnames} --contigs {input.reference} --search_engine diamond \
        --write_bin_evals 1 --create_plots 1 --threads {threads} --debug &> {log};
        mv {params.output_prefix}_DASTool_scaffolds2bin.txt {output.cluster_attribution} &>> {log}
        """

rule das_tool_euk:
    input:
        binfiles = expand("binning/{bins}/{assembly_group}/{bins}.euk.clusters.tab", bins=BINS, assembly_group=getAssemblyBaseName(config["assemblies"])),
        reference = "assembly/{assembly_group}.fa"
    output:
        expand("binning/DASTool/{{assembly_group}}.euk{postfix}",
               postfix=["_DASTool_summary.txt", "_DASTool_hqBins.pdf", "_DASTool_scores.pdf"]),
        expand("binning/DASTool/{{assembly_group}}.euk_{bins}.eval",
               bins= BINS),
        cluster_attribution = "binning/DASTool/{assembly_group}.euk_cluster_attribution.tsv"
    threads: 10
    log:
        config['logdir'] + "/dastool.{assembly_group}.euk.log"
    conda:
        "../envs/dastool.yaml"
    params:
        binnames = ",".join(BINS),
        scaffolds2bin = lambda wildcards, input: ",".join(input.binfiles),
        output_prefix = "binning/DASTool/{assembly_group}.euk"
    shell:
        """
        module load usearch/11.0.667
        echo {params.output_prefix} {params.scaffolds2bin} {params.binnames}
        DAS_Tool --outputbasename {params.output_prefix} --bins {params.scaffolds2bin} \
        --labels {params.binnames} --contigs {input.reference} --search_engine diamond \
        --write_bin_evals 1 --create_plots 1 --threads {threads} --debug &> {log};
        mv {params.output_prefix}_DASTool_scaffolds2bin.txt {output.cluster_attribution} &>> {log}
        """

rule mag_generation:
    input:
        cluster = "binning/DASTool/{assembly_group}.{king}_cluster_attribution.tsv",
        reference = "assembly/{assembly_group}.fa"
    output:
        files = dynamic("mags/{assembly_group}/{id}.fa")
    params:
        king = "{king}",
        dir = "mags/{assembly_group}"
    conda:
        "../envs/concoct.yaml"
    script:
        "../scripts/magGeneration.py"

def getIds():
    assembly, ids, = glob_wildcards("mags/{assembly_group}/{id}.fa")
    return ids

rule mag_counts:
    input:
        files = expand("mags/{assembly_group}/{id}.fa", assembly_group=getAssemblyBaseName(config["assemblies"]), id=getIds())
        cluster = "binning/DASTool/{assembly_group}.{king}_cluster_attribution.tsv"
    output:
        counts = "binning/DASTool/{assembly_group}.{king}_cluster_counts.tab"
    run:
        from collections import defaultdict
        with open(input.cluster, 'r') as clust, open(output.counts, 'w') as out:
            bins = defaultdict(int)
            for l in clust:
                s = l.rstrip().split()
                bins[s[1]] += 1
            for b, c in bins.items():
                out.write(f'{b}\t{c}\n')
