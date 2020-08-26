import re

rule all_taxified:
    input:
        "FinishedTaxify"

rule diamond:
    input:
        "assembly/{assembly_group}.fa"
    output:
        "blobtools/{assembly_group}.diamondout.tsv"
    threads: 16
    params:
        db = config['diamonddb']
    conda:
        "../envs/blobtools.yaml"
    shell:
        """
        diamond blastx --query {input} --db {params.db} --threads {threads} --outfmt 6 --sensitive --max-target-seqs 1 --evalue 1e-25 -o {output}
        """

rule blobtools_taxify:
    input:
        "blobtools/{assembly_group}.diamondout.tsv"
    output:
        "blobtools/taxify.{assembly_group}.diamondout.tsv.taxified.out"
    threads: 2
    params:
        tax = config['diamondtaxid'],
        blobtools = config['blobtools'],
        outbase = "blobtools/taxify"
    log:
    conda:
        "../envs/blobtools.yaml"
    shell:
        """
        {params.blobtools} taxify -f {input} -m {params.tax} -s 0 -t 2 -o {params.outbase}
        """

rule blobtools_cov:
    input:
        bams = "mapping/{assembly_group}/{sample}.bam",
        bais = "mapping/{assembly_group}/{sample}.bam.bai",
        fasta = "assembly/{assembly_group}.fa"
    output:
        cov = "blobtools/{assembly_group}.{sample}.bam.cov"
    params:
        outbase = "blobtools/{assembly_group}",
        blobtools = config['blobtools']
    conda:
        "../envs/blobtools.yaml"
    shell:
        """
        {params.blobtools} map2cov -i {input.fasta} -b {input.bams} -o {params.outbase}
        """

def getCovStr(covs):
    return "-c " + " -c ".join(covs)

rule blobtools_create:
    input:
        contigs = "assembly/{assembly_group}.fa",
        covs = expand("blobtools/{assembly_group}.{sample}.bam.cov", assembly_group=getAssemblyBaseName(config["assemblies"]), sample=config["samples"]),
        tax = "blobtools/taxify.{assembly_group}.diamondout.tsv.taxified.out"
    output:
        blobdb = "blobtools/{assembly_group}.blobDB.json"
    params:
        cstr = getCovStr(expand("blobtools/{assembly_group}.{sample}.bam.cov", assembly_group=getAssemblyBaseName(config["assemblies"]), sample=config["samples"])),
        outpre = "blobtools/{assembly_group}",
        db = config["ncbidb"],
        blobtools = config['blobtools']
    log:
    conda:
        "../envs/blobtools.yaml"
    shell:
        """
        echo "using {params.cstr}"
        {params.blobtools} create -i {input.contigs} {params.cstr} -t {input.tax} -o {params.outpre} --db {params.db}
        """

rule blobtools_viewplot:
    input:
        blobdb = "blobtools/{assembly_group}.blobDB.json"
    output:
        supplot = "blobtools/supkingdom.{assembly_group}.blobDB.json.bestsum.superkingdom.p8.span.100.blobplot.covsum.pdf",
        phylumplot = "blobtools/phylum.{assembly_group}.blobDB.json.bestsum.phylum.p8.span.100.blobplot.covsum.pdf",
        table = "blobtools/table.{assembly_group}.blobDB.table.txt"
    conda:
        "../envs/blobtools.yaml"
    params:
        blobtools = config['blobtools']
    shell:
        """
        {params.blobtools} plot -i {input.blobdb} --notitle --format pdf -r superkingdom -o blobtools/supkingdom
        {params.blobtools} plot -i {input.blobdb} --notitle --format pdf -r phylum -o blobtools/phylum

        {params.blobtools} view -i {input.blobdb} -o blobtools/table -r all
        """

rule blobtab_condense:
    input:
        table = "blobtools/table.{assembly_group}.blobDB.table.txt"
    output:
        crop = "blobtools/table.{assembly_group}.lens.tab"
    params:
        script = workflow.basedir + "/scripts/trimBlobTable.py"
    shell:
        """
        {params.script} {input.table} {output.crop}
        """

def interlace(input):
    v = list()
    for i in input.crops:
        m = re.search(r'blobtools/table.(.+).lens.tab', i)
        v.append(m.group(1))
        v.append(i)
    return ' '.join(v)

rule blobview_plot:
    input:
        crops = expand("blobtools/table.{assembly_group}.lens.tab", assembly_group=getAssemblyBaseName(config["assemblies"]))
    output:
        plot = "blobtools/summary_taxonomic_plot.pdf"
    conda:
        "../envs/rlibs.yaml"
    params:
        interlaced = lambda wildcards, input : interlace(input),
        script = workflow.basedir + "/scripts/gcSupKingPlot.R"
    shell:
        """
        Rscript {params.script} {params.interlaced}
        """


rule blobtools_complete:
    input:
        expand("blobtools/table.{assembly_group}.blobDB.table.txt", assembly_group=getAssemblyBaseName(config["assemblies"])),
        "blobtools/summary_taxonomic_plot.pdf"
    output:
        temp(touch("FinishedTaxify"))
