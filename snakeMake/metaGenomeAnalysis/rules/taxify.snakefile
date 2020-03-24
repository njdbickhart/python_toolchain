

rule diamond:
    input:
        "assembly/{assembly_group}.fa"
    output:
        "blobtools/{assembly_group}.diamondout.tsv"
    threads: 16
    params:
        db = config['diamonddb']
    conda:
        "envs/blobtools.yaml"
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
        tax = config['diamondtaxid']
    log:
    conda:
        "envs/blobtools.yaml"
    shell:
        """
        blobtools taxify -f {input} -m {params.tax} -s 0 -t 2 -o {output}
        """

rule blobtools_cov:
    input:
        bams = "/mapping/{assembly_group}/{sample}.bam",
        fasta = "assembly/{assembly_group}.fa"
    output:
        cov = "blobtools/{assembly_group}_{sample}.cov"
    params:
        outbase = "blobtools/{assembly_group}_{sample}"
    conda:
        "envs/blobtools.yaml"
    shell:
        """
        blobtools map2cov -i {input.fasta} -b {input.bams} -o {params.outbase}
        """

rule blobtools_create:
    input:
        contigs = "assembly/{assembly_group}.fa",
        covs = expand("blobtools/{assembly_group}_{sample}.cov", assembly_group=wildcards.assembly_group, sample=wildcards.sample),
        tax = "blobtools/taxify.{assembly_group}.diamondout.tsv.taxified.out"
    output:
        blobdb = "blobtools/{assembly_group}.blobDB.json"
    params:
        cstr = "-c " + " -c ".join(input.covs),
        outpre = "blobtools/{assembly_group}",
        db = config["ncbidb"]
    log:
    conda:
        "envs/blobtools.yaml"
    shell:
        """
        echo "using {params.cstr}"
        blobtools create -i {input.contigs} {params.cstr} -t {input.tax} -o {params.outpre} --db {params.db}
        """

rule blobtools_viewplot:
    input:
        blobdb = "blobtools/{assembly_group}.blobDB.json"
    output:
        supplot = "blobtools/supkingdom.{assembly_group}.blobDB.json.bestsum.superkingdom.p8.span.100.blobplot.covsum.pdf",
        phylumplot = "blobtools/phylum.{assembly_group}.blobDB.json.bestsum.phylum.p8.span.100.blobplot.covsum.pdf",
        table = "blobtools/table.{assembly_group}.blobDB.table.txt"
    params:
    log:
    conda:
        "envs/blobtools.yaml"
    shell:
    """
    blobtools plot -i {input.blobdb} --notitle --format pdf -r superkingdom -o blobtools/supkingdom
    blobtools plot -i {input.blobdb} --notitle --format pdf -r phylum -o blobtools/phylum

    blobtools view -i {input.blobdb} -o blobtools/table -r all
    """
