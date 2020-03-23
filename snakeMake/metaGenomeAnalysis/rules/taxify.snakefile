# TODO
rule diamond:
    input:
    output:
    params:
    conda:
        "envs/blobtools.yaml"
    shell:
        """
        diamond blastx --query --db --threads --outfmt 6 --sensitive --max-target-seqs 1 --evalue 1e-25 -o
        """

#TODO
rule blobtools_taxify:
    input:
    output:
    params:
    log:
    conda:
        "envs/blobtools.yaml"
    shell:
        """
        blobtools taxify -f -m -s 0 -t 2 -o
        """

#TODO
rule blobtools_create:
    input:
    output:
    params:
    log:
    conda:
        "envs/blobtools.yaml"
    shell:
        """
        blobtools create -i -b -t -o --db
        """

#TODO
rule blobtools_viewplot:
    input:
    output:
    params:
    log:
    conda:
        "envs/blobtools.yaml"
    shell:
    """
    blobtools plot -i --notitle --format pdf -r superkingdom -o
    blobtools plot -i --notitle --format pdf -r phylum -o

    blobtools view -i -o -r all
    """
