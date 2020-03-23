rule prodigal:
    input:
        assembly = "assembly/{assembly_group}.fa",
    output:
        proteins = "/prodigal/{assembly_group}/proteins.faa",
        genes = "/prodigal/{assembly_group}/genes.gff"
    conda:
        "envs/prodigal.yaml"
    shell:
        """
        prodigal -i {input.assembly} -f gff -o {output.genes} -a {output.proteins} -p meta
        """
