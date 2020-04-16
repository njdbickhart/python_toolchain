# Note: many of these rules were refactored from MagPy
# Credit to: https://github.com/WatsonLab/MAGpy
import os

#TODO: refactor for my pipeline

localrules: sourmash_report, diamond_bin_summary



rule all_stats:
    input:
        "FinishedStats"

rule stats_completion:
    input:
        expand("stats/{assembly_group}/total.checkm_plus.txt", assembly_group=getAssemblyBaseName(config["assemblies"])),
        "stats/diamond_bin_report_plus.tsv",
        'stats/sourmash/sourmash_report.csv'
    output:
        temp(touch("FinishedStats"))

rule checkm:
    input:
        mags = expand("mags/{assembly_group}/{id}.fa", assembly_group=getAssemblyBaseName(config["assemblies"]), id=getIds()),
        fini = "FinishedBinning"
    output: "stats/{assembly_group}/total.checkm.txt"
    threads: 16
    conda:
        "../envs/checkm.yaml"
    params:
        cdr=config["checkm_dataroot"],
        dir="mags/{assembly_group}"
    shell:
        """
        checkm_db={params.cdr}
        echo ${{checkm_db}} | checkm data setRoot ${{checkm_db}}
        checkm lineage_wf -f checkm.txt --reduced_tree -t {threads} -x fa {params.dir} ./checkm
        """

rule checkm_plus:
    input: "stats/{assembly_group}/total.checkm.txt"
    output: "stats/{assembly_group}/total.checkm_plus.txt"
    threads: 1
    conda: "../envs/ete3.yaml"
    script: "../scripts/add_tax.py"

rule prodigal_stats:
    input:
        mags = 'mags/{assembly_group}/{id}.fa',
        fini = "FinishedBinning"
    output:
        faa='stats/{assembly_group}/proteins/{id}.faa',
        gff='stats/{assembly_group}/proteins/{id}_prodigal.gff'
    conda: "../envs/prodigal.yaml"
    shell: 'prodigal -p meta -a {output.faa} -q -i {input.mags} -f gff -o {output.gff}'

rule diamond_prer:
    input: 'stats/{assembly_group}/proteins/{id}.faa'
    output: 'stats/{assembly_group}/diamond/{id}.diamond.tsv'
    threads: 16
    params:
        db=config["diamonddb"],
        of="6 qseqid sseqid stitle pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore"
    conda: "../envs/blobtools.yaml"
    shell: "diamond blastp --threads {threads} --max-target-seqs 10 --db {params.db} --query {input} --outfmt {params.of} --out {output}"

rule diamond_report:
    input:
        tsv='stats/{assembly_group}/diamond/{id}.diamond.tsv',
        faa='stats/{assembly_group}/proteins/{id}.faa'
    output:
        'stats/{assembly_group}/diamond_report/bin.{id}.tsv',
        'stats/{assembly_group}/diamond_report/con.{id}.tsv'
    params:
        outdir="stats/{assembly_group}/diamond_report",
        dir = workflow.basedir + '/scripts'
    conda: "../envs/bioperl.yaml"
    shell: "{params.dir}/diamond_report.pl {input.tsv} {input.faa} {params.outdir}"

rule diamond_bin_summary:
    input: expand("stats/{assembly_group}/diamond_report/bin.{id}.tsv", assembly_group=getAssemblyBaseName(config["assemblies"]), id=getIds())
    output: temp("stats/diamond_bin_report.tsv")
    shell:
        """
        echo -e 'name\tnprots\tnhits\tnfull\tgenus\tngenus\tspecies\tnspecies\tavgpid' >> {output}
        find stats/*/*/diamond_report/ -name "bin*.tsv" | xargs -I {{}} cat {{}} >> {output}
    	"""

rule diamond_bin_summary_plus:
    input: "stats/diamond_bin_report.tsv"
    output: "stats/diamond_bin_report_plus.tsv"
    conda: "../envs/ete3.yaml"
    script: "../scripts/add_tax_diamond.py"


rule sourmash_sig:
    input:
        mags = 'mags/{assembly_group}/{id}.fa',
        fini = "FinishedBinning"
    output: 'stats/sourmash/{assembly_group}/{id}.sig'
    conda: "../envs/sourmash.yaml"
    shell: "sourmash compute --scaled 1000 -k 31 -o {output} {input.mags}"

rule sourmash_gather:
    input: 'stats/sourmash/{assembly_group}/{id}.sig'
    output:
        csv='stats/sourmash/{assembly_group}/{id}.csv',
        out='stats/sourmash/{assembly_group}/{id}.sm'
    params:
        gb=config["sourmash_gbk"]
    conda: "../envs/sourmash.yaml"
    shell: "sourmash gather -k 31 {input} {params.gb} -o {output.csv} > {output.out}"

rule sourmash_report:
    input:
        csvs = expand("stats/sourmash/{assembly_group}/{id}.csv", assembly_group=getAssemblyBaseName(config["assemblies"]), id=getIds())
    output:
        file = 'stats/sourmash/sourmash_report.csv'
    run:
        import os
        import re
        with open(output.file, 'w') as out:
            out.write("assembly,file,intersect_bp,f_orig_query,f_match,f_unique_to_query,f_unique_weighted,average_abund,median_abund,std_abund,name,filename,md5\n")
            for i in input.csvs:
                fname = re.sub('\.csv', '', os.path.basename(i))
                assembly = os.path.abspath(i).split('/')
                with open(i, 'r') as text:
                    text.readline()
                    for l in text:
                        out.write(f'{assembly[-2]},{fname},' + l)
