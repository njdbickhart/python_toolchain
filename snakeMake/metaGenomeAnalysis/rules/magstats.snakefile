# Note: many of these rules were refactored from MagPy
# Credit to: https://github.com/WatsonLab/MAGpy
import os

#TODO: refactor for my pipeline

rule checkm:
	input: "mags"
	output: "checkm.txt"
	threads: 16
	conda: "envs/checkm.yaml"
	params:
		cdr=config["checkm_dataroot"]
	shell:
		'''
		checkm_db={params.cdr}
		echo ${{checkm_db}} | checkm data setRoot ${{checkm_db}}
		checkm lineage_wf -f checkm.txt --reduced_tree -t {threads} -x fa {input} ./checkm
		'''

rule checkm_plus:
	input: "checkm.txt"
	output: "checkm_plus.txt"
	threads: 1
	conda: "envs/ete3.yaml"
	shell: "scripts/add_tax.py {input} > {output}"

rule prodigal:
	input: 'mags/{id}.fa'
	output:
		faa='proteins/{id}.faa',
		gff='proteins/{id}_prodigal.gff'
	conda: "envs/prodigal.yaml"
	shell: 'prodigal -p meta -a {output.faa} -q -i {input} -f gff -o {output.gff}'

rule diamond:
        input: 'proteins/{id}.faa'
        output: 'diamond/{id}.diamond.tsv'
        threads: 16
        params:
                db=config["uniprot_sprot"],
                of="6 qseqid sseqid stitle pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore"
        conda: "envs/diamond.yaml"
	shell: "diamond blastp --threads {threads} --max-target-seqs 10 --db {params.db} --query {input} --outfmt {params.of} --out {output}"

rule diamond_report:
	input:
		tsv='diamond/{id}.diamond.tsv',
		faa='proteins/{id}.faa'
	output: 'diamond_report/bin.{id}.tsv', 'diamond_report/con.{id}.tsv'
	params:
		outdir="diamond_report"
	conda: "envs/bioperl.yaml"
	shell: "scripts/diamond_report.pl {input.tsv} {input.faa} {params.outdir}"

rule diamond_bin_summary:
        input: expand("diamond_report/bin.{sample}.tsv", sample=IDS)
        output: "diamond_bin_report.tsv"
	shell:
                '''
		echo -e 'name\tnprots\tnhits\tnfull\tgenus\tngenus\tspecies\tnspecies\tavgpid' >> {output}
            	find diamond_report/ -name "bin*.tsv" | xargs -I {{}} cat {{}} >> {output}
		'''

rule diamond_bin_summary_plus:
        input: "diamond_bin_report.tsv"
        output: "diamond_bin_report_plus.tsv"
	conda: "envs/ete3.yaml"
        shell: "scripts/add_tax_diamond.py {input} > {output}"


rule sourmash_sig:
        input: 'mags/{id}.fa'
        output: 'sourmash/{id}.sig'
	conda: "envs/sourmash.yaml"
        shell: "sourmash compute --scaled 1000 -k 31 -o {output} {input}"

rule sourmash_gather:
        input: 'sourmash/{id}.sig'
        output:
                csv='sourmash/{id}.csv',
                out='sourmash/{id}.sm'
	params:
		gb=config["sourmash_gbk"]
	conda: "envs/sourmash.yaml"
	shell: "sourmash gather -k 31 {input} {params.gb} -o {output.csv} > {output.out}"

rule sourmash_report:
	input: expand("sourmash/{sample}.csv", sample=IDS)
	output: 'sourmash_report.csv'
	shell: "echo 'intersect_bp,f_orig_query,f_match,f_unique_to_query,f_unique_weighted,average_abund,median_abund,std_abund,name,filename,md5' >> {output} && scripts/sourmash_report.pl {input} >> {output}"
