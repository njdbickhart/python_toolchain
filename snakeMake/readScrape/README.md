# ReadScrape
---
A pipeline to gather unmapped reads from a reference genome alignment and convert them into contigs for additional analysis. 

## Requirements

* Python 3.6+
* [SPAdes](http://cab.spbu.ru/software/spades/) v3.12+
* [Centrifuge](https://ccb.jhu.edu/software/centrifuge/manual.shtml) 
* Bedtools v2.5+
* Samtools v1.6+
* BWA v1.7+
* Minimap2
* [cd-hit](http://weizhongli-lab.org/cd-hit/)

## Installation

First, clone my python_toolchain github repository:

```bash
git clone https://github.com/njdbickhart/python_toolchain.git
```

Afterwards, there are several key steps to make this pipeline suited to your cluster environment. It is HIGHLY recommended that you run this on a cluster, as this is a highly parallel process workflow!

### 1. Install all programs on your path

Install all of the above required software on your path or make it available on your path before you submit the job. The snakemake workflow assumes that everything is available on your path, and that there are no collisions with job names! 

### 2. Prepare the run parameter JSON file

When you run the pipeline, it expects a **default.json** file to be in the current directory. This file contains a list of all of the locations for key files that will be used in the rules of the pipeline. It should also contain a list of all of the locations for the bam files that will be processed. To help with creating this file, I've written a script that will help start the file:

> scripts/generateJSONForPipeline.py

Here is the usage statement:

```bash
Usage: python3 generateJSONForPipeline.py <space delimited list of folders that contain bams>
```

The script will search one additional subdirectory down for the bam files, so make sure that you specify the top-level directory on the command line. Also, the script expects that the sample name is in the bam file name and that this name is delimited by either a "." or a "_". Here are some example file names that it recognizes:

```bash
sample1.realign.dedup.sort.bam  -> is converted to: {sample1 : sample1.realign.dedup.sort.bam}
sample2_myworkflow_is_cool.bam  -> is converted to: {sample2 : sample2_myworkflow_is_cool.bam}
```

#### Note: you need to manually edit the output file!

After the script runs successfully, you will need to edit the default.json file to add additional file locations. Please use full paths when possible! Here are the files that need to be manually edited to make the pipeline work:

* repeatmasker:  (you can ignore this for now, it's not needed)
* centrifugedb: This is the location of the centrifuge database you will use to screen the scaffolds
* taxids: This is the location of all NCBI taxonomic IDs that you want to keep after centrifuge screening
* logdir: This is the location of the log directory that you will use to monitor the tasks
* samples: (this should be populated by the script, but manually verify that the file locations and contents are correct before running the pipeline)

### 3. Prepare the cluster JSON file

Because your cluster parameters are likely different than mine, you may need to edit the cluster.json parameter file. It is located in the base directory of this folder, but you can copy it to your working directory and edit it as needed. 

If you are using SLURM, most of the JSON categories should make sense. In that case,  you just need to edit the "partition" and add in any other specifications that your cluster admin requires. Otherwise, here are the important parameters to check and confirm:

* mem : This is the amount of RAM in megabytes
* ntasks-per-node : This is the number of tasks per node (ie. cpus)
* partition : We run on a partition system that requires us to specify this -- remove if needed!

After you prepare all of this, you should be ready to run!

## Running the pipeline.

Snakemake requires a bit of tinkering to work on a cluster (for more information, [see this](https://snakemake.readthedocs.io/en/v5.1.4/executable.html#cluster-execution)). In order to run this pipeline on slurm, you would execute the following command:

```bash
sbatch --nodes=1 --mem=5000 --ntasks-per-node=2 -p msn -q msn  \			# submit as a normal batch job
snakemake \																	# The snakemake executable
--cluster-config ~/python_toolchain/snakeMake/readScrape/cluster.json   \ 	# Your cluster configuration file
--cluster "sbatch --nodes={cluster.nodes} --ntasks-per-node={cluster.ntasks-per-node} --mem={cluster.mem} --partition={cluster.partition} -q msn -o {cluster.stdout}"    \				# Cluster command for each job. The suffices correspond to the cluster.json file
--jobs 999 \																# max number of jobs that the cluster can maintain per user
-s ~/python_toolchain/snakeMake/readScrape/readScrape						# The snakemake file
```

You will need to modify the cluster submission text to accommodate your cluster conditions and job submission executable.

## Directory and file structure

If the pipeline runs successfully, it will generate the following folders in the current directory. The beginning portion of the pipeline will have separate sample-level subfolders that contain intermediary files. 

```bash
   .
   |-rawreads			<- folder containing raw reads from "extract_reads"
   |---SAMPLE
   |
   |-assembly			<- folder with assembled scaffolds from "run_spades"
   |---SAMPLE
   |-----corrected		<- Indexed, renamed assembly files
   |
   |-classification		<- folder with scaffold-level centrifuge classification data
   |---SAMPLE
   |
   |-filtered			<- scaffolds that were selected by "filter_centrifuge"
   |---SAMPLE
   |
   |-linkage			<- list of reads and metadata that link scaffolds to chromosomes
   |---SAMPLE
   |
   |-collation			<- Clustered scaffolds for use in genotyping
   |
   |-genotypes			<- Final folder, containing genotype and location information on each scaffold
```

### Final files

The important files that summarize the output of the pipeline can be found in the "collation" and "genotypes" folders. I will describe them in order of their generation:

* Collation
	* **main_c_file.cdhit	**	<-	Fasta file containing the best representative sequence for each clustered novel sequence
	* **main_c_file.cdhit.clstr** <- Cd-Hit cluster file showing which scaffolds were found to cluster in each cluster designation

* Genotypes
	* **genotype_summary.tab**	<- (has header) lists cluster statistics and sample presence
	* **association_counts.tab** <- (limited header) variable column length file that shows linkage counts for chromosome assignment of each cluster
	* **(sample).genotype.tab** <- (has header; printed for each sample) has list of filtered novel sequence scaffolds that belong to that sample, along with possible chromosome assignment