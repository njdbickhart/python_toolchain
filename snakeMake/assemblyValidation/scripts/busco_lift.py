"""Snakemake wrapper for BUSCO assessment"""

__author__ = "Tessa Pierce"
__copyright__ = "Copyright 2018, Tessa Pierce"
__email__ = "ntpierce@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell
from os import path

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")
mode = snakemake.params.get("mode")
assert mode is not None, "please input a run mode: genome, transcriptome or proteins"
lineage = snakemake.params.get("lineage_path")
assert lineage is not None, "please input the path to a lineage for busco assessment"

# busco does not allow you to direct output location: handle this by moving output
outdir = path.dirname(snakemake.output[0])
out_name = "btemp_" + snakemake.params.get("asm", "t")

# note: --force allows snakemake to handle rewriting files as necessary
# without needing to specify *all* busco outputs as snakemake outputs
shell(
    "busco --in {snakemake.input[1]} --out {out_name} --force "
    " --cpu {snakemake.threads} --mode {mode} --lineage {lineage} "
    " {extra} {log}"
)

# move to intended location
shell("cp {out_name}/short_summary*.txt {outdir}/busco_summary.txt")
shell("rm -rf {out_name}")
