#decompress fastq
#Samuel Ahuno
#April 1st, 2021
#Paz Polak lab
#Mount Sinai,NY

"""
#command to run snakemake (remove -np at end when done validating):
snakemake -s compressFastq.snakefile --latency-wait 60 --restart-times 2 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 100 -np
"""

configfile: "config/samples.yaml"

rule all:
    input:
        expand("results/{samples}/{samples}.gz", samples=config["samples"])

rule compressfastq:
    input:
        fastq1=lambda wildcards: config["samples"][wildcards.samples]
    output:
         fastq_gz="results/{samples}/{samples}.gz"
    log:
        "logs/compressfastq/{samples}_compressfastq.txt"
    shell:
        "gzip {input.fastq1} > {output.fastq_gz} 2> {log}"