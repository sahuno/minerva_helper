#demultiplex_fastq
#Samuel Ahuno
#April 1st, 2021
#Paz Polak lab
#Mount Sinai,NY

"""
#command to run snakemake (remove -np at end when done validating):
snakemake -s compressFastq.snakefile --latency-wait 60 --restart-times 2 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 100 -np
"""

configfile: "config/input_paths_fastq.yaml"

rule all:
    input:
        expand("results/{samples}/{samples}.fastq", samples=config["samples"])

rule demultiplex_fastq:
    input:
        fastq1=lambda wildcards: config["samples"][wildcards.samples]
    output:
         fastq_gz="results/{samples}/{samples}.fastq"
    log:
        "logs/demultiplex_fastq/{samples}_demultiplexed_fastq.txt"
    shell:
        """
       zcat {input.fastq1} | 
        awk -v sname="{wildcards.samples}" '
            BEGIN {FS = ":"} 
            {
                lane=$4;
                fileName = sname".lane."lane".fastq" 
                print > fileName
                for (i = 1; i <= 3; i++) {
                    getline
                    print > fileName
                }
            }' 2> {log}
        """
