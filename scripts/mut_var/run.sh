conda activate apes
snakemake -p --snakefile mutvar.snake --profile slurm --restart-times 3 --latency-wait 30
