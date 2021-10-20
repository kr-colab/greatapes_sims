conda activate apes
snakemake -p --snakefile ga_sims.snake --profile slurm --restart-times 3 --latency-wait 30 --configfile chr12_strongbgs.yaml
