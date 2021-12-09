conda activate apes
snakemake -p --snakefile ga_sims.snake --profile slurm --restart-times 6 --latency-wait 30 --configfile partialsweeps.yaml --nolock
snakemake -p --snakefile ga_sims.snake --profile slurm --restart-times 6 --latency-wait 30 --configfile chr12_strongbgs.yaml --nolock
