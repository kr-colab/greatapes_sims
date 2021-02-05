conda activate apes
module use /projects/apps/shared/modulefiles/
module load SLiM bedtools
snakemake -p --snakefile ga_sims.snake --profile slurm --restart-times 3 --latency-wait 30
