conda activate greatapes
module use /projects/apps/shared/modulefiles/
module load SLiM bedtools
snakemake -p --snakefile ga_sims.snake --profile slurm --jobs 1000 --cluster-config talapas.yaml
