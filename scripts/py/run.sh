conda activate greatapes
module use /projects/apps/shared/modulefiles/
module load SLiM/dev bedtools
snakemake -p --snakefile chr12.snake --profile slurm --jobs 2000
