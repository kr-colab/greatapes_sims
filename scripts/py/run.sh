conda activate greatapes
module use /projects/apps/shared/modulefiles/
module load SLiM bedtools
snakemake -p --snakefile chr.snake --profile slurm --jobs 100 --cluster-config talapas.yaml
snakemake -p --snakefile windowed.snake --profile slurm --jobs 100 --cluster-config talapas.yaml
