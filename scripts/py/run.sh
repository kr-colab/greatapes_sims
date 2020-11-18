#srun -p kern -A kernlab --mem 8GB --ntasks 2 --time 25-00:00:00 --pty bash
conda activate apes
module use /projects/apps/shared/modulefiles/
module load SLiM bedtools
snakemake -p --snakefile ga_sims.snake --profile slurm --restart-times 3
