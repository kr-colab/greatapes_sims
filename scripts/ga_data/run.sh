srun -p kern -A kernlab --mem 2GB --ntasks 1 --time 15-00:00:00 --pty bash
conda activate apes
snakemake -p --snakefile pidxy.snake --profile slurm --configfile stats.yaml --latency-wait 30
