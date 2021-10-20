import pandas as pd
import subprocess
import argparse

parser = argparse.ArgumentParser(
    description='Splits species level VCFs into subspecies.',
)

parser.add_argument('metadata_path', type=str, help='Path to a .tsv file with at least two columns (species and subspecies).')
parser.add_argument("vcf_path", type=str, help='Path where species VCFs reside.')
parser.add_argument("out_path", type=str, help='Path where splitted VCFs will be saved')
parser.add_argument("--species-col", "-spp", type=str, default="spp", help='Species column name.')
parser.add_argument("--subspecies-col", "-subspp", type=str, default="subspp", help='Subspecies column name. Note: rows encoded with other as subspecies are ignored.')
parser.add_argument("--sample-col", type=str, default="sample_name", help='Subspecies column name.')

args = parser.parse_args()

meta = pd.read_csv(args.metadata_path, sep="\t", encoding = "ISO-8859-1")

for subspp in meta[args.subspecies_col].unique():
    if subspp == 'other':
        continue
    print(f'Preparing VCF for subspecies {subspp}.')
    subset = meta[meta[args.subspecies_col]==subspp]
    origin = args.vcf_path+subset[args.species_col].iloc[0]+".vcf.gz"
    samples = ','.join(subset[args.sample_col])
    subprocess.Popen(f'srun --account kernlab --partition kern --mem 16G --time 24:00 --cpus-per-task=4 bcftools view -s {samples} {origin} -O z -o {args.out_path}{subspp}.vcf.gz --threads 4', shell=True)

