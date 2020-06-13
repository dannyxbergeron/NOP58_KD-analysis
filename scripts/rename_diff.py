import os
import pandas as pd
from snakemake.shell import shell

genes_dir = snakemake.params.genes_dir
transcripts_dir = snakemake.params.transcripts_dir

design = snakemake.params.samples

for dir in [genes_dir, transcripts_dir]:

    diff_files = [x for x in os.listdir(dir) if x.endswith('.csv')]

    df = pd.read_csv(design, sep=' ')
    df['sample'] = df['sample'].str[:-1]
    df = df.drop_duplicates()

    eq = dict(zip(df['condition'], df['sample']))
    print(eq)

    for n in diff_files:
        n_list = n.replace('.csv', '').split('-')
        print(n_list[0], n_list[1])
        new_name = '{}-{}.csv'.format(eq[n_list[0]], eq[n_list[1]])
        print(n, new_name)
        os.rename(os.path.join(dir, n), os.path.join(dir, new_name))

shell("touch {}".format(snakemake.output.tok))
