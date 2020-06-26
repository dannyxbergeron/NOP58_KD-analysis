from collections import Counter

import pandas as pd

from snakemake.shell import shell

kallisto_input = snakemake.input.kallisto
coco_input = snakemake.input.coco

def load_dfs():
    dfs = {}

    # Load dfs
    kallisto_df = pd.read_csv(kallisto_input, sep='\t')
    coco_df = pd.read_csv(coco_input, sep='\t')

    # Edit dfs
    kallisto_df.columns = [x.replace('gene', 'gene_name') for x in kallisto_df.columns]
    coco_df.drop(columns=['gene_id'], inplace=True)

    # To remove de snoRNAs with duplicates names
    counter = Counter(list(coco_df.gene_name))
    dups = [x for (x, y) in counter.most_common() if y > 1]
    # print(dups)

    kallisto_df = kallisto_df.loc[~(kallisto_df.gene_name.isin(dups))]
    coco_df = coco_df.loc[~(coco_df.gene_name.isin(dups))]

    # sort the dfs by gene_name
    kallisto_df.sort_values('gene_name', inplace=True)
    coco_df.sort_values('gene_name', inplace=True)

    # print(len(dfs['kallisto']))
    # print(len(dfs['coco']))

    dfs['kallisto'] = kallisto_df
    dfs['coco'] = coco_df
    return dfs


def main():

    dfs = load_dfs()

    for name, df in dfs.items():
        print(name)
        print(df.loc[df.gene_name == 'NOP58'])

    shell('touch {}'.format(snakemake.output.tok))



if __name__ == '__main__':
    main()
