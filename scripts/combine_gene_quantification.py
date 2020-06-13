import pandas as pd

# Creating tpm and est_counts datasets
matrix = dict()
quants = snakemake.output.keys()
for quant in quants:
    matrix[quant] = pd.read_csv(snakemake.input.map, sep='\t', names=['transcript', 'gene'])


# For every datasets
for dataset in snakemake.input.datasets:
    id = dataset.split('/')[-2].split('/')[0]

    data = pd.read_csv(dataset, sep='\t')
    data.set_index('target_id', inplace=True)

    for quant in quants:
        quant_type = ""
        if 'tpm' in quant:
            quant_type = 'tpm'
        else:
            quant_type = 'est_counts'

        _dict = data[quant_type].to_dict()
        matrix[quant][id] = matrix[quant_type]['transcript'].map(_dict)

# Simplyfing to gene quantification
for quant in quants:
    if 'transcript' not in quant:
        matrix[quant].drop('transcript', axis=1, inplace=True)
        matrix[quant] = matrix[quant].groupby('gene').sum()
        matrix[quant].reset_index(inplace=True)
    else:
        matrix[quant].drop('gene', axis=1, inplace=True)

    # Write to file
    matrix[quant].to_csv(snakemake.output[quant], sep='\t', index=False)
