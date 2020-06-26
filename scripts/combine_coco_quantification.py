import pandas as pd

input_files = snakemake.input.tpm_files
output_file = snakemake.output.merged

for i, file in enumerate(input_files):
    sample_name = file.replace('.tsv', '').split('/')[-1]

    if i == 0:
        master_df = pd.read_csv(file, sep='\t')
        master_df.drop(columns=['count', 'cpm'], inplace=True)
        master_df.columns = ['gene_id', 'gene_name', sample_name]
    else:
        df = pd.read_csv(file, sep='\t')
        sample_dict = dict(zip(df.gene_id, df.tpm))

        master_df[sample_name] = master_df.gene_id.map(sample_dict)

master_df.to_csv(output_file, sep='\t', index=False)
