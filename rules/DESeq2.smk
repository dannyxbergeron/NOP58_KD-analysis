rule DESeq2_genes:
    """ Differential expression for the different conditions """
    input:
        counts = 'results/kallisto/est_counts.tsv',
        samples = "data/design.tsv"
    output:
        results = directory("results/DESeq2/genes"),
    log:
        "logs/DESeq2/genes.log"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/DESeq2_genes.R"

rule DESeq2_transcripts:
    """ Differential expression for the different conditions """
    input:
        counts = expand("results/kallisto/{id}/abundance.h5",
                            id=simple_id),
        samples = "data/design.tsv"
    output:
        results = directory("results/DESeq2/transcripts")
    params:
        names = expand('{id}', id=simple_id)
    log:
        "logs/DESeq2/transcripts.log"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/DESeq2_transcripts.R"

rule rename:
    input:
        genes = "logs/DESeq2/genes.log",
        transcripts = "logs/DESeq2/transcripts.log"
    output:
        tok = "logs/DESeq2/rename.tok"
    params:
        genes_dir = "results/DESeq2/genes/",
        transcripts_dir = "results/DESeq2/transcripts/",
        samples = "data/design.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/rename_diff.py"
