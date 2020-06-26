rule heatmap_analysis:
    input:
        kallisto = "results/kallisto/tpm.tsv",
        coco = "results/coco/merged/tpm.tsv"
    output:
        tok = "tok/heatmap_analysis.tok"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/coco_vs_kallisto.py"
