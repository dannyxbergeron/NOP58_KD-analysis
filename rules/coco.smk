rule create_coco_annotation:
    input:
        gtf = config["path"]['annotation']
    output:
        correct_gtf = config["path"]['correct_annotation']
    conda:
        "../envs/coco.yaml"
    shell:
        "coco correct_annotation {input.gtf}"

rule coco_cc:
    input:
        annotation = config["path"]['correct_annotation'],
        bam_file = "results/STAR/{id}/Aligned.sortedByCoord.out.bam"
    output:
        out_file = "results/coco/{id}.tsv"
    log:
        "logs/coco/{id}.log"
    threads:
        16
    conda:
        "../envs/coco.yaml"
    shell:
        "coco cc -s 1 "
        "-t {threads} "
        "-p {input.annotation} "
        "{input.bam_file} "
        "{output.out_file}"

rule coco_merge:
    input:
        tpm_files = expand("results/coco/{id}.tsv", id=simple_id)
    output:
        merged = "results/coco/merged/tpm.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/combine_coco_quantification.py"
