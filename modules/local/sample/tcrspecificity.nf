process VDJDB_GET {
    label 'process_low'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    output:
    path("vdjdb-2025-02-21/"), emit: ref_db

    script:
    """
    wget https://github.com/antigenomics/vdjdb-db/releases/download/pyvdjdb-2025-02-21/vdjdb-2025-02-21.zip
    unzip vdjdb-2025-02-21.zip
    """
}

process VDJDB_VDJMATCH {
    tag "${sample_meta.sample}"
    label 'process_medium'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    tuple val(sample_meta), path(count_table)
    path(ref_db)

    output:
    path("${sample_meta.sample}.vdjmatch.txt")
    path("${sample_meta.sample}.annot.summary.txt")
    path "logs/${sample_meta.sample}.vdjmatch.log"

    script:
    def memGb = (task.memory.toMega() * 0.8 / 1024).intValue()
    
    """
    python - <<EOF
    import pandas as pd

    df = pd.read_csv("${count_table}", sep="\t")
    rename_map = {
            "duplicate_count": "cloneCount",
            "duplicate_frequency_percent": "cloneFraction",
            "sequence": "cdr3nt",
            "junction_aa": "cdr3aa",
            "v_call": "v",
            "d_call": "d",
            "j_call": "j"
        }
    df.rename(columns=rename_map, inplace=True)
    final_columns = ['cloneCount', 'cloneFraction', 'cdr3nt', 'cdr3aa', 'v', 'd', 'j']
    df = df[final_columns]

    df.to_csv("vdjmatch.tsv", sep="\t", index=False, na_rep='NA')
    EOF

    mkdir -p logs

    java -Xmx${memGb}G -jar /usr/local/bin/vdjmatch.jar match \
        -S human -R TRB --search-scope 2,2,2 --vdjdb-conf 0 --database "${ref_db}/vdjdb" \
        "vdjmatch.tsv" "${sample_meta.sample}" \
        > logs/${sample_meta.sample}.vdjmatch.log 2>&1
    """
}