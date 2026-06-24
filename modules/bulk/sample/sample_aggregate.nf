process SAMPLE_AGGREGATE {
    tag "${output_file}"
    label 'process_low'
    publishDir "${params.outdir}/sample", mode: 'copy', overwrite: true

    input:
    path csv_files
    val output_file

    output:
    path output_file, emit: aggregated_csv

    script:
    """
    cat > aggregate.py <<EOF
    import sys
    import pandas as pd

    input_files = sys.argv[1:]
    dfs = [pd.read_csv(f) for f in input_files]
    merged = pd.concat(dfs, axis=0, ignore_index=True)
    merged.to_csv("${output_file}", index=False)
    EOF

    python3 aggregate.py ${csv_files}
    """
}