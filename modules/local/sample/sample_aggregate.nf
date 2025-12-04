process SAMPLE_AGGREGATE {
    tag "${output_file}"
    label 'process_low'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    path csv_files
    val output_file

    output:
    path output_file, emit: aggregated_csv

    script:
    """
    python3 <<EOF
    import pandas as pd

    input_files = [${csv_files.collect { '"' + it.getName() + '"' }.join(', ')}]
    dfs = [pd.read_csv(f) for f in input_files]
    merged = pd.concat(dfs, axis=0, ignore_index=True)
    merged.to_csv("${output_file}", index=False)
    EOF
    """
}