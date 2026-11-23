process CONVERGENCE {
    tag "${sample_meta.sample}"
    label 'process_low'

    input:
    tuple val(sample_meta), path(count_table)

    output:
    path "${count_table.baseName}_tcr_convergence.tsv", emit: "convergence_output"

    script:
    """
    # Extract vector of cdr3 aa, dropping null values
    python - <<EOF
    import pandas as pd

    # Load your TCR data (make sure the file has 'cdr3_aa' and 'cdr3_nt' columns)
    df = pd.read_csv("${count_table}", sep="\t", usecols=["junction_aa", "sequence"])
    df = df.dropna(subset=["junction_aa"])

    # Group by amino acid sequence and count unique nucleotide sequences (convergence)
    convergence_df = (
        df.groupby("junction_aa")["sequence"]
        .nunique()
        .reset_index(name="convergence")
    )

    # Sort by convergence count, descending
    convergence_df = convergence_df.sort_values(by="convergence", ascending=False)

    # Export
    convergence_df.to_csv("${count_table.baseName}_tcr_convergence.tsv", sep="\t", index=False)
    EOF
    """
}
