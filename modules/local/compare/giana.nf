process GIANA_CALC {
    tag "${patient}"
    label 'process_medium'

    input:
    tuple val(patient), path(concat_cdr3)
    val threshold
    val threshold_score
    val threshold_vgene

    output:
    path "${patient}_VgeneScores.txt", emit: 'vgene_scores'
    path "${patient}_giana.txt", emit: 'giana_output'
    // path "giana_EncodingMatrix.txt"

    script:
    def dedup_clonotypes = (params.giana_dedup_clonotypes ?: false) ? 'True' : 'False'
  
    """
    python3 - <<EOF
    import pandas as pd
    df = pd.read_csv("${concat_cdr3}", sep="\t")
    df = df.rename(columns={"junction_aa": "CDR3b", "v_call": "TRBV"})

    # PATIENT_CONCATENATE pools a patient's samples by stacking rows, so one clonotype seen
    # in N samples arrives as N identical (CDR3b, TRBV, j_call) rows and GIANA reports those
    # duplicates as clusters. On an 8-sample single-cell test set that was 76 of 77
    # "clusters", and patients with a single sample produced no output at all.
    #
    # Off by default so existing bulk results are unchanged; the single-cell route enables it
    # via params.giana_dedup_clonotypes. Bulk can opt in after review.
    if ${dedup_clonotypes}:
        key = [c for c in ("CDR3b", "TRBV", "j_call") if c in df.columns]
        if key:
            n_before = len(df)
            agg = {}
            if "duplicate_count" in df.columns:
                agg["duplicate_count"] = "sum"
            if "sample" in df.columns:
                agg["sample"] = lambda s: ";".join(sorted(set(s.astype(str))))
            for c in df.columns:
                if c not in key and c not in agg:
                    agg[c] = "first"
            cols = list(df.columns)
            df = df.groupby(key, as_index=False, sort=False).agg(agg)
            if "duplicate_frequency_percent" in df.columns and "duplicate_count" in df.columns:
                total = df["duplicate_count"].sum()
                if total:
                    df["duplicate_frequency_percent"] = 100.0 * df["duplicate_count"] / total
            df = df[cols]
            print(f"[GIANA] collapsed {n_before} rows -> {len(df)} unique clonotypes", flush=True)

    df.to_csv("concat_cdr3_renamed.tsv", sep="\t", index=False)
    EOF
    
    GIANA --file "concat_cdr3_renamed.tsv" \
        --output . \
        --outfile giana_RotationEncodingBL62.txt \
        --EncodingMatrix true \
        --threshold ${threshold} \
        --threshold_score ${threshold_score} \
        --threshold_vgene ${threshold_vgene} \
        --NumberOfThreads ${task.cpus} \
        --Verbose

    # Insert header after GIANA comments
    insert=\$(head -n 1 "concat_cdr3_renamed.tsv")
    insert=\$(echo "\$insert" | awk -F'\t' 'BEGIN{OFS="\t"} {
        out = \$1 OFS "cluster"
        for (i=2; i<=NF; i++) {
            out = out OFS \$i
        }
        print out
    }')

    awk -v insert="\$insert" '
    /^##/ { print; next }
    !inserted { print insert; inserted=1 }
    { print }
    ' giana_RotationEncodingBL62.txt > ${patient}_giana.txt

    # mv giana_RotationEncodingBL62.txt_EncodingMatrix.txt giana_EncodingMatrix.txt
    mv VgeneScores.txt ${patient}_VgeneScores.txt
    """
}
