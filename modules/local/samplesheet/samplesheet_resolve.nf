process SAMPLESHEET_RESOLVE {
    label 'process_single'

    input:
    path samplesheet_utf8
    val(resolved_rows)     // List of tab-separated strings
    val(resolved_header)   // Comma-separated header line

    output:
    path "samplesheet_resolved.csv", emit: samplesheet_resolved

    script:
    """
# Write resolved rows to a temp file
cat << 'EOF' > resolved.tmp
${resolved_rows.join('\n')}
EOF

# Emit header
cat << 'EOF' > samplesheet_resolved.csv
${resolved_header}
EOF

# Two-pass awk:
#  - pass 1: read original samplesheet, store sample order
#  - pass 2: read resolved rows, store rows by sample
awk -F',' '
    NR==FNR {
        if (FNR > 1) order[++n] = \$1
        next
    }
    {
        resolved[\$1] = \$0
    }
    END {
        for (i = 1; i <= n; i++) {
            s = order[i]
            if (!(s in resolved)) {
                printf "ERROR: missing resolved row for %s\\n", s > "/dev/stderr"
                exit 1
            }
            print resolved[s]
        }
    }
' "${samplesheet_utf8}" resolved.tmp >> samplesheet_resolved.csv
    """
}