#!/usr/bin/env python3

import sys
import json
from pathlib import Path

# The script now takes two simple arguments:
# 1. The input JSON file
# 2. The output CSV filename
meta_json_file = sys.argv[1]
output_csv = sys.argv[2]

# Load the list from the JSON file
with open(meta_json_file, 'r') as f_in:
    # This will be a list, e.g.,
    # [ [ {'sample': 'Patient01_Base_CD4', ...}, 'file1.tsv' ],
    #   [ {'sample': 'Patient01_Base_CD8', ...}, 'file2.tsv' ] ]
    all_new_meta = json.load(f_in)

with open(output_csv, 'w') as f_out:
    
    # Write the new header
    f_out.write('sample,subject_id,timepoint,origin,phenotype,file\n')

    # Iterate over each item in the list
    for item in all_new_meta:
        meta = item                   # 'item' is the meta map
        file_path = str(meta['file']) # Get the file path from the map

        f_out.write(
            f"{meta['sample']},"
            f"{meta['subject_id']},"
            f"{meta['timepoint']},"
            f"{meta['origin']},"
            f"{meta['phenotype']},"
            f"{file_path}\n"
        )
