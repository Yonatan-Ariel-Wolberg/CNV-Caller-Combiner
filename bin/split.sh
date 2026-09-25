#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Detailed documentation: split.sh
# Purpose: This pipeline helper script splits a combined multi-sample read count matrix
#   generated for CANOES (`canoes.reads_new.txt`) into sample-specific text files.
#   It uses an awk script to iterate over the columns of each row, creating output
#   files named after column indices containing coordinates and sample counts.
#
# Inputs:
#   - canoes.reads_new.txt: A whitespace-separated file where:
#       - Column 1: Chromosome
#       - Column 2: Start position
#       - Column 3: End position
#       - Columns 4 to end: Read counts for each individual sample (typically spaced or structured).
#
# Outputs:
#   - <column_index>.txt: Sample-specific text files containing target regions and their
#                         respective counts for that sample (e.g. `4.txt`, `14.txt`, etc.).
# =============================================================================
awk '
    {
        # Loop through columns starting from column 4 (first sample count)
        # stepping by 10 columns (or as defined by workflow sample layout).
        for( i = 4; i <= NF; i=i+10 )
        {
        print(i)
        # Create output filename based on column index
        fn = i ".txt";
        # Write Chromosome ($1), Start ($2), End ($3) and count for sample i ($(i)) to the file
        print($1, $2, $3, $(i))>>fn;
        # Close file to avoid "too many open files" errors during long runs
        close( fn );
        }
    }
' canoes.reads_new.txt

# fn = i ".txt";
# printf( "%s %s %s\n", $1, $2, $(i) ) >>fn;
# close( fn );
##printf( "%s %s %s %s\n", $1, $2, $3, $(i) ) >>fn;
