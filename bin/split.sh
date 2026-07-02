# =============================================================================
# Detailed documentation: split.sh
# Purpose: This pipeline helper script supports CNV calling and/or post-processing tasks in
#   this project, with parameters and I/O conventions aligned to repository workflows.
# Scope: General utility execution.
# Inputs: CLI arguments, repository config/params files, and sample-level CNV data.
# Outputs: Generated intermediate artifacts, final result files, and/or run metadata.
# Operational notes: Keep parameter names synchronized with related workflow/process
#   definitions and update this file alongside upstream interface changes.
# =============================================================================
awk '
    {
        for( i = 4; i <= NF; i=i+10 )
        {
        print(i)
        fn = i ".txt";
        print($1, $2, $3, $(i))>>fn;
        close( fn );
        }
    }
' canoes.reads_new.txt

# fn = i ".txt";
# printf( "%s %s %s\n", $1, $2, $(i) ) >>fn;
# close( fn );
##printf( "%s %s %s %s\n", $1, $2, $3, $(i) ) >>fn;
