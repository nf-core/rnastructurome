#!/usr/bin/awk -f
# Author: Victoria Begley (begley@ebi.ac.uk)
# License: MIT
# Parse the rf-count / rf-count-genome console summary table into a one-row TSV.
#
# Variables (pass with -v):
#   sample      sample prefix to match on the summary's first column
#   is_map      "1" for MaP (mutation) mode, "0" for RT-stop; selects the column layout
#   match_mode  "exact"  -> row where $1 == sample, emit $1 as the sample id (rf-count)
#               "prefix" -> row where $1 starts with sample, emit sample as the id
#                           (rf-count-genome, whose $1 carries the staged BAM's suffix)
#
# Emits a header row plus the last matching data row (rf-count prints one summary line
# per sample; the last match is the completed one).

BEGIN {
    if (is_map == "1")
        print "sample\tcovered\tmutated_alignments\tpct_mutated\tpct_a_muts\tpct_c_muts\tpct_g_muts\tpct_u_muts"
    else
        print "sample\tcovered\tpct_a_muts\tpct_c_muts\tpct_g_muts\tpct_u_muts"
}

{
    if (match_mode == "prefix") { matched = (index($1, sample) == 1); id = sample }
    else                        { matched = ($1 == sample);          id = $1 }
    if (!matched) next

    if (is_map == "1") {
        if (index($3, "/") > 0) {
            # MaP: "<mutated>/<total> (<pct>%)" — $3=mutated/total, $4=(pct%), $5-$8=%A/C/G/U
            pct_mut = $4; gsub("[()%]", "", pct_mut)
            row = id "\t" $2 "\t" $3 "\t" pct_mut "\t" $5 "\t" $6 "\t" $7 "\t" $8
        } else {
            # MaP with zero counted alignments: rf-count prints "-" ($3), bases at $4-$7
            row = id "\t" $2 "\t" $3 "\t" "NA" "\t" $4 "\t" $5 "\t" $6 "\t" $7
        }
    } else {
        # RT-stop: no Mutated-alignments column; $3-$6 = per-base stop percentages
        row = id "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6
    }
}

END { if (row != "") print row }
