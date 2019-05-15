### Bash cell from python notebook

        %%bash
        # comments:
        # skip this cell, included as a reference only
        # depending on the workflow, you can either use a single sample or concatenate multiple samples
        # but you MUST keep read pairs separate
        cat 39526_Mamu-exon*_matches_r1* > 39526_Mamu-exonAll_allR1.txt
        cat 39526_Mamu-exon*_matches_r2* > 39526_Mamu-exonAll_allR2.txt
