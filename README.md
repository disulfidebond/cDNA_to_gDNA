# cDNA_to_gDNA
Work on reconstructing gDNA from cDNA Macca mulatta MHC sequences

Current Workflow is:
1) preprocess IPD file and create multifasta file of allele:exon
2) import cdna2gdna.py
3) run combine_output_using_bash.md
4) run parse_check_output.md
5) collect sequences where both pairs of reads mapped, run _de novo_ assembly
6) use draft consensus from step 5 with unpaired sequences, then re-run _de novo_ assembly
7) test results by mapping against gDNA, both published Mamu10 and SRA sequences.  Accuracy should be > 90% in both cases.

Parallel tasks:
* retrieve exome sequences from SRA, done by scanning pre-downloaded SRR IDs:

        # not pretty, but it does the job, parses based on file size
        # run from directory of files with the format SRRXXXXX_1.fastq
        ls -lah | grep '\.[[:digit:]]G*' | rev | cut -d. -f2 | cut -d\  -f1 | rev | cut -d_ -f1 | sort -n | uniq > listOfSRIDs.txt
        
* These are genomic reads, and can be used for subsequent downstream steps to find gDNA. Retrieve from SRA using curl or similar tool.
