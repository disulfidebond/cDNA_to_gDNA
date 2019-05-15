# cDNA_to_gDNA
Work on reconstructing gDNA from cDNA Macca mulatta MHC sequences

Current Workflow is:
1) preprocess IPD file and create multifasta file of allele:exon
2) import cdna2gdna.py
3) run combine_output_using_bash.md
4) run parse_check_output.md

Parallel tasks:
* retrieve exome sequences from SRA, done by scanning pre-downloaded SRR IDs:

        # not pretty, but it does the job, parses based on file size
        # run from directory of files with the format SRRXXXXX_1.fastq
        ls -lah | grep '\.[[:digit:]]G*' | rev | cut -d. -f2 | cut -d\  -f1 | rev | cut -d_ -f1 | sort -n | uniq > listOfSRIDs.txt
        
* These are genomic reads, and can be used for subsequent downstream steps to find gDNA. Retrieve from SRA using curl or similar tool.
