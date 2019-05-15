## Overview
The steps involved are:
1) quality trim fq files, and convert fastq to fasta files
2) Split the number of reads for each converted fasta file roughly in half


        # step 1
        ARR=($(ls | grep gz | cut -d_ -f1 | sort -n | uniq))
        for i in "${ARR[@]}" ; do 
          PVAR='/Volumes/archiveGSFS2/Baylor_10_11_all/' ; 
          OUTP='/Volumes/MgSpace/seqC_workspace/' ;
          V=($(ls /Volumes/archiveGSFS2/Baylor_10_11_all | grep "$i")) ; 
          V1=$(echo "${V[0]}") ; 
          V2=$(echo "${V[1]}") ; 
          R1="${PVAR}${V1}" ; 
          R2="${PVAR}${V2}" ; 
          FA1="${OUTP}${V1}" ; 
          FA2="${OUTP}${V2}" ; 
          FAOUT1=$(echo "${FA1}" | cut -d. -f1) ; 
          FAOUT2=$(echo "${FA2}" | cut -d. -f1) ; 
          FAOUT1=$(echo "${FAOUT1}.qtrim.1.fasta") ; 
          FAOUT2=$(echo "${FAOUT2}.qtrim.2.fasta") ; 
          reformat.sh trimq=20 qtrim=t fastawrap=160 in=$R1 in2=$R2 out=$FAOUT1 out2=$FAOUT2
          sleep 1
          echo "Step 1 of parsing $i finished"
        done
        
        # step 2: length trim
        # repeat 2x for length trim, set fastawrap=160 reads=20000000, then fastawrap=160 skipreads=20000000
        
* These steps can be completed in a HTC environment with the following DAG. Start at **step a**. Stop indicates completion of task, or an error.


Stop <-- step b: forward trim <-- **step a: quality trim** --> step c: reverse trim --> Stop

                                      
