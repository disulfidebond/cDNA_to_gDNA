# cDNA_to_gDNA
Work on reconstructing gDNA from cDNA Macca mulatta MHC sequences

## Preprocessing
1) Split the fastq files using:

        # create bash array of identifiers
        ls | grep fastq.gz | cut -d_ -f1 | sort -n | uniq > idList
        ARR=($(<idList))
        
        # filter only high quality reads
        # note that combining these for loops seems to make java unhappy
        for i in "${ARRV[@]}" ; do 
          V=($(ls ${PWD} | grep "$i")) ; 
          V1=$(echo "${V[0]}") ; 
          V2=$(echo "${V[1]}") ; 
          R1="${PWD}/${V1}" ; 
          R2="${PWD}/${V2}" ; 
          FAOUT1=$(echo "${R1}" | cut -d. -f1-2) ; 
          FAOUT2=$(echo "${R2}" | cut -d. -f1-2) ; 
          FAOUT1=$(echo "${FAOUT1}.fasta.gz")
          FAOUT2=$(echo "${FAOUT2}.fasta.gz")
          reformat.sh fastawrap=160 qtrim=t in=$R1 in2=$R2 out=$FAOUT1 out2=$FAOUT2 2>> errlog.txt
        done
        
        # split fastq file, reads=20000000 is about 1/2 of the reads for this dataset
        for i in "${ARRV[@]}" ; do 
          V=($(ls ${PWD} | grep "$i")) ; 
          V1=$(echo "${V[0]}") ; 
          V2=$(echo "${V[1]}") ; 
          R1="${PWD}/${V1}" ; 
          R2="${PWD}/${V2}" ; 
          FAOUT1=$(echo "${R1}" | cut -d. -f1-2) ; 
          FAOUT2=$(echo "${R2}" | cut -d. -f1-2) ;      
          FAOUT1R1=$(echo "${FAOUT1}.1.fasta.gz") ; 
          FAOUT1R2=$(echo "${FAOUT2}.1.fasta.gz") ;  
          reformat.sh fastawrap=160 qtrim=t reads=20000000 in=$R1 in2=$R2 out=$FAOUT2R1 out2=$FAOUT2R2 2>> errlog2.txt ;
          echo "part 1 parsing $i finished" ;
        done
        
        # split fastq file, skipreads=20000000 is the remainder
        for i in "${ARRV[@]}" ; do 
          V=($(ls ${PWD} | grep "$i")) ; 
          V1=$(echo "${V[0]}") ; 
          V2=$(echo "${V[1]}") ; 
          R1="${PWD}/${V1}" ; 
          R2="${PWD}/${V2}" ; 
          FAOUT1=$(echo "${R1}" | cut -d. -f1-2) ; 
          FAOUT2=$(echo "${R2}" | cut -d. -f1-2) ;      
          FAOUT2R1=$(echo "${FAOUT1}.2.fasta.gz") ; 
          FAOUT2R2=$(echo "${FAOUT2}.2.fasta.gz") ;  
          reformat.sh fastawrap=160 qtrim=t skipreads=20000000 in=$R1 in2=$R2 out=$FAOUT2R1 out2=$FAOUT2R2 2>> errlog3.txt ;
          echo "part 1 parsing $i finished" ;
        done

2) Use [parseIPD_file.py](), which requires [pyParseIPD.py](https://github.com/disulfidebond/cDNA_to_gDNA/blob/master/preprocessing/pyParseIPD.py) and [parseIPD](https://github.com/disulfidebond/cDNA_to_gDNA/blob/master/preprocessing/parseIPD) to parse out each exon, and then create parsed file of 5' and 3' anchored sequences. A length of 125 basepairs seemed to be a good tradeoff between sensitivity and specificity. The file [cdha2gdna.py](https://github.com/disulfidebond/cDNA_to_gDNA/blob/master/cdna2gdna.py) has been created for this purpose.
3) run the python script [seqC_workflow.py](https://github.com/disulfidebond/cDNA_to_gDNA/blob/master/preprocessing/seqC_workflow.py)
4) collect sequences where either pair of reads mapped. This can be accomplished with something similar to:

        parsedFile = ParseFasta('39771_Mamu-exonAll_allR1.txt').parseFastaOutputFile()
        parsedFileObject_R1 = ParseOutputFile()
        parsedFileObject_R1.set_parseFile(parsedFile)
        parsedFileAsList_R1 = parsedFileObject_R1.get_parseFile()
        parsedDictOfReads_R1 = dict()
        for i in parsedFileAsList_R1:
            k = i[1].get_parseLine()[0]
            v = i[0].get_parseLine()
            if k not in parsedDictOfReads_R1:
                parsedDictOfReads_R1[k] = v
            else:
                tmp = parsedDictOfReads_R1[k]
                tmp.extend(v)
                tmpSet = set(tmp)
                parsedDictOfReads_R1[k] = list(tmpSet)
        print('R1\n')
        quickCheck(parsedDictOfReads_R1)

        parsedFile = ParseFasta('39771_Mamu-exonAll_allR2.txt').parseFastaOutputFile()
        parsedFileObject_R2 = ParseOutputFile()
        parsedFileObject_R2.set_parseFile(parsedFile)
        parsedFileAsList_R2 = parsedFileObject_R2.get_parseFile()
        parsedDictOfReads_R2 = dict()
        for i in parsedFileAsList_R2:
            k = i[1].get_parseLine()[0]
            v = []
            v = i[0].get_parseLine()
            if k not in parsedDictOfReads_R2:
                parsedDictOfReads_R2[k] = v
            else:
                tmp = parsedDictOfReads_R2[k]
                tmp.extend(v)
                tmpSet = set(tmp)
                parsedDictOfReads_R2[k] = list(tmpSet)
        print('\nR2\n')
        
        listOfR1_keys = parsedDictOfReads_R1.keys()
        listOfR2_keys = parsedDictOfReads_R2.keys()
        listOfR1R2_keys = []
        for k in listOfR1_keys:
            if k in listOfR2_keys:
                listOfR1R2_keys.append(k)
        for k in listOfR1R2_keys:
            check_R1List = parsedDictOfReads_R1[k]
            check_R2List = parsedDictOfReads_R2[k]
            check_R1List.sort()
            check_R2List.sort()
            if check_R1List == check_R2List:
                print(str(k) + ' matches to')
                print(parsedDictOfReads_R1[k])
            else:
                # print('# Warning!  Check:')
                # print(check_R1List)
                # print(check_R2List)
                pass
6) Use Geneious to run _de novo_ assembly of output fasta sequences.  Set stringency high (75%).
7) filter out any sequences that did not fully encompass the exon, i.e.

          for s in seqs:
            if s == '?':
              return None
            else:
              l.append(s)
          return l    
8) Map the filtered exon extensions to the draft gDNA assembly
