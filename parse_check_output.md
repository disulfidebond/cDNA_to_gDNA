* The next step is to verify the output, which should be 
fastq_header [allele1, allele2,...,alleleN]

* Preliminary work verified that there should not be alleles from different types grouped together, i.e. a Mamu-B should not be in the same groups as a Mamu-A
* Preliminary work also verified that the forward and reverse complement of the allele sequences were not grouped together
* Example output:

        A00433:14:H3W3NDSXX:1:1112:17951:18427 ['Mamu-A1*008:01:01', 'Mamu-A1*027:01', 'Mamu-A1*008:01:04', 'Mamu-A1*006:03', 'Mamu-A1*090:01', 'Mamu-A1*008:03', 'Mamu-A1*003:12', 'Mamu-A1*008:02', 'Mamu-A1*006:01', 'Mamu-A1*006:02', 'Mamu-A1*003:13', 'Mamu-A1*003:04', 'Mamu-A1*065:02']

* Code:
        
        def quickCheck(d):
            ct = 0
            if type(d) is dict:
                for k,v in d.items():
                    ct += 1
                    if ct < 10:
                        print(k,v)
                    else:
                        break
            else:
                try:
                    for i in d:
                        ct += 1
                        if ct < 10:
                            print(i)
                        else:
                            break
                except TypeError:
                    print('Not a list or set, cannot show preview')

        parsedFile = ParseFasta('39526_Mamu-exonAll_allR1.txt').parseFastaOutputFile()
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

        parsedFile = ParseFasta('39526_Mamu-exonAll_allR2.txt').parseFastaOutputFile()
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
        quickCheck(parsedDictOfReads_R2)

* Then, the lists of R1 and R2 reads are compared to determine if any reads did not match to their pair, and if any did. In the case of the former, these are collected for the next step. In the case of the latter, these reads are discarded, and/or flagged for review

* Code

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
                print('Warning!  Check:')
                print(check_R1List)
                print(check_R2List)

### Comments
Note that this workflow is dependent on sequential stages of assembly. If only one pair of the reads matched, this still could be a legitimate mapping; for example, the unmatched pair could have sequenced a non-coding region, and therefore would not report a mapping to a known cDNA region. These reads are not completely discarded, but are instead set aside until reads with a known basis for comparison are used for an initial assembly. After this initial draft assembly, reads with only a single pair that successfully mapped will be used for an expanded assembly, because the draft assembly will provide anchors for mapping.
