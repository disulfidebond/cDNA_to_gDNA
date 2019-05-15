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
