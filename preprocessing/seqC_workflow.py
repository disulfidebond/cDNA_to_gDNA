#!/Users/jrcaskey/anaconda3/bin/python

import re
import hashlib
import ahocorasick
import argparse
import time

def parseOutExonSeq(parsedList, exonN):
    return list(filter(lambda x: x[0] == exonN, parsedList))

def revComp_nucleotides(s):
    s_rev = s[::-1]
    # s_rev = s
    s_rc = ''
    for i in s_rev:
        if i == 'A':
            s_rc += 'T'
        elif i == 'T':
            s_rc += 'A'
        elif i == 'C':
            s_rc += 'G'
        elif i == 'G':
            s_rc += 'C'
        else:
            print('Error, unfamiliar base')
    return s_rc

def parseFasta(f):
    l = []
    with open(f, 'r') as fOpen:
        h = ''
        for i in fOpen:
            i = i.rstrip('\r\n')
            if i[0] == '>':
                h = i
            else:
                l.append((h,i))
    return l

def parseFastaWithRC(f):
    l = []
    with open(f, 'r') as fOpen:
        h = ''
        for i in fOpen:
            i = i.rstrip('\r\n')
            if i[0] == '>':
                h = i
            else:
                h_f = h + '_f'
                h_rc = h + '_rc'
                l.append((h_f,i))
                i_rc = revComp_nucleotides(i)
                l.append((h_rc, i_rc))
    return l

def parseFastaIntoExons(l):
    parsedFastaList = []
    currentSeq = l[0][0].split(':')
    currentSeq = ':'.join(currentSeq[:-1])
    seqTuple = []
    for i in l:
        checkSeq = i[0].split(':')
        checkSeq = ':'.join(checkSeq[:-1])
        if checkSeq != currentSeq:
            parsedFastaList.append((currentSeq, seqTuple))
            currentSeq = checkSeq
            seqTuple = []
            h = i[0].split(':')[-1]
            seqTuple.append((h, i[1]))
        else:
            h = i[0].split(':')[-1]
            seqTuple.append((h, i[1]))
    parsedFastaList.append((currentSeq, seqTuple))
    return parsedFastaList

def writeExonToFastaFile(l, exonN):
    fastaName_exonN = 'Mamu_' + exonN + '.fasta.txt'
    for i in l:
        list2Parse = i[1]
        res = parseOutExonSeq(list2Parse, exonN)
        if res:
            # assumes unique exon identifiers within each sublist
            hString = str(i[0]) + ':' + res[0][0]
            hSeq = res[0][1]
            with open(fastaName_exonN, 'a') as fWrite:
                fWrite.write(hString + '\n' + hSeq + '\n')

def parseAndSplitSeqs(l):
    r = []
    for i in l:
        if len(i[1]) < 60:
            continue
        elif len(i[1]) > 59 and len(i[1]) < 105:
            h = i[0] + '_1'
            parsedSeq_1 = i[1]
            r.append((h, parsedSeq_1))
            h = i[0] + '_1rc'
            parsedSeq_1rc = revComp_nucleotides(parsedSeq_1)
            r.append((h, parsedSeq_1rc))
        else:
            unparsedSeq = i[1]
            h = i[0] + '_1'
            parsedSeq_1 = unparsedSeq[0:106]
            r.append((h, parsedSeq_1))
            h = i[0] + '_1rc'
            parsedSeq_1rc = revComp_nucleotides(parsedSeq_1)
            r.append((h, parsedSeq_1rc))
            h = i[0] + '_2'
            parsedSeq_2 = unparsedSeq[-106:]
            r.append((h, parsedSeq_2))
            h = i[0] + '_2rc'
            parsedSeq_2rc = revComp_nucleotides(parsedSeq_2)
            r.append((h, parsedSeq_2rc))
    return r

# merge/mark duplicated sequences
def createParsedDictOfSeqs(l):
    parsedDict = dict()
    for x,y in l:
        if y not in parsedDict:
            parsedDict[y] = [x]
        else:
            tmp = parsedDict[y]
            tmp.append(x)
            parsedDict[y] = tmp
    return parsedDict

# convert to list of tuples to preserve indexing
def createListOfTuples(d):
    lOfTuples = []
    for k,v in d.items():
        t1 = v
        t2 = k
        unhashedString = t1[0]
        hashedString = hashlib.md5(unhashedString.encode('utf-8')).hexdigest()
        lOfTuples.append((hashedString, (t1, t2)))
    return lOfTuples

def matchSeqsToReads(l, seqL):
    mamu_exonN = ahocorasick.Automaton()
    for idx, val in enumerate(l):
        mamu_exonN.add_word(val[1][1], (idx, val[1][1]))
    mamu_exonN.make_automaton()
    res_matches = []
    for s in seqL:
        queryLen = len(s[1]) - 1
        for end_idx, (insert_order, origVal) in mamu_exonN.iter(s[1]):
            start_idx = end_idx - len(origVal) + 1
            if start_idx == 0 or (queryLen - end_idx) == 0: # anchored at either 5' or 3'
                # print((s[0], s[1], start_idx, end_idx, (insert_order, origVal)))
                matchedAlleles = l[insert_order][1][0]
                parsedFastaReadLine = s[0] + '\n' + s[1]
                res_matches.append((parsedFastaReadLine, matchedAlleles))
    return res_matches

def parseSequences(seqs2p, ofp, rp, ifp, tupleListAll):
    for i in range(0, 5):
        exonNameInt = i + 1
        exonNameString = str(exonNameInt)
        outputFile = ofp + 'Mamu-exon' + exonNameString + '_matches_r' + rp + '_' + ifp + '.txt'
        res_matches = matchSeqsToReads(tupleListAll[i], seqs2p)
        for itm in res_matches:
            with open(outputFile, 'a') as fWrite:
                fWrite.write('# ' + ','.join(itm[1]) + '\n')
                fWrite.write(itm[0] + '\n')
        print('Exon ' + str(exonNameString) + ' completed!')

def main():
    parser = argparse.ArgumentParser(description='Map anchored sequences to IPD cDNA sequences.')
    parser.add_argument('-f', '--filePrefix', help='File prefix for sample sequence fasta file. Example: 39771')
    parser.add_argument('-t', '--trimPrefix', choices=['1', '2'], help='Integer identifying if the trimmed sample is 1 or 2')
    parser.add_argument('-r', '--readNum', choices=['1', '2'], help='Integer indicating whether the sample fasta file is read 1 or read 2')
    parser.add_argument('-i', '--inputFasta', help='File name for the input multifasta file from IPD')
    args = parser.parse_args()
    # ### Required Variables
    outputFilePrefix = args.filePrefix
    outputFilePrefix += '_'
    inputFilePrefix = args.trimPrefix
    fastaName = args.inputFasta
    readPairPrefix = []

    if not args.readNum:
        readPairPrefix = ['1', '2']
    else:
        readPairPrefix = args.trimPrefix
    # for i in range(1,6):
    #    fastaListAll = parseFasta(fastaName)
    #    parsedFastaList = parseFastaIntoExons(fastaListAll)
    #    exonN_str = 'exon' + str(i)
    #    writeExonToFastaFile(parsedFastaList, exonN_str)
    # print('parsed fasta files')
    # time.sleep(5)

    exonL_1 = parseFasta('Mamu_exon1.fasta.txt')
    exonL_2 = parseFasta('Mamu_exon2.fasta.txt')
    exonL_3 = parseFasta('Mamu_exon3.fasta.txt')
    exonL_4 = parseFasta('Mamu_exon4.fasta.txt')
    exonL_5 = parseFasta('Mamu_exon5.fasta.txt')

    exonL_all = [exonL_1, exonL_2, exonL_3, exonL_4, exonL_5]
    exon_tupleListAll = []
    for i in range(0, 5):
        exonN_parsedList = parseAndSplitSeqs(exonL_all[i])
        exonN_parsedDict = createParsedDictOfSeqs(exonN_parsedList)
        exonN_tupleList = createListOfTuples(exonN_parsedDict)
        exon_tupleListAll.append(exonN_tupleList)

    if type(readPairPrefix) is list:
        for i in readPairPrefix:
            fastaFileString = outputFilePrefix + 'R' + i + '_' + '001.qtrim.' + inputFilePrefix + '.fasta'
            seqs_to_parse = parseFasta(fastaFileString)
            parseSequences(seqs_to_parse, outputFilePrefix, readPairPrefix, inputFilePrefix, exon_tupleListAll)
    else:
        fastaFileString = outputFilePrefix + 'R' + readPairPrefix + '_' + '001.qtrim.' + inputFilePrefix + '.fasta'
        seqs_to_parse = parseFasta(fastaFileString)
        parseSequences(seqs_to_parse, outputFilePrefix, readPairPrefix, inputFilePrefix, exon_tupleListAll)
if __name__ == '__main__':
    main()
