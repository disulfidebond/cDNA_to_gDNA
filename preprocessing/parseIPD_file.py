#!/usr/bin/python
import argparse


def main():
    parser = argparse.ArgumentParser(description='Parse IPD database text file.')
    parser.add_argument('-t', '--typeOfOutputFile', choices=['exonList', 'mergedFastaList', 'fastaExonList'], help='The format of the output file.')
    parser.add_argument('-o', '--outputFile', help='The name of the outputFile.')
    args = parser.parse_args()
    # output section
    wroteToFile = False
    remDuplicates_list = []
    if args.typeOfOutputFile == 'exonList':
        res = parsedExonOutput(alleleList) # FIX THIS! need to import class
        for i in res:
            print(i)
    else:
        if args.typeOfOutputFile == 'mergedFastaList':
            for i in alleleList:
                hString = '>' + str(i[0]) + ' ' + str(i[1])
                seqString = str(i[-1])
                print(hString + '\n' + seqString)
        elif args.typeOfOutputFile == 'fastaExonList':
            res = parsedExonOutput(alleleList)
            for i in res:
                iSplit = i.split(',')
                parseExonNames = iSplit[0].split(':')
                parsedExonName = parseExonNames[-1]
                hString = '>' + str(iSplit[0])
                seqString = str(iSplit[1])
                with open(args.outputFile, 'a') as fWrite:
                    fWrite.write(hString + '\n' + seqString + '\n')
            wroteToFile = True
        else:
            print('Warning!  Unrecognized outputType.  Please double check what you entered, and rerun this notebook from the beginning.')
if wroteToFile:
    print('output Finished!')
    print('The output is in file ' + str(args.outputFile))
