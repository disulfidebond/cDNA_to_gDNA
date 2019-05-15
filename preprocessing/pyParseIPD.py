#!/user/bin/python

class IPDfile:
    def __init__(self):
        self.l = []
class ImportIPD:
    def __init__(self):
        self.l = []
        self.f = ''
    def setFileName(self, fname):
        self.f = fname
    def getFileName(self):
        return self.f
    def setFileData(self):
        # requires file with entries delimited by '//', ended by line with '//'
        with open(self.f, 'r') as fOpen:
            currentItem = []
            for i in fOpen:
                i = i.rstrip('\r\n')
                splitItem = i.split(' ')
                pItem = list(filter(lambda x: x != '', splitItem))
                pItemString = ' '.join(pItem)
                if i[0:3] == '//':
                    self.l.append(currentItem)
                    currentItem = []
                else:
                    currentItem.append(pItemString)
    def getFileData(self):
        return self.l
class ParseIPDexonLine:
    def __init__(self):
        self.s = []
    def getParsedLine(self):
        return self.s
    def setParsedLine(self, l):
        for i in l:
            listOfValues = i
            exonSeq = i[-1]
            exon_ct = []
            exon_ct_list = []
            mxLen = len(listOfValues[3])
            for j in listOfValues[3]:
                exon_ct.append(j[0])
                exon_ct_list.append((j[1][0], j[1][1]))
            indexTransition = 0
            checkForIntron = 0
            for x in range(0, mxLen, 1):
                exonStartAndStop = exon_ct_list[x]
                exonStart = exonStartAndStop[0]
                if x == 0:
                    exonStart = exonStartAndStop[0] - 1
                exonStop = exonStartAndStop[1] + 1
                exonSequence = exonSeq[exonStart:exonStop]
                indexTransition = exonStartAndStop[1]+1
                r = str(listOfValues[0]) + ':' + str(exon_ct[x]) + ',' + exonSequence
                self.s.append(r)

class ParseIPD:
    def __init__(self):
        self.l = []
    def getIPDList(self, lName):
        return self.l
    def setIPDList(self, iFile):
        s_id = ''
        s_alleleID = ''
        s_exonList = []
        alleleExonRange = []
        alleleExonSeq = []
        exonSeqs = ''
        exonCt = 1
        for i in iFile:
            seqFlag = False
            for x in range(0, len(i)):
                e_entry = i[x]
                sLine = e_entry.split(' ')
                sLineList = list(filter(lambda x: x != '', sLine)) # remove?
                if seqFlag:
                    s = ''.join(sLineList[:-1])
                    exonSeqs += s
                    continue
                if sLineList[0] == 'KW':
                    s_id = sLineList[1]
                elif sLineList[0] == 'ID':
                    s_alleleID = sLineList[1]
                elif sLineList[0] == 'FT':
                    if sLineList[1] == 'allele':
                        alleleExonRange = parsedAndCastRange(sLineList[2])
                    elif sLineList[1] == 'exon':
                        exonSeqRange = parsedAndCastRange(sLineList[2])
                        # exonN = 'exon' + str(exonCt)
                        exonNum = i[x+1] # will never be greater than len(i)
                        exonNum_list = exonNum.split('=')
                        exonN = -1
                        try:
                            exonN = int(exonNum_list[1])
                        except ValueError:
                            print('Warning! Unable to parse exon number.  Check value:')
                            print(e_entry)
                            exonN = exonCt
                        exonIntAsString = 'exon' + str(exonN)
                        alleleExonSeq.append((exonIntAsString, exonSeqRange))
                        exonCt += 1
                    else:
                        continue
                elif sLineList[0] == 'SQ':
                    seqFlag = True
                else:
                    continue
            self.l.append([s_id, s_alleleID, alleleExonRange, alleleExonSeq, exonSeqs])
            alleleExonSeq = []
            exonCt = 1
            exonSeqs = ''
