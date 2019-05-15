class DNAstring:
    def __init__(self, s):
        self.s = s
    def fwdStrand(self):
        return self.s
    def revStrand(self):
        revDNA = self.s
        return revDNA[::-1]
    def compStrand(self):
        compDNA = self.s
        compDNA_string = ''
        for i in compDNA:
            if i == 'A':
                compDNA_string += 'T'
            elif i == 'T':
                compDNA_string += 'A'
            elif i == 'C':
                compDNA_string += 'G'
            elif i == 'G':
                compDNA_string += 'C'
            else:
                print('Error, unfamiliar base')
        return compDNA_string
    def revCompStrand(self):
        s_rev = self.s[::-1]
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
class ParseFasta:
    def __init__(self, f):
        self.f = f
    def parseFastaFile(self):
        fastaFile = self.f
        l = []
        with open(fastaFile, 'r') as fOpen:
            h = ''
            for i in fOpen:
                i = i.rstrip('\r\n')
                if i == '':
                    continue
                if i[0] == '>':
                    h = i
                else:
                    l.append((h,i))
        return l
    def parseFastaFileWithRC(self):
        fastaFile = self.f
        l = []
        with open(fastaFile, 'r') as fOpen:
            h = ''
            for i in fOpen:
                i = i.rstrip('\r\n')
                if i == '':
                    continue
                if i[0] == '>':
                    h = i
                else:
                    h_f = h + '_f'
                    h_rc = h + '_rc'
                    l.append((h_f,i))
                    i_rc = DNAstring(i).revCompStrand()
                    l.append((h_rc, i_rc))
        return l
    def parseFastaOutputFile(self):
        fastaFile = self.f
        l = []
        with open(fastaFile, 'r') as fOpen:
            for i in fOpen:
                i = i.rstrip('\r\n')
                l.append(i)
        return l
class OutputFasta:
    def __init__(self):
        self.l = []
    def createOutputList(self, c):
        self.l = c
    def outputFastaToFile(self, fName):
        if not self.l:
            print('Error, call createOutputList() function first!')
            return None
        for i in self.l:
            with open(fName, 'a') as fWrite:
                fWrite.write(i[0] + '\n')
                fWrite.write(i[1] + '\n')
class ParseOutputLine:
    def __init__(self):
        self.l = []
    def get_parseLine(self):
        return self.l
    def set_parseLine(self, b, s):
        if b:
            parsedLineForExon = []
            parse1 = s[2:] # remove '# ' from line
            leadingCharsRemoved = parse1.split(',')
            leadingCharRemovedFromAll = [x[1:] for x in leadingCharsRemoved] # list is now ['seq1:exon1', 'seq2:exon2',...]
            parsedLineAsLoL = [x.split(':')[:-1] for x in leadingCharRemovedFromAll]
            self.l = [':'.join(x) for x in parsedLineAsLoL] # output is ['seq1:01:01', 'seq2:01',...]
        else:
            parse1 = s[1:]
            splitParse = parse1.split(' ')
            self.l = parse1.split(' ')
class ParseOutputFile:
    def __init__(self):
        self.o = []
    def get_parseFile(self):
        return self.o
    def set_parseFile(self, f):
        returnedList = []
        parsedAlleleLine = ParseOutputLine()
        parsedFqHeaderLine = ParseOutputLine()
        if not f:
            return None
        for i in f:
            if i[0] == '#':
                parsedAlleleLine.set_parseLine(True, i)
            elif i[0] == '>':
                parsedFqHeaderLine.set_parseLine(False, i)
            else:
                returnedList.append((parsedAlleleLine, parsedFqHeaderLine, [i]))
                parsedFqHeaderLine = ParseOutputLine()
                parsedAlleleLine = ParseOutputLine()
        # return returnedList
        self.o = returnedList
