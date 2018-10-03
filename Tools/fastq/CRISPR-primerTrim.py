import opt, sys, re
from fastq.FastQ import FastQFile

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["inputFastq=", "outputFastq=" "CRISPR-primers="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
        
crisprPrimersFilename = None
inputFilename = None
outputFilename = None

for o, a in opts:
    if o=="--inputFastq":
        inputFilename = a
    elif o=="--outputFastq":
        outputFilename = a
    elif o=="--CRISPR-primers":
        crisprPrimersFilename = a
        
def getCRISPRPrimers(filename):
    """
    The primers file is a tab delimited file containing the PRIMER sequence in the first column and each
    individual primer on a new line.
    """
    primers = []
    with open(filename, "r") as input:
        lines = input.readlines()
        for line in lines:
            primers.append(line[0].strip())
    return primers

primers = getCRISPRPrimers(crisprPrimersFilename)

fastqFile = FastQFile(inputFilename)

with open(outputFilename, "w") as outputFile:
    
    for fastqBlock in fastqFile.readlines():
        try:
            assert fastqBlock.header[0] == r"@"
        except:
            print "Fastq file read failed.  Check file consistency!"
            break
        
        for primer in primers:
            
            fastqBlock.sequence
        