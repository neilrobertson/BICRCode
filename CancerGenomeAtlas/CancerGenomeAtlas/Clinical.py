'''
Created on 29 Aug 2014

@author: neilrobertson
'''

class Clinical(object):
    def __init__(self):
        pass
    
    
    
    @staticmethod
    def obtainSpecificClinicalInformation_fromSampleID(inputFilename, outputFilename, headerLineCount = int(3), delimiter = "\t"):
        assert int(headerLineCount)
        assert inputFilename
        
        from DataFiles.DataFiles import DataFiles
        patientIds = DataFiles.getClinicalPatientIDs()
        
        patientDetails = ""
        with open(inputFilename, "r") as inputFile:
            print "Starting to gather clinical patient details for %s patients..." % (str(len(patientIds)))
            headerLines = ""
            for i, row in enumerate(inputFile):
                if i < headerLineCount: #headerLines
                    headerLines += row 
                else:
                    rowID = row.split(delimiter)[0]
                    if rowID in patientIds:
                        print "Patient found: %s" % rowID
                        patientDetails += row

        with open(outputFilename, "w") as outputFile:
            outputFile.write(headerLines)
            outputFile.write(patientDetails)
            print "Output file completed."
        return outputFilename


Clinical.obtainSpecificClinicalInformation_fromSampleID("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Clinical/Biotab/nationwidechildrens.org_clinical_patient_coad.txt", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Outputs/Clinical/TCGA_ColonAdenocarcinoma_ClinicalData_Non-CIMP_Details.tsv")