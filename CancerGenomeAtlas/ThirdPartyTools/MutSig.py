'''
Created on 21 Oct 2014

MutSig is an algorithm that allows you to access a .maf file of somatic mutations signatures for significance.  It requires the 2013a MatLAB MRC compiler only tool.

@ http://www.broadinstitute.org/cancer/cga/mutsig_run
@ http://www.mathworks.co.uk/products/compiler/mcr

/mnt/50tb/repository/3rdparty/MutSig/MutSigCV_1.4/run_MutSigCV.sh \
/usr/local/MATLAB/MATLAB_Compiler_Runtime/v81 \
/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Somatic_Mutations/BCM__IlluminaGA_DNASeq/Level_2/hgsc.bcm.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf \
/mnt/50tb/repository/3rdparty/MutSig/MutSigCV_1.4/PublicData/exome_full192.coverage.txt \
/mnt/50tb/repository/3rdparty/MutSig/MutSigCV_1.4/PublicData/gene.covariates.txt \
/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Somatic_Mutations/BCM__IlluminaGA_DNASeq/Level_2/Illumina_Genome_Analyzer_DNA_Sequencing_level2.MutSig \
/mnt/50tb/repository/3rdparty/MutSig/MutSigCV_1.4/PublicData/mutation_type_dictionary_file.txt \
/mnt/50tb/repository/3rdparty/MutSig/MutSigCV_1.4/PublicData/chr_files_hg19

@author: neilrobertson
'''


class MutSig(object):
    def __init__(self):
        pass
    
    @staticmethod
    def run_MutSig(mafFile, outputDir, mutSigDirectory = None, matLabMRCDirectory = None):
        import subprocess
        from ThirdPartyTools import ThirdPartyTools
        if mutSigDirectory == None: mutSigDirectory = ThirdPartyTools.getToolsRepositoryPath() + r"/MutSig/MutSigCV_1.4"
        print "MutSig Dir: %s"% (mutSigDirectory)
        if matLabMRCDirectory == None: matLabMRCDirectory= r"/usr/local/MATLAB/MATLAB_Compiler_Runtime/v81"
        
        command = [mutSigDirectory + r"/run_MutSigCV.sh", matLabMRCDirectory, mafFile, mutSigDirectory + r"/PublicData/exome_full192.coverage.txt", mutSigDirectory + r"/PublicData/gene.covariates.txt", outputDir , mutSigDirectory + r"/PublicData/mutation_type_dictionary_file.txt", mutSigDirectory + r"/PublicData/chr_files_hg19"]
        print " ".join(command)
        try:
            subprocess.call(command)
            print " *** Running MutSig Processes Internally ***"
            print "MutSig COMPLETED."
        except:
            print "MutSig Process Failed."
        
#MutSig.run_MutSig("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Somatic_Mutations/BCM__IlluminaGA_DNASeq/Level_2/hgsc.bcm.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf", "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_16_09_2014/Somatic_Mutations/BCM__IlluminaGA_DNASeq/Level_2/MutSig.Results")