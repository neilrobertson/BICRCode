'''
Created on 21 Oct 2014

./gp_gistic2_from_seg -b $basedir -seg $segfile -mk $markersfile -refgene $refgenefile -alf $alf -cnv $cnvfile -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme

@author: neilrobertson
'''

class GISTIC(object):
    
    def __init__(self):
        pass
    
    @staticmethod
    def run_GISTIC(outputDir, segmentationFile, sampleListFile, genomeBuild = "hg19"):
        import subprocess
        from ThirdPartyTools import ThirdPartyTools
        
        GISTIC_DIR = "%s/GISTIC/GISTIC_2_0_22" % (ThirdPartyTools.getToolsRepositoryPath())
        
        markersFile = None
        if genomeBuild == "hg18": markersFile = "/mnt/50tb/publicdata/TCGA/TCGA_DataFiles/TCGA.DCC.GenomeWideSNP6.marker.na30.lst"
        else: markersFile = "/mnt/50tb/publicdata/TCGA/TCGA_DataFiles/TCGA.DCC.GenomeWideSNP6.marker.na31.lst"
        
        refGeneFile = None
        if genomeBuild == "hg18": refGeneFile = "%s/refgenefiles/hg18.mat" % (GISTIC_DIR)
        else: refGeneFile = "%s/refgenefiles/hg19.mat" % (GISTIC_DIR)
        
        command = [GISTIC_DIR + r"/run_gistic", outputDir, segmentationFile, markersFile, refGeneFile, sampleListFile]
        print " ".join(command)
        #try:
        subprocess.call(command)
        #    print " *** Running GISTIC Processes Internally ***"
        #    print "GISTIC PROCESS COMPLETED."
        #except:
        #    print "GISTIC Process Failed."