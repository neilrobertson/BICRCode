''' File Containing Code for Douglas Thesis '''

from CancerGenomeAtlas.CancerGenomeAtlas import CancerGenomeAtlas
from CancerGenomeAtlas.DNA_Methylation import DNA_Methylation
from CancerGenomeAtlas.RNASeqV2 import RNASeqV2
from CancerGenomeAtlas.Somatic_Mutations import Somatic_Mutations
from CancerGenomeAtlas.CNV_Low_Pass_DNASeq import CNV_SNP_Array
from CancerGenomeAtlas.Expression_Genes import Expression_Genes
from CancerGenomeAtlas.Expression_Protein import Expression_Protein

from Analysis.Filter import Filter

#Filter.filterfiles_forSharedSampleIDs("/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_11_02_15/Somatic_Mutations/BCM__SOLiD_DNASeq/Level_2/hgsc.bcm.edu__ABI_SOLiD_DNA_Sequencing_level2.maf", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNMT3B-RNA_BRAFv600E/UNC__IlluminaGA_RNASeqV2.rsem.genes.results.output_matrix.DNMT3B.tsv", "/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/thesis/TCGA/DNMT3B-RNA_BRAFv600E/BCM__IlluminaGA_DNASeq.mutationSignatureBySample.BRAF.SharedIDs.tsv", delimiter = "\t")

homeDir = "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_11_02_15/"

tcga = CancerGenomeAtlas(homeDir)
#tcga.buildRNAMatrix()


for rnaSeqType in RNASeqV2.getRNASequencers():
    
    
    rnaSeqDir = homeDir + r"RNASeqV2/"+ rnaSeqType + r"/Level_3"
    print rnaSeqDir
    print CancerGenomeAtlas.checkDirectory(rnaSeqDir)
    
    if CancerGenomeAtlas.checkDirectory(rnaSeqDir) != None:
        
        for fileType in (r".rsem.genes.results", r".rsem.genes.normalized_results"):
            
            print fileType, rnaSeqDir
            
            fileName = rnaSeqType + fileType + ".output_matrix.tsv"

            outputDir = tcga.getOutputPath() + r"/" + fileName
            
            rnaSeqMatrix = RNASeqV2(rnaSeqDir, outputDir, fileType, tcga, delimiter = "\t")
            rnaSeqMatrix.buildMatrix()

            outputFile = rnaSeqMatrix.writeToOutput()

            rnaSeqMatrix.dispose()

            if fileType == ".rsem.genes.normalized_results":
                RNASeqV2.foldChange_RNASeqMatrix(outputFile, 1, 0.001)

