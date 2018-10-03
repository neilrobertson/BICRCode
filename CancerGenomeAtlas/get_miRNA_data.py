
from os import listdir
from os.path import isfile, join

miRNA_directory = "/mnt/50tb/publicdata/TCGA/Colon_Adeno_Carcenoma_11_02_15/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3/"

outputFile = "/mnt/50tb/privatedata/Douglas/BRAF_DNMT3B_RNA-Seq/analysis/TCGA/miRNA/COAD-TCGA.miRNA-expression.HISeq.csv"

files = [f for f in listdir(miRNA_directory) if isfile(join(miRNA_directory, f))]
hg_files = [f for f in files if f.find(".hg19.mirna.quantification.txt") != -1]

print "{0} number of hg19 miRNA files in directory.".format(str(len(hg_files)))

sample_ids = []

matrix_map = {}

for filename in hg_files:
	id = "-".join(filename.strip().split(".")[0].strip().split("-")[0:4])[:-1]
	sample_ids.append(id)

	filename = miRNA_directory + "/" + filename
	with open(filename, "r") as file:

		for i, line in enumerate(file):
			if i == 0:
				pass
			else:
				line_parts = line.strip().split("\t")
				miRNA_name = line_parts[0]
				miRNA_reads = line_parts[1]
				miRNA_FPKM = line_parts[2]

				try:
					miRNA_values = matrix_map[miRNA_name]
					miRNA_values.append(miRNA_FPKM)
					matrix_map[miRNA_name] = miRNA_values
				except:
					matrix_map[miRNA_name] = [miRNA_FPKM]


print len(sample_ids)

with open(outputFile, "w") as output:

	header_line = "miRNA_name" + "\t" + "\t".join(sample_ids) + "\n"
	output.write(header_line)

	for key in matrix_map.keys():

		line = key + "\t" + "\t".join(matrix_map[key]) + "\n"
		output.write(line)

print "Complete"