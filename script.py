'''
Script for obtaining gene IDs from gene groups obtained in scoary output.

input:  gene_presence_absence.csv file from roary
		scoary output csv file for the genes of interest, otherwise it will perfom 
		analysis for all group genes in roary csv
output: txt file with geneId to be used in DAVID
'''

def parseGeneCSVFile(filename):
	geneGroup={}
	with open(gene_filename, 'r') as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			geneGroup[row[0]]={}  #initialize dictionary of dictionaries
	return geneGroup





def main():

pa_file='/home/ines/Dropbox/Tese/roary/roary_ines_n61/gene_presence_absence.csv'

dic=parseGeneCSVFile(pa_file)
print dic

if __name__ == "__main__":
    main()