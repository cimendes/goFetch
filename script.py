'''
Script for obtaining gene IDs from gene groups obtained in scoary output.

input:  gene_presence_absence.csv file from roary
		scoary output csv file for the genes of interest, otherwise it will perfom 
		analysis for all group genes in roary csv
output: txt file with geneId to be used in DAVID
'''

import csv

def parseGeneCSVFile(filename):
	geneGroup={}
	with open(filename, 'r') as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			geneGroup[row[0]]={}  #initialize dictionary of dictionaries
	return geneGroup




def main():

	pa_file='/home/ines/Dropbox/Tese/roary/roary_ines_n61/gene_presence_absence.csv'

	scoary_file = '/Users/cimendes/Dropbox/Tese/Scoary_1.1.2/Horse_05_05_2016_1144.csv exclusive_present.csv'
	dic_allgenes=parseGeneCSVFile(pa_file)
	print len(dic)
	dic_scoary_humans=parseGeneCSVFile(pa_file)
	print dic_scoary_humans

	

if __name__ == "__main__":
    main()