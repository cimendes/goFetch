'''
Script for obtaining gene IDs from gene groups obtained in scoary output.

input:  gene_presence_absence.csv file from roary
		scoary output csv file for the genes of interest, otherwise it will perfom 
		analysis for all group genes in roary csv
output: txt file with geneId to be used in DAVID
'''

import csv

def detectDelimiter(csvFile):
	#detect file delimiter, ',' or ';'
    with open(csvFile, 'r') as myCsvfile:
        header=myCsvfile.readline()
        if header.find(";")!=-1:
            return ";"
        if header.find(",")!=-1:
            return ","
    return ";"

def parseGeneCSVFile(filename):
	#Reads csv file with gene information, wither from Roary or Scoary. Gene group id must be in the first col.
	#Saves gene group ids as key in dictionary.
	geneGroup={}
	Delimiter=detectDelimiter(filename)

	with open(filename, 'r') as csvfile:
		reader = csv.reader(csvfile, delimiter=Delimiter)
		reader.next() #skip header line
		for row in reader:
			#print row[0]
			geneGroup[row[0]]=[] #initialize dictionary containing empty list
	return geneGroup

def parsePAGeneFile(filename, dic):
	with open(filename, 'r') as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			#print row[0]
			#print dic.keys()
			if row[0] in dic.keys():
				#geneID starts in the 15 camp
				#print row[14:]
				temp_dic={}
				for item in row[14:]:
					#check if item is empty, only saves non-empty ids
					if item == '':
						pass
					else:
						temp_dic[item]=[]
				dic[row[0]]=temp_dic
	return dic

def main():

	pa_file='/home/ines/Dropbox/Tese/roary/roary_ines_n61/gene_presence_absence.csv'

	scoary_file = '/home/ines/Dropbox/Tese/Scoary_1.1.2/Horse_05_05_2016_1144.csv exclusive_present.csv'
	dic_allgenes=parseGeneCSVFile(pa_file)
	#print len(dic_allgenes)
	#print dic_allgenes
	dic_scoary_set=parseGeneCSVFile(scoary_file)
	#print len(dic_scoary_set)
	#print dic_scoary_set

	set_genes=parsePAGeneFile(pa_file,dic_scoary_set)

	all_genes=parsePAGeneFile(pa_file,dic_allgenes)
	'''
	for key, value in all_genes.items():
		print "chave: " + str(key)
		print len(value)'''

	#print len(all_genes)


	#TODO-parse GFF file (where to save filenames?)


if __name__ == "__main__":
    main()