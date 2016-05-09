'''
Script for obtaining gene IDs from gene groups obtained in scoary output.

input:  gene_presence_absence.csv file from roary
		scoary output csv file for the genes of interest, otherwise it will perfom 
		analysis for all group genes in roary csv
output: txt file with geneId to be used in DAVID
'''

import csv
import os

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
		reader.next()
		for row in reader:
			geneGroup[row[0]]=[] #initialize dictionary containing empty list
	return geneGroup


def parsePAGeneFile(filename, dic):
	with open(filename, 'r') as csvfile:
		reader = csv.reader(csvfile)
		filenames=reader.next()[14:] #save filenames
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
	return dic, filenames

def cleanFilenames(filenames, dic):
	#not sure if necessary
	fnames=[]
	for key,value in dic.items():
		name =value.keys()#
		for item in name:
			fname=item.split('_')
			fname='_'.join(fname[:-1])
			#print fname
			if fname not in fnames:
				fnames.append(fname)
	#print fnames
	finalList=[]
	for lala in filenames:
		realName=item.split('_')
		realName='_'.join(realName[:-1])
		for name in fnames:
			if fname == realName and lala not in finalList:
				finalList.append(lala)
	#print len(finalList)

	return finalList

def removeGeneID(dic, listToRemove):

	#remove geneIDs, if groupID is empty, remove that too
	for k, value in dic.items():
		for key in value.keys():
			if key in listToRemove:
				#value.pop(key)
				del value[key]
	for k,v in dic.items():
		if len(v)==0:
			del dic[k]
	#print dic
	return dic

def parseGFFs(gffdir, filenames, dic):

	geneIDtoRemove=[]

	for item in filenames:
		with open(gffdir+item+'.gff', 'r') as gffFile:
			for line in gffFile:
				if line.startswith('gnl'):
					line=line.split('\t')
					line=line[-1].split(';')
					locusID=str(line[0].split('=')[1])
					if '_gene' in locusID:
						#print locusID
						pass
					else:
						try:
							line=line[2].split('|')
							#print line
							#print line
							gi= line[1]
							ref=line[3]

							for key, value in dic.items():
								if locusID in value.keys():
									value[locusID]=(gi, ref)

						except:
							#no id found! remove from dic!
							geneIDtoRemove.append(locusID)


	'''
	for key, value in dic.items():
		print "group ID: " + str(key)
		for k,v in value.items():
			print "Gene ID: " + str(k)
			print v'''
	print len(dic)
	newDic=removeGeneID(dic, geneIDtoRemove)
	print len(newDic)

	for key, value in newDic.items():
		print "group ID: " + str(key)
		for k,v in value.items():
			print "Gene ID: " + str(k)
			print v
	#return dic, geneIDtoRemove


def main():

	pa_file='/home/ines/Dropbox/Tese/roary/roary_ines_n61/gene_presence_absence.csv'

	scoary_file = '/home/ines/Dropbox/Tese/Scoary_1.1.2/Horse_05_05_2016_1144.csv exclusive_present.csv'
	
	gffdir='/home/ines/Dropbox/Tese/roary/gff_roary/'

	dic_allgenes=parseGeneCSVFile(pa_file)
	#print len(dic_allgenes)
	#print dic_allgenes
	dic_scoary_set=parseGeneCSVFile(scoary_file)
	#print len(dic_scoary_set)
	#print dic_scoary_set

	set_genes, filenames =parsePAGeneFile(pa_file,dic_scoary_set)

	all_genes, filenames=parsePAGeneFile(pa_file,dic_allgenes)
	'''
	for key, value in all_genes.items():
		print "chave: " + str(key)
		print len(value)'''

	#print len(all_genes)

	#print len(filenames)
	#TODO-parse GFF file

	#cleanFilenames(filenames,set_genes)

	parseGFFs(gffdir, filenames, set_genes)




if __name__ == "__main__":
    main()