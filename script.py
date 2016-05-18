'''
Script for obtaining gene IDs from gene groups obtained in scoary output.

input:  gene_presence_absence.csv file from roary
		scoary output csv file for the genes of interest, otherwise it will perfom 
		analysis for all group genes in roary csv
output: txt file with geneId to be used in DAVID

requires: bioservices, internet connection

For each genegroup it only saves unique ids

Dictionary Structure:
	{GeneGroup:
		-{GeneID:
			-GI Number
			-Ref Number
			-{Uniprot ID:
				-{Process:
					-[(GOid, GOName)]}
				-{Component:
					-[(GOid, GOName)]}
				-{Function:
					-[(GOid, GOName)]}
			}
		}
	}
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

def parseInputCSVFile(filename):
	#Reads csv file with gene information, wither from Roary or Scoary. Gene group id must be in the first col.
	#Saves gene group ids as key in dictionary.
	geneGroup={}
	Delimiter=detectDelimiter(filename) #detects file delimiter
	with open(filename, 'r') as csvfile:
		reader = csv.reader(csvfile, delimiter=Delimiter)
		reader.next()
		for row in reader:
			geneGroupID=row[0]
			annotation=row[1]
			otherAnnotation=row[2]
			geneGroup[geneGroupID]=[annotation,otherAnnotation] #initialize dictionary containing annotations
	return geneGroup


def parsePAGeneFile(filename, dic):
	#parser for Roary's gene_presence_absence.csv file
	with open(filename, 'r') as csvfile:
		reader = csv.reader(csvfile)
		filenames=reader.next()[14:] #save all filenames
		for row in reader:
			if row[0] in dic.keys():
				#geneID starts in the 15th col
				temp_dic={}
				for item in row[14:]:
					#check if item is empty, only saves non-empty ids
					if item == '':
						pass
					else:
						temp_dic[item]=[]
				dic[row[0]].append(temp_dic)
	return dic, filenames

def cleanFilenames(filenames, dic):
	#TODO: not sure if necessary
	#Cleans not needed filenames
	fnames=[]
	for key,value in dic.items():
		name =value.keys()
		for item in name:
			fname=item.split('_')
			fname='_'.join(fname[:-1])
			if fname not in fnames:
				fnames.append(fname)
	finalList=[]
	for lala in filenames:
		realName=item.split('_')
		realName='_'.join(realName[:-1])
		for name in fnames:
			if fname == realName and lala not in finalList:
				finalList.append(lala)
	return finalList

def removeGeneID(dic, listToRemove):
	#TODO! Verify is it's ok!
	#remove geneIDs, if groupID is empty, removes that too
	for k, value in dic.items():
		#print value
		for key in value[2].keys():
			if key in listToRemove:
				del value[2][key]
	for k,v in dic.items():
		if len(v[2])==0:
			del dic[k]
	#only leave unique IDs
	for k,v in dic.items():
		#print v
		first_key=v[2].keys()[0] #save first geneID in genegroup
		#print first_key
		previous_value=v[2][first_key] 
		#print values
		for key,value in v[2].items():
			if key != first_key:
				if value == previous_value:
					del v[2][key]
	'''
	for key,value in dic.items():
		print "Group ID: " + str(key)
		for k,v in value.items():
			print "Gene ID: " + str(k)
			print " Gi, ref: " + str(v)'''

	return dic

def parseGFFs(gffdir, filenames, dic):
	#parse gff files. Requires directory containing the gff files
	geneIDtoRemove=[]
	for item in filenames:
		with open(gffdir+item+'.gff', 'r') as gffFile:
			for line in gffFile:
				if line.startswith('gnl'):
					line=line.split('\t')
					line=line[-1].split(';')
					locusID=str(line[0].split('=')[1])
					if '_gene' in locusID:
						pass
					else:
						try:
							line=line[2].split('|')
							gi= line[1]
							ref=line[3]
							for key, value in dic.items():
								if locusID in value[2].keys():
									value[2][locusID]=[gi, ref]
						except:
							#no id found! remove from dic!
							geneIDtoRemove.append(locusID)
	
	newDic=removeGeneID(dic, geneIDtoRemove) #cleanup!
	
	return newDic


def printFile(dic):
	#TODO
	#printing 3 files: . gi numbers, ref with version and ref w/o version

	#printing GI:
	toPrint_Gi=[]
	for key, value in dic.items():
		for k,v in value.items():
			if v[0] not in toPrint_Gi:
				toPrint_Gi.append(v[0])
	with open("gi_ids.txt", "w") as fileGI:
		for item in toPrint_Gi:
			fileGI.write(str(item) + '\n')

	#printing Ref with version
	toPrint_Ref=[]
	for key, value in dic.items():
		for k,v in value.items():
			if v[1] not in toPrint_Ref:
				toPrint_Ref.append(v[1])
	with open("ref_wVersion.txt", "w") as fileRef_V:
		for item in toPrint_Ref:
			fileRef_V.write(str(item)+ '\n')

	#printing Ref without version:
	toPrint_Ref_noVersion=[]
	for key, value in dic.items():
		for k,v in value.items():
			toprint=v[1].split('.')
			toprint=toprint[0]
			if toprint not in toPrint_Ref_noVersion:
				toPrint_Ref_noVersion.append(toprint)
	with open("ref_noVersion.txt", "w") as fileRef:
		for item in toPrint_Ref_noVersion:
			fileRef.write(str(item) + '\n')

def retrieveUniprot(dic):
	#TODO: improve!
	#source: http://www.uniprot.org/help/programmatic_access
	import urllib,urllib2

	toPrint_Gi=[]
	for key, value in dic.items():
		for k,v in value[2].items():
			lala=v[1].split('.')[0]
			#print lala
			if lala not in toPrint_Gi: #v[0]
				toPrint_Gi.append(v[0]) #0 - GI

	url = 'http://www.uniprot.org/mapping/'

	#retrieve uniparc IDs from GI
	params = {
	'from':'P_GI', #P_GI P_REFSEQ_AC
	'to':'ACC',
	'format':'tab',
	'query':'	'.join(toPrint_Gi)
	}

	data = urllib.urlencode(params)
	request = urllib2.Request(url, data)
	contact = "" # Please set your email address here to help us debug in case of problems.
	request.add_header('User-Agent', 'Python %s' % contact)
	response = urllib2.urlopen(request)
	page = response.read(200000)

	gi_to_uniparc={}
	toPrint_uniparc=[]
	line_uniparc=page.split('\n')[1:-1]
	for item in line_uniparc:
		item=item.split('\t')
		gi_to_uniparc[item[0]]=item[1]
		toPrint_uniparc.append(item[1])

	#Retrieve uniprot IDs from uniparc
	params2 = {
	'from':'UPARC',
	'to':'ACC',
	'format':'tab',
	'query':'	'.join(toPrint_uniparc)
	}

	data2 = urllib.urlencode(params2)
	request2 = urllib2.Request(url, data2)
	contact2 = "" # Please set your email address here to help us debug in case of problems.
	request2.add_header('User-Agent', 'Python %s' % contact2)
	response2 = urllib2.urlopen(request2)
	page2 = response2.read(200000)
	
	uniparc_to_uniprot={}
	line_uniprot=page2.split('\n')[1:-1]
	for item in line_uniprot:
		item=item.split('\t')
		uniparc_to_uniprot[item[0]]=item[1]

	for key, value in dic.items(): #genegroup
		for k,v in value[2].items(): #geneID
			gi_ID=v[0]
			try:
				uniprot_dic={}
				uniparc_ID=gi_to_uniparc[gi_ID]
				uniprot_ID=uniparc_to_uniprot[uniparc_ID]
				uniprot_dic[uniprot_ID]=[]
				v.append(uniprot_dic)
			except:
				print "No UniParc or UniProt ID found for gene %s" % (k)

	return dic

def getGOnumbers(dic):
	#from the uniprot id list, recovers GOids. returns dict with the list of
	#GO ids for each uniprot id. 

	from bioservices import QuickGO
	import unicodedata

	s = QuickGO() #init QuickGO

	for geneGroup,value in dic.items():
		for geneID, values in value[2].items():
			for uniprotID, v in values[2].items():
				#check frmt as dict - saves as unicode
				result=s.Annotation(protein=str(uniprotID),source="UniProt",frmt='tsv',col="goID,aspect,goName")
				result= result.encode('ascii','ignore').strip() #convert unicode string to regular str
				#print result
				goTerm=result.split('\n')[2:]
				if len(goTerm)>0:
					listOfTerms={}
					for terms in goTerm:
						terms=terms.split('\t')
						if terms[1] not in listOfTerms:
							listOfTerms[terms[1]]=[(terms[0], terms[2])]
						else:
							values=listOfTerms[terms[1]]
							values.append((terms[0], terms[2]))
							listOfTerms[terms[1]]=values
				v.append(listOfTerms)
	return dic

def printReport(dic):
	#TODO - print final report

	with open("report.tsv",'w') as outfile:
		outfile.write("Gene Group\tNon-unique gene name\tAnnotation\tGene ID\tGI Number\tRef Number\tUniProt ID\t" + \
		"Cellular Component\tBiological Process\tMolecular Function\n")

		for GeneGroupID, general in dic.items():
			geneGroupID=GeneGroupID
			annotation = general[0]
			otherAnnotation = general[1]
			for GeneID, geneInfo in general[2].items():
				geneID=GeneID
				giNumber=geneInfo[0]
				refNumber=geneInfo[1]
				for UniprotID, geneOntology in geneInfo[2].items():
					uniprotID=UniprotID
					component=''
					process=''
					function=''
					for domain in geneOntology:
						#print domain
						for domainName, term in domain.items():
							#print domainName
							#print term
							if domainName == 'Component':
								for touple in term:
									component+=str(touple[0])+'='+str(touple[1])
									if touple != term[-1]:
										component+=';'
							elif domainName == 'Process':
								for touple in term:
									process+=str(touple[0])+'='+str(touple[1])
									if touple != term[-1]:
										process+=';'
							elif domainName == 'Function':
								for touple in term:
									function+=str(touple[0])+'='+str(touple[1])
									if touple != term[-1]:
										function+=';'


					#toPrint=str(domainName)+':'+'"'+str(touple[0])+','+str(touple[1])+'"'
							#print toPrint
							#GOID=term[0]
							#GOName=term[1]
					outfile.write(geneGroupID+'\t'+annotation+'\t'+otherAnnotation+'\t' + \
					geneID+'\t'+giNumber+'\t'+refNumber+'\t'+uniprotID+'\t'+component+'\t'+process+'\t'+function+'\n')
					
					#print uniprotID
					#print geneOntology
					'''
					for domain, terms in geneOntology:

						outfile.write(geneGroupID+'\t'+annotation+'\t'+otherAnnotation+'\t' + \
						geneID+'\t'+giNumber+'\t'+refNumber+'\t'+uniprotID+'\t'+domain+'\t'+terms)
					'''
		

def main():

	#TODO: Add argparse!

	import argparse

	VERSION=0.9

	parser = argparse.ArgumentParser(description='Gene Ontology fetcher for Roary and Scoary outputs.', epilog='by Catarina Mendes (cimendes@medicina.ulisboa.pt)')
	parser.add_argument('-g', '--genes', help='Input gene presence/absence table (comma-separated-values) from Roary (https:/sanger-pathogens.github.io/Roary)')
	parser.add_argument('-i', '--input', help='Input interest gene presence/avsence table (comma or semicolon-separated-values from Roary or Scoary')
	parser.add_argument('-d', '--dir', help='Path to directory containing all gff files used in the Roary analysis.')
	#parser.add_argument('--delimiter', help='The delimiter between cells in the gene presence/absence and trait files. NOTE: Even though commas are the default they might mess with the annotation column, and it is therefore recommended to save your files using semicolon or tab ("\t") instead. SCOARY will output files delimited by semicolon', default=',', type=str)
	#parser.add_argument('--version', help='Display version, and exit.', default=False,action='store_true')

	args = parser.parse_args()



	#input files
	pa_file='/home/ines/Dropbox/Tese/roary/roary_ines_n61/gene_presence_absence.csv'

	scoary_file = '/home/ines/Dropbox/Tese/Scoary_1.1.2/Horse_05_05_2016_1144.csv exclusive_present.csv'
	
	gffdir='/home/ines/Dropbox/Tese/roary/gff_roary/'

	#parsing 
	#dic_allgenes=parseInputCSVFile(pa_file)
	#print len(dic_allgenes)
	#print dic_allgenes
	dic_scoary_set=parseInputCSVFile(scoary_file)
	#print len(dic_scoary_set)
	#print dic_scoary_set

	set_genes, filenames =parsePAGeneFile(pa_file,dic_scoary_set)

	#all_genes, filenames=parsePAGeneFile(pa_file,dic_allgenes)
	'''
	for key, value in all_genes.items():
		print "chave: " + str(key)
		print len(value)'''

	#print len(all_genes)

	#print len(filenames)
	#TODO-parse GFF file

	#cleanFilenames(filenames,set_genes)

	#parse GFFs
	set_genes=parseGFFs(gffdir, filenames, set_genes)

	#printing files
	#printFile(set_genes) TODO!

	#retrieve uniprotKB ids
	set_genes_uniprot=retrieveUniprot(set_genes)

	'''
	for key, value in set_genes_uniprot.items():
		print "group ID: " + str(key)
		print "Annotation: " + str(value[0])
		print "Other Annotation: " + str(value[1])
		for k,v in value[2].items():
			print "Gene ID: " + str(k)
			print v'''

	#retrieve GO terms
	last_dic=getGOnumbers(set_genes_uniprot)
	
	for key, value in set_genes_uniprot.items():
		print "group ID: " + str(key)
		print "Annotation: " + str(value[0])
		print "Other Annotation: " + str(value[1])
		for k,v in value[2].items():
			print "Gene ID: " + str(k)
			print 'GI: ' + str(v[0])
			print 'Ref: ' + str(v[1])
			#print 'Uniprot ID: ' + str(v[2])
			for uniprotID, GOterms in v[2].items():
				print 'Uniprot ID: ' + str(uniprotID)
				print GOterms

	printReport(last_dic)


if __name__ == "__main__":
    main()