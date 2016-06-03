'''
Gene Ontology fetcher for Roary and Scoary outputs
Idea and implementation: Catarina Mendes (cimendes@medicina.ulisboa.pt)

input:  gene_presence_absence.csv file from roary
		scoary output csv file for the genes of interest or gene_presence_absence.csv file from roary
		Path to directory containing all gff files used in the Roary analysis
output: tsv file containing gene information, GI number, RefSeq Protein Number, Uniprot ID and Gene
		ontology terms found for Cellular Component, Biological Process and Molecular Function
		For each genegroup it only saves unique IDs

requires: bioservices, internet connection


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

import csv, os, argparse, sys, unicodedata, time
from bioservices import QuickGO

# Redirect stdout
# Modified from:
# http://stackoverflow.com/questions/14906764/how-to-redirect-stdout-to-both-file-and-console-with-scripting
class Logger(object):
	def __init__(self, out_directory):
		self.logfile = os.path.join(out_directory, "run.log")
		if os.path.isfile(self.logfile):
			print "Logfile already exists! It will be overwritten..." + "\n"
		self.terminal = sys.stdout
		self.log = open(self.logfile, "w")
	def write(self, message):
		self.terminal.write(message)
		self.log.write(message)
		self.log.flush()
	def flush(self):
		pass


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
						continue
					elif '\t' in item: #in case of paralogs
						item=item.split('\t')
						for geneID in item:
							temp_dic[geneID]=[]
					else:
						temp_dic[item]=[]
				dic[row[0]].append(temp_dic)
	return dic, filenames

def cleanFilenames(filenames, dic):
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
	for pos in filenames:
		realName=pos.split('_')
		realName='_'.join(realName[:-1])
		for name in fnames:
			if fname == realName and pos not in finalList:
				finalList.append(pos)
	return finalList

def removeGeneID(dic, listToRemove):
	#remove geneIDs, if groupID is empty, removes that too
	for k, value in dic.items():
		for key in value[2].keys():
			if key in listToRemove:
				del value[2][key]
	for k,v in dic.items():
		if len(v[2])==0:
			print "No IDs found. Group ID removed: " + str(k) 
			del dic[k]
	#only leave unique IDs
	for k,v in dic.items():
		first_key=v[2].keys()[0] #save first geneID in genegroup
		previous_value=[v[2][first_key]]
		for key,value in v[2].items():
			if key != first_key:
				if value in previous_value:
					del v[2][key]
				else:
					previous_value.append(value)

	return dic

def parseGFFs(gffdir, filenames, dic):
	#parse gff files. Requires directory containing the gff files
	geneIDtoRemove=[]
	for item in filenames:
		with open(gffdir+item+'.gff', 'r') as gffFile:
			for line in gffFile:
				if not line.startswith('##'):
					if line.startswith('gnl') or line.startswith('NZ') or line.startswith('gi'):
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
								#print 'No Id found! Removing.. ' +str(locusID)

	newDic=removeGeneID(dic, geneIDtoRemove) #cleanup!
	
	return newDic


def printFile(dic):
	#Not implemented
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
	#source: http://www.uniprot.org/help/programmatic_access
	import urllib,urllib2

	badID_genes=[]
	toPrint_Gi=[]
	for key, value in dic.items():
		for k,v in value[2].items():
			GI_ID=v[1].split('.')[0]
			if GI_ID not in toPrint_Gi: #v[0]
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
	request2.add_header('User-Agent', 'Python %s' % contact)
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
				print "No UniParc or UniProt ID found for gene %s. Being removed from analysis." % (k)
				badID_genes.append(k)
				del value[2][k] #removing
				if len(value[2])==0:
					print "No IDs found. Removed Gene Group ID: %s" % (key) 

	return dic,badID_genes

def getGOnumbers(dic):
	#from the uniprot id list, recovers GOids. returns dict with the list of
	#GO ids for each uniprot id. 

	s = QuickGO() #init QuickGO
	noGO=[]

	for geneGroup,value in dic.items():
		for geneID, values in value[2].items():
			for uniprotID, v in values[2].items():
				#check frmt as dict - saves as unicode
				result=s.Annotation(protein=str(uniprotID),source="UniProt",frmt='tsv',col="goID,aspect,goName")
				result= result.encode('ascii','ignore').strip() #convert unicode string to regular string
				goTerm=result.split('\n')[2:]
				listOfTerms={}
				if len(goTerm)>0:
					for terms in goTerm:
						terms=terms.split('\t')
						if terms[1] not in listOfTerms:
							listOfTerms[terms[1]]=[(terms[0], terms[2])]
						else:
							values=listOfTerms[terms[1]]
							values.append((terms[0], terms[2]))
							listOfTerms[terms[1]]=values
				else:
					print "No Gene Ontology terms found for %s" % geneID
					noGO.append(geneID)
				v.append(listOfTerms)
	return dic,noGO

def badResultsReport(badID,badGO):

	with open("failReport"+time.strftime("_%d_%m_%Y_%H%M")+".txt",'w') as failReport:
		failReport.write("No UniParc or UniProt ID found:\n")
		for geneID1 in badID:
			failReport.write(str(geneID1)+'\n')
		failReport.write('\n')
		failReport.write("No Gene Ontology terms found:\n")
		for geneID2 in badGO:
			failReport.write(str(geneID2)+'\n')

def printReport(dic):
	#Method to print the tsv report file.

	with open("report"+time.strftime("_%d_%m_%Y_%H%M")+".tsv",'w') as outfile:
		outfile.write("Gene Group\tNon-unique gene name\tAnnotation\tUnique Gene ID\tUnique GI Number\tUnique RefSeq Protein Number\tUnique UniProt ID\t" + \
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
						for domainName, term in domain.items():
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

					outfile.write(geneGroupID+'\t'+annotation+'\t'+otherAnnotation+'\t' + \
					geneID+'\t'+giNumber+'\t'+refNumber+'\t'+uniprotID+'\t'+component+'\t'+process+'\t'+function+'\n')
		

def main():

	version=1.1

	parser = argparse.ArgumentParser(description='Gene Ontology fetcher for Roary and Scoary outputs.', epilog='by Catarina Mendes (cimendes@medicina.ulisboa.pt)')
	parser.add_argument('-g', '--genes', help='Input gene presence/absence table (comma-separated-values) from Roary (https:/sanger-pathogens.github.io/Roary)')
	parser.add_argument('-i', '--input', help='Input interest gene presence/absence table (comma or semicolon-separated-values from Roary or Scoary')
	parser.add_argument('-d', '--gffdir', help='Path to directory containing all gff files used in the Roary analysis.')
	parser.add_argument('--version', help='Display version, and exit.', default=False,action='store_true')

	args = parser.parse_args()

	#version
	if args.version:
		print sys.stdout, "Current version: %s" %(version)
		sys.exit(0)

	#check if something is missing
	if args.genes is None or args.input is None or args.gffdir is None:
		parser.print_usage()
		if args.genes is None:
			print "error: argument -g/--genes is required"
		if args.input is None:
			print "error: argument -i/--input is required"
		if args.gffdir is None:
			print "error: argument -d/--gffdir is required"
		sys.exit(1)

	start_time = time.time()

	sys.stdout = Logger("./")

	#Reading files
	print "Reading input file."
	input_file=parseInputCSVFile(args.input)

	print "Reading gene file."
	input_genes, filenames=parsePAGeneFile(args.genes,input_file)

	print "Parsing GFF files."
	input_genes=parseGFFs(args.gffdir,filenames,input_genes)

	print "Retrieving Uniprot IDs."
	input_genes,badID=retrieveUniprot(input_genes)

	print "Retrieving Gene Ontology terms."
	input_genes,noGO=getGOnumbers(input_genes)

	print "Printing final report."
	printReport(input_genes)
	badResultsReport(badID,noGO)

	end_time = time.time()
	time_taken = end_time - start_time
	hours, rest = divmod(time_taken,3600)
	minutes, seconds = divmod(rest, 60)
	print "Runtime :" + str(hours) + "h:" + str(minutes) + "m:" + str(round(seconds, 2)) + "s" + "\n"

	print "Finished"
	sys.exit(0)

if __name__ == "__main__":
    main()