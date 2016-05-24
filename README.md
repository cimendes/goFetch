# goFetch
Gene Ontology fetcher for Roary and Scoary outputs.

It's designed to take the gene_presence_absence.csv file from [Roary] (https://sanger-pathogens.github.io/Roary/), a gene input csv file obtained from Roary or [Scoary] (https://github.com/AdmiralenOla/Scoary) and the path to the directory containing the GFF files used in the Roary analysis, reporting a list of the input genes with information on IDs and Gene Ontology terms found.

## Updates:

Current release - v1.1 

- Fixed bug related to the presence of paralog genes in the gene_presence_absence.csv file from Roary.
- Printing a report for the geneIDs for which the script couldn't retrieve UniprotID or	any Gene Ontology terms.
- Printing log file.

## Usage

    goFetch.py -g gene_presence_absence.csv -i genes_file.csv -d '/path/to/GFF/files/'

    usage: goFetch.py [-h] [-g GENES] [-i INPUT] [-d GFFDIR] [--version]

	Gene Ontology fetcher for Roary and Scoary outputs.

	optional arguments:
		-h, --help            show this help message and exit
		-g GENES, --genes GENES
						Input gene presence/absence table (comma-separated-
						values) from Roary (https:/sanger-
						pathogens.github.io/Roary)
		-i INPUT, --input INPUT
						Input interest gene presence/absence table (comma or
						semicolon-separated-values from Roary or Scoary
		-d GFFDIR, --gffdir GFFDIR
						Path to directory containing all gff files used in the
						Roary analysis.
		--version  		Display version, and exit.

## Dependencies

- Python (2.7.x)
- [bioservices] (https://pythonhosted.org/bioservices/)
- Internet connection

## Installation

goFetch is a standalone python script and does not require any installation. Simply clone the git repository:

    git clone https://github.com/cimendes/goFetch

## Input
goFetch requires two input files, the gene_presence_absence.csv file from [Roary] (https://sanger-pathogens.github.io/Roary/) and a gene input csv file obtained from [Roary] (https://sanger-pathogens.github.io/Roary/) or [Scoary] (https://github.com/AdmiralenOla/Scoary), and the path to the directory containing the GFF files used in the Roary analysis.

## Output
goFetch outputs a single tsv file contaning the gene group and annotation. For each group with unique IDs it will print the GI Number, RefSeq Protein Number, UniprotID and GO terms found (Cellular Component, Biological Process and Molecular Function).

## License
goFetch is freely available under a GPLv3 license.

## Contact
Catarina Mendes (cimendes@medicina.ulisboa.pt)