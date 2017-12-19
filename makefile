.PHONY : all
.DEFAULT: all
PACKAGE      = package
VERSION      = ` date "+%Y.%m%d%" `
RELEASE_DIR  = ..
RELEASE_FILE = $(PACKAGE)-$(VERSION)
SILENT 	     = yes
DATA_DIR     = data/datasets

all: generate_files generate_figures	

# !!!Download the datasets from https://doi.org/10.6084/m9.figshare.5715880 and place them in the data/datasets folder

#Generates all the files necessary for the analyses
generate_files: 
	perl generate_files.pl;
	perl add_orthology_info.pl
	#Generates all the figures from the analyses described in the article: TISSUES: An integrative web resource on mammalian tissue expression
generate_figures:
	Rscript analyses.R




