#!/bin/bash
mkdir -p GSE15434_RAW/stash
mkdir -p GSE6891_RAW/stash
mkdir -p GSE13204_RAW/stash
mkdir -p GSE1159_RAW/stash
mkdir -p GSE22845_RAW/stash
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE15434&format=file \
	-O GSE15434_RAW/GSE15434_RAW.tar 
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE6891&format=file \
	-O GSE6891_RAW/GSE6891_RAW.tar
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE13204&format=file \
	-O GSE13204_RAW/GSE13204_RAW.tar
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE1159&format=file \
	-O GSE1159_RAW/GSE1159_RAW.tar
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE22845&format=file \
	-O GSE22845_RAW/GSE22845_RAW.tar
