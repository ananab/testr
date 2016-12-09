#!/usr/bin/make -f

# tag all. testdata.pdf is dependent on all other outputs and scripts its the only one needed here
all: results/testdata.pdf

# give clean tag to remove all results in results folder as well as data downloaded
clean: 
	rm -f results/*

# cleaned lkup table is dependent on testdata.Rmd as well as the original lkup
results/testdata.pdf: src/launch.R data/testData.csv data/uniprot-all.fasta
	Rscript src/launch.R ../data/testData.csv ../data/uniprot-all.fasta
