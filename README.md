# PSSM_features
A set of codes written in Python and R which help you to calculate PSSM based features for protein sequences.

## Pre-Requisites 
0. Necessary libraries and dependencies (subprocess, Biopython, Numpy, Pandas, os, glob, shutil, csv, etc)
1. You need to have a locally installed database (consturcted using the blast makedb .... command), and you need to specify the prefix for the same in the bigram pssm and psuedo pssm codes, in the db_prefix variable.

## RUNBLAST 
### Generating PSSM Matrix from FASTA Sequence
The RUNBLAST.py code helps you to generate PSSM matrices for a protein, given its fasta sequence. The function inputs a fasta file, run PSI-BLAST locally and generates a *.pssm file containing the PSSM Matrix, along with a log_file.txt.
You can specify which database to run the PSI-BLAST depending upon your locally installed database, by modifying the 'db_prefix' variable.

## Bi-Gram PSSM
Generates PSSM scores for bigrams i.e amino acids at 'i' and 'i+1' position in the inputted protein sequence
Saves the result in a *.csv file

## Psuedo PSSM (PsePSSM)
Generates a 20 + 20 x lambda feature values for an input pssm matrix. The first 20 values are normalised Amino Acid Composition derived from PSSM matrices, and the next 20 x lambda feature values are PsePSSM feature values calculated using the input pssm matrix and the input lag parammeter (lambda)

## PSSM Standard Features (PSSM_Completepipe)
This code written in R utilises a R package namely PSSMCOOL and generates various PSSM based features at a single go and saves all the results in a *.csv file with appropriate headers. This piepline enhances the current functionality of the PSSMCOOL package, imporves the readability of the result and can be used for high throughput result generation.


