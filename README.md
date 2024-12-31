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



