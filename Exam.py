#!/usr/bin/python3

import sys, subprocess, shutil, os
import string, re
import collections
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Bio import Entrez


def getInput():
	choice1 = input("Please indicate nucleotide or protein to be searched?, nucleotide/protein\n")
	if choice1 == "nucleotide":
		nucle = input("Please enter the gene name\n")
		taxon_gp = input("Please enter the taxonomic group name\n")
		es1 = 'esearch -db nucleotide -query \" '+ taxon_gp +'[orgn] AND '+ nucle + '[Gene Name]\" | efetch -db protein -format fasta > nucleotide_seq.fa '
	if choice1 == "protein":
		protein = input("Please enter the protein name\n")
		taxon_gp = input("Please enter the taxonomic group name\n")
		print("Thanks, you have chosen " + protein)
		es1 = 'esearch -db protein -query \" '+ taxon_gp +'[orgn] AND '+ protein + '[Protein Name]\" | efetch -db protein -format fasta > protein_seq.fa '
	print("This is what will be run: " + es1)
#Run the esearch and efetch command to get the dataset from NCBI
	return subprocess.call(es1,shell=True)

getInput()
print("The sequences have been downloaded and stored")

#Find how many sequences are there in the dataset user have chosen
def findSeq():
	seqnumber= subprocess.getoutput("grep -c \">\" seq.fa")
	seq_number = int(seqnumber)
	return seq_number
seq_Num = findSeq()


#Test if the the number of the sequences are appropriate and give the user options to decide if they want to continue
def sequence_check(seq_Num):
	while seq_Num <= 1 or seq_Num >= 1000:
		print("Sorry,there is something wrong with your input - no result or over 1000 sequence. Please input again!")
		if seq_Num >= 1 and seq_Num <= 1000:
                        break
		getInput()
		seq_Num = findSeq()
		return seq_Num
sequence_check(seq_Num)

#Tell the user the sequence number and ask them if they are willing to continue
print("There are " + str(seq_Num) + " sequences in the dataset you have chosen")
choice1 = input("Do you want to continue?,Y/N\n")
while choice1 == "N":
	if choice1 == "Y":
		break
	print("Thanks, please then input what you want again")
	getInput()
	seq_Num = findSeq()
	print("There are " + str(seq_Num) + " sequences in the dataset you have chosen")
	choice1 = input("Do you want to continue?,Y/N\n")
	
#find the number of species in this dataset
def findSpec():
	get_header = "grep \">\" seq.fa > seq_header.fa"
	subprocess.call(get_header,shell=True)
	seq = open("seq_header.fa")
	spe_num = []
	for i in seq:
		l1 = i.find('[')
		l2 = i.find(']')
		spe_num.append(i[l1:l2])
	spe_n = len(set(spe_num))
	return spe_n
spe_N = findSpec()
#Tell the user the number of species in this dataset and ask them if you want to continue
print("There are " + str(spe_N) + " species in this dataset")
choice2 = input("Do you want to continue with the current dataset? Y/N\n")
while choice2 == "N":
	if choice2 == "Y":
		break
	print("Thanks, please then input what you want again")
	getInput()
	seq_Num = findSeq()
	spe_N = findSpec()
	print("There are " + str(spe_N) + " species in this dataset")
	choice2 = input("Do you want to continue?,Y/N\n")

mdb = "makeblastdb -in seq.fa -dbtype prot -out " + taxon_gp              #make BLAST database
#print(mdb)
subprocess.call(mdb,shell=True)
