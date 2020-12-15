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


print("OK, Let's make BLAST database first")
if choice1 == "nucleotide":
	mdb = "makeblastdb -in nucleotide_seq.fa -dbtype prot -out " + taxon_gp              #make BLAST database
if choice1 == "protein":
	mdb = "makeblastdb -in protein_seq.fa -dbtype prot -out " + taxon_gp
#print(mdb)
subprocess.call(mdb,shell=True)
#asking the user to determine the query sequence for BLAST analysis
print("You may have 2 options to determine the query sequence. 1: enter a Seq ID; 2: enter the whole sequence")

def check_ID(ID):     
	#when you do esearch it print out ENTREZ_DIRECT, where you can find found data count, I try to save that value and check
	if choice1 == "nucleotide":
		count = subprocess.check_output('esearch -db nucleotide -query "{0}" | xtract -pattern ENTREZ_DIRECT -element Count'.format(ID),shell=True) 
	if choice1 == "protein":
		count = subprocess.check_output('esearch -db protein -query "{0}" | xtract -pattern ENTREZ_DIRECT -element Count'.format(ID),shell=True)
	if int(count) == 0:      #if there is no data found thenfunctionss return 0 further it gaves error message to the user
			return 0
	return 1         #if there is at least one data found regarding to the user input we continue p

choice4 = input("Please enter your choice. 1/2\n")
if choice4 == "1":
	ID = input("Please enter a Seq ID that you would like to do BLAST\n")
	while check_ID(ID) != 1:
		ID = input("Sorry, your Seq ID may not be correct. Please enter the correct one!")
	print("Thanks, the sequence you have chosen will be downloaded and served as the test sequence for BLAST.")
	if choice1 == "nucleotide":
		os.system("esearch -db nucleotide -query " + ID + " | efetch -db nucleotide -format fasta > test_seq.fa")
		blt = "blastp -db " + taxon_gp + " -query test_seq.fa -outfmt 6 > blastoutput.out"
		#print(blt)
		subprocess.call(blt,shell=True)
		print("BLAST analysis has been successfully done!")
	if choice1 == "protein":
		os.system("esearch -db protein -query " + ID + " | efetch -db protein -format fasta > test_seq.fa")
		blt = "blastp -db " + taxon_gp + " -query test_seq.fa -outfmt 6 > blastoutput.out"
		#print(blt)
		subprocess.call(blt,shell=True)
		print("BLAST analysis has been successfully done!")
if choice4 == "2":
	iseq = input("Please enter the whole sequence")
	print("Thanks, we will use this sequence to do BLAST!")
	seq_of_interest = open('seq_of_interest.fa',"w")                 #creat a new file for storing the motifs information, make it writable
	seq_of_interest.write(">Sequence_of_Interest\n")
	seq_of_interest.write('{0}'.format(iseq))  
	seq_of_interest.close()
	blt = "blastp -db " + taxon_gp + " -query seq_of_interest.fa -outfmt 6 > blastoutput.out"
	#print(blt)
	subprocess.call(blt,shell=True)
	print("BLAST analysis has been successfully done!")
