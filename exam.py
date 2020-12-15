#!/usr/bin/python3

import sys, subprocess, shutil, os
import string, re
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Bio import Entrez

taxon_gp = ""
choice1 = ""

def getInput():
	global taxon_gp           #set it as a global varibale
	global choice1            #set it as global variable
	choice1 = input("Please indicate nucleotide or protein to be searched?, n/p\n")
	if choice1 == "n":
		nucle = input("Please enter the gene name\n")
		taxon_gp = input("Please enter the taxonomic group name\n")
		print("Just a reminder, there may have some partial sequences in this dataset.")
		pseq = input("Do you want to include those partial sequences when searching?, Y/N\n")
		if pseq == "Y":
			es1 = 'esearch -db nucleotide -query \" '+ taxon_gp +'[orgn] AND '+ nucle + '[Gene Name]\" | efetch -db nucleotide -format fasta > seq.fa '
		if pseq == "N":
			es1 = 'esearch -db nucleotide -query \" '+ taxon_gp +'[orgn] AND '+ nucle + '[Gene Name] Not partial\" | efetch -db nucleotide -format fasta > seq.fa '
		print("Thanks, you have chosen " + nucle + " in " + taxon_gp + "\n")
	if choice1 == "p":
		protein = input("Please enter the protein name\n")
		taxon_gp = input("Please enter the taxonomic group name\n")
		print("Just a reminder, there may have some partial sequences in this dataset.")
		pseq = input("Do you want to include those partial sequences when searching?, Y/N\n")
		if pseq == "Y":
			es1 = 'esearch -db protein -query \" '+ taxon_gp +'[orgn] AND '+ protein + '[Protein Name]\" | efetch -db protein -format fasta > seq.fa '
		if pseq == "N":
			es1 = 'esearch -db protein -query \" '+ taxon_gp +'[orgn] AND '+ protein + '[Protein Name] Not partial\" | efetch -db protein -format fasta > seq.fa '
		print("Thanks, you have chosen " + protein + " in " + taxon_gp + "\n")
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
	while seq_Num <= 1:
		print("Sorry,there is something wrong with your input - no result. Please input again!")
		if seq_Num >= 1:
			break
		getInput()
		seq_Num = findSeq()
		return seq_Num
sequence_check(seq_Num)

#Tell the user the sequence number and ask them if they are willing to continue
print("There are " + str(seq_Num) + " sequences in the dataset you have chosen")
choice2 = input("Do you want to continue?,Y/N\n")
while choice2 == "N":
	if choice2 == "Y":
		break
	print("Thanks, please then input what you want again")
	getInput()
	seq_Num = findSeq()
	sequence_check(seq_Num)
	print("There are " + str(seq_Num) + " sequences in the dataset you have chosen")
	choice2 = input("Do you want to continue?,Y/N\n")
	
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
choice3 = input("Do you want to continue with the current dataset? Y/N\n")
while choice3 == "N":
	if choice3 == "Y":
		break
	print("Thanks, please then input what you want again")
	getInput()
	seq_Num = findSeq()
	spe_N = findSpec()
	sequence_check(seq_Num)
	print("There are " + str(seq_Num) + " sequences in the dataset you have chosen")
	print("There are " + str(spe_N) + " species in this dataset")
	choice3 = input("Do you want to continue?,Y/N\n")

#move on to the BLAST analysis
print("OK, Let's make BLAST database first")
if choice1 == "n":
	mdb = "makeblastdb -in seq.fa -dbtype nucl -out " + taxon_gp              #make BLAST database
if choice1 == "p":
	mdb = "makeblastdb -in seq.fa -dbtype prot -out " + taxon_gp
#print(mdb)
subprocess.call(mdb,shell=True)
#asking the user to determine the query sequence for BLAST analysis
print("You may have 2 options to determine the query sequence. 1: enter a Seq ID; 2: enter the whole sequence")



##define a function to check the Seq ID that the user has inputted. This is for the next step
def check_ID(ID):     
	#when you do esearch it print out ENTREZ_DIRECT, where data count can be found, check the value
	if choice1 == "n":
		value = subprocess.check_output('esearch -db nucleotide -query "+ ID + " | xtract -pattern ENTREZ_DIRECT -element Count',shell=True) 
	if choice1 == "p":
		value = subprocess.check_output('esearch -db protein -query " + ID+ " | xtract -pattern ENTREZ_DIRECT -element Count',shell=True)
	if int(value) == 0:      #if there is no data found thenfunctionss return 0
		return 0
	return 1         #the ID is correct, we can continue
##let the user to make the decision
choice4 = input("Please enter your choice. 1/2\n")
if choice4 == "1":
	ID = input("Please enter a Seq ID that you would like to do BLAST\n")
	while check_ID(ID) != 1:
		ID = input("Sorry, your Seq ID may not be correct. Please enter the correct one!\n")
	print("Thanks, the sequence you have chosen will be downloaded and served as the test sequence for BLAST.")
	if choice1 == "n":
		os.system("esearch -db nucleotide -query " + ID + " | efetch -db nucleotide -format fasta > test_seq.fa")
		blt = "blastn -db " + taxon_gp + " -query test_seq.fa -outfmt 6 > blastoutput.out"
		#print(blt)
		subprocess.call(blt,shell=True)
		print("BLAST analysis has been successfully done!")
	if choice1 == "p":
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
	if choice1 == "n":
		os.system("esearch -db nucleotide -query " + ID + " | efetch -db nucleotide -format fasta > test_seq.fa")
		blt = "blastn -db " + taxon_gp + " -query test_seq.fa -outfmt 6 > blastoutput.out"
		#print(blt)
		subprocess.call(blt,shell=True)
		print("BLAST analysis has been successfully done!")
	if choice1 == "p":
		os.system("esearch -db protein -query " + ID + " | efetch -db protein -format fasta > test_seq.fa")
		blt = "blastp -db " + taxon_gp + " -query test_seq.fa -outfmt 6 > blastoutput.out"
		#print(blt)
		subprocess.call(blt,shell=True)
		print("BLAST analysis has been successfully done!")

#ask the user if they want to see the results
choice5 = input("Now BLAST results are ready. Do you want to display it?, Y/N\n")
if choice5 == "Y":
	dis = "cat blastoutput.out"
	subprocess.call(dis,shell = True)
if choice5 == "N":
	print("OK, you can check the results later.")	


print("Thanks for using this programme!")
