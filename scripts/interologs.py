"""
	Step  1 in the pipeline - s.cere vs s.pombe
	
	
"""

import os
import pickle
from collections import defaultdict,namedtuple
import csv
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
from tqdm import tqdm
from  ast import literal_eval
from modeled_PPIs import PPI,Protein,Interface,Residue
# import sys
# import zipfile
# import subprocess
# # import matplotlib.pyplot as plt
# # import numpy as np
# # import seaborn as sns


def manually_curated_orthologs(orthologues_file,spom_PPIs,scer_PPIs):
	"""
	Get the ortholog mapping between S. pombe and S. cerevisiae ORF names from manually curated file on PomBase
	Format mapping into two dictionaries (keys going both directions: pombe -> cerevisiae and cerevisiae -> pombe)
	Args:
		orthologues_file (str): path to the file with manual orthology mapping from PomBase
		spom_PPIs/scer_PPIs (dict): dictionaries with the PPIs gathered in each species. Format: {('ORF1','ORF2'):int(number of reports)} 	with ('ORF1','ORF2') sorted alphabetically
		multiple_reports (bool): whether or not we want to print reports only for PPIs reported more than once (different pumed IDs)
	"""
	print ("-------- Get ortholog mapping (manually curated) --------")
	orthologues_dictionary = {"S_cere_keys":defaultdict(dict),"S_pombe_keys":defaultdict(dict)} #Format: {"S_cere_keys":{Scer_orf:set(Spom_orf,Spom_orf...)} ,"S_pombe_keys":{Spom_orf:set(Scer_orf,Scer_orf...)}}

	S_pombe_orfs=set([orf for ppi in spom_PPIs.keys() for orf in ppi])
	S_cere_orfs=set([orf for ppi in scer_PPIs.keys() for orf in ppi])
	
	with open(orthologues_file, 'r') as f:
		reader = csv.reader(f,delimiter = '\t', )
		next(reader)
		for row in reader:
			spomORF = row[0]
			# format the scer orf list: split on "|", remove text in parenthesis ((N), (C),(FUSION-N), (FUSION-C))
			scerORFs= [re.sub(r'\(.*\)', '', orf) for scerORF in [row[1].split('|')] for orf in scerORF ]

			if (scerORFs) == ["NONE"]:
				scerORFs = []
			# Only add to the ortholog dictionary if the ORF is involved in a PPI
			if spomORF in S_pombe_orfs:
				orthologues_dictionary["S_pombe_keys"].setdefault(spomORF,set()).update(scerORFs)
			else:
				continue
			# Add the reverse mapping
			for scerORF in scerORFs:
				if scerORF in S_cere_orfs:
					orthologues_dictionary["S_cere_keys"].setdefault(scerORF,set()).update([spomORF])

		
	print("Ortholog mapping info for",len(orthologues_dictionary["S_pombe_keys"].keys()),"/",len(S_pombe_orfs),"S. pombe proteins")
	print("Ortholog mapping info for",len(orthologues_dictionary["S_cere_keys"].keys()),"/",len(S_cere_orfs),"S. cerevisiae proteins")

	return orthologues_dictionary

def reciprocal_BLAST_orthologs(spom_PPIs,scer_PPIs,inDir):
	"""
	Use reciprocal best hit BLAST to find orthologs between S. cerevisiae and S. pombe orf in PPIs
	Format mapping into two dictionaries (keys going both directions: pombe -> cerevisiae and cerevisiae -> pombe)
	Args:
		spom_PPIs/scer_PPIs (dict): dictionaries with the PPIs gathered in each species. Format: {('ORF1','ORF2'):int(number of reports)} 	with ('ORF1','ORF2') sorted alphabetically
		inDir (str): path to directory with sequence info for the two species
		multiple_reports (bool): whether or not we want to print reports only for PPIs reported more than once (different pumed IDs)
	"""
	print ("-------- Get ortholog mapping (reciprocal BLAST) --------")
	Scer_aa,Spom_aa = write_sequences(spom_PPIs,scer_PPIs,inDir)
	eValue = '1e-5'
	
	ofiles = ['Scer_reduced_aa-Spom_reduced_aa-blast.best.txt','Spom_reduced_aa-Scer_reduced_aa-blast.best.txt']

	# Run Blast both ways
	for species1,species2,ofile  in zip([Scer_aa,Spom_aa],[Spom_aa,Scer_aa],ofiles):
		opath = os.path.join(inDir,ofile)
		run_blast(species1, species2, opath, eValue, singleBest=True)

	fwd_results = pd.read_csv(os.path.join(inDir,ofiles[0]), sep="\t", header=None)
	rev_results = pd.read_csv(os.path.join(inDir,ofiles[1]), sep="\t", header=None)
	headers = ["query", "subject", "identity", "coverage","qlength", "slength", "alength","bitscore", "E-value"]
	fwd_results.columns = headers
	rev_results.columns = headers

	# Get reciprocal best hits (RBH)
	rbh = pd.merge(fwd_results, rev_results[['query', 'subject']],left_on='subject', right_on='query',how='outer') # Merge forward and reverse results
	rbh = rbh.loc[rbh.query_x == rbh.subject_y] # Discard rows that are not RBH
	rbh = rbh.groupby(['query_x', 'subject_x']).max() # Group duplicate RBH rows, taking the maximum value in each column

	orthologues_dictionary = {"S_cere_keys":defaultdict(dict),"S_pombe_keys":defaultdict(dict)} #Format: {"S_cere_keys":{Scer_orf:set(Spom_orf,Spom_orf...)} ,"S_pombe_keys":{Spom_orf:set(Scer_orf,Scer_orf...)}}
	for index, row in rbh.iterrows():
		scerORF = row['query_y']
		spomORF = row['subject_y']
		orthologues_dictionary["S_cere_keys"].setdefault(scerORF,set()).update([spomORF])
		orthologues_dictionary["S_pombe_keys"].setdefault(spomORF,set()).update([scerORF])

	
	# Check how many of the orf in both species are not in this ortholog dictionary
	S_pombe_orfs=set([orf for ppi in spom_PPIs.keys() for orf in ppi])
	S_cere_orfs=set([orf for ppi in scer_PPIs.keys() for orf in ppi])
	has_mapping_spom = [orf for orf in S_pombe_orfs if orf in orthologues_dictionary["S_pombe_keys"].keys()]
	has_mapping_scer = [orf for orf in S_cere_orfs if orf in orthologues_dictionary["S_cere_keys"].keys()]
	print("Ortholog mapping info for",len(has_mapping_spom),"/",len(S_pombe_orfs),"S. pombe proteins")
	print("Ortholog mapping info for",len(has_mapping_scer),"/",len(S_cere_orfs),"S. cerevisiae proteins")
	

	return orthologues_dictionary

def write_sequences(spom_PPIs,scer_PPIs,inDir):
	"""
	Write fasta files for the ORF involved in S. cerevisiae and S.pombe PPIs (subset of the full sequence file for each species). 
	Args:
		spom_PPIs/scer_PPIs (dict): dictionaries with the PPIs gathered in each species. Format: {('ORF1','ORF2'):int(number of reports)} 	with ('ORF1','ORF2') sorted alphabetically
		inDir (str): path to directory with sequence info for the two species
	"""

	if not os.path.exists(inDir):
		os.makedirs(inDir)

	if not os.path.exists(os.path.join(inDir,"Scer_reduced.aa")):

		S_pombe_orfs=set([orf for ppi in spom_PPIs.keys() for orf in ppi])
		S_cere_orfs=set([orf for ppi in scer_PPIs.keys() for orf in ppi])
		
		scer_aa_file = os.path.join(inDir,"Scer.aa")
		spom_aa_file = os.path.join(inDir,"Spom.aa")

		for inFile,orfs_to_keep,charsToStrip  in zip([scer_aa_file,spom_aa_file],[S_cere_orfs,S_pombe_orfs],[5,6]):
			outFile = os.path.join(inDir,os.path.basename(inFile)[:-3]+"_reduced.aa") 	# Path for new file with only the sequences for ORFs involved in PPIs
			recordsToKeep = []										 					# Sequences for ORFs involved in PPIs
			for record in SeqIO.parse(inFile, "fasta"):
				orf = record.id[:-charsToStrip]
				sequence = record.seq
				if(orf in orfs_to_keep):
					newRecord = SeqRecord(sequence, orf,'','')
					recordsToKeep.append(newRecord)
					
			SeqIO.write(recordsToKeep,outFile, "fasta")
	
	return os.path.join(inDir,"Scer_reduced.aa"),os.path.join(inDir,"Spom_reduced.aa")

def run_blast(databaseFile, queryFile, opath, eValue, singleBest = False):
	"""
	Run a blast protein alignment to find either 
	(1) the single best match between database and query (if singleBest == True)
	(2) all matches between database and query (if singleBest == False)

	Args:
		databaseFile (str): path to file containing the protein aa sequences used as the target database
		queryFile (str): path to file containing the protein aa sequences used as the query 
		opath (str): path to output file containing results of blast.
		eValue (str): e-value cutoff for results.
	"""
	if not os.path.exists(databaseFile+'.pin'):
		os.system('makeblastdb -in ' + databaseFile + ' -dbtype prot')
	if os.path.exists(opath):
		print ('Blast matching statistics already exist. Skipping.')
		return
	
	fmt = "\'6 qseqid sseqid pident qcovs qlen slen length bitscore evalue\'"
	if(singleBest):
		maxResNum = ' -max_target_seqs 1'
	else:
		maxResNum = ''
	os.system('blastp -query ' + queryFile +
			  ' -db ' + databaseFile +
			  ' -evalue ' + eValue +
			  ' -outfmt ' + fmt +
			  ' -num_threads 4 > ' + opath +
			  maxResNum +
			  ' 2>/dev/null') # Redirect stderr to /dev/null to supress warnings
	return

def find_homologs_BLAST(spom_PPIs,scer_PPIs,inDir):
	"""
	Run BLAST alignements to identify orfs without any homologs in the other species
	If an ORF has no homolog in the other species, then we can conclude that it has no ortholog
	Args:
		spom_PPIs/scer_PPIs (dict): dictionaries with the PPIs gathered in each species. Format: {('ORF1','ORF2'):int(number of reports)} 	with ('ORF1','ORF2') sorted alphabetically
		inDir (str): path to directory with sequence info for the two species
		multiple_reports (bool): whether or not we want to print reports only for PPIs reported more than once (different pumed IDs)

	"""
	print ("-------- Get homolog mapping (BLAST) --------")
	Scer_aa,Spom_aa = write_sequences(spom_PPIs,scer_PPIs,inDir)
	eValue = '1e-5'
	
	ofiles = ['Scer_reduced_aa-Spom_reduced_aa-blast.txt','Spom_reduced_aa-Scer_reduced_aa-blast.txt']
	# Run Blast both ways
	for species1,species2,ofile in zip([Scer_aa,Spom_aa],[Spom_aa,Scer_aa],ofiles):
		opath = os.path.join(inDir,ofile)
		run_blast(species1, species2, opath, eValue,singleBest = False)

	homologs_dictionary = {"S_cere_keys":defaultdict(dict),"S_pombe_keys":defaultdict(dict)} #Format: {"S_cere_keys":{Scer_orf:set(Spom_orf,Spom_orf...)} ,"S_pombe_keys":{Spom_orf:set(Scer_orf,Scer_orf...)}}
	
	for result,key in zip(ofiles,["S_cere_keys","S_pombe_keys"]):
		with open(os.path.join(inDir,result)) as f:
			for items in csv.reader(f, dialect="excel-tab" ):
				orf = items[1]
				homologue = items[0]
				homologs_dictionary[key].setdefault(orf,set()).update([homologue])
	
	S_pombe_orfs=set([orf for ppi in spom_PPIs.keys() for orf in ppi])
	S_cere_orfs=set([orf for ppi in scer_PPIs.keys() for orf in ppi])
	has_noHomolog_spom = [orf for orf in S_pombe_orfs if orf not in homologs_dictionary["S_pombe_keys"].keys()]
	has_noHomolog_scer = [orf for orf in S_cere_orfs if orf not in homologs_dictionary["S_cere_keys"].keys()]
	no_homologs_dictionary = {"S_cere_keys":set(has_noHomolog_scer),"S_pombe_keys":set(has_noHomolog_spom)}  #Format: {"S_cere_keys":set(Scere_orf,Scere_orf...)} ,"S_pombe_keys":{set(Spom_orf,Spom_orf...)}} with Scere_orf and Spom_orf orfs with no homologs
		
	
	# Check how many of the orf in both species are not in the homologs dictionary
	S_pombe_orfs=set([orf for ppi in spom_PPIs.keys() for orf in ppi])
	S_cere_orfs=set([orf for ppi in scer_PPIs.keys() for orf in ppi])
	has_noHomolog_spom = [orf for orf in S_pombe_orfs if orf in no_homologs_dictionary["S_pombe_keys"]]
	has_noHomolog_scer = [orf for orf in S_cere_orfs if orf in no_homologs_dictionary["S_cere_keys"]]
	print("No homologs for",len(has_noHomolog_spom),"/",len(S_pombe_orfs),"S. pombe proteins")
	print("No homologs for",len(has_noHomolog_scer),"/",len(S_cere_orfs),"S. cerevisiae proteins")
	
	return no_homologs_dictionary

def get_ortholog_mapping(spom_PPIs,scer_PPIs,orthologs_manual,orthologs_BLAST,no_homologues):
	"""
	Combine the orthologs obtained from the manually curated file and the orthologs obtained from reciprocal BLAST 
	Format mapping into two dictionaries (keys going both directions: pombe -> cerevisiae and cerevisiae -> pombe)
	Args:
		spom_PPIs/scer_PPIs (dict): dictionaries with the PPIs gathered in each species. Format: {('ORF1','ORF2'):int(number of reports)} 	with ('ORF1','ORF2') sorted alphabetically
		orthologs_manual (dict): dictionary with manually curated orthologs. Format: #Format: {"S_cere_keys":{Scer_orf:set(Spom_orf,Spom_orf...)} ,"S_pombe_keys":{Spom_orf:set(Scer_orf,Scer_orf...)}}
		orthologs_BLAST (dict): dictionary orthologs from reciprocal BLAST. Format: #Format: {"S_cere_keys":{Scer_orf:set(Spom_orf,Spom_orf...)} ,"S_pombe_keys":{Spom_orf:set(Scer_orf,Scer_orf...)}}
	"""

	print ("----------- Combine ortholog mappings -----------")
	orthologues_dictionary = {"S_cere_keys":defaultdict(dict),"S_pombe_keys":defaultdict(dict)} # Format: #Format: {"S_cere_keys":{Scer_orf:set(Spom_orf,Spom_orf...)} ,"S_pombe_keys":{Spom_orf:set(Scer_orf,Scer_orf...)}}
	
	for dictionary in [orthologs_manual,orthologs_BLAST]:
		for key in ["S_cere_keys","S_pombe_keys"]:
			for orf in dictionary[key].keys():
				orthologues_dictionary[key].setdefault(orf,set()).update(dictionary[key][orf])
	
	for key in ["S_cere_keys","S_pombe_keys"]:
		for orf in no_homologues[key]:
			orthologues_dictionary[key].setdefault(orf,set())
	
	#  Check how many of the orf in both species are not in this ortholog dictionary
	S_pombe_orfs=set([orf for ppi in spom_PPIs.keys() for orf in ppi])
	S_cere_orfs=set([orf for ppi in scer_PPIs.keys() for orf in ppi])
	has_mapping_spom = [orf for orf in S_pombe_orfs if orf in orthologues_dictionary["S_pombe_keys"].keys()]
	has_mapping_scer = [orf for orf in S_cere_orfs if orf in orthologues_dictionary["S_cere_keys"].keys()]
	
	print("Ortholog mapping info for",len(has_mapping_spom),"/",len(S_pombe_orfs),"S. pombe proteins")
	print("Ortholog mapping info for",len(has_mapping_scer),"/",len(S_cere_orfs),"S. cerevisiae proteins")

	return orthologues_dictionary

def	edit_ortholog_mapping(spom_PPIs,scer_PPIs,orthologs):
	"""
	Set ORFs with no ortholog mapping (either manual or reciprocal BLAST) to ortholog == NONE
	Args:
		spom_PPIs/scer_PPIs (dict): dictionaries with the PPIs gathered in each species. Format: {('ORF1','ORF2'):int(number of reports)} 	with ('ORF1','ORF2') sorted alphabetically
		orthologs (dict): dictionary with manually curated + reciprocal BLAST orthologs. Format: #Format: {"S_cere_keys":{Scer_orf:set(Spom_orf,Spom_orf...)} ,"S_pombe_keys":{Spom_orf:set(Scer_orf,Scer_orf...)}}
	"""
	print ("----------- Edit ortholog mappings -----------")
	previous_scer, previous_spom = len(orthologs["S_cere_keys"].keys()),len(orthologs["S_pombe_keys"].keys())
	
	new_orthologs = dict(orthologs)

	S_pombe_orfs=set([orf for ppi in spom_PPIs.keys() for orf in ppi])
	S_cere_orfs=set([orf for ppi in scer_PPIs.keys() for orf in ppi])
	
	for orfs, key in zip([S_pombe_orfs,S_cere_orfs],["S_pombe_keys","S_cere_keys"]):
		for orf in orfs:
			if orf not in new_orthologs[key].keys():
				new_orthologs[key][orf] = set()
				
	current_scer, current_spom = len(new_orthologs["S_cere_keys"].keys()),len(new_orthologs["S_pombe_keys"].keys())
	
	print("Set",current_spom-previous_spom,"S. pombe proteins and",current_scer-previous_scer, "S. cerevisiae proteins as having no ortholog in the other species",)
	print(current_spom,"/",len(S_pombe_orfs),"S. pombe proteins and",current_scer,"/",len(S_cere_orfs), "S. cerevisiae proteins have ortholog mapping information",)

	return orthologs

def print_summary(orthologs,summary_spom,summary_scer):

	print ("----------- Ortholog summary files -----------")

	if not os.path.exists(summary_spom):
		with open(summary_spom, 'w') as f:

			for ORF in orthologs["S_pombe_keys"]:
				
				orthos=list(orthologs["S_pombe_keys"][ORF])

				if len(orthos)==0:
					f.write(ORF+"\t"+"NONE"+"\n")
				else:
					f.write(ORF+"\t"+'|'.join([ortho for ortho in orthologs["S_pombe_keys"][ORF]])+"\n")

	if not os.path.exists(summary_scer):
		with open(summary_scer, 'w') as f:

			for ORF in orthologs["S_cere_keys"]:
				
				orthos=list(orthologs["S_cere_keys"][ORF])

				if len(orthos)==0:
					f.write(ORF+"\t"+"NONE"+"\n")
				else:
					f.write(ORF+"\t"+'|'.join([ortho for ortho in orthologs["S_cere_keys"][ORF]])+"\n")

	return

def compare_PPIs_between_species(spom_PPIs,scer_PPIs,orthologs,ppi_type_summary_file,re_run=False):
	"""
	Compare the list of PPIs between the two species to classify them as 
	- preserved: (A and B) & (A' and B') interact (with X and X' orthologs between S. cerevisiae and S. pombe)
	- rewired: (A and B) interact & (A' and B') do not interact (with X and X' orthologs between S. cerevisiae and S. pombe)
	- lost: (A and B) interact & (A' or B') does not exist (with X and X' orthologs between S. cerevisiae and S. pombe)
	Args:
		spom_PPIs/scer_PPIs (dict): dictionaries with the PPIs gathered in each species. Format: {('ORF1','ORF2'):int(number of reports)} 	with ('ORF1','ORF2') sorted alphabetically
		orthologs (dict): dictionary with orthologs between S. cerevisiae and S. pombe. Format: #Format: {"S_cere_keys":{Scer_orf:set(Spom_orf,Spom_orf...)} ,"S_pombe_keys":{Spom_orf:set(Scer_orf,Scer_orf...)}}
		ppi_type_summary_file (str): path to the desired location for the output file
		re_run (bool): if True, re-run the comparison, overwritting any previously writted PPI type summary files
	"""
	print ("----------- Compare PPIs between species -----------")
	PPI_comparison = {"preserved":[],"lost":[],"rewired":[]}
	missingORFs = 0
	lost_count=0
	preserved_count=0
	rewired_count=0
	# Not already done
	if ((not os.path.exists(ppi_type_summary_file)) or re_run):

		# Iterate through S. pombe and S. cerevisiae PPIs
		for PPIs_in_species,PPIs_in_ortho_species, key, text in zip([scer_PPIs,spom_PPIs],[spom_PPIs,scer_PPIs], ["S_cere_keys","S_pombe_keys"],["S. cerevisiae","S. pombe"]):
			print("Checking orthologs for",text,"PPIs")
			for (orf1,orf2) in tqdm(PPIs_in_species.keys()):
				# print(orf1,orf2)
				# print(key)
				
				PPI = namedtuple('PPI', ['SpomPPI', 'ScerPPI',"reportsInSpom","reportsInScer","ppiType"])
				# (1) One of the two partners is not in the orthologs dict (cRNS, snRNA, snoRNA and orf that changed names are excluded here)
				if ((orf1 not in orthologs[key].keys()) or (orf2 not in orthologs[key].keys())):
					missingORFs = missingORFs+1
					continue
				# (2) Either or both partner(s) do(es) not have an ortholog in the other species (lost PPI)
				ortho1 = None if orthologs[key][orf1]==set() else tuple(orthologs[key][orf1])
				ortho2 = None if orthologs[key][orf2]==set() else tuple(orthologs[key][orf2])
				# print(ortho1)
				# print(ortho2)
				# print(PPIs_in_species[(orf1,orf2)])
				if(ortho1 == None or ortho2 == None):
					if key == "S_pombe_keys":
						ppi = PPI((tuple([orf1]),tuple([orf2])),(ortho1,ortho2),str(PPIs_in_species[(orf1,orf2)].numberOfReports),0,"lost")
						
					else:
						ppi = PPI((ortho1,ortho2),(tuple([orf1]),tuple([orf2])),0,str(PPIs_in_species[(orf1,orf2)].numberOfReports),"lost")
						
					
					lost_count = lost_count+1
					PPI_comparison["lost"].append(ppi)
					
				else: # Both partners have an ortholog in the other species
					ortho_PPIs = [(x,y) for x in ortho1 for y in ortho2]
					# (3) The two orthologous partners do not interact in the other species (rewired PPI)
					if not any(x in PPIs_in_ortho_species.keys() for x in ortho_PPIs):
						if key == "S_pombe_keys":
							ppi = PPI((tuple([orf1]),tuple([orf2])),(ortho1,ortho2),str(PPIs_in_species[(orf1,orf2)].numberOfReports),0,"rewired")
							
						else:
							ppi = PPI((ortho1,ortho2),(tuple([orf1]),tuple([orf2])),0,str(PPIs_in_species[(orf1,orf2)].numberOfReports),"rewired")
							

						rewired_count = rewired_count+1
						PPI_comparison["rewired"].append(ppi)
					# (4) The two orthologous partners interact in the other species (preserved PPI)
					else:
						preserved_count=preserved_count+1
						for ortho_PPI in ortho_PPIs:
							if ortho_PPI in PPIs_in_ortho_species.keys():
								if key == "S_pombe_keys":
									ppi = PPI((tuple([orf1]),tuple([orf2])),(tuple([ortho_PPI[0]]),tuple([ortho_PPI[1]])),str(PPIs_in_species[(orf1,orf2)].numberOfReports),str(PPIs_in_ortho_species[ortho_PPI].numberOfReports),"preserved")
									
								else:
									ppi = PPI((tuple([ortho_PPI[0]]),tuple([ortho_PPI[1]])),(tuple([orf1]),tuple([orf2])),str(PPIs_in_ortho_species[ortho_PPI].numberOfReports),str(PPIs_in_species[(orf1,orf2)].numberOfReports),"preserved")
									
								PPI_comparison["preserved"].append(ppi)

		# Removing duplicates
		PPI_comparison["lost"]=list(set([i for i in PPI_comparison["lost"]]))
		PPI_comparison["preserved"]=list(set([i for i in PPI_comparison["preserved"]]))
		PPI_comparison["rewired"]=list(set([i for i in PPI_comparison["rewired"]]))
		
		# Write summary file
		write_PPI_type_summary_files(PPI_comparison,ppi_type_summary_file)
		
	# Already done
	else:
		with open(ppi_type_summary_file) as f:
			PPI = namedtuple('PPI', ['SpomPPI', 'ScerPPI',"reportsInSpom","reportsInScer","ppiType"])
			for items in csv.reader(f, dialect="excel-tab" ):
				if (items[0]=="SpomPPI"): # Skip the header row
					continue
				SpomPPI,ScerPPI=literal_eval(items[0]),literal_eval(items[1])
				reportsInSpom,reportsInScer=str(items[2]),str(items[3])
				ppiType=items[4]
				ppi = PPI(SpomPPI,ScerPPI,reportsInSpom,reportsInScer,ppiType)
				
				PPI_comparison[ppiType].append(ppi)

	lost,preserved,rewired = len(PPI_comparison["lost"]),len(PPI_comparison["preserved"]),len(PPI_comparison["rewired"])	
	total = len([key for key in scer_PPIs.keys()])+len([key for key in spom_PPIs.keys()])
	duplicates= total - (lost+preserved+rewired+missingORFs)
	print(lost,"/",total,"PPIs are lost")
	print(preserved,"/",total,"PPIs are preserved")
	print(rewired,"/",total,"PPIs are rewired")
	print(missingORFs,"/",total,"PPIs were skipped (orfs not in orthologs dictionary)")
	print(duplicates,"/",total,"PPIs were skipped (duplicates)")
	print("Total :",lost,"+",preserved,"+",rewired,"+",missingORFs,"+",duplicates,"=",lost+preserved+rewired+missingORFs+duplicates,"/",total)
	
	return PPI_comparison

def write_PPI_type_summary_files(PPIs, opath):
	""" 
	Write summary files for the different PPI types in the analysis 
	- One file with all combined
	- One file for each of lost, preserved and rewired PPIs
	Args:
		PPIs (dict): dictionary of PPIs in the species. Keys= sorted tuples of interacting ORF names, Values= number of unique pubmed ids reporting the interaction
		opath (str): path to the desired location for the filtered output file (all PPIs)
	"""
	opath_lost = opath[:-4]+'_lost.txt'
	opath_preserved = opath[:-4]+'_preserved.txt'
	opath_rewired = opath[:-4]+'_rewired.txt'

	with open(opath, 'w') as fh0:
		with open(opath_lost, 'w') as fh1:
			with open(opath_preserved, 'w') as fh2:
				with open(opath_rewired, 'w') as fh3:
					header = '\t'.join(['SpomPPI','ScerPPI','reportsInSpom','reportsInScer','ppiType'])+'\n'
					fh0.write(header)
					fh1.write(header)
					fh2.write(header)
					fh3.write(header)

					for type in PPIs.keys():
						for ppi in PPIs[type]:
							fh0.write('\t'.join([str(ppi.SpomPPI),str(ppi.ScerPPI),str(ppi.reportsInSpom),str(ppi.reportsInScer),type])+'\n')
							if(type=="lost"): 			# Lost PPI
								fh1.write('\t'.join([str(ppi.SpomPPI),str(ppi.ScerPPI),str(ppi.reportsInSpom),str(ppi.reportsInScer),type])+'\n')
							elif(type=="preserved"):	# Preserved PPI
								fh2.write('\t'.join([str(ppi.SpomPPI),str(ppi.ScerPPI),str(ppi.reportsInSpom),str(ppi.reportsInScer),type])+'\n')
							else:						# Rewired PPI
								fh3.write('\t'.join([str(ppi.SpomPPI),str(ppi.ScerPPI),str(ppi.reportsInSpom),str(ppi.reportsInScer),type])+'\n')
	return

def readSummary(modeledPPIs_pombe,modeledPPIs_cere,ppi_type_summary_file):
	"""
		Creates the interologs dictionary with the key and ppi type
	"""

	print ("----------- PPI type summary -----------")
	scerPPIs = {}
	spomPPIs = {}
	

	# interologs = defaultdict(set)
	# Read in the PPI type file and write out only the PPIs with structural mapping
	for items in csv.reader(open(ppi_type_summary_file), dialect="excel-tab" ):
		# Skip header
		if (items[0]== "SpomPPI"): continue
			
		# Read in basic info on the PPI
		ppiType = items[4]
		

		# 1) The interaction is preserved, check if structure available for 2 PPIs (1 in S. cerevisiae, 1 in S. pombe)
		if ppiType == "preserved": 
				
			spom_ppi = (literal_eval(items[0])[0][0],literal_eval(items[0])[1][0])
			scer_ppi = (literal_eval(items[1])[0][0],literal_eval(items[1])[1][0])
			
			spomPPIs.update({tuple(sorted(spom_ppi)):"preserved"})	
			scerPPIs.update({tuple(sorted(scer_ppi)):"preserved"})		

		# 2) The interaction is lost, check if structure available for 1 PPI (the one that is not lost) and add 1 single protein
		elif ppiType == "lost":	

			spom_ppi_temp = (literal_eval(items[0])[0], literal_eval(items[0])[1])
			scer_ppi_temp = (literal_eval(items[1])[0],literal_eval(items[1])[1])

			
			# Ortholog lost in S. pombe
			if None in spom_ppi_temp:
				scer_ppi = (scer_ppi_temp[0][0]),(scer_ppi_temp[1][0])
				scerPPIs.update({tuple(sorted(scer_ppi)):"lost"})

				# if spom_ppi_temp[0] == None and spom_ppi_temp[1] == None:
				# 	scerPPIs.update({tuple(sorted(scer_ppi)):"lost"})

				# elif spom_ppi_temp[0] == None:
				# 	spom_ppis = [(None,singleProt) for singleProt in spom_ppi_temp[1]]
				# 	scerPPIs.setdefault("lost partner",[]).append(tuple(sorted(scer_ppi)))
				# 	modeledPPIs_cere[scer_ppi].add_ppiType('lost partner')  
				# else : #spom_ppi_temp[1] == None
				# 	spom_ppis = [(singleProt,None) for singleProt in spom_ppi_temp[0]]
				# 	scerPPIs.setdefault("lost partner",[]).append(tuple(sorted(scer_ppi)))
				# 	modeledPPIs_cere[scer_ppi].add_ppiType('lost partner')  
					
			# ortholog lost in S. cere
			else: # None in scer_ppi_temp:
				spom_ppi = (spom_ppi_temp[0][0]),(spom_ppi_temp[1][0]) 
				spomPPIs.update({tuple(sorted(spom_ppi)):"lost"})
				
				# if scer_ppi_temp[0] == None and scer_ppi_temp[1] == None:
				# 	scer_ppis = [(None,None)]
				# 	spomPPIs.setdefault("lost interolog",[]).append(tuple(sorted(spom_ppi)))
				# 	modeledPPIs_pombe[spom_ppi].add_ppiType('lost interolog')
				# elif scer_ppi_temp[0] == None:
				# 	scer_ppis = [(None,singleProt) for singleProt in scer_ppi_temp[1]]
				# 	spomPPIs.setdefault("lost partner",[]).append(tuple(sorted(spom_ppi)))
				# 	modeledPPIs_pombe[spom_ppi].add_ppiType('lost partner') 
				# else : #scer_ppi_temp[1] == None
				# 	scer_ppis = [(singleProt,None) for singleProt in scer_ppi_temp[0]]
				# 	spomPPIs.setdefault("lost partner",[]).append(tuple(sorted(spom_ppi)))
				# 	modeledPPIs_pombe[spom_ppi].add_ppiType('lost partner') 

				
				
		# 3) The interaction is rewired
		elif ppiType == "rewired":	
			spom_ppi = (literal_eval(items[0])[0][0],literal_eval(items[0])[1][0])
			scer_ppi = (literal_eval(items[1])[0][0],literal_eval(items[1])[1][0])

			keyspom = tuple(sorted(spom_ppi))
			keyscer = tuple(sorted(scer_ppi))

			if tuple(sorted(keyspom)) not in spomPPIs:
				spomPPIs.update({tuple(sorted(keyspom)):"rewired"})
			if tuple(sorted(keyscer)) not in scerPPIs:
				scerPPIs.update({tuple(sorted(keyscer)):"rewired"})
			

		else:
			continue
	

	PPIs_summary = {"Scer":scerPPIs, "Spom":spomPPIs}
	
	# print("S pombe")
	# print("lost partner",len([ppi for ppi in modeledPPIs_pombe.values() if ppi.ppiType == ["lost partner"]]))
	# print("lost interolog",len([ppi for ppi in modeledPPIs_pombe.values() if ppi.ppiType == ["lost interolog"]]))
	# print("preserved",len([ppi for ppi in modeledPPIs_pombe.values() if ppi.ppiType == ["preserved"]]))
	# print("rewired",len([ppi for ppi in modeledPPIs_pombe.values() if ppi.ppiType == ["rewired"]]))
	# # print("total=",len(spomPPIs['lost partner'])+len(spomPPIs['lost interolog'])+len(spomPPIs['preserved'])+len(spomPPIs['rewired']))
	
	# print("\nS cere")
	# print("lost partner",len([ppi for ppi in modeledPPIs_cere.values() if ppi.ppiType == ["lost partner"]]))
	# print("lost interolog",len([ppi for ppi in modeledPPIs_cere.values() if ppi.ppiType == ["lost interolog"]]))
	# print("preserved",len([ppi for ppi in modeledPPIs_cere.values() if ppi.ppiType == ["preserved"]]))
	# print("rewired",len([ppi for ppi in modeledPPIs_cere.values() if ppi.ppiType == ["rewired"]]))
	# # print("total=",len(scerPPIs['lost partner'])+len(scerPPIs['lost interolog'])+len(scerPPIs['preserved'])+len(scerPPIs['rewired']))
	
	return PPIs_summary

def add_to_ppi_models(modeledPPIs_pombe,modeledPPIs_cere,pickled_path_cere_new,pickled_path_pombe_new,PPIs_summary):
	# ----
	# update models to add PPiType
	
	new_modeledPPIs_cere={}
	for ppi in modeledPPIs_cere.values():

		key= tuple(sorted([ppi.orf1,ppi.orf2]))
		
		newPPI=PPI(ppi.orf1,ppi.orf2,ppi.numberOfReports)
		newPPI.add_ppiType(PPIs_summary["Scer"][key])

		for prot in ppi.proteins:
			newPPI.add_protein(prot)
		new_modeledPPIs_cere.update({key:newPPI})
	pickle.dump(new_modeledPPIs_cere, open(pickled_path_cere_new, "wb"))

	# update models to add PPiType
	# pombe
	new_modeledPPIs_pombe={}
	for ppi in modeledPPIs_pombe.values():
		key= tuple(sorted([ppi.orf1,ppi.orf2]))
		
		newPPI=PPI(ppi.orf1,ppi.orf2,ppi.numberOfReports)
		newPPI.add_ppiType(PPIs_summary["Spom"][key])
		for prot in ppi.proteins:
			newPPI.add_protein(prot)
		new_modeledPPIs_pombe.update({key:newPPI})
	pickle.dump(new_modeledPPIs_pombe, open(pickled_path_pombe_new, "wb"))


	return

def main():
	pombe_data = '../data/pombe'
	cere_data = '../data/cerevisiae'
	orthologs_data='../data/orthologs'
	ppi_type_data='../data/PPI_type'
	# blastToPDBData = os.path.join(externalData, 'blastToPDB')

	print ("-------- Loading PPI models --------")
	pickled_path_pombe =os.path.join(pombe_data,'modeled_ppis.pkl')
	summary_path_pombe =os.path.join(pombe_data,'model_summary.txt')

	pickled_path_cere =os.path.join(cere_data,'modeled_ppis.pkl')
	summary_path_cere =os.path.join(cere_data,'model_summary.txt')

	modeledPPIs_pombe = pickle.load(open(pickled_path_pombe,"rb"))
	modeledPPIs_cere = pickle.load(open(pickled_path_cere,"rb"))
	print ("-------- PPI models loaded --------")


	# 1) Get orthologues between S.cere and S.pombe from orthologs file (manually curated from https://www.pombase.org/data/orthologs/cerevisiae-orthologs.txt)
	path = os.path.join(orthologs_data,"orthologs.txt")
	orthologs_manual= manually_curated_orthologs(path,modeledPPIs_pombe,modeledPPIs_cere)
	
	# 2) Get orthologues between S.cere and S.pombe from BLAST reciprocal best hits
	path = os.path.join(orthologs_data,"RBH/")
	orthologs_BLAST=reciprocal_BLAST_orthologs(modeledPPIs_pombe,modeledPPIs_cere,path)

	# 3) Get ORFs in S.cere and S.pombe that have no homologs in the other species (no homolog = no ortholog)
	path = os.path.join(orthologs_data,"homologs/")
	no_homologues = find_homologs_BLAST(modeledPPIs_pombe,modeledPPIs_cere,path)
	
	# 4) Combine the two ortholog mapping dictionaries, 
	orthologs = get_ortholog_mapping(modeledPPIs_pombe,modeledPPIs_cere,orthologs_manual,orthologs_BLAST,no_homologues)
	
	# 5) Set ORFs with no ortholog mapping (either manual or reciprocal BLAST) to ortholog == NONE
	orthologs = edit_ortholog_mapping(modeledPPIs_pombe,modeledPPIs_cere,orthologs)
	
	# 6) Print summary file of orthology mapping
	path_spom = os.path.join(orthologs_data,"orthologs_summary_spom.txt")
	path_scer = os.path.join(orthologs_data,"orthologs_summary_scer.txt")
	print_summary(orthologs,path_spom,path_scer)


	# 6) Get preserved, lost, rewired PPIs between the two species
	ppi_type_summary_file = os.path.join(ppi_type_data, 'PPI_type.txt')
	PPI_comparison= compare_PPIs_between_species(modeledPPIs_pombe,modeledPPIs_cere,orthologs,ppi_type_summary_file)

	
	
	ppi_type_summary_file = os.path.join(ppi_type_data, 'PPI_type.txt')
	PPIs_summary = readSummary(modeledPPIs_pombe,modeledPPIs_cere,ppi_type_summary_file)
	
	# 7) add to models?


	pickled_path_cere_new =os.path.join(cere_data,'modeled_ppis_new.pkl')
	pickled_path_pombe_new =os.path.join(pombe_data,'modeled_ppis_new.pkl')
	add_to_ppi_models(modeledPPIs_pombe,modeledPPIs_cere,pickled_path_cere_new,pickled_path_pombe_new,PPIs_summary)
	
if __name__ == '__main__':
	main()

