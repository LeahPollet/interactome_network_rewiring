"""
Script to read and format results from the analysis for plotting

Author Leah Pollet
leah.pollet@mail.mcgill.ca
"""
import os,re, sys
import pandas as pd
import pyarrow
import pyarrow.feather as feather
import pickle
from modeled_PPIs import PPI,Protein,Interface,Residue
from tqdm import tqdm
import csv
import numpy as np
from scipy import optimize as opt
import numpy as np
# from collections import Counter
#feather errors useful link: https://bobbyhadz.com/blog/python-no-module-named-pyarrow



def bin_residues(modeledPPIs, bin):
	"""
	Returns a list of all residues in a model according to a gived binning category
	Args:
		model (str): path to the pickled model
		bin (str): interf, noninterf, interf_preserved, interf_lost,interf_rewired,noninterf_preserved,noninterf_lost,noninterf_rewired,buried,exposed,buried_preserved,exposed_preserved,buried_lost,exposed_lost,buried_rewired, exposed_rewired


	"""
	
	if bin == "interf":
		# Interfacial, all residues
		return [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == True]
	elif bin == "noninterf":
		# Non-interfacial, all residues
		return [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == False]
	elif bin == "interf_preserved":
		# Interfacial residues in preserved PPIs
		return [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == True if ppi.ppiType == ["preserved"]]
	elif bin == "interf_lost":
		# Interfacial residues in lost PPIs
		return [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == True if ppi.ppiType == ["lost"]]
	elif bin == "interf_rewired":	
		# Interfacial residues in rewired PPIs
		return [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == True if ppi.ppiType == ["rewired"]]
	elif bin == "noninterf_preserved":
		# Non-interfacial residues in preserved PPIs
		return [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == False if ppi.ppiType == ["preserved"]]
	elif bin == "noninterf_lost":	
		# Non-interfacial residues in lost PPIs
		return [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == False if ppi.ppiType == ["lost"]]
	elif bin == "noninterf_rewired":	
		# Non-interfacial residues in rewired PPIs
		return [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == False if ppi.ppiType == ["rewired"]]
	elif bin == "buried":
		# Buried, non-interfacial residues
		return [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == False if res.rsaComplex<0.25]
	elif bin == "exposed":
		# Exposed, non-interfacial residues
		return [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == False if res.rsaComplex>=0.25]
	elif bin == "buried_preserved":
		# Buried, non-interfacial residues in preserved PPIs
		return [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == False if res.rsaComplex<0.25 if ppi.ppiType == ["preserved"]]
	elif bin == "exposed_preserved":
		# Exposed, non-interfacial residues in preserved PPIs
		return [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == False if res.rsaComplex>=0.25 if ppi.ppiType == ["preserved"]]
	elif bin == "buried_lost":
		# Buried, non-interfacial residues in lost PPIs
		return [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == False if res.rsaComplex<0.25 if ppi.ppiType == ["lost"]]
	elif bin == "exposed_lost":
		# Exposed, non-interfacial residues in lost PPIs
		return [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == False if res.rsaComplex>=0.25 if ppi.ppiType == ["lost"]]
	elif bin == "buried_rewired":
		# Buried, non-interfacial residues in rewired PPIs
		return [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == False if res.rsaComplex<0.25 if ppi.ppiType == ["rewired"]]
	elif bin == "exposed_rewired":
		# Exposed, non-interfacial residues in rewired PPIs
		return [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == False if res.rsaComplex>=0.25 if ppi.ppiType == ["rewired"]]
	else :
		print("Error, wrong bin type")
		return []

def write_summary_file(model,summary_file):

	with open(model, 'rb') as f:
		modeledPPIs = pickle.load(f)


	df = pd.DataFrame(columns=['consurfRate4Site_score', 'consurfRate4Site_std', 'consurfRate4Site_resNum','consurfDB_score', 
								   'consurfDB_std','consurfDB_resNum'])

	bins = ["interf", "noninterf", "interf_preserved", "interf_lost","interf_rewired","noninterf_preserved","noninterf_lost","noninterf_rewired","buried","exposed","buried_preserved","exposed_preserved","buried_lost","exposed_lost","buried_rewired", "exposed_rewired"]
		
	for bin in bins:
		residues = bin_residues(modeledPPIs,bin)
		row ={}
		row["consurfRate4Site_score"] = np.nanmean([r.consurf_score for r in residues])
		row["consurfRate4Site_std"] = np.sqrt(np.nanmean([pow(r.consurf_std,2) for r in residues])) # average the variance and take sqrt
		row["consurfRate4Site_resNum"] = sum(~np.isnan([r.consurf_score for r in residues]))
		row["consurfDB_score"] = np.nanmean([r.consurfDB_score for r in residues])
		row["consurfDB_std"] = np.nanstd([r.consurfDB_score for r in residues])
		row["consurfDB_resNum"] = sum(~np.isnan([r.consurfDB_score for r in residues]))

		row = pd.DataFrame([row])
		df = pd.concat([df, row],axis = 0,ignore_index=True)
		

	
	df.index = bins
	print(df.to_string())
	
	with open(summary_file, 'w') as f:
	    df_string = df.to_string()
	    f.write(df_string)

	return


def main():

	# 1) Set up
	SPOM_PPI_MODELS_PATH = '../data/pombe/modeled_ppis_new.pkl'
	SCER_PPI_MODELS_PATH = '../data/cerevisiae/modeled_ppis_new.pkl'
	plotDir= "../plots"
	summary_file_spom = os.path.join(plotDir, 'spom_summary.txt')
	summary_file_scer =os.path.join(plotDir, 'scer_summary.txt')

	# Write summary for S.pombe
	write_summary_file(SPOM_PPI_MODELS_PATH,summary_file_spom)

	# Write summary for S.cere
	write_summary_file(SCER_PPI_MODELS_PATH,summary_file_scer)
	
	
	
	return


if __name__ == '__main__':
	main()







