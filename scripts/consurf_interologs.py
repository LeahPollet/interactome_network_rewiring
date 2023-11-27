"""
Manually compute consurf scores using the rate4site program

"""
import os, sys, stat
from modeled_PPIs import PPI,Protein,Interface,Residue
import pickle
from collections import Counter
import re
from tqdm import tqdm
import zipfile
import subprocess
from shutil import copyfile
from Bio import SeqIO
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
sns.set(style="whitegrid")

def download_rate4sites(): 
	# useful if rate4site is not compling properly:
	# https://www.tau.ac.il/~itaymay/cp/rate4site.html
	oPath = "../data"
	url = "https://www.tau.ac.il/~itaymay/cp/rate4site.3.2.source.zip"

	# Dowload the source code
	zipped = oPath+"/"+os.path.basename(url)
	path_to_rate4site = oPath+"/"+os.path.basename(url)[:-4]
	if not os.path.exists(path_to_rate4site):
		print("Downloading rate4site")
		os.system('wget -P '+ oPath +' '+url)
		with zipfile.ZipFile(zipped, 'r') as zip_ref:
			zip_ref.extractall(oPath)
		os.remove(zipped)

	# Compiling the source code
	path_to_rate4site_exec = oPath+"/"+os.path.basename(url)[:-4]+"/sourceMar09/rate4site"
	if not os.path.exists(path_to_rate4site_exec):	
		path = path_to_rate4site_exec[:-10]
		subprocess.call('make', shell=True,cwd=path)
	
	print(path_to_rate4site_exec)
	return path_to_rate4site+"/sourceMar09"

def set_up_for_rate4sites(model_path,alignment_path,IN_DIR,OUT_DIR,treefile,txt): 
	'''Create the consurf directory structure and the input files for consurf analysis (MSA and phylo tree for each protein)
	Returns: 
		inputs (list): list of ORF names with input available for running rate4site
		proteins (dict): PPI models from the build_protein_models.py script
	'''
	print("------ Set up for rate4Site ------ ")
	# Set up the consurf directory
	
	if not os.path.exists(OUT_DIR):
		os.makedirs(OUT_DIR)
	if not os.path.exists(IN_DIR):
		os.makedirs(IN_DIR)
	
	# Load protein models
	if not os.path.exists(model_path):
		raise UserWarning('Protein models not created')
	with open(model_path, 'rb') as f:
		model = pickle.load(f)
		print ("Loaded protein models")

	

	# Create 1 input dir per ORF in the model
	# Uses same tree and alignment file  as in dnds calculations
	# Alignment file seq.id format slightly modifed to match the tree format
	
	inputs = list(set([protein.orf for ppi in model.values() for protein in ppi.proteins]))
	proteins = [protein for ppi in model.values() for protein in ppi.proteins]
	
	for protein in tqdm(proteins):
		ORF_folder = IN_DIR+"/"+protein.orf
		if not os.path.exists(ORF_folder):
			os.makedirs(ORF_folder)

			# Copy the tree file into the new directory
			copyfile(treefile, ORF_folder+"/yeasts.tree")
			alignment_file_original = alignment_path+"/"+protein.orf+"/"+protein.orf+".aa.aln"
			
			# Copy the alignment file into the new directory, modifying its sequence names to match the tree file
			alignment_file_corrected = ORF_folder+"/"+protein.orf+".aln"
			with open(alignment_file_original) as original, open(alignment_file_corrected, 'w') as corrected:
					records = SeqIO.parse(original, 'fasta')
					for record in records:
						record.id = record.description.split()[1][:4]
						record.description = record.id 
						SeqIO.write(record, corrected, 'fasta')
	return inputs, model

def run_rate4sites(path_to_rate4site,input_list,in_dir,out_dir):
	print("-------- Run rate4Site -------- ")
	
	for orf in tqdm(input_list):
		
		tree_file = in_dir+"/"+orf+"/yeasts.tree"
		aln_file = in_dir+"/"+orf+"/"+orf+".aln"
		outfile = out_dir+"/"+orf+".consurf.txt"

		if not os.path.exists(outfile):
			
			cmd = "./rate4site -s " + aln_file+" -t "+tree_file+" -o "+outfile
			# subprocess.Popen(["./autogen.sh"], stdout=subprocess.PIPE, cwd=path_to_dssp)
			# stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,cwd=path_to_rate4site,
			subprocess.call("../data/rate4site.3.2.source/sourceMar09/rate4site -s "+aln_file+" -t "+tree_file+" -o "+outfile,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,shell=True)
			# subprocess.call(cmd, shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,cwd=path_to_rate4site)
			
	return

def add_to_protein_models(inputs,modeledPPIs,OUTFILE_DIR_pombe,model_path):
	
	'''
	Add the consurf results to proteins in our protein models and prints a log file
	Args:
		inputs (list): list of yeast ORFs with ConSurf scores computed
		modeledPPIs (dict):  PPI models from the build_protein_models.py script
		TODO, add logfile?
	'''
	print("-------- Add to models -------- ")
	# --- For all proteins modeled
	noConsurfScore = 0
	proteins = [protein for ppi in modeledPPIs.values() for protein in ppi.proteins]
	
	for protein in tqdm(proteins):
		if protein.orf in inputs:
			# print(">" ,protein.orf )
			consurf_file = OUTFILE_DIR_pombe+"/"+protein.orf+".consurf.txt"
			consurf_summary = get_results_from_consurf_file(consurf_file)
			for res in protein.residues:
				if (res.indexInSeq in consurf_summary):
					if (res.aa.letter != consurf_summary[res.indexInSeq][0]):
						print("Problem, Sequence mismatch in the mapping of Consurf scores")

					res.consurf_score = consurf_summary[res.indexInSeq][1]
					res.conf_int_consurf_score=consurf_summary[res.indexInSeq][2]
					res.consurf_std =consurf_summary[res.indexInSeq][3]
					res.conf_MSA_data =consurf_summary[res.indexInSeq][4]
				else:
					noConsurfScore += 1

		# --- Check counts in the models:
	
	total_res 		= len([r for ppi in modeledPPIs.values() for p in ppi.proteins for r in p.residues])
	score_for_res 	= len([r for ppi in modeledPPIs.values() for p in ppi.proteins for r in p.residues if not np.isnan(r.consurf_score)]) 
	print(score_for_res,"/",total_res," residues in the PPI models have a rate4site consurf score.")
	
	# --- Save the new protein models updated with consurf (rate4site) scores
	pickle.dump(modeledPPIs, open(model_path, 'wb'))

	return

def get_results_from_consurf_file(consurf_file):
	# Read out the results from a consurf file to a dictionary

	consurf_summary = {} # Dictionary to save the consurf results

	with open(consurf_file, 'r') as h:
		rows = list((line.split() for line in h))
		rows = rows[11:-2] # skip the header and footer

		for row in rows:
			if len(row )== 5:
				row=[row[0],row[1],row[2],row[3].split("]")[0]+"]",row[3].split("]")[1],row[4]]
				
			index = float(row[0])-1
			residue = row[1]
			score = float(row[2])
			confidence_interval= eval(row[3])
			std = float(row[4])
			data = eval(row[5])


			consurf_summary.update({index: (residue,score,confidence_interval,std,data)})
					
	return consurf_summary

# def statistical_tests(a,b,score,label1,label2):

# 	print("> Test for difference in",score, "between",label1,"and",label2)
# 	difference = round(np.mean(a)-np.mean(b),3)
# 	print("\tavg(",label1,") - avg(",label2,")=",difference)
# 	t1, p1 = stats.ttest_ind(a, b)
# 	t2, p2 = stats.mannwhitneyu(a, b)
# 	print("\tT-test: p-value =",p1)
# 	print("\tMannWhitneYu-test: p-value =",p2)

# 	significance = ""
# 	if (p1<0.05 and p2<0.05):
# 		significance = "*"
# 	if (p1<0.001 and p2<0.001):
# 		significance = "**"

# 	if(significance =="*" or significance =="**"):
# 		if (difference>0):
# 			print(label2,"are more conserved than",label1,". Significance =",significance)
# 		elif (difference<0):
# 			print(label1,"are more conserved than",label2,". Significance =",significance)
# 		else:
# 			print("Error")

# 	return significance

# def boxplots(data,title,ylabel,xlabels,colors,legend,saveFig,scale):
	
# 	# Figure basics
# 	fig,ax = plt.subplots()
# 	ax.set_title(title)
# 	ax.set_ylabel(ylabel)
	
# 	box = ax.boxplot(data,sym="",patch_artist=True,showmeans=True,meanprops={"markerfacecolor":"darkred", "markeredgecolor":"darkred"})
# 	plt.xticks([1.5,3.5,5.5],xlabels)
# 	for patch, color in zip(box['boxes'], colors):
# 	    patch.set_facecolor(color)

# 	# Statistical annotations
# 	signif = statistical_tests(data[0],data[1],ylabel,legend[0]+" residues in "+xlabels[0],legend[1]+" residues in "+xlabels[0])
# 	x1, x2 = 1, 2   
# 	y, h, col = scale, 0, 'k'
# 	plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
# 	plt.text((x1+x2)*.5, y+h, signif, ha='center', va='bottom', color=col)

# 	signif = statistical_tests(data[2],data[3],ylabel,legend[0]+" residues in "+xlabels[1],legend[1]+" residues in "+xlabels[1])
# 	x1, x2 = 3, 4   
# 	y, h, col = scale, 0, 'k'
# 	plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
# 	plt.text((x1+x2)*.5, y+h, signif, ha='center', va='bottom', color=col)

# 	signif = statistical_tests(data[4],data[5],ylabel,legend[0]+" residues in "+xlabels[2],legend[1]+" residues in "+xlabels[2])
# 	x1, x2 = 5, 6   
# 	y, h, col = scale, 0, 'k'
# 	plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
# 	plt.text((x1+x2)*.5, y+h, signif, ha='center', va='bottom', color=col)

# 	# Legend
# 	ax.legend([box["boxes"][0], box["boxes"][1]], [legend[0], legend[1]], loc=(0.87,1.02))

# 	fig.set_size_inches(6, 6)
# 	plt.tight_layout()
# 	plt.savefig(saveFig,dpi=200)

# 	return
	

# def test_plots(modeledPPIs):
# 	# Start by plotting even with duplicates
# 	# If not as expected, we can start thinking about removing proteins that are modeled multiple times, using the same PDB file, using the same chain, etc.

# 	if not os.path.exists(PLOT_DIR):
# 		os.makedirs(PLOT_DIR)

# 	# # Print general summary 
# 	# total = [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True]
# 	# interfacial = [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if res.interfacial == True]
# 	# consurfDBscore = [res for res in total if not np.isnan(res.consurfDB_score)]
# 	# consurfDBscore_interf = [res for res in interfacial if not np.isnan(res.consurfDB_score)]
# 	# consurfscore = [res for res in total if not np.isnan(res.consurf_score)]
# 	# consurfscore_interf = [res for res in interfacial if not np.isnan(res.consurf_score)]
	
# 	# print("Total number of res in models:",len(total))
# 	# print("Total number of interfacial res in models:",len(interfacial))
# 	# print(len(consurfDBscore),"/",len(total),"=",round(len(consurfDBscore)/len(total)*100,2),"%","of all residues have a ConsurfDB score")
# 	# print(len(consurfDBscore_interf),"/",len(interfacial),"=",round(len(consurfDBscore_interf)/len(interfacial)*100,2),"%","of all interfacial residues have a ConsurfDB score")
# 	# print(len(consurfscore),"/",len(total),"=",round(len(consurfscore)/len(total)*100,2),"%","of all residues have a Consurf rate4site score")
# 	# print(len(consurfscore_interf),"/",len(interfacial),"=",round(len(consurfscore_interf)/len(interfacial)*100,2),"%","of all interfacial residues have a Consurf rate4site score")

# 	# # Print Consurf scores summary
# 	# consurf = [res.consurf_score for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True]
# 	# consurfDB = [res.consurfDB_score for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True]
# 	# combined = [(res.consurf_score,res.consurfDB_score) for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True]
# 	# consurfs = pd.DataFrame({'Consurf DB': consurfDB,'Consurf rate4site': consurf})
# 	# consurfs = consurfs[consurfs['Consurf DB'].notna()]
# 	# consurfs = consurfs[consurfs['Consurf rate4site'].notna()]
	
# 	# print(len([score for score in consurf if not np.isnan(score)]),"/",len(consurf),"residue have consurf (rate4site) score calculated")
# 	# print(len([score for score in consurfDB if not np.isnan(score)]),"/",len(consurfDB),"residue have consurf (DB) score calculated")
# 	# print(len(consurfs),"/",len(consurfDB),"residue have both consurf (DB) score and consurf (rate4site) score calculated")
# 	# print("Pearson correlation:",consurfs["Consurf DB"].corr(consurfs["Consurf rate4site"]))
# 	# print("\n\n")

# 	# Boxplots
# 	# -------------------
# 	# 1) INTERFACIAL VS NON INTERFACIAL, allPPIs/ScerPPis/SpomPPIs
# 	print("--- interfacialResidues vs nonInterfacialResidues ---")
# 	interf= [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if res.interfacial == True]
# 	noninterf =[res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if res.interfacial == False]
# 	interfScer= [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if res.interfacial == True if ppi.species=="Scer"]
# 	noninterfScer =[res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if res.interfacial == False if ppi.species=="Scer"]
# 	interfSpom= [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if res.interfacial == True if ppi.species=="Spom"]
# 	noninterfSpom =[res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if res.interfacial == False if ppi.species=="Spom"]
	
# 	# Consurf DB
# 	a=[res.consurfDB_score for res in interf if not np.isnan(res.consurfDB_score)]
# 	b=[res.consurfDB_score for res in noninterf if not np.isnan(res.consurfDB_score)]
# 	c=[res.consurfDB_score for res in interfScer if not np.isnan(res.consurfDB_score)]
# 	d=[res.consurfDB_score for res in noninterfScer if not np.isnan(res.consurfDB_score)]
# 	e=[res.consurfDB_score for res in interfSpom if not np.isnan(res.consurfDB_score)]
# 	f=[res.consurfDB_score for res in noninterfSpom if not np.isnan(res.consurfDB_score)]
	
# 	# print("Scer PPIs")
# 	# print("Interfacial: avgConSurfDB:",np.nanmean(c),"resNum:",len(c),"std:",np.nanstd(c))
# 	# print("Non-Interfacial: avgConSurfDB:",np.nanmean(d),"resNum:",len(d),"std:",np.nanstd(d))
# 	# print("Spom PPIs")
# 	# print("Interfacial: avgConSurfDB:",np.nanmean(e),"resNum:",len(e),"std:",np.nanstd(e))
# 	# print("Non-Interfacial: avgConSurfDB:",np.nanmean(f),"resNum:",len(f),"std:",np.nanstd(f))

# 	xlabels = ["All PPIs","Scer PPIs","Spom PPIs"]
# 	colors = ['sandybrown', 'palegreen', 'sandybrown', 'palegreen','sandybrown', 'palegreen']
# 	legend=["Interfacial","Non-Interfacial"]
# 	boxplots([a,b,c,d,e,f],"Interfacial vs. Non-Interfacial residues","ConSurf score (DB)",xlabels,colors,legend,PLOT_DIR+"/interf_vs_nonInterf_consurfDB",3)

# 	# Consurf rate4site
# 	a=[res.consurf_score for res in interf if not np.isnan(res.consurf_score)]
# 	b=[res.consurf_score for res in noninterf if not np.isnan(res.consurf_score)]
# 	c=[res.consurf_score for res in interfScer if not np.isnan(res.consurf_score)]
# 	d=[res.consurf_score for res in noninterfScer if not np.isnan(res.consurf_score)]
# 	e=[res.consurf_score for res in interfSpom if not np.isnan(res.consurf_score)]
# 	f=[res.consurf_score for res in noninterfSpom if not np.isnan(res.consurf_score)]
	
# 	# print("Scer PPIs")
# 	# print("Interfacial: avgConSurfrate4site:",np.nanmean(c),"resNum:",len(c),"std:",np.nanstd(c))
# 	# print("Non-Interfacial: avgConSurfrate4site:",np.nanmean(d),"resNum:",len(d),"std:",np.nanstd(d))
# 	# print("Spom PPIs")
# 	# print("Interfacial: avgConSurfrate4site:",np.nanmean(e),"resNum:",len(e),"std:",np.nanstd(e))
# 	# print("Non-Interfacial: avgConSurfrate4site:",np.nanmean(f),"resNum:",len(f),"std:",np.nanstd(f))

# 	xlabels = ["All PPIs","Scer PPIs","Spom PPIs"]
# 	colors = ['sandybrown', 'palegreen', 'sandybrown', 'palegreen','sandybrown', 'palegreen']
# 	legend=["Interfacial","Non-Interfacial"]
# 	boxplots([a,b,c,d,e,f],"Interfacial vs. Non-Interfacial residues","ConSurf score (rate4site)",xlabels,colors,legend,PLOT_DIR+"/interf_vs_nonInterf_consurfrate4site",0.5)
# 	print("\n")
	
# 	# -------------------
# 	# 2) INTERFACIAL VS NON INTERFACIAL, allPPIs/lostPPIs/preservedPPIs
# 	print("--- Interfacial vs non interfacial in lost and preserved PPIs ---")
# 	interf= [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if res.interfacial == True]
# 	noninterf =[res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if res.interfacial == False]
# 	interfacial_lost= [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if ppi.type == 'lost' if res.interfacial == True]
# 	non_interfacial_lost=[res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if ppi.type == 'lost' if res.interfacial == False]
# 	interfacial_preserved =[res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if ppi.type == 'preserved' if res.interfacial == True]
# 	non_interfacial_preserved = [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if ppi.type == 'preserved'  if res.interfacial == False]
	
# 	# Consurf DB
# 	a=[res.consurfDB_score for res in interf if not np.isnan(res.consurfDB_score)]
# 	b=[res.consurfDB_score for res in noninterf if not np.isnan(res.consurfDB_score)]
# 	c=[res.consurfDB_score for res in interfacial_lost if not np.isnan(res.consurfDB_score)]
# 	d=[res.consurfDB_score for res in non_interfacial_lost if not np.isnan(res.consurfDB_score)]
# 	e=[res.consurfDB_score for res in interfacial_preserved if not np.isnan(res.consurfDB_score)]
# 	f=[res.consurfDB_score for res in non_interfacial_preserved if not np.isnan(res.consurfDB_score)]
# 	xlabels = ["All PPIs","Lost PPIs","Preserved PPIs"]
# 	colors = ['sandybrown', 'palegreen', 'sandybrown', 'palegreen','sandybrown', 'palegreen']
# 	legend=["Interfacial","Non-Interfacial"]
# 	boxplots([a,b,c,d,e,f],"Interfacial vs. non-interfacial residues","ConSurf score (DB)",xlabels,colors,legend,PLOT_DIR+"/interf_vs_nonInterf_LostPreservedconsurfDB",3)

# 	# Consurf rate4site
# 	a=[res.consurf_score for res in interf if not np.isnan(res.consurf_score)]
# 	b=[res.consurf_score for res in noninterf if not np.isnan(res.consurf_score)]
# 	c=[res.consurf_score for res in interfacial_lost if not np.isnan(res.consurf_score)]
# 	d=[res.consurf_score for res in non_interfacial_lost if not np.isnan(res.consurf_score)]
# 	e=[res.consurf_score for res in interfacial_preserved if not np.isnan(res.consurf_score)]
# 	f=[res.consurf_score for res in non_interfacial_preserved if not np.isnan(res.consurf_score)]
# 	xlabels = ["All PPIs","Lost PPIs","Preserved PPIs"]
# 	colors = ['sandybrown', 'palegreen', 'sandybrown', 'palegreen','sandybrown', 'palegreen']
# 	legend=["Interfacial","Non-Interfacial"]
# 	boxplots([a,b,c,d,e,f],"Interfacial vs. non-interfacial residues","ConSurf score (rate4site)",xlabels,colors,legend,PLOT_DIR+"/interf_vs_nonInterf_LostPreservedconsurfrate4site",0.5)
# 	print("\n")

# 	# -------------------
# 	# 3) LOST VS PRESERVED, allPPIs/Scer PPIs/Spom PPIs
# 	print("--- lost vs preserved ---")
# 	preserved= [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if ppi.type == 'preserved']
# 	lost =[res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if ppi.type == 'lost']
# 	preservedScer= [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if ppi.type == 'preserved' if ppi.species=="Scer"]
# 	lostScer =[res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if ppi.type == 'lost' if ppi.species=="Scer"]
# 	preservedSpom= [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if ppi.type == 'preserved' if ppi.species=="Spom"]
# 	lostSpom =[res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if ppi.type == 'lost' if ppi.species=="Spom"]
	
# 	# Consurf DB
# 	a=[res.consurfDB_score for res in preserved if not np.isnan(res.consurfDB_score)]
# 	b=[res.consurfDB_score for res in lost if not np.isnan(res.consurfDB_score)]
# 	c=[res.consurfDB_score for res in preservedScer if not np.isnan(res.consurfDB_score)]
# 	d=[res.consurfDB_score for res in lostScer if not np.isnan(res.consurfDB_score)]
# 	e=[res.consurfDB_score for res in preservedSpom if not np.isnan(res.consurfDB_score)]
# 	f=[res.consurfDB_score for res in lostSpom if not np.isnan(res.consurfDB_score)]
	
# 	# print("Scer PPIs")
# 	# print("Preserved: avgConSurfDB:",np.nanmean(c),"resNum:",len(c),"std:",np.nanstd(c))
# 	# print("Lost: avgConSurfDB:",np.nanmean(d),"resNum:",len(d),"std:",np.nanstd(d))
# 	# print("Spom PPIs")
# 	# print("Preserved: avgConSurfDB:",np.nanmean(e),"resNum:",len(e),"std:",np.nanstd(e))
# 	# print("Lost: avgConSurfDB:",np.nanmean(f),"resNum:",len(f),"std:",np.nanstd(f))

# 	xlabels = ["All PPIs","Scer PPIs","Spom PPIs"]
# 	colors = ['powderblue', 'lightcoral', 'powderblue', 'lightcoral','powderblue', 'lightcoral']
# 	legend=["Preserved","Lost"]
# 	boxplots([a,b,c,d,e,f],"Preserved vs. lost residues","ConSurf score (DB)",xlabels,colors,legend,PLOT_DIR+"/preserved_vs_lost_consurfDB",3)

# 	# Consurf rate4site
# 	a=[res.consurf_score for res in preserved if not np.isnan(res.consurf_score)]
# 	b=[res.consurf_score for res in lost if not np.isnan(res.consurf_score)]
# 	c=[res.consurf_score for res in preservedScer if not np.isnan(res.consurf_score)]
# 	d=[res.consurf_score for res in lostScer if not np.isnan(res.consurf_score)]
# 	e=[res.consurf_score for res in preservedSpom if not np.isnan(res.consurf_score)]
# 	f=[res.consurf_score for res in lostSpom if not np.isnan(res.consurf_score)]
	
# 	# print("\nScer PPIs")
# 	# print("Preserved: avgConSurfrate4site:",np.nanmean(c),"resNum:",len(c),"std:",np.nanstd(c))
# 	# print("Lost: avgConSurfrate4site:",np.nanmean(d),"resNum:",len(d),"std:",np.nanstd(d))
# 	# print("Spom PPIs")
# 	# print("Preserved: avgConSurfrate4site:",np.nanmean(e),"resNum:",len(e),"std:",np.nanstd(e))
# 	# print("Lost: avgConSurfrate4site:",np.nanmean(f),"resNum:",len(f),"std:",np.nanstd(f))

# 	xlabels = ["All PPIs","Scer PPIs","Spom PPIs"]
# 	colors = ['powderblue', 'lightcoral', 'powderblue', 'lightcoral','powderblue', 'lightcoral']
# 	legend=["Preserved","Lost"]
# 	boxplots([a,b,c,d,e,f],"Preserved vs. lost residues","ConSurf score (rate4site)",xlabels,colors,legend,PLOT_DIR+"/preserved_vs_lost_consurfrate4site",1)
# 	print("\n")
	
# 	# -------------------
# 	# 3) LOST VS PRESERVED INTERFACIAL, allPPIs/Scer PPIs/Spom PPIs
# 	print("--- lost vs preserved interfacial ---")
# 	preserved= [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if ppi.type == 'preserved' if res.interfacial == True]
# 	lost =[res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if ppi.type == 'lost' if res.interfacial == True]
# 	preservedScer= [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if ppi.type == 'preserved' if ppi.species=="Scer" if res.interfacial == True]
# 	lostScer =[res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if ppi.type == 'lost' if ppi.species=="Scer" if res.interfacial == True]
# 	preservedSpom= [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if ppi.type == 'preserved' if ppi.species=="Spom" if res.interfacial == True]
# 	lostSpom =[res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if ppi.multipleReports==True if ppi.type == 'lost' if ppi.species=="Spom" if res.interfacial == True]
	
# 	# Consurf DB
# 	a=[res.consurfDB_score for res in preserved if not np.isnan(res.consurfDB_score)]
# 	b=[res.consurfDB_score for res in lost if not np.isnan(res.consurfDB_score)]
# 	c=[res.consurfDB_score for res in preservedScer if not np.isnan(res.consurfDB_score)]
# 	d=[res.consurfDB_score for res in lostScer if not np.isnan(res.consurfDB_score)]
# 	e=[res.consurfDB_score for res in preservedSpom if not np.isnan(res.consurfDB_score)]
# 	f=[res.consurfDB_score for res in lostSpom if not np.isnan(res.consurfDB_score)]
# 	# print("Scer PPIs")
# 	# print("Preserved: avgConSurfDB:",np.nanmean(c),"resNum:",len(c),"std:",np.nanstd(c))
# 	# print("Lost: avgConSurfDB:",np.nanmean(d),"resNum:",len(d),"std:",np.nanstd(d))
# 	# print("Spom PPIs")
# 	# print("Preserved: avgConSurfDB:",np.nanmean(e),"resNum:",len(e),"std:",np.nanstd(e))
# 	# print("Lost: avgConSurfDB:",np.nanmean(f),"resNum:",len(f),"std:",np.nanstd(f))

# 	xlabels = ["All PPIs","Scer PPIs","Spom PPIs"]
# 	colors = ['powderblue', 'lightcoral', 'powderblue', 'lightcoral','powderblue', 'lightcoral']
# 	legend=["Preserved","Lost"]
# 	boxplots([a,b,c,d,e,f],"Preserved vs. lost residues","ConSurf score (DB)",xlabels,colors,legend,PLOT_DIR+"/preserved_vs_lost_consurfDB_interfacial",3)

# 	# Consurf rate4site
# 	a=[res.consurf_score for res in preserved if not np.isnan(res.consurf_score)]
# 	b=[res.consurf_score for res in lost if not np.isnan(res.consurf_score)]
# 	c=[res.consurf_score for res in preservedScer if not np.isnan(res.consurf_score)]
# 	d=[res.consurf_score for res in lostScer if not np.isnan(res.consurf_score)]
# 	e=[res.consurf_score for res in preservedSpom if not np.isnan(res.consurf_score)]
# 	f=[res.consurf_score for res in lostSpom if not np.isnan(res.consurf_score)]
# 	# print("\nScer PPIs")
# 	# print("Preserved: avgConSurfrate4site:",np.nanmean(c),"resNum:",len(c),"std:",np.nanstd(c))
# 	# print("Lost: avgConSurfrate4site:",np.nanmean(d),"resNum:",len(d),"std:",np.nanstd(d))
# 	# print("Spom PPIs")
# 	# print("Preserved: avgConSurfrate4site:",np.nanmean(e),"resNum:",len(e),"std:",np.nanstd(e))
# 	# print("Lost: avgConSurfrate4site:",np.nanmean(f),"resNum:",len(f),"std:",np.nanstd(f))
	
# 	xlabels = ["All PPIs","Scer PPIs","Spom PPIs"]
# 	colors = ['powderblue', 'lightcoral', 'powderblue', 'lightcoral','powderblue', 'lightcoral']
# 	legend=["Preserved","Lost"]
# 	boxplots([a,b,c,d,e,f],"Preserved vs. lost residues","ConSurf score (rate4site)",xlabels,colors,legend,PLOT_DIR+"/preserved_vs_lost_consurfrate4site_interfacial",1)

	
# 	return
def msa_run_clustalw(fasta_path,path_to_clustalo):
	""" 
	Run CLUSTALW to generate an MSA from a fasta file
	Uses local version of CLUSTALW for mac.
	Downloaded in download_clustal 
	
	"""
	print ("----------- Run MSA -----------")
	test_ORF = os.listdir(fasta_path)[0]
	if not os.path.exists(os.path.join(fasta_path,test_ORF,test_ORF+".nt.aln")) :
		for ORF in tqdm(os.listdir(fasta_path)):
			if ORF == ".DS_Store":
					continue
			nt_fasta = os.path.join(fasta_path,ORF,ORF+".nt.fa")

			if not os.path.exists(nt_fasta[:-3]+'.aln'):
				os.system(path_to_clustalo+" -i "+ nt_fasta + ' -o '+nt_fasta[:-3]+'.aln')
				# os.system("../../data/clustal/clustalo -i "+ nt_fasta + ' -o '+nt_fasta[:-3]+'.aln')

	if not os.path.exists(os.path.join(fasta_path,test_ORF,test_ORF+".aa.aln")) :
		for ORF in tqdm(os.listdir(fasta_path)):
			if ORF == ".DS_Store":
					continue
			aa_fasta = os.path.join(fasta_path,ORF,ORF+".aa.fa")

			if not os.path.exists(aa_fasta[:-3]+'.aln'):
				os.system(path_to_clustalo+" -i "+ aa_fasta + ' -o '+aa_fasta[:-3]+'.aln')
				# os.system("../../data/clustal/clustalo -i "+ nt_fasta + ' -o '+nt_fasta[:-3]+'.aln')

	if os.path.exists(os.path.join(fasta_path,test_ORF,test_ORF+".aa.aln")) and os.path.exists(os.path.join(fasta_path,test_ORF,test_ORF+".nt.aln")) :
		print ("CLUSTALW MSA already performed. Skipping.")

	exit()
		
	return 

def main():
	SPOM_PPI_MODELS_PATH = '../data/pombe/modeled_ppis_new.pkl'
	SCER_PPI_MODELS_PATH = '../data/cerevisiae/modeled_ppis_new.pkl'
	tree_file_spom = "../data/pombe/Spom.tree"
	tree_file_scer = "../data/cerevisiae/Scer.tree"
	alignment_folder_spom = "../data/pombe/MSA"
	alignment_folder_scer = "../data/cerevisiae/MSA"
	# PLOT_DIR = '../plots'
	rate4site_dir = '../data/rate4site' 
	INFILES_DIR_pombe = os.path.join(rate4site_dir, 'pombe/input_files')
	OUTFILE_DIR_pombe = os.path.join(rate4site_dir, 'pombe/output_files')
	INFILES_DIR_cere = os.path.join(rate4site_dir, 'cere/input_files')
	OUTFILE_DIR_cere = os.path.join(rate4site_dir, 'cere/output_files')
	
	
	path_to_rate4site_exec= download_rate4sites()

	# Spom
	inputs,modeledPPIs = set_up_for_rate4sites(SPOM_PPI_MODELS_PATH,alignment_folder_spom,INFILES_DIR_pombe,OUTFILE_DIR_pombe,tree_file_spom,"pombe")

	run_rate4sites(path_to_rate4site_exec,inputs,INFILES_DIR_pombe,OUTFILE_DIR_pombe)
	add_to_protein_models(inputs,modeledPPIs,OUTFILE_DIR_pombe,SPOM_PPI_MODELS_PATH)

	
	# Scer
	inputs,modeledPPIs = set_up_for_rate4sites(SCER_PPI_MODELS_PATH,alignment_folder_scer,INFILES_DIR_cere,OUTFILE_DIR_cere,tree_file_scer,"scer")
	run_rate4sites(path_to_rate4site_exec,inputs,INFILES_DIR_cere,OUTFILE_DIR_cere)
	add_to_protein_models(inputs,modeledPPIs,OUTFILE_DIR_cere,SCER_PPI_MODELS_PATH)

if __name__ == '__main__':
    main()




