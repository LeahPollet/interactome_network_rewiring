Scripts needed to run the analyis
-----------
Requirements:
-----------
Software and packages:
----------
(Unless otherwise specified all software and packages used in the pipeline were downloaded and compiled as per website instructions, and the path to the final executable needed was added to the PATH environment variable)

Python 3.7.3
Protein-Protein BLAST 2.7.1+
Clustal Omega 1.2.3 - Downloaded and compiled locally in the process_data.py script
DSSP 2.0.4 - Downloaded and compiled locally in the build_protein_models.py script
paml4.9j - Downloaded and compiled locally in the calc_dnds.py script

Libraries:
----------
Standard Python 3.7.3 libraries 
feather 0.4.0
scipy 1.4.1
json 2.0.9
tqdm 4.42.0
Bio 1.76
numpy 1.18.1
pandas 1.0.1
matplotlib 3.1.2
Databases:
----------
http://thebiogrid.org
https://pdb101.rcsb.org
https://www.ebi.ac.uk
http://ensemblgenomes.org
https://yeastgenome.org
https://consurfdb.tau.ac.il

Running the pipeline will download and generate files into a data and a plot folder

Any additional information can be requested by email at leah.pollet@mail.mcgill.ca.
