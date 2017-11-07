### IMPORT LIBRARIES ###
from os.path import abspath
from os import makedirs
from os.path import isdir
import re
from glob import glob
import pandas as pd
import BBH_functions as BF

### GET PARAMETERS ###
args = BF.get_comp_parms()

### IF INPUT IS GB FILE, UNPACK TO FASTA; OTHERWISE ADD INDEX  ###
#For synteny analysis it will be necessary to have some sort of index value in the FASTA header.
#Also, for separate contigs, it will be necessary to do SOMETHING to the index to avoid cross-contig artifacts.
ingroup_proteomes  = BF.index_proteomes(args.g1)
outgroup_proteomes = BF.index_proteomes(args.g2)

### MAKE BLAST/DIAMOND INDICES ###
BF.makeblastdb(ingroup_proteomes)
BF.makeblastdb(outgroup_proteomes)

### RUN BLAST/DIAMOND SEARCHES ###
if args.R :
    reference = abspath(args.g1 + args.R.split(".")[0]+".faa")
else :
    reference = ingroup_proteomes[0]
#Make a folder for the BLAST results.
blastoutputdir = args.o + "/BLAST_data"
if not isdir(blastoutputdir) :
    makedirs(abspath(blastoutputdir))

#The hitlist is a list of pandas dataframes with the results
in_hitlist, out_hitlist = [], []

#Loop thru each proteome; run BLAST; parse results
for i in (ingroup_proteomes) :
    if i != reference :
        in_hitdict = BF.runblast(reference, i, blastoutputdir, evalue = args.E)
        in_hitlist.append(pd.DataFrame(hitdict).transpose())

for i in (outgroup_proteomes) :
    out_hitdict = BF.runblast(reference, i, blastoutputdir, evalue = args.E)
    out_hitlist.append(pd.DataFrame(hitdict).transpose())

### PARSE RESULTS TO OBTAIN SYNTENY SCORES ###
results = BF.parse_results(hitlist, window = args.W)

### PRINT RESULT TABLES ###
Verbose or regular.

### INTERACTIVE RESULTS ###
Plot ingroup genome #1 on X axis, ingroup genome #2 on Y axis, color points to indicate homology scores.
Show protein ID and annotation with mouseover.
