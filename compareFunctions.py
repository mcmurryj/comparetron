#!/usr/bin/python
"""Functions for the comparison of microbial genomes."""
import argparse
from os.path import isdir
from os.path import isfile
from os.path import abspath
from os import makedirs
import subprocess
import re
from glob import glob
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
from matplotlib.pyplot import figure, show
import numpy as npy
import pandas as pd

def get_comp_parms():
    """A function to pull parameters for running comparative genomic analysis."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-g1",
                        help="A directory with fasta files for proteomes of ingroup.")
    parser.add_argument("-g2",
                        help="A directory with fasta files for proteomes of outgroup.")
    parser.add_argument("-o",
                        help="the name of the directory where you wish to store output.")
    parser.add_argument("--V",
                        default=False,
                        help="If true, output more detailed tables rather than summary tables.")
    parser.add_argument("--E",
                        default=1E-4,
                        help="Change the default expecation value cutoff for determining BBHes.")
    parser.add_argument("--R",
                        default=False,
                        help="Select an ingroup proteome to use as the reference.")
    parser.add_argument("--W",
                        default=20000,
                        help="Window in bp for synteny scoring.")
    args = parser.parse_args()
    return args

def skin_gb(gbDir):
    """A function that renders a genbank file into a multi-sequence amino acid fasta file.
    Important contextual information is stored in the name of each protein.
    This allows one to run BLAST or similar while retaining genomic context information.
    The information currently includes the following:
    organism, contig (acession), protein ID, annotation, nucleic acid location.
    Usage:
    skin_gb(path/to/directory/full_of/genbank/files/forexample/foo.gb)
    Returns:
    There is no return value, it just prints files.
    The format is kludgy-improvements in compatibility with DIAMOND and better
    sanitization of input strings would help.  But it does work."""

    #Initialize misc variables and read the files
    gbList = glob(gbDir + "*.gb*")
    delim = "|"
    faaList = []

    #Iterate through the list of genbank files.
    for gbfile in gbList:
        gb = SeqIO.parse(gbfile, "genbank")		#parse needed for multirecord
        faa = abspath(gbDir + gbfile.split(".")[0]+".faa") #the name of a corresponding FASTA file
        faaList.append(faa) #stick the file name in a list
        out = open(faa, 'a') #open a file handle corresponding to the faa file name.
        first = True #Assign variable to tell script if it is the first run thru a genbank file.
        #For each record in the genbank, extract info and write it to a FASTA file.
        for rec in gb:
            #stuff you only have to do for first record, stays the same afterwards.
            if first is True:
                first = False
                #Retrieve organism name; sanitize illegal chars; fill spaces with _
                organism = re.sub('[$%^&*#|\[\]]',
                                  "",
                                  rec.annotations['organism'].replace(" ", "_"))
                spec = "_".join(organism.split("_")[:2])    #Species is first 2x chunks of the organism name

            #Check that the contig is of a reasonable length. If it <20,000 bp, skip the contig.
            if len(rec.seq) < 20000:
                continue
            #Pull the contig name.  Need to do for every record.
            accession = re.sub('[$%^&*#|]', "", rec.id)
            ###note: feature parsing works for refseq records from ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/
            for feat in rec.features:
                if 'protein_id' in feat.qualifiers:	###Not funtional with EG JGI genomes, yet.
                    wp_id = re.sub('[$%^&*#|]', "", feat.qualifiers['protein_id'][0])
                    if 'product' in feat.qualifiers:
                        annotation = re.sub('[$%^&*#|]', "", feat.qualifiers['product'][0].replace(" ", "_"))
                    else:
                        annotation = "No annotation."
                    ss = str(feat.location.start) + "-" + str(feat.location.end)
                    #I forget, may need to add a newline as I've switched from print to write.
                    out.write(">" + delim.join([organism, accession, wp_id, annotation, ss]) + "\n")
                    out.write((feat.extract(rec.seq).translate()) + "\n")
    #return a list of amino acid fasta file names to use.
    return faaList

def runblast(db1, db2, blastoutputdir, evalue=1E-30):
    """A wrapper function that runs all vs. all blast between two genomes, and parses the results.
    The returned object is a dictionary containing important attributes of the results.
    For example, BLAST bit score, query ID, hit ID, and so on.
    Input:
    db1 and db2 are paths to faa files (proteomes) that have been configured
    as BLAST databases.
    output_dir is a blank directory for putting the BLAST results.
    evalue is the threshold value to use for BLAST searches.
    Returns:
    Dictionary.  Keys are BLAST query protein IDS.  Values are dictionaries with
    information about the query and the top hit.
    """
    #make output directory
    if not isdir(blastoutputdir):
        makedirs(blastoutputdir)
    #make an identifier for the BLAST output.
    outputID = "/" + re.split('[\\\/.]+', db1)[-2] + "_vs_" + re.split('[\\\/.]+', db2)[-2] + ".XML"
    blastoutputfile1 = abspath(blastoutputdir + outputID)
    # run the BLAST search.  NOTE-not configured for BBH-ing right now.
    cline = NcbiblastpCommandline(query=db1,
                                  db=db2,
    							  evalue=evalue,
    							  outfmt=5,
    							  out=blastoutputfile1)
    cline()
    #Parse the BLAST
    result_handle = open(blastoutputfile1, "r")
    BLASTrecs = NCBIXML.parse(result_handle)
    idDict = {}
    #For each record, pull all the info.
    for B in BLASTrecs:
        qspecies, qContig, qId, qAnnotation, qLocation = B.query.split("|")
        qSTART, qSTOP = qLocation.split("-")
        if B.alignments:
            ali = B.alignments[0]
            species, hContig, hId, hAnnotation, hLocation = ali.hit_def.split("|")
            hStart, hStop = hLocation.split("-")
            #THis is only the score of the top HSP, not for the whole thing
            score = ali.hsps[0].score
            idDict[qId] = {"h_ID_" + qspecies      : hId,
                             "q_contig_" + qspecies  : qContig,
                             "h_contig_" + qspecies  : hContig,
                             "q_location_" + qspecies: qLocation,
                             "h_location_" + qspecies: (int(hStart) + int(hStop))/2,
                             "q_annot_" + qspecies   : qAnnotation,
                             "h_annot_" + qspecies   : hAnnotation,
                             "score_" + qspecies     : score}
    return idDict

def colstart(cols, starter):
    """Helper function for slicing dataframes with varied column names.
    Returns a boolean vector; Trues are for columns that begin with the
    string specified by variable starter."""
    return [c.startswith(starter) for c in cols]

def colstartName(cols, starter):
    """Helper function for slicing dataframes with varied column names.
    Returns string corresponding to the column that begins with
    the string specified by variable starter. May give unexpected results if
    more than one column starts with this string."""
    return [c for c in cols if c.startswith(starter)][0]

def get_syn(hitlist, window):
    """Retrieve synteny information from dict of all v. all BLAST results.
    Performs calculations, then returns a pandas dataframe.
    Input:
    hitlist is a list of pandas dataframes with BLAST results in tabular form.
    window is the # of nucleotides that defines if 2 proteins are adjacent for synteny purposes.
    Output:
    Returns a single dataframe combining the results of multiple BLAST searches.
    """
    for hf in hitlist:
        hf['syn_score'] = pd.Series([0] * len(hf.index))
        hLocIndex = colstartName(hf.columns, "h_location_")
        hContigIndex = colstartName(hf.columns, "h_contig_")
        for i in hf.index:
            ####must fix for compatibility with species names
            h_loc = hf.loc[i, hLocIndex]
            h_contig = hf.loc[i, hContigIndex]
            hf.loc[i, 'syn_score'] = sum((abs(hf.loc[:, hLocIndex] - h_loc) <  window)
                                     and
                                          (hf.loc[:, hContigIndex] == h_contig))
        frame = pd.merge(frame, hf, right_index=True, left_index=True, how="left", sort=True)
        #FILL IN NA VALUES
        frame = frame.fillna(value=0)
        return frame

def summarize(frame):
    """Get summary stats from BLAST scores and synteny scores."""
    #In the future multi-indexing might be better but this is easier.
    score_cols = colstart(frame.columns, "score_")
    frame['mean_score'] = frame.loc[:, score_cols].mean(axis=1)
    frame['sd_score'] = frame.loc[:, score_cols].std(axis=1)
    syn_cols = colstart(frame.columns, "syn_score_")
    frame['mean_syn'] = frame.loc[:, syn_cols].mean(axis=1)
    frame['sd_syn'] = frame.loc[:, syn_cols].std(axis=1)
    return frame

def parse_results(in_hitlist, out_hitlist, window=args.W):
    """Wrapper script.  Takes ingroup data and outgroup data.
    Processes the data using the get_syn and summarize functions.
    Return a merged frame with all the good stuff."""
    in_frame = get_syn(in_hitlist,  window)
    in_frame = summarize(in_frame)
    out_frame = get_syn(out_hitlist, window)
    out_frame = summarize(out_frame)
    merg_frame = pd.merge(in_frame, out_frame, right_index=True, left_index=True, how="left", sort=True)
    return merg_frame


def makeblastdb(db):
    """Simple function to make protein BLAST db from faa file.
       Input: path to a faa file.
       Output: no return value, just runs makeblastdb command."""
    #Does this file exist already?
    if isfile(db + ".psq"):
        pass
    #If it doesn't make us a protein blast DB.
    else:
        make2ndBLASTdbcmd = "makeblastdb -in " + db + " -input_type fasta -dbtype prot"
        subprocess.call(make2ndBLASTdbcmd, shell=True)

def plotdata(hitframe):
    """Plot the results.
       X index is the nucleotide coords of the reference strain.
       Y index is the coords of some other strain.
       Color indicates how good the blast score is.
       Point size indicates how good the synteny is.
       Not yet resolved: how to handle multiple contigs.  Maybe only plot one contig at a time.
       Very incomplete. """
    #calculations.  the _x and _y is from merging dfs with identical column names.
    syn_score = ((hitframe['mean_syn_x'] - hitframe['mean_syn_y'])/
                (hitframe['sd_syn_x']**2 + hitframe['sd_syn_y']**2)**.5)
    hit_score = ((hitframe['mean_score_x'] - hitframe['mean_score_y'])/
                (hitframe['sd_score_x']**2 + hitframe['sd_score_y']**2)**.5)
    # Matplotlib voodoo related to making data show on mouseover click. Code from:
    #https://stackoverflow.com/questions/7908636/possible-to-make-labels-appear-when-hovering-over-a-point-in-matplotlib
    def onpick3(event):
        """Mystery function for interactive plotting..."""
        ind = event.ind
        print('onpick3 scatter:', ind, npy.take(x, ind), npy.take(y, ind))
    #initialize figue
    fig = figure()
    #Make a white square so you don't look like excel 2003
    rect = fig.patch
    rect.set_facecolor('white')
    # 111 means 1 x 1 grid, 1st subplot.
    # for add_subplot(xyz) x = x dimension, y = y dimension, z = zth subplot
    ax1 = fig.add_subplot(111)
    # other settable parms of note: 'marker' and 'linewidths'
    col = ax1.scatter(hitlist["q_location_XXX"],
                      hitlist["h_location_XXX"],
                      s=syn_score, #s sets the size of the point
                      c=hit_score, #c sets the color of the point; low is blue high is red
                      picker=True)
    fig.canvas.mpl_connect('pick_event', onpick3)
    #Show the plot....
    show()
