#!/usr/bin/python
#in: folder with in group faas, folder with out group faas.
#out: For each in in group, table with WP no, annotation?, [score, top homologue ID] for every species
from glob import glob
from os.path import abspath
from itertools import combinations
import BBH_functions as BF
import argparse
import re
import pickle

parser = argparse.ArgumentParser()
parser.add_argument("-g1",
                    help = "A directory with fasta files for proteomes of ingroup")
parser.add_argument("-g2",
                    help = "A directory with fasta files for proteomes of outgroup")
parser.add_argument("-o",
                    help = "the name of the directory where you wish to store output")
args = parser.parse_args()

ingroup    = glob(abspath(args.g1 + "/*.faa"))
outgroup    = glob(abspath(args.g2 + "/*.faa"))


import pandas as pd
import pickle

for i in ingroup :
    query_spec  = re.split('[\\\/.]+', i)[-2]
    hitlist = []
    for x in (ingroup + outgroup) :
        hitdict  = BF.diamond_index(i, x, args.o)
        hitframe = pd.DataFrame(hitdict).transpose()
        hitlist.append(hitframe)
    frame = hitlist[0]
    for h in hitlist[1:]:
        frame   = pd.merge(frame, h, right_index = True, left_index = True, how = "left", sort = True)
    #replace empties for ones where there is no hit with a zero.
    pickle.dump(frame, open(abspath(args.o + "/" + "dump"), "wb"))
    frame = frame.fillna(value = 0)
    frame.to_csv(abspath(args.o + "/" + query_spec + ".csv"))
