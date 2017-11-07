## Get_BBH

#### Description
Get_BBH is a tool for retrieving Best Bilateral BLAST hits (BBHs) between two fasta databases.

#### Usage
*get_BBH -db1 your_first_faa_file.faa -db2 your_2nd_FAA_file.faa -out /some/path/youroutputdir --evalue 1E-50*

#### Arguments

#### Output


## comparetron

#### Description

The comparetron is a special-purpose tool for the comparison of multiple microbial genomes.  This tool is useful in cases where at least two species are known to have a particular molecular phenotype of unknown genetic origin.  For example, they produce a natural product with unknown biosynthesis.  

#### Usage
*comparetron -ingroup_folder -outgroup_folder*

#### Arguments

#### Output

Right now, prints the IDs of best bilateral BLAST hit pairs.  Pretty messy.  

#### Example

### TODO:
- Finish the get_homologues script, which will retrieve sequences with the same best hit from 2x reference sets.  For assembling MSA of large complex type stuff.
