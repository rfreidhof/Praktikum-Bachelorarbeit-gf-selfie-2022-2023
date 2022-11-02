import sys
from argparse import ArgumentParser
import xml.etree.ElementTree as ET

#I don't know whether this checks whether a file is actually in the proper format. 
#This script is designed for XML created by NCBI BLAST as of 16.9.2022

zerteiler = ArgumentParser()

zerteiler.add_argument("--xml", dest="xml", required=True, help="Amount of sequences in the smaller fasta files.")
zerteiler.add_argument("--out", dest="out", required=False, help="Output file. Defaults to ./descriptions.tsv.")


args = zerteiler.parse_args()


if args.out :
    out_file = open(args.out, 'w')
else: 
    out_file = open("./descriptions.tsv", 'w')

tree = ET.parse(args.xml)
root = tree.getroot() #<BlastOutput>
out_file.write("Gene" + "\t" + "Query length" + "\t" + "E-value" + "\t" + "Percent Identity" + "\t" + "Description"  + "\t" + "Accession Number" + "\n")
for iteration in root[8]:
    gene = iteration[2].text
    query_length = iteration[3].text
    if len(iteration[4]) > 0: 
        description = iteration[4][0][2].text
        accession = iteration[4][0][3].text
        measurement = iteration[4][0][5][0][3].text #Here e-value
        identity = iteration[4][0][5][0][10].text
        alignment_length = iteration[4][0][5][0][13].text
        percent_identity = str((float(identity)/float(alignment_length))*100)
    
        out_file.write(gene + "\t" + query_length + "\t" + measurement + "\t" + percent_identity + "\t" + description  + "\t" + accession + "\n")
    else:
        out_file.write(gene + "\t" + query_length + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "No match" + "\n")