#!/usr/bin/env python3 
# This ^ suddenly became necessary? Why?
import sys
from argparse import ArgumentParser

zerteiler = ArgumentParser()

zerteiler.add_argument("--gff-file", dest="gff", required=True, help="gff file from which genes are to be extracted.")
zerteiler.add_argument("--candidates", dest="candidates", required=True, help="Candidate genes .tsv file. Will determine exon count and total cds length.")
zerteiler.add_argument("--out", dest="out", required=False, help="Output file. Defaults to ./exons_cds.tsv.")


args = zerteiler.parse_args()

gff_file = open(args.gff, 'r')

c_file = open(args.candidates, 'r')
names={}
c_file.readline() #skip headers
for line in c_file:
    names[line.split("\t")[0]] = (0,0)
c_file.close()

for line in gff_file:
    components = line.split("\t")
    if len(components)>5:
        chromosome=components[0]
        designation=components[2]
        start=int(components[3])
        stop=int(components[4])
        name=components[8]

        if chromosome == 'Lp_chr6_0' and (designation == 'exon' or designation == 'CDS')  and start > 242000000 and stop <255000000: #Hardcoded for my purposes
            name=name.split(";")[1]
            name=name.split("=")[1][:-3] #get parent name. First element at the end ist '\n', second is '1', third is '.' all need to be removed for matching.
            print(name)
            if name in names:
                
                (exons, cds_length) = names[name]
                if designation == 'exon':
                    exons += 1
                    names[name] = (exons, cds_length)
                elif designation == 'CDS':
                    cds_length += (stop - start)
                    names[name] = (exons, cds_length)
gff_file.close()

if args.out :
    out_file = open(args.out, 'w')
else: 
    out_file = open("./exons_cds.tsv", 'w')


out_file.write("Gene \t Number of exons \t Total CDS length \n") 
for gene, (exons, cds_length) in names.items(): #TODO
    out_file.write(gene + "\t" + str(exons) + "\t" + str(cds_length) + "\n")

out_file.close()
print("Found exons and cds lengths for " + str(len(names)) + " genes.")
print("All done.")

