import sys
from argparse import ArgumentParser

#This script is designed for InterProScan 5.52.
#I did something wrong in preprocessing, the genes still have a 1 at the end of their names. This script removes them.

zerteiler = ArgumentParser()

zerteiler.add_argument("--genes-file", dest="genes", required=True, help=".txt file with gene names as lines.")
zerteiler.add_argument("--ips-file", dest="ips", required=True, help="File created by InterProScan in .tsv format.")
zerteiler.add_argument("--out", dest="out", required=False, help="Output file. Defaults to ./output.tsv.")


args = zerteiler.parse_args()

genes_file = open(args.genes, 'r')
names = []
dict = {}

for line in genes_file:
    gene = line.strip()
    names.append(gene)
    dict[gene] = "NA" 

genes_file.close()

ips_file = open(args.ips, 'r')
if args.out :
    out_file = open(args.out, 'w')
else: 
    out_file = open("./output.tsv", 'w')



for line in ips_file:
    gene = line.split("\t")[0][:-1]
    if gene in names:
        dict[gene] = line

for gene in names:
    if dict[gene] is not "NA":
        components = dict[gene].split("\t")
        components[0] = gene

        out_file.write("\t".join(components))
    else:
        out_file.write(gene + "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")

ips_file.close()
out_file.close()
print("All done.")

