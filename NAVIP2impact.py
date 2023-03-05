import sys
from argparse import ArgumentParser

#This program is designed to work with files created by NAVIP. I discontinued working on this script, it is not relevant to my bachelor in 2023.

zerteiler = ArgumentParser()

zerteiler.add_argument("--in1", dest="in1", required=True, help=".vcf file created by NAVIP. Assumed 50 percent fertility.")
zerteiler.add_argument("--in2", dest="in2", required=True, help=".vcf file created by NAVIP. Assumed 100 percent fertility.")
zerteiler.add_argument("--genes", dest="genes", required=True, help="A .txt files containing relevant genes.")
zerteiler.add_argument("--out1", dest="out1", required=False, help="Output file with Variations, their position, impact and probable effect. Defaults to ./NAVIP_Variants_impact.tsv.")
zerteiler.add_argument("--out2", dest="out2", required=False, help="Output file with only genes. Formatted to fit an excel sheet. Defaults to ./NAVIP_Genes_impact.tsv.")



args = zerteiler.parse_args()


genes_file = open(args.genes, 'r')
top_genes =[]
fifty_genes = {}
hundred_genes = {}
for line in genes_file:
    top_genes.append(line.strip())
    fifty_genes[line.strip()] = []
    hundred_genes[line.strip()] = []
genes_file.close()

if args.out1 :
    out1_file = open(args.out1, 'w')
else: 
    out1_file = open("./NAVIP_Variants_impact.tsv", 'w')
if args.out2 :
    out2_file = open(args.out2, 'w')
else: 
    out2_file = open("./NAVIP_Genes_impact.tsv", 'w')

out2_file.write("Gene\tNAVIP_50_Max_Impact\tNAVIP_100_Max_Impact\n")

fifty_file = open(args.in1, 'r')

hundred_file = open(args.in2, 'r')

chromosome = "Lp_chr6_0"

for line in fifty_file:
    if line[0] is not '#': #Skip comments
        components = line.split("\t")
        position = int(components[1])
        info = components[7]
        NAV2 = info.split("NAV2=")[1].split("NAV1=")[0].split('|')
        NAV1 = info.split("NAV1=")[1].split('|')
        impact = NAV2[2]
        gene = NAV2[6]
        effect = NAV2[1]

    for i in range(len(genes)):
        if '-' not in genes[i] and genes[i] in fifty_genes: #If it is not intergenic (Might still be downstream or something, but those are alwazs MODIFIER) | if gene is relevant.
            fifty_genes[genes[i]].append(impacts[i])

for line in hundred_file:
    components = line.split("\t")               
    genes = components[1].split(",")
    impacts = components[3].split(",")

    for i in range(len(genes)):
        if '-' not in genes[i] and genes[i] in hundred_genes: #If it is not intergenic (Might still be downstream or something, but those are alwazs MODIFIER) | if gene is relevant.
            hundred_genes[genes[i]].append(impacts[i])
       

for gene in top_genes:
    fifty ="NA"
    hundred = "NA"
    if "HIGH" in fifty_genes[gene]:
        fifty = "HIGH"
    elif "MODERATE" in fifty_genes[gene]:
        fifty = "MODERATE"
    elif "LOW" in fifty_genes[gene]:
        fifty = "LOW"
    else:
        fifty = "MODIFIER"

    if "HIGH" in hundred_genes[gene]:
        hundred = "HIGH"
    elif "MODERATE" in hundred_genes[gene]:
        hundred = "MODERATE"
    elif "LOW" in hundred_genes[gene]:
        hundred = "LOW"
    else:
        hundred = "MODIFIER"

    out_file.write(gene + "\t" + fifty + "\t" + hundred + "\n")



fifty_file.close()
hundred_file.close()


out_file.close()
print("Tables created.")
print("All done.")