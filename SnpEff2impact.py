import sys
from argparse import ArgumentParser

#This program is designed to work with files created by SnpEff as of 7.10.2022. The SnpSift files have to have the extracted fields: CHROM ANN[*].GENEID POS ANN[*].IMPACT ANN[*].EFFECT ANN[*].BIOTYPE 

zerteiler = ArgumentParser()

zerteiler.add_argument("--in1", dest="in1", required=True, help=".ann.vcf.tsv file created bz SnpSift with fields CHROM ANN[*].GENEID POS ANN[*].IMPACT ANN[*].EFFECT ANN[*].BIOTYPE. Assumed 50 percent fertility.")
zerteiler.add_argument("--in2", dest="in2", required=True, help=".ann.vcf file created bz SnpSift with fields CHROM ANN[*].GENEID POS ANN[*].IMPACT ANN[*].EFFECT ANN[*].BIOTYPE. Assumed 100 percent fertility.")
zerteiler.add_argument("--genes", dest="genes", required=True, help="A .tsv files containing relevant genes and their start and stop positions.")
zerteiler.add_argument("--out", dest="out", required=False, help="Output file with only genes. Formatted to fit an excel sheet. Defaults to ./SNPEFF_Genes_impact.tsv.")



args = zerteiler.parse_args()


genes_file = open(args.genes, 'r')
genes_file.readline()
top_genes =[]
fifty_genes = {}
hundred_genes = {}
for line in genes_file:
    components = line.split("\t")
    gene = components[0]
    start = float(components[1])
    stop = float(components[2].strip())
    top_genes.append((gene,start,stop))
    fifty_genes[gene] = []
    hundred_genes[gene] = []
genes_file.close()


if args.out :
    out_file = open(args.out, 'w')
else: 
    out_file = open("./SNPEFF_Genes_impact.tsv", 'w')

out_file.write("Gene\tSNPEff_50_Max_Impact\tSNPEff_50_Max_Impact_Variation_Type\tSNPEff_50_Max_Impact_Relative_Position\tSNPEff_100_Max_Impact\tSNPEff_100_Max_Impact_Variation_Type\tSNPEff_100_Max_Impact_Relative_Position\n")

fifty_file = open(args.in1, 'r')

hundred_file = open(args.in2, 'r')

fifty_file.readline() #skip headers
hundred_file.readline()

for line in fifty_file:
    components = line.split("\t")               
    genes = components[1].split(",")
    position = components[2]
    impacts = components[3].split(",")
    effects = components[4].split(",")

    for i in range(len(genes)):
        if '-' not in genes[i] and genes[i] in fifty_genes: #If it is not intergenic (Might still be downstream, but those are always MODIFIER) | if gene is relevant.
            fifty_genes[genes[i]].append((impacts[i],float(position),effects[i]))

for line in hundred_file:
    components = line.split("\t")               
    genes = components[1].split(",")
    position = components[2]
    impacts = components[3].split(",")
    effects = components[4].split(",")

    for i in range(len(genes)):
        if '-' not in genes[i] and genes[i] in hundred_genes: #If it is not intergenic (Might still be downstream, but those are always MODIFIER) | if gene is relevant.
            hundred_genes[genes[i]].append((impacts[i],float(position),effects[i]))
       

for (gene,start,stop) in top_genes:
    fifty ="NA\tNA\tNA"
    hundred = "NA\tNA\tNA"
    high = [(str(round((x[1]-start)/(stop-start)*100,1))+"%",x[2]) for x in fifty_genes[gene] if x[0] == "HIGH" and (x[1]-start)/(stop-start)*100 >=0 and (x[1]-start)/(stop-start)*100 <=100 and ("stop_gained" in x[2] or "frameshift_variant" in x[2])]
    moderate = [(str(round((x[1]-start)/(stop-start)*100,1))+"%",x[2]) for x in fifty_genes[gene] if x[0] == "MODERATE" and (x[1]-start)/(stop-start)*100 >=0 and (x[1]-start)/(stop-start)*100 <=100 and ("stop_gained" in x[2] or "frameshift_variant" in x[2])]
    low = [(str(round((x[1]-start)/(stop-start)*100,1))+"%",x[2]) for x in fifty_genes[gene] if x[0] == "LOW" and (x[1]-start)/(stop-start)*100 >=0 and (x[1]-start)/(stop-start)*100 <=100 and ("stop_gained" in x[2] or "frameshift_variant" in x[2])]
    modifier = [(str(round((x[1]-start)/(stop-start)*100,1))+"%",x[2]) for x in fifty_genes[gene] if x[0] == "MODIFIER" and (x[1]-start)/(stop-start)*100 >=0 and (x[1]-start)/(stop-start)*100 <=100 and ("stop_gained" in x[2] or "frameshift_variant" in x[2])]

    if len(high) > 0:
        percentile, effect = zip(*high)
        fifty = "HIGH\t" + high[0][0] + "\t" + high[0][1] #This selects only the first variant. The following lines are technically unneeded since no stop_gained/frameshift variant will have lower than high impact. They're kept here for future inspiration.
    elif len(moderate) > 0:
        percentile, effect = zip(*moderate)
        fifty = "MODERATE\t" + ",".join(effect) + "\t" + ",".join(percentile)
    elif len(low) > 0:
        percentile, effect = zip(*low)
        fifty = "LOW\t" + ",".join(effect) + "\t" + ",".join(percentile)
    elif len(modifier) > 0:
        percentile, effect = zip(*modifier)
        fifty = "MODIFIER\t" + ",".join(effect) + "\t" + ",".join(percentile)
    
    high = [(str(round((x[1]-start)/(stop-start)*100,1))+"%",x[2]) for x in hundred_genes[gene] if x[0] == "HIGH" and (x[1]-start)/(stop-start)*100 >=0 and (x[1]-start)/(stop-start)*100 <=100 and ("stop_gained" in x[2] or "frameshift_variant" in x[2])]
    moderate = [(str(round((x[1]-start)/(stop-start)*100,1))+"%",x[2]) for x in hundred_genes[gene] if x[0] == "MODERATE" and (x[1]-start)/(stop-start)*100 >=0 and (x[1]-start)/(stop-start)*100 <=100 and ("stop_gained" in x[2] or "frameshift_variant" in x[2])]
    low = [(str(round((x[1]-start)/(stop-start)*100,1))+"%",x[2]) for x in hundred_genes[gene] if x[0] == "LOW" and (x[1]-start)/(stop-start)*100 >=0 and (x[1]-start)/(stop-start)*100 <=100 and ("stop_gained" in x[2] or "frameshift_variant" in x[2])]
    modifier = [(str(round((x[1]-start)/(stop-start)*100,1))+"%",x[2]) for x in hundred_genes[gene] if x[0] == "MODIFIER" and (x[1]-start)/(stop-start)*100 >=0 and (x[1]-start)/(stop-start)*100 <=100 and ("stop_gained" in x[2] or "frameshift_variant" in x[2])]

    if len(high) > 0:
        percentile, effect = zip(*high)
        hundred = "HIGH\t" + high[0][0] + "\t" + high[0][1] #This selects only the first variant. The following lines are technically unneeded since no stop_gained/frameshift variant will have lower than high impact. They're kept here for future inspiration.
    elif len(moderate) > 0:
        percentile, effect = zip(*moderate)
        hundred = "MODERATE\t" + ",".join(effect) + "\t" + ",".join(percentile)
    elif len(low) > 0:
        percentile, effect = zip(*low)
        hundred = "LOW\t" + ",".join(effect) + "\t" + ",".join(percentile)
    elif len(modifier) > 0:
        percentile, effect = zip(*modifier)
        hundred = "MODIFIER\t" + ",".join(effect) + "\t" + ",".join(percentile)

    out_file.write(gene + "\t" + fifty + "\t" + hundred + "\n")



fifty_file.close()
hundred_file.close()


out_file.close()
print("Tables created.")
print("All done.")