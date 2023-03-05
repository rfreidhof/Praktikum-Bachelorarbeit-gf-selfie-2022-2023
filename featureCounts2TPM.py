import sys
from argparse import ArgumentParser

#This programm is specifically for featureCounts output files as of 20.09.2022 

zerteiler = ArgumentParser()

zerteiler.add_argument("--in", dest="input", required=True, help="featureCounts output file.")
zerteiler.add_argument("--candidates", dest="candidates", required=True, help="Candidate gene file. Used to shrink the output to a useful size.")
zerteiler.add_argument("--out", dest="out", required=False, help="Output file. Defaults to ./tpms.tsv.")



args = zerteiler.parse_args()

in_file = open(args.input, 'r')

c_file = open(args.candidates, 'r')

if args.out :
    out_file = open(args.out, 'w')
else: 
    out_file = open("./tpms.tsv", 'w')


names=[]
c_file.readline() #skip headers
for line in c_file:
    names.append(line.split("\t")[0])
c_file.close()


in_file.readline() #skip comment line
headers = in_file.readline()
geneHeader = headers.split("\t")[0]
sample = headers.split("\t")[6].strip()

rpksum = float(0)
rpks = []
geneids = []

for line in in_file:
    components = line.split("\t")
    gene = components[0]
    gene = gene.split(".")[0] + "." + gene.split(".")[1] #Removes transcript numbers.
    length = float(components[5])
    reads = float(components[6])
    rpk = reads/(length/1000)
    rpksum += rpk
    rpks.append(rpk)
    geneids.append(gene)

per_million_scaling_factor = rpksum/1000000

if per_million_scaling_factor == 0:
    print("No reads were mapped. This might be ok, but probably points to some error.")
    print(args.input)

tpms = [item/per_million_scaling_factor for item in rpks]

out_file.write(geneHeader + "\t" + "tpm(" + sample + ")" + "\n") 

for i in range(len(geneids)):
    if geneids[i]in names:
        out_file.write(geneids[i] + "\t" + str(tpms[i]) + "\n") 

in_file.close()
out_file.close()
print("Table created.")
print("All done.")

