import sys
from argparse import ArgumentParser

zerteiler = ArgumentParser()

zerteiler.add_argument("--gff-file", dest="gff", required=True, help="gff file from which genes are to be extracted.")
zerteiler.add_argument("--fasta-file", dest="fasta", required=True, help="fasta file which contains the gene sequences.")
zerteiler.add_argument("--out", dest="out", required=False, help="Output file. Defaults to ./output.fa.")


args = zerteiler.parse_args()

gff_file = open(args.gff, 'r')
names = []

for line in gff_file:
    components = line.split("\t")
    if len(components)>5:
        chromosome=components[0]
        designation=components[2]
        start=int(components[3])
        stop=int(components[4])
        name=components[8]

        if chromosome == 'Lp_chr6_0' and designation == 'gene' and start > 242000000 and stop < 255000000: #Hardcoded for my purposes
            name=name.split(";")[0]
            name=name.split("=")[1]
            names.append(name)
gff_file.close()

fasta_file = open(args.fasta, 'r')
if args.out :
    out_file = open(args.out, 'w')
else: 
    out_file = open("./output.fa", 'w')

currentnname=""
currentsequence=""
inscope = False
count=0
for line in fasta_file:
    if line[0] == '>':
        name = line.split("gene")[1]
        name = name.split(" ")[0]
        name = name[1:].strip()
        if name in names:
            if count >= 1:
                out_file.write(">"+currentnname+"\n"+currentsequence)  
            currentnname=name
            currentsequence=""
            inscope=True
            count+=1
        else:
            inscope=False
    else:
        if inscope:
            currentsequence+=line
out_file.write(">"+currentnname+"\n"+currentsequence) 
fasta_file.close()
out_file.close()
print(str(count)+" sequences extracted.")
print("All done.")

