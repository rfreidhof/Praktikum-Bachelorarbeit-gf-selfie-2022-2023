import sys
from argparse import ArgumentParser

#This script was designed to show a flaw in a vcf file. It is not used in the final version of my bachelor of 2023.

zerteiler = ArgumentParser()

zerteiler.add_argument("--vcf", dest="vcf", required=True, help="vcf file. This needs to be sorted for the program to provide correct output. if it isn\'t then a quicksort like algorithm would be more efficient for determining variants\' positions.")
zerteiler.add_argument("--gff", dest="gff", required=True, help="gff3 file. This needs to be sorted for the program to provide correct output. if it isn\'t then a quicksort like algorithm would be more efficient for determining variants\' positions.")
zerteiler.add_argument("--out", dest="out", required=False, help="Output file. Defaults to ./combined.tsv.")



args = zerteiler.parse_args()

vcf_file = open(args.vcf, 'r')

gff_file = open(args.gff, 'r')

if args.out :
    out_file = open(args.out, 'w')
else: 
    out_file = open("./combined.tsv", 'w')


genes=[]

for line in gff_file:
    components = line.split("\t")
    if len(components)>5:
        chromosome=components[0]
        designation=components[2]
        start=int(components[3])
        stop=int(components[4])
        name=components[8]

        if chromosome == 'Lp_chr6_0' and designation == 'gene' and start > 242000000 and stop < 255000000: #Change this if you want to adapt the programm.
            name=name.split(";")[0]
            name=name.split("=")[1]
            genes.append([name,start,stop])
gff_file.close()


out_file.write("Chromosome\tGene\tGene Position\tPosition\tReference\tAlternative\tFormat\tData_set\n") 



def func(genes):
    vcf_file.readline() #skip headers
    for line in vcf_file:
        components=line.split("\t")
        chromosome = components[0]
        position = int(components[1])
        if chromosome == 'Lp_chr6_0' and position > 242000000 and position < 255000000: #The vcf and gff are sorted, thus you can do all of this in one pass.
            if position > genes[0][1]: #If the variation might be in the gene
                while True:
                    if position < genes[0][2]: #If the variation is in the gene
                        reference_sequence = components[3]
                        alternative_sequence = components[4]
                        format = components[8].split(":")
                        format = format[0]+":"+format[1]+":"+format[2]
                        data_set = components[9].split(":")
                        data_set = data_set[0]+":"+data_set[1]+":"+data_set[2]
                        out_file.write(chromosome + "\t" + str(genes[0][0]) + "\t" + str(genes[0][1])+"-"+str(genes[0][2]) + "\t" + str(position) + "\t" + reference_sequence + "\t" + alternative_sequence + "\t" + format + "\t" + data_set + "\n")
                        break
                    else: #If the variation is further along than the gene.
                        genes=genes[1:] #Remove the head (The position has moved past it.)
                        if len(genes) ==0:
                            return
func(genes=genes)
vcf_file.close()



out_file.close()
print("Table created.")
print("All done.")
