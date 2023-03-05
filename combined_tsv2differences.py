import sys
from argparse import ArgumentParser

#This program is designed to work with files created by vcfgff_combiner.py, also written by me.

zerteiler = ArgumentParser()

zerteiler.add_argument("--in1", dest="tsv1", required=True, help="tsv file created by vcfgff_combiner.py. Assumed 50 percent fertility.")
zerteiler.add_argument("--in2", dest="tsv2", required=True, help="tsv file created by vcfgff_combiner.py. Assumed 100 percent fertility. ")
zerteiler.add_argument("--gff", dest="gff", required=True, help="gff3 file. This needs to be sorted for the program to provide correct output. If it isn\'t then the program won\'t work through the data properly.")
zerteiler.add_argument("--out1", dest="out1", required=False, help="Output file with detailed output. Defaults to ./Variations_marked.tsv.")
zerteiler.add_argument("--out2", dest="out2", required=False, help="Output file with only genes. Formatted to fit an excel sheet. Defaults to ./Genes_marked.tsv.")



args = zerteiler.parse_args()

gff_file = open(args.gff, 'r')

if args.out1 :
    out1_file = open(args.out1, 'w')
else: 
    out1_file = open("./Variations_marked.tsv", 'w')

if args.out2 :
    out2_file = open(args.out2, 'w')
else: 
    out2_file = open("./Genes_marked.tsv", 'w')

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
            genes.append(name)
gff_file.close()


out1_file.write("Chromosome\tGene\tPosition\tStatus\n") 
out2_file.write("Gene\tStatus\n")

fifty_file = open(args.tsv1, 'r')

hundred_file = open(args.tsv2, 'r')

fifty_file.readline() #skip headers
hundred_file.readline() #skip headers

fifty_line = fifty_file.readline()
fifty_components = fifty_line.split("\t")
chromosome = fifty_components[0] #Only assigned once because of my use case.
fifty_gene = fifty_components[1]
fifty_position = int(fifty_components[3])

hundred_line = hundred_file.readline()
hundred_components = hundred_line.split("\t")
hundred_gene = hundred_components[1]
hundred_position = int(hundred_components[3])

fifty = False
hundred = False

#def fifty_iteration(): Put the repetitive code into a function if you want the code to look cleaner.




for gene in genes:
    if fifty_gene == gene:
        if hundred_gene == gene:
            fifty = True
            hundred = True
            while fifty_gene == gene or hundred_gene == gene:
                print("position " + str(fifty_position) + " or " + str(hundred_position))
                if hundred_position == fifty_position:
                    out1_file.write(chromosome+"\t"+gene+"\t"+str(fifty_position)+"\t2\n")
                    
                    fifty_line = fifty_file.readline()
                    if len(fifty_line) == 0:
                        fifty_gene = "-1" 
                    else:
                        fifty_components = fifty_line.split("\t")
                        fifty_gene = fifty_components[1]
                        fifty_position = fifty_components[3]

                    hundred_line = hundred_file.readline()
                    if len(hundred_line) == 0:
                        hundred_gene = "-2"
                    else:
                        hundred_components = hundred_line.split("\t")
                        hundred_gene = hundred_components[1]
                        hundred_position = hundred_components[3]
                elif hundred_position > fifty_position:
                    out1_file.write(chromosome+"\t"+gene+"\t"+str(fifty_position)+"\t50\n")
                    fifty_line = fifty_file.readline()
                    if len(fifty_line) == 0:
                        fifty_gene = "-1" 
                    else:
                        fifty_components = fifty_line.split("\t")
                        fifty_gene = fifty_components[1]
                        fifty_position = fifty_components[3]
                else:
                    out1_file.write(chromosome+"\t"+gene+"\t"+str(hundred_position)+"\t100\n")
                    hundred_line = hundred_file.readline()
                    if len(hundred_line) == 0:
                        hundred_gene = "-2"
                    else:
                        hundred_components = hundred_line.split("\t")
                        hundred_gene = hundred_components[1]
                        hundred_position = hundred_components[3]

        else:
            fifty = True
            out1_file.write(chromosome+"\t"+gene+"\t"+str(fifty_position)+"\t50\n")
            fifty_line = fifty_file.readline()
            if len(fifty_line) == 0:
                    fifty_gene = "-1" 
            else:
                fifty_components = fifty_line.split("\t")
                fifty_gene = fifty_components[1]
                fifty_position = fifty_components[3]
            while fifty_gene == gene:
                out1_file.write(chromosome+"\t"+gene+"\t"+str(fifty_position)+"\t50\n")
                fifty_line = fifty_file.readline()
                if len(fifty_line) == 0:
                    fifty_gene = "-1" 
                else:
                    fifty_components = fifty_line.split("\t")
                    fifty_gene = fifty_components[1]
                    fifty_position = fifty_components[3]
    elif hundred_gene == gene:
        hundred = True
        out1_file.write(chromosome+"\t"+gene+"\t"+str(hundred_position)+"\t100\n")
        hundred_line = hundred_file.readline()
        if len(hundred_line) == 0:
                hundred_gene = "-2"
        else:
            hundred_components = hundred_line.split("\t")
            hundred_gene = hundred_components[1]
            hundred_position = hundred_components[3]
        while hundred_gene == gene:
            out1_file.write(chromosome+"\t"+gene+"\t"+str(hundred_position)+"\t100\n")
            hundred_line = hundred_file.readline()
            if len(hundred_line) == 0:
                hundred_gene = "-2"
            else:
                hundred_components = hundred_line.split("\t")
                hundred_gene = hundred_components[1]
                hundred_position = hundred_components[3]
            
    if fifty and hundred:
        out2_file.write(gene+"\t2\n")
    elif fifty:
        out2_file.write(gene+"\t50\n")
    elif hundred:
        out2_file.write(gene+"\t100\n")
    else:
        out1_file.write(chromosome+"\t"+gene+"\tNA\t0\n")
        out2_file.write(gene+"\t0\n")
    fifty = False
    hundred = False
    print(gene + " done")


fifty_file.close()
hundred_file.close()


out1_file.close()
out2_file.close()
print("Tables created.")
print("All done.")
