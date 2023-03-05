import sys
from argparse import ArgumentParser

#This is purely a helper script to decrease the computation time of future scripts.

zerteiler = ArgumentParser()

zerteiler.add_argument("--vcf", dest="vcf", required=True, help="vcf file. This needs to be sorted for the program to provide correct output. if it isn\'t then a quicksort like algorithm would be more efficient for determining variants\' positions.")
zerteiler.add_argument("--out", dest="out", required=False, help="Output file. Defaults to ./smaller.vcf.")



args = zerteiler.parse_args()

vcf_file = open(args.vcf, 'r')


if args.out :
    out_file = open(args.out, 'w')
else: 
    out_file = open("./smaller.vcf", 'w')


header = vcf_file.readline()

out_file.write(header) 


last_position = -1

first = True
last = True

first_line = header

for line in vcf_file:
    components=line.split("\t")
    chromosome = components[0]
    position = int(components[1])
    if chromosome == 'Lp_chr6_0' and position > 242000000 and position < 255000000 and position != last_position: #Change this to adapt the program. | Also removes duplicates, if something went wrong in the vcf.
            if first: #Make sure the first variation doesn't start in a gene.
                first = False
                if first_line != header:
                    out_file.write(first_line)
            last_position = position
            out_file.write(line)
    first_line = line
    if chromosome == 'Lp_chr6_0' and position > 255000000 and position != last_position: #This might need to be more tweaked if the last variation of the peak is the last of the chromosome
        out_file.write(line)#Make sure the last variation doesn't end in a gene.
        break            
                
vcf_file.close()
out_file.close()

print("File cut.")
print("All done.")
