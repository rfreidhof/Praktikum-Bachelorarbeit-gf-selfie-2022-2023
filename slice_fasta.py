import sys
from argparse import ArgumentParser

#I was too lazy to implement checks whether a file is actually in the fasta format.

zerteiler = ArgumentParser()

zerteiler.add_argument("--size", dest="size", required=True, help="Amount of sequences in the smaller fasta files.")
zerteiler.add_argument("--fasta", dest="fasta", required=True, help="fasta file to be sliced.")
zerteiler.add_argument("--out", dest="out", required=False, help="Output files. Defaults to ./[fasta_file]_sliced[i].fa. Where [fasta_file] is the name of the fasta file and [i] is the index of files created.")


args = zerteiler.parse_args()

fasta_file = open(args.fasta, 'r')


slices=0

if args.out :
    out_file = open(args.out+str(slices), 'w')
else: 
    out_file = open("./" + args.fasta + "_sliced" + str(slices) + ".fa", 'w')

currentnname=""
currentsequence=""
count = 0
for line in fasta_file:
    if line[0] == '>':
            out_file.write(currentnname + currentsequence)       
            currentnname=line
            currentsequence=""
            count += 1
            if count > int(args.size):
                out_file.close()
                slices += 1
                if args.out :
                    out_file = open(args.out+str(slices), 'w')
                else: 
                    out_file = open("./" + args.fasta + "_sliced"+str(slices) + ".fa", 'w')
                count = 1
    else:
        currentsequence += line
out_file.write(currentnname + currentsequence) 
fasta_file.close()
out_file.close()
print(str(slices+1) + " slices created.")
print("All done.")

