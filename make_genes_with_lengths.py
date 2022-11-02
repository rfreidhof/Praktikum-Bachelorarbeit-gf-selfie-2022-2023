import sys
from argparse import ArgumentParser

zerteiler = ArgumentParser()

zerteiler.add_argument("--gff-file", dest="gff", required=True, help="gff file from which genes are to be extracted.")
zerteiler.add_argument("--out", dest="out", required=False, help="Output file. Defaults to ./genes_with_lengths.tsv.")


args = zerteiler.parse_args()

gff_file = open(args.gff, 'r')

if args.out :
    out_file = open(args.out, 'w')
else: 
    out_file = open("./genes_with_lengths.tsv", 'w')

out_file.write("Gene\tStart\tStop\n")

for line in gff_file:
    components = line.split("\t")
    if len(components)>5:
        chromosome=components[0]
        designation=components[2]
        start=int(components[3])
        stop=int(components[4])
        name=components[8]

        if chromosome == 'Lp_chr6_0' and designation == 'gene' and start > 242000000 and stop <255000000:
            name=name.split(";")[0]
            name=name.split("=")[1]
            out_file.write(name + "\t" + str(start) + "\t" + str(stop) + "\n" )
gff_file.close()



out_file.close()
print("All done.")

