import sys
from argparse import ArgumentParser

#This program takes the positions of genes and gives out their sequence in a reference, and their modified sequence in a 50% and 100% fertile variant.


zerteiler = ArgumentParser()

zerteiler.add_argument("--chr6", dest="chr", required=True, help="A fasta file containing the sequence of the 6th chromosome of L. perenne.")
zerteiler.add_argument("--genes", dest="genes", required=True, help="A .tsv file containing the gene names, their start-, and their stop-positions, ordered by position")
zerteiler.add_argument("--half", dest="half", required=True, help="A .vcf file containing the SNPs of half fertile plants.")
zerteiler.add_argument("--full", dest="full", required=True, help="A .vcf file containing the SNPs of fullz fertile plants.")
zerteiler.add_argument("--out", dest="out", required=False, help="Output file. Defaults to ./sequence_variations.txt.")

#There is nothing implemented to check for the proper format, make sure you use the right files.

args = zerteiler.parse_args()






fasta_file = open(args.chr, 'r')
sequence=""
fasta_file.readline() #Skip header
for line in fasta_file:
    sequence += line.strip() #Accumulate the reference sequence.
fasta_file.close()

genes = []
genes_file = open(args.genes, 'r')
for line in genes_file:
    components = line.split("\t")
    geneName = components[0]
    geneStart = int(components[1])
    geneEnd = int(components[2].strip())
    geneSequence = sequence[geneStart-1:geneEnd-1] #These sequences are not identical to the ones given in the transcript files.
    genes.append((geneName, geneStart, geneEnd, geneSequence))
genes_file.close()



if args.out :
    out_file = open(args.out, 'w')
else: 
    out_file = open("./sequence_variations.txt", 'w')

half_file = open(args.half, "r")
full_file = open(args.full, "r")
half_file.readline() #Skip headers
full_file.readline()

half_components = half_file.readline().split('\t')
half_position = int(half_components[1])
half_reference = half_components[3]
half_alternative = half_components[4]

full_components = full_file.readline().split('\t')
full_position = int(full_components[1])
full_reference = full_components[3]
full_alternative = full_components[4]



for (geneName, geneStart, geneEnd, geneSequence) in genes:
    half_sequence = geneSequence    
    half_difference = 0
    while half_position < geneStart: #Skip to relevant SNPs
        half_components = half_file.readline().split('\t')
        half_position = int(half_components[1])
        half_reference = half_components[3]
        half_alternative = half_components[4]
    while half_position >= geneStart and half_position <= geneEnd:
        if len(half_reference) == len(half_alternative):
            half_sequence = half_sequence[:(half_position-1)-(geneStart-1)-half_difference] + half_alternative + half_sequence[(half_position-1)-(geneStart-1)+1-half_difference:]
        if len(half_reference) > len(half_alternative): #If a deletion occured:
            difference = len(half_reference) - len(half_alternative) #A deletion always leaves the first part of the reference. 
            half_sequence = half_sequence[:(half_position-1)-(geneStart-1)-half_difference] + half_alternative + half_sequence[(half_position-1)-(geneStart-1)+1+difference-half_difference:]
            half_difference += difference #Scope moving happens afterwards as to not modify how much of the head and tail is taken.
        if len(half_reference) < len(half_alternative): #If an insertion occured:
            difference = len(half_reference) - len(half_alternative)
            half_sequence = half_sequence[:(half_position-1)-(geneStart-1)-half_difference] + half_alternative + half_sequence[(half_position-1)-(geneStart-1)+1-half_difference:] 
            half_difference += difference #Scope moving happens afterwards as to not modify how much of the head and tail is taken.
        half_components = half_file.readline().split('\t')
        half_position = int(half_components[1])
        half_reference = half_components[3]
        half_alternative = half_components[4]

    full_sequence = geneSequence #This could probably be put into a function instead of doing basically the same thing twice.
    full_difference = 0
    while full_position < geneStart: #Skip to relevant SNPs
        full_components = full_file.readline().split('\t')
        full_position = int(full_components[1])
        full_reference = full_components[3]
        full_alternative = full_components[4]
    while full_position >= geneStart and full_position <= geneEnd:
        if len(full_reference) == len(full_alternative):
            full_sequence = full_sequence[:(full_position-1)-(geneStart-1)-full_difference] + full_alternative + full_sequence[(full_position-1)-(geneStart-1)+1-full_difference:]
        if len(full_reference) > len(full_alternative): #If a deletion occured:
            difference = len(full_reference) - len(full_alternative) #A deletion always leaves the first part of the reference. 
            full_sequence = full_sequence[:(full_position-1)-(geneStart-1)-full_difference] + full_alternative + full_sequence[(full_position-1)-(geneStart-1)+1+difference-full_difference:]
            full_difference += difference #Scope moving happens afterwards as to not modify how much of the head and tail is taken.
        if len(full_reference) < len(full_alternative): #If an insertion occured:
            difference = len(full_reference) - len(full_alternative)
            full_sequence = full_sequence[:(full_position-1)-(geneStart-1)-full_difference] + full_alternative + full_sequence[(full_position-1)-(geneStart-1)+1-full_difference:] 
            full_difference += difference #Scope moving happens afterwards as to not modify how much of the head and tail is taken.
        full_components = full_file.readline().split('\t')
        full_position = int(full_components[1])
        full_reference = full_components[3]
        full_alternative = full_components[4]

    out_file.write(geneName + "\n")
    out_file.write("Reference: \t" + geneSequence + "\n")
    out_file.write("Half fertile: \t" + full_sequence + "\n")
    out_file.write("Fully fertile: \t" + full_sequence + "\n\n")






out_file.close()
print("All done.")
