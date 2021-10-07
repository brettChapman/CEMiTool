#!/usr/bin/python3

import os, sys
import math
import pandas as pd

from optparse import OptionParser

def main():

    parser = OptionParser()
    parser.add_option("-i", "--GFF", help = "GFF file input (required)")
    parser.add_option("-o", "--GMT", help = "GMT file output [default: gene_set.gmt]")
    (options, args) = parser.parse_args()

    if not options.GFF:
        print("No GFF file found!")
        parser.print_help()
    else:
        gff_file = options.GFF
        gmt_file = options.GMT

        gff_file = pd.read_csv(gff_file, header=None, sep='\t')

        gene_descriptions = list(gff_file[gff_file[8].str.contains(' Description=')][8].reset_index(drop=True))

        gene_set = {}
        for gene in gene_descriptions:
            gene_name = gene.split(';')[0].replace('ID=', '')
            description = gene.split(';AHRD')[1].replace(' Description=', '')
            
            if not gene_set:
                gene_set[description] = ['GFF']
                gene_set[description].append(gene_name)
            else:
                if description not in gene_set.keys():
                    gene_set[description] = ['GFF']
                    gene_set[description].append(gene_name)
                else:
                    gene_set[description].append(gene_name)
        
        df = pd.DataFrame.from_dict(gene_set, orient='index')
        df = df.reset_index(drop=False).replace([None], [''], regex=True)
        df = df.reset_index(drop=True).replace(' ', '_', regex=True)

        if gmt_file:
            df.to_csv(gmt_file, sep='\t', index=False, header=False)
        else:
            df.to_csv('gene_set.gmt', sep='\t', index=False, header=False)

if __name__ == "__main__":
    main()
