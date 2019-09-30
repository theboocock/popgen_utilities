# Copyright (c) 2019 Boocock James <james.boocock@otago.ac.nz>
# Author: Boocock James <james.boocock@otago.ac.nz>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import argparse

from popgen_utils_bio import *

# Make sure these are installed
import pandas
import gffpandas.gffpandas as gffpd 
# Load the scer GFF file. 

from multiprocessing import Pool
from functools import partial

def get_ssearch36(ssearch36_in):
    pandas_in = pandas.read_csv(ssearch36_in,sep="\t", header=None)
    (pandas_in.columns) = ["gene","strain","percentage_match","3","4","5","ref_start","ref_end","start","end","E","11"]
    return(pandas_in) 

def main():
    parser = argparse.ArgumentParser(description="Generate outputs from ssearch36 output files")
    parser.add_argument("-s","--ssearch36-in", dest="ssearch_in")
    parser.add_argument("-g","--genome-in", dest="genome_in")
    parser.add_argument("-o","--output-dir",dest="output_dir")
    parser.add_argument("-g","--gff-yeast", dest="gff_yeast", default="/media/theboocock/data/Dropbox/PHDTHESIS/projects/gal_final_github/data/annotations/saccharomyces_cerevisiae.gff")
    args = parser.parse_args()
    ssearch36_in = args.ssearch_in 
    annot = args.gff_yeast
    genome_in = {os.path.basename(args.genome_in),args.genome_in}
    output_dir = args.output_dir
    out_prefix = os.path.basename(genome_in[0])
    out_new_all = out_prefix.split(".ss")[0]
    
    ssearch_df = get_ssearch36(args.ssearch_in)

    gff_yeast = gffpd.read_gff3(annot)
    gff_yeast_annot = (gff_yeast.attributes_to_columns())
    # Extract orfs that are only Verified and genes.
    gff_yeast_annot = gff_yeast_annot[(gff_yeast_annot["orf_classification"] == "Verified") & (gff_yeast_annot["type"] == "gene")]
    lengths = gff_yeast_annot["end"]- gff_yeast_annot["start"] + 1
    gff_yeast_annot["lengths"] = (lengths)

    with open(os.path.join(output_dir,out_new_all +".strict.fasta"),"w") as out_strict:
        with open(os.path.join(output_dir,out_new_all +".permissive.fasta"),"w") as out_f:
            for gene in pd_df_names[0]:
                gene = str(gene)
                if(any(gff_yeast_annot["Name"].isin([gene]))):
                    gene_length = int(gff_yeast_annot[gff_yeast_annot["Name"] == gene]["lengths"])
                    f_gene = partial(get_gene_from_ssearch36_files, gene, genome_in, gene_length)
                    ###
                    df = (genome_in.keys()[0], ssearch_df) 
                    gene_out = [f_gene(df)]
    


if __name__ == "__main__":
    main()
