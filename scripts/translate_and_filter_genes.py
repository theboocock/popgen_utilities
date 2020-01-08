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

from popgen_utils_bio import analysis_functions_introgressions

import pyfasta
from Bio.SeqRecord import SeqRecord 
from Bio.Alphabet import IUPAC   
from Bio import Seq 

import random

def filter_genes(input_fasta, reference_hash):
    bad_genes = set()
    in_fasta = pyfasta.Fasta(input_fasta)
    for keys in in_fasta.keys():
        tmp_gene = str(in_fasta[keys]).upper()
        try:
            tmp_gene = analysis_functions_introgressions.extend_ambiguous_dna(tmp_gene)
        except:
            print(tmp_gene)
            sys.exit(1)
        ref_gene = reference_hash[keys]
        ref_gene_translate= SeqRecord(Seq.Seq(str(ref_gene).replace("-",""), alphabet=IUPAC.IUPACUnambiguousDNA()), id="REF").seq.translate(to_stop=True)
        if ((len(ref_gene_translate)+1) !=  len(ref_gene)/3):
            bad_genes.add(keys)
        for options in (tmp_gene):
            tmp_gene_translate= SeqRecord(Seq.Seq(str(options).replace("-",""), alphabet=IUPAC.IUPACUnambiguousDNA()), id="REF").seq.translate(to_stop=True)

            if ("strict" in input_fasta):
                if ((len(tmp_gene_translate)+1) !=  len(options)/3):
                    bad_genes.add(keys)
            else:
                if ((len(tmp_gene_translate)+1) != int(len(options))/3):
                    bad_genes.add(keys) 
            if(len(tmp_gene_translate) <= len(ref_gene_translate) *.9 or len(tmp_gene_translate) * .9 >= len(ref_gene_translate)):
                bad_genes.add(keys)
    with open(input_fasta + ".filt_genes","w") as output_filt:
        for bad_gene in bad_genes:
            output_filt.write(bad_gene + "\n")

def main():
    parser = argparse.ArgumentParser(description="Get lists of genes to remove because they contain early stop codons")
    parser.add_argument("input_fasta")
    args = parser.parse_args()
    input_fasta = args.input_fasta
    reference_hash = analysis_functions_introgressions._get_protein_hash()
    print(len(reference_hash.keys()))
    filter_genes(input_fasta, reference_hash)


if __name__ == "__main__":
    main()
