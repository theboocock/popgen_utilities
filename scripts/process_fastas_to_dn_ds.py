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
from popgen_utils_bio.analysis_functions_introgressions import get_dn_ds_from_alignment
import pyfasta
from collections import OrderedDict

import os

def get_dn_ds_from_fasta(input_fasta, output_prefix, window, step):
    try:
        os.mkdir(output_prefix)
    except:
        pass

    window = int(window)
    step = int(step)

    fasta_in = pyfasta.Fasta(input_fasta)
    genes = list(fasta_in.keys())
    output_dn_ds = OrderedDict()
    if os.path.basename(input_fasta).startswith("N"):
        if "permissive" in input_fasta:
            output_file = os.path.join(output_prefix, os.path.basename(input_fasta).split(".permissive.fasta")[0] + ".permissive.dn_ds")
        else:
            output_file = os.path.join(output_prefix, os.path.basename(input_fasta).split(".strict.fasta")[0] + ".strict.dn_ds")
    else:
        output_file = os.path.join(output_prefix, os.path.basename(input_fasta).split(".fasta")[0] + ".dn_ds")
    if os.path.exists(output_file): 
        with open(output_file) as out_f:
            for line in out_f:
                line_s = line.split("\t")
                last_gene = line_s[0]
        idx = genes.index(last_gene)
    else:
    # Do the whole thing
        idx = 0
    with open(output_file, "w") as out_f:
        #with open(output_file) as out_f_old:
        #    for line in out_f_old:
        #        out_f.write(line)
        for gene in genes:
            out_ds =  get_dn_ds_from_alignment(input_fasta,these_samples=[gene],do_window=True,gene_name=gene,cbs_reference=False,window=window,step=step, hoffman=True)
            if out_ds is not None:
                rows = out_ds 
                out_f.write(str(gene) + "\tOVERALL\t" + str(rows[0][0]) + "\t" + str(rows[0][1]) + "\n")
                for row in rows[1][gene]:
                    out_f.write(str(gene) + "\tWINDOW\t" + str(row[0]) + "\t" + str(row[1]) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Calculate dn/ds from ssearch36 outputs")
    parser.add_argument("-o","--output-prefix", dest="output_prefix", help="Output prefix")
    parser.add_argument("-w","--window", dest="Window", default=200, help="Window size for dn ds")
    parser.add_argument("-s","--step", dest="step", default=10, help="Step size for dn ds")
    parser.add_argument("input_fasta", help="input fasta file")
    args = parser.parse_args()
    input_fasta_file = args.input_fasta
    get_dn_ds_from_fasta(args.input_fasta, args.output_prefix, args.window, args.step)

if __name__=="__main__":
    main()
