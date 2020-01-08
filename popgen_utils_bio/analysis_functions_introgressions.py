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

def _genes_to_names():
    """
        Convert yeast geneids to names

        file:
        geneid genename
    """
    gene_names = {}
    with open("/media/theboocock/data//PHDTHESIS/Single-CellRNASEQ/eqtls/ref_data/yeast/genes/gene_names.txt") as gene_input:
        for line in gene_input:
            line = line.strip().split(" ")
            gene = line[0]
            alias = line[1]
            gene_names[gene] = alias
    return(gene_names)

import os
import pandas

def get_gene_from_ssearch36_perm_files(gene,genomes_to_ssearch, gene_length, dfs,debug=False):
        get_gene_from_ssearch36_perm_files(gene,genomes_to_ssearch,gene_length,dfs,permissive=True)


def get_gene_from_ssearch36_script(fasta_file, ssearch36_in, gene, gene_length):
    """
        @author James Boocock
        @date 30 sept 2019
    """
    row_gal_gene = (ssearch36_in[ssearch36_in["gene"] == gene])
    fasta_file = fasta_file
    genome = os.path.basename(fasta_file)
    gene_strict = _get_gene_from_ssearch36_files(genome, row_gal_gene,fasta_file,gene_length,gene)
    gene_permissive = _get_gene_from_ssearch36_files(genome, row_gal_gene,fasta_file,gene_length,gene, permissive=True)
    return(gene_strict, gene_permissive)

def get_gene_from_ssearch36_files(gene,genomes_to_ssearch, gene_length, dfs,debug=False,permissive=False):
        genome = dfs[0]
        df = dfs[1]
        GENOMES_DIR= "/media/theboocock/data/PHDTHESIS/projects/gal_final_github/data/genomes//"
        try:
            (dfs[2])
            fasta_file= dfs[2] 
        except:
            fasta_file = os.path.join(GENOMES_DIR, os.path.basename(genomes_to_ssearch[genome]).split(".ss")[0])
        row_gal_gene = (df[df["gene"] == gene])
        old_gal_gene = row_gal_gene
        if debug:
            print(row_gal_gene)
        #print(row_gal_gene)
        #print(row_gal_gene)
        _get_gene_from_ssearch36_files(genome, row_gal_gene,fasta_file,gene_length,gene)

def _get_gene_from_ssearch36_files(genome, row_gal_gene,fasta_file,gene_length,gene,debug=False,permissive=False, alignment_fraction=0.9):
        rev_strand =False
        if (row_gal_gene.shape[0] == 0):
            return(genome, None)
        row_gal_gene = row_gal_gene.sort_values(by="11", ascending=False).head(n=1)
        #length_first = row_gal_gene["end"].iloc[0]- row_gal_gene["start"].loc[0] + 1
        #if length_first < gene_length * 0.9:
        #    row_gal_gene = row_gal_gene.sort_values(by="3", ascending=False).head(n=1)
        #else:
        #    row_gal_gene = row_gal_gene.head(n=1)
        start = int(row_gal_gene["start"])
        end = int(row_gal_gene["end"])
        ref_start = int(row_gal_gene["ref_start"])
        ref_end = int(row_gal_gene["ref_end"])
        length = gene_length
        #length = ref_end - ref_start  + 1
        # Make sure that the length is really 
        # TODO: make sure that the length is really 
        # abs(end-start) + 1 it might not be .
        #print(row_gal_gene)
        contig = str(row_gal_gene["strain"].astype(str)).split()[1]
        length_alignment = (abs(end -start) + 1)
        #print(length_alignment)
        if (float(length) * alignment_fraction) <= length_alignment:
            old_start = start
            start_offset = ref_start - 1
            # exclusive
            end_offset = length - ref_end
            #print(start_offset)
            #print(end_offset)
            if end > start and ref_end < ref_start:
                rev_strand = True
                tmp = ref_end
                ref_end = ref_start 
                ref_start = tmp
                start_offset = ref_start - 1
                # exclusive
                end_offset = length - ref_end
        #        print(end_offset)
                start = start -end_offset 
                end = end + start_offset+ 1
                length = abs(start - end) 
            else:
             #   print(ref_start)
                start_offset = ref_start - 1
             #   print(start_offset)
                # exclusive
                end_offset = length - ref_end
                start = start - start_offset
                end = end + end_offset  + 1
                length = abs(start - end) 
          #      print(start)
         #   #print(gene)
            #print(length)
            # Get this weird alignment working...
            if gene == "YKL127W" and length == 1715:
                if rev_strand:
                    seq = (extract_gene_from_alignment(fasta_file, contig,old_start,length-2)).decode("utf-8").strip()
                else:
                    seq = (extract_gene_from_alignment(fasta_file, contig,old_start-1,length-2)).decode("utf-8").strip()
                #print( (extract_gene_from_alignment(fasta_file, contig,start,length-2)).decode("utf-8").strip())
                ### NEED to fix the start and end because something fucked happened.
                ### NEED to fix the start and end because something fucked happened.
            else:
                seq = (extract_gene_from_alignment(fasta_file, contig,start,length)).decode("utf-8").strip()
            #seq = (extract_gene_from_alignment(fasta_file, contig,start,length)).decode("utf-8").strip()
            os.sync()
            if(debug):
                print(start_offset)
                print(end_offset)
                print(start)
                print(row_gal_gene)
                print(" ",(seq))
            if rev_strand:
                seq = revcomp(seq)
            if ("N" in seq):
                return(genome, None)
            if (seq.startswith("ATG") and (seq.endswith("TAA") or seq.endswith("TGA") or seq.endswith("TAG"))):
                return(genome, seq)
            elif permissive:
                return(genome, seq)
            else:
                return(genome, None)
        return(genome, None)


import subprocess
def extract_gene_from_alignment(fasta_query,contig,start,length):
    cmd="""cat {0} | /u/home/s/smilefre/bin/bioawk -c fastx '{{if($name=="{1}"){{print substr($seq,{2},{3}) }}}}' """.format(fasta_query,contig, start,length) 
    output_gene = subprocess.check_output(cmd, shell=True)
    return(output_gene) 
def get_gene_alias(gene, gene_names):
    """
        Check if alias is present if not return the geneid.
    
    """
    # Return alias if I can identify it. If I can't just return
    try:
        return(gene_names[gene])
    except:
        return(gene)

def revcomp(string):
    """
        Reverse complement string.
    """
    new_string = []
    #print(string)
    for char in list(string):
       # print(char)
        if(char == "A"):
            new_string.append("T")
        elif(char == "G"):
            new_string.append("C")
        elif(char == "T"):
            new_string.append("A")
        elif(char == "C"):
            new_string.append("G")
        else:
            new_string.append(char)
    return("".join(new_string[::-1]))       
                    
import tempfile
import pysam
import os
import pandas

def extract_gene_name(string):
    try:
        tmp_string = string.split(" ")
    except:
        return("NA")
    gene_name = [tmp_s for tmp_s in tmp_string if tmp_s.startswith("Y")]
    try:
        gene_name = gene_name[0]
    except:
        gene_name = "NA"
    return(gene_name.upper())

def get_ygob_genes(og_name, genes, input_fasta,input_tab):

    pd_tab = pandas.read_csv(input_tab, sep="\t",header=None)
    pd_tab["GENE_NAME"] = pd_tab[8].apply(extract_gene_name) 
    fast_input = pyfasta.Fasta(input_fasta)
    with open("data/cgla/"+ og_name,"w") as f:
        for gene in genes:
            if gene == "YBR020W":
                lookup = "YDR009W"
                pdtab = pd_tab["GENE_NAME"].str.contains(lookup)
            else:
                pdtab = pd_tab["GENE_NAME"].str.contains(gene)
            try:
                get_gene_row = pd_tab[pdtab].head(n=1)
            except:
                continue
            df_length = (len(get_gene_row.index))
            if df_length == 0:
                continue
            chrom = get_gene_row.iloc[0][5]
            key_fasta = "knag_Chr_" + str(chrom)
            for key in  (list(fast_input.keys())):
                if key_fasta in key:
                    key_chrom = key
            start = get_gene_row.iloc[0][2]
            end = get_gene_row.iloc[0][3]
            strand = get_gene_row.iloc[0][1]
            if strand == 1:
                sequence =  (fast_input[key_chrom][(start-1):end])
            else:
                sequence = (revcomp(fast_input[key_chrom][(start-1):end]))

            
            f.write(">" + gene + "\n" + sequence + "\n")
def get_ygob_genes(og_name, genes, input_fasta,input_tab):

    pd_tab = pandas.read_csv(input_tab, sep="\t",header=None)
    pd_tab["GENE_NAME"] = pd_tab[8].apply(extract_gene_name) 
    fast_input = pyfasta.Fasta(input_fasta)
    with open("data/cgla/"+ og_name,"w") as f:
        for gene in genes:
            if gene == "YBR020W":
                lookup = "YDR009W"
                pdtab = pd_tab["GENE_NAME"].str.contains(lookup)
            else:
                pdtab = pd_tab["GENE_NAME"].str.contains(gene)
            try:
                get_gene_row = pd_tab[pdtab].head(n=1)
            except:
                continue
            df_length = (len(get_gene_row.index))
            if df_length == 0:
                continue
            chrom = get_gene_row.iloc[0][5]
            key_fasta = "knag_Chr_" + str(chrom)
            for key in  (list(fast_input.keys())):
                if key_fasta in key:
                    key_chrom = key
            start = get_gene_row.iloc[0][2]
            end = get_gene_row.iloc[0][3]
            strand = get_gene_row.iloc[0][1]
            if strand == 1:
                sequence =  (fast_input[key_chrom][(start-1):end])
            else:
                sequence = (revcomp(fast_input[key_chrom][(start-1):end]))

            
            f.write(">" + gene + "\n" + sequence + "\n")
        
        #pdtab = pd_tab.query("GENE_NAME == {0}".format(gene))

def get_gene_fasta_alignment(genes, input_strains, output_dir = "data/fasta_output", window =0,
                            promoter = False):
        tmp_out_gene = tempfile.NamedTemporaryFile(delete=False)
        cmd="""cat {0} | /u/home/s/smilefre/bin/bioawk -c fastx '{{if($name=="{1}"){{print ">"$name"\\n"$seq}}}}' """.format(fasta_query,gene)
        output_gene = subprocess.check_output(cmd, shell=True)
        tmp_out_gene.write(output_gene)
        tmp_out_gene.close()
        blast_search = """blastn -outfmt "6 qseqid sseqid pident qlen length mismatch gapope evalue bitscore sstart send sstrand" -db scaffold_blast -query {0} -num_alignments 15 -word_size 10""".format(tmp_out_gene.name)
        blast_output = (subprocess.check_output(blast_search, shell=True))
        species_list = []
        with open(output_dir+ "/" + gene + ".fasta","w") as gene_fasta:
            for line in blast_output.splitlines():
                line = line.decode()
                line_s = line.split("\t")
                #print(line_s[1])
                species_match = line_s[1].split("_")[1].split(".")[0]
                kafr = line_s[1].split("_")[0]
                if "CBS288A" in line_s[1]:
                    #print("HERE")
                    species_match = "CBS288A"
                    #print("")
                if "kafr" in kafr:
                    species_match = "kafr"
                #print(line)
                #print(species_match)
                if species_match in input_strains and (species_match not in species_list):
                    #print(line_s)
                    subject_reference=line_s[1]
                    #print(all_scaffolds_fasta.references)
                    subject_start = int(line_s[8])
                    subject_end = int(line_s[9])
                    #print(len(all_scaffolds_fasta.fetch(subject_reference)))
                    #print(all_scaf_seq)
                    species_list.append(species_match)
                    #print(len(all_scaf_seq))
                    #subject_reference = subject_reference.decode()
                    gene_fasta.write(">"+gene+species_match+"\n")
                    if line_s[10] == "plus":
                        all_scaf_seq = all_scaffolds_fasta.fetch(subject_reference)[(subject_start-1 - window):(subject_end + window)]
                        gene_fasta.write(all_scaf_seq+"\n")
                    else:
                        all_scaf_seq = all_scaffolds_fasta.fetch(subject_reference)[(subject_end-1 - window):(subject_start + window)]
                        gene_fasta.write(revcomp(all_scaf_seq)+"\n")

def create_fasta_outputs_from_annotations(genes = ["YLR081W","YBR018C","YBR019C","YBR020W"], 
                        input_strains=["Smik","Spar","Port","Sbay","Scer","CBS288A"], 
                        output_dir = "data/fasta_output", out_group_fasta = None,
                        fasta_amino="/media/theboocock/data//PHDTHESIS/Single-CellRNASEQ/eqtls/ref_data/yeast/intro/yeast_senso/all_genes_with_gal.fasta"):
    
    try:
        os.makedirs(output_dir)
    except:
        pass
    all_scaffolds_fasta = pyfasta.Fasta(fasta_amino)

    all_scaffolds_hash = {} 
    if out_group_fasta is not None:
        out_group_fasta_file = pyfasta.Fasta(out_group_fasta)
    for key, seq in all_scaffolds_fasta.items():
        if key.startswith("Y"):
            # must be CBS
            all_scaffolds_hash[key.strip()] = str(seq) 
        scer_column = key.split(" ")
        sac_cer_column = [scer_c for scer_c in scer_column if "SCER:" in scer_c]
        try:
            gene =  sac_cer_column[0].split(":")[1].split(";")[0]
            species = scer_column[0]
            #species = (species.split("_")[0])
            all_scaffolds_hash[species] = str(seq) 
        except:
            continue
    OUT_BLAST="blast_out"
    blast_cmd="makeblastdb -in {0} -dbtype 'nucl' -out {1}".format(fasta_amino, OUT_BLAST)
    subprocess.check_call(blast_cmd,shell=True)
    for gene in genes:
        gene_og = None
        if out_group_fasta is not None:
            try:
                gene_og = out_group_fasta_file[gene]
                print(gene_og)
            except:
                continue
        with open(output_dir+ "/" + gene + ".fasta","w") as gene_fasta:
            if gene_og is not None:
                gene_fasta.write(">" + gene +"_og\n" + str(gene_og)+ "\n")
            with open("blast_search.tmp","w") as blast_search:
                tmp_fasta = all_scaffolds_hash[gene]
                blast_search.write(">" + gene+"\n")
                blast_search.write(tmp_fasta + "\n")


            blast_search = """blastn -outfmt "6 qseqid sseqid pident qlen length mismatch gapope evalue bitscore sstart send sstrand" -db {0} -query {1} -num_alignments 15 -word_size 10""".format(OUT_BLAST,"blast_search.tmp")
            
            blast_output = (subprocess.check_output(blast_search, shell=True))
            print(gene)
            species_not_found = []
            for line in blast_output.splitlines():
                line = line.decode()
                line_s = line.split("\t")
                species = line_s[1]
                species = (species.split("_")[0])
                #print(line)
                if species.startswith("Y"):
                    species = "CBS288A"
                if species not in species_not_found and species in input_strains:
                    species_not_found.append(species)
                    sequence = all_scaffolds_hash[line_s[1]]
                    gene_fasta.write(">" + species + "_" + gene +"\n")
                    gene_fasta.write(sequence+"\n")


def create_fasta_outputs(genes = ["YLR081W","YBR018C","YBR019C","YBR020W"], 
                        input_strains=["Smik","Spar","Port","Sbay","Scer","CBS288A"], 
                        output_dir = "data/fasta_output",
                        window=0, 
                        promoter=False, out_group_fasta=None, 
                        fasta_search="/media/theboocock/data//PHDTHESIS/Single-CellRNASEQ/eqtls/ref_data/yeast/intro/sense_pgm1/all_scaffolds.fa", 
                        fasta_query="/media/theboocock/data//PHDTHESIS/Single-CellRNASEQ/eqtls/ref_data/yeast/genes/orf_coding_all.fasta",
                        check_start_and_end=False):
    OUT_BLAST="scaffold_blast"
    blast_cmd="makeblastdb -in {0} -dbtype 'nucl' -out scaffold_blast".format(fasta_search)
    
    try:
        os.makedirs(output_dir)
    except:
        pass
    all_scaffolds_fasta = pysam.FastaFile(fasta_search)
    subprocess.check_call(blast_cmd,shell=True)
    if out_group_fasta is not None:
        out_group_fasta_file = pyfasta.Fasta(out_group_fasta)
    for gene in genes:
        tmp_out_gene = tempfile.NamedTemporaryFile(delete=False)
        cmd="""cat {0} | /u/home/s/smilefre/bioawk -c fastx '{{if($name=="{1}"){{print ">"$name"\\n"$seq}}}}' """.format(fasta_query,gene)
        output_gene = subprocess.check_output(cmd, shell=True)
        tmp_out_gene.write(output_gene)
        tmp_out_gene.close()
        blast_search = """blastn -outfmt "6 qseqid sseqid pident qlen length mismatch gapope evalue bitscore sstart send sstrand" -db scaffold_blast -query {0} -num_alignments 15 -word_size 10""".format(tmp_out_gene.name)
        blast_output = (subprocess.check_output(blast_search, shell=True))
        species_list = []
        gene_og = None
        if out_group_fasta is not None:
            try:
                gene_og = out_group_fasta_file[gene]
            except:
                continue
        print("here")
        print(gene)
        with open(output_dir+ "/" + gene + ".fasta","w") as gene_fasta:
            if gene_og is not None:
                gene_fasta.write(">" + gene +"_og\n" + str(gene_og)+ "\n")
            for line in blast_output.splitlines():
                line = line.decode()
                line_s = line.split("\t")
                #print(line_s[1])
                species_match = line_s[1].split("_")[1].split(".")[0]
                kafr = line_s[1].split("_")[0]
                if "CBS288A" in line_s[1]:
                    #print("HERE")
                    species_match = "CBS288A"
                    #print("")
                if "kafr" in kafr:
                    species_match = "kafr"
                if "cgla" in kafr:
                    species_match = "cgla"
                #print(line)

                if species_match in input_strains and (species_match not in species_list):
                    #print(line_s)
                    subject_reference=line_s[1]
                    #print(all_scaffolds_fasta.references)
                    subject_start = int(line_s[8])
                    subject_end = int(line_s[9])
                    #print(len(all_scaffolds_fasta.fetch(subject_reference)))
                    #print(all_scaf_seq)
                    species_list.append(species_match)
                    #print(len(all_scaf_seq))
                    #subject_reference = subject_reference.decode()
                    gene_fasta.write(">"+gene+"_"+species_match+"\n")
                    if line_s[10] == "plus":
                        all_scaf_seq = all_scaffolds_fasta.fetch(subject_reference)[(subject_start-1 - window):(subject_end + window)]
                        gene_fasta.write(all_scaf_seq+"\n")
                    else:
                        all_scaf_seq = all_scaffolds_fasta.fetch(subject_reference)[(subject_end-1 - window):(subject_start + window)]
                        gene_fasta.write(revcomp(all_scaf_seq)+"\n")

import pyfasta

from Bio.Alphabet import IUPAC
from Bio.Data.IUPACData  import ambiguous_dna_values 
#from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import codonalign
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Align.Applications import MuscleCommandline
import subprocess
from Bio import SeqIO, AlignIO
import pyfasta
from Bio.Phylo.TreeConstruction import DistanceCalculator

from multiprocessing import Pool
from functools import partial


from Bio import Seq
from itertools import product
import numpy as np

def extend_ambiguous_dna(seq):
   """return list of all possible sequences given an ambiguous DNA input"""
   d = Seq.IUPAC.IUPACData.ambiguous_dna_values
   return  list(map("".join, product(*map(d.get, seq)))) 

import random

def get_dn_ds_from_seq_no_align(genes,seqs, bootstrap,bam_genes_cbs=None,method="NG86",all_genes=False):
    #gal_genes_cbs = pyfasta.Fasta(fasta_input)
    #if(all_genes):
    #    genes = gal_genes_cbs.keys()
    dn_ds = {}
    #name = seq[0]
    #seq= seq[1]
    s288c_protein_hash = _get_protein_hash()
    #for gene in genes:
    prot1_l = []
    prot2_l = []
    dna1_l = []
    dna2_l =[]
     
    for gene in genes:
        name = gene
        #print(gene)
        seq =str(seqs[gene])
        if bam_genes_cbs is not None:
#            print("HERE")
            by_sequence = bam_genes_cbs[gene]
        else:
            by_sequence = s288c_protein_hash[gene]
        #print(len(by_sequence))
        seq_list = (extend_ambiguous_dna(seq))
        seq = random.choice(seq_list)
        seq_cbs = SeqRecord(Seq.Seq(seq, alphabet=IUPAC.IUPACUnambiguousDNA()), id=name)
        seq_by = SeqRecord(Seq.Seq(by_sequence,alphabet=IUPAC.IUPACUnambiguousDNA()),id='BY')
        prot_cbs = seq_cbs.seq.translate(to_stop=True)
        prot_by = seq_by.seq.translate(to_stop=True)
        in_file_prot = "tmp3/"+ name + ".in.prot.fasta"
        with open(in_file_prot, "w") as tmp_fasta:
            tmp_fasta.write(">"+name+"\n")
            tmp_fasta.write(str(prot_cbs) + "\n")
            tmp_fasta.write(">BY\n")
            tmp_fasta.write(str(prot_by))
        out_file_prot= "tmp3/" + name + ".out.prot.fasta" 
        muscle_cline = ClustalOmegaCommandline(infile=in_file_prot, outfile=out_file_prot,force=True,outfmt="clu")
        subprocess.check_call(str(muscle_cline),shell=True)
        proteins = AlignIO.read(out_file_prot,'clustal',alphabet=IUPAC.protein)
        prot_start1 = str(proteins[0].seq)
        prot_start2 = str(proteins[1].seq)
        dna_start1 = str(seq_cbs.seq)
        dna_start2 = str(seq_by.seq)
        prot_start1_idx = [ i for i, x in enumerate(prot_start1) if x == "-" ]
        #print(prot_start1_idx)
        prot_start2_idx = [ i for i, x in enumerate(prot_start2) if x == "-" ] 
        #print(prot_start2_idx)
        # Identify the indels....
        prot_start1_idx = [ i for i, x in enumerate(prot_start1) if x != "-" ]
        prot_start2_idx = [ i for i, x in enumerate(prot_start2) if x != "-" ] 
        prot_start = [ i for i, (x, y) in enumerate(zip(prot_start1,prot_start2)) if x != "-" and y !="-" ] 
        prot_start1 = "".join([ prot_start1[i] for i in prot_start])
        prot_start2 = "".join([ prot_start2[i] for i in prot_start])
        dna_start1 = "".join([dna_start1[(x*3):(x*3 + 3)] for x in prot_start2_idx]) 
        dna_start2 = "".join([dna_start2[(x*3):(x*3 + 3)] for x in prot_start1_idx])
        if dna_start1.endswith("TAG") or dna_start1.endswith("TAA") or dna_start1.endswith("TGA"):
            dna_start1 = dna_start1[:(len(dna_start1)-3)]
        if dna_start2.endswith("TAG") or dna_start2.endswith("TAA") or dna_start2.endswith("TGA"):
            dna_start2 = dna_start2[:(len(dna_start2)-3)]
        prot1_l.append(prot_start1)
        prot2_l.append(prot_start2)
        dna1_l.append(dna_start1)
        #print("DNA")
        #print(dna_start1)
        #print(len(dna_start1))
        #print(len(prot_start1))
        dna2_l.append(dna_start2)
        #in_file = "tmp3/"+ name + ".in.fasta"
        #with open(in_file, "w") as tmp_fasta:
        #    tmp_fasta.write(">"+name+"\n")
        #    tmp_fasta.write(str(seq_cbs.seq) + "\n")
        #    tmp_fasta.write(">BY\n")
        #    tmp_fasta.write(str(seq_by.seq))
        #out_file = "tmp3/" +name +".out.fasta" 
        #muscle_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file,force=True,outfmt="clu")
        #subprocess.check_call(str(muscle_cline),shell=True)
        #dna = AlignIO.read(out_file,'clustal',alphabet=IUPAC.IUPACUnambiguousDNA())
    #print(len(prot_start1))
    prot_start1 = "".join(prot1_l)
    prot_start2 = "".join(prot2_l)
    #print(prot_start1)
    #print(len(prot_start1))
    dna_start1= "".join(dna1_l)
    #print(dna_start1)
    dna_start2= "".join(dna2_l)
    #print(len(dna_start1))
    seq_cbs.seq = Seq.Seq(dna_start1, alphabet=IUPAC.IUPACUnambiguousDNA()) 
    #print#(seq_cbs)
    seq_by.seq = Seq.Seq(dna_start2, alphabet=IUPAC.IUPACUnambiguousDNA()) 
    if(bootstrap):
        max_dn_ds = 3.23
        dn_array = [] 
        #prot_start1 = str(proteins[0].seq)
        #prot_start2 = str(proteins[1].seq)
        #dna_start1 = str(dna[0].seq)
        #dna_start2 = str(dna[1].seq)
        # Identify the indels....
        #prot_start = [ i for i, (x, y) in enumerate(zip(prot_start1,prot_start2)) if x != "-" and y !="-" ] 
        #prot_start1 = "".join([ prot_start1[i] for i in prot_start])
        #prot_start2 = "".join([ prot_start2[i] for i in prot_start])
        prot_length = len(prot_start1) 
        #print(len(prot_start1))
        #if ("-" in prot_start1):
        #    print(prot_start.index("-"))
        #dna_start1 =  "".join([ dna_start1[(i*3):((i*3)+3)] for i in prot_start])
        #print(dna_start1)
        #dna_start2 =  "".join([ dna_start2[(i*3):((i*3)+3)] for i in prot_start])
        #print(len(dna_start2))
        #prot_start = [  for x, y in zip(prot_start1,prot_start2) if x != "-" and y !="-" ] 
        for i in range(200):
            idxs = np.random.choice(prot_length, prot_length)
            #proteins[0].seq[idxs]
            #print(idxs)
            prot1 = "".join([prot_start1[i] for i in idxs])
            prot2 = "".join([prot_start2[i] for i in idxs])
            proteins[0].seq = Seq.Seq(prot1, alphabet=IUPAC.IUPACProtein()) 
            proteins[1].seq = Seq.Seq(prot2, alphabet=IUPAC.IUPACProtein()) 
            # bootstrap DNA.
            seq_by.seq = Seq.Seq( "".join([ dna_start2[(i*3):((i*3)+3)] for i in idxs ] ),alphabet=IUPAC.IUPACUnambiguousDNA()) 
            seq_cbs.seq = Seq.Seq( "".join([ dna_start1[(i*3):((i*3)+3)] for i in idxs ]), alphabet=IUPAC.IUPACUnambiguousDNA()) 
            #print(proteins)
            codon_aln = codonalign.build(proteins,[seq_cbs, seq_by], alphabet=codonalign.default_codon_alphabet)
            #dna_alignment = AlignIO.read(ouput_file,"clustal")
            #calculator = DistanceCalculator('blosum62')
            #dm = calculator.get_distance(dna_alignment)
            dn_ds= (codonalign.codonseq.cal_dn_ds(codon_aln[0],codon_aln[1]))
            #print(dn_ds)
            if(dn_ds[1] == -1):
                dn_array.append((max_dn_ds))
            else:
                dn_array.append((dn_ds[1])) 
        #print(np.var(dn_array))
        #print(np.median(dn_array))
        return(dn_array)
    else:
        proteins[0].seq = Seq.Seq(prot_start1, alphabet=IUPAC.IUPACProtein()) 
        proteins[1].seq = Seq.Seq(prot_start2, alphabet=IUPAC.IUPACProtein()) 
        #ouput_file = out_file 
        codon_aln = codonalign.build(proteins,[seq_cbs, seq_by], alphabet=codonalign.default_codon_alphabet)
        #seq_by.seq = Seq.Seq( "".join([ dna_start2[(,alphabet=IUPAC.IUPACUnambiguousDNA()) 
        #seq_cbs.seq = Seq.Seq( "".join([ dna_start1[(i*3):((i*3)+3)] for i in idxs ]), alphabet=IUPAC.IUPACUnambiguousDNA()) 
        #dna_alignment = AlignIO.read(ouput_file,"clustal")
        calculator = DistanceCalculator('blosum62')
        #dm = calculator.get_distance(dna_alignment)
        dn_ds= (codonalign.codonseq.cal_dn_ds(codon_aln[0],codon_aln[1]))
        return([name, dn_ds[0], dn_ds[1]])
#def get_dn_ds_rom_seq(gene,bootstrap,seq,method="NG86",all_genes=False, skip_alignment=False):
#    #gal_genes_cbs = pyfasta.Fasta(fasta_input)
#    #if(all_genes):
#    #    genes = gal_genes_cbs.keys()
#    dn_ds = {}
#    name = seq[0]
#    seq= seq[1]
#    s288c_protein_hash = _get_protein_hash()
#    #for gene in genes:
#    #name = gene
#    seq =str(seq)
#    by_sequence = s288c_protein_hash[gene]
#    seq_list = (extend_ambiguous_dna(seq))
#    seq = random.choice(seq_list)
#    seq_cbs = SeqRecord(Seq.Seq(seq, alphabet=IUPAC.IUPACUnambiguousDNA()), id=name)
#    seq_by = SeqRecord(Seq.Seq(by_sequence,alphabet=IUPAC.IUPACUnambiguousDNA()),id='BY')
#    prot_cbs = seq_cbs.seq.translate(to_stop=True)
#    prot_by = seq_by.seq.translate(to_stop=True)
#    print(prot_cbs)
#    #print(prot_by)
#    in_file_prot = "tmp3/"+ name + ".in.prot.fasta"
#    with open(in_file_prot, "w") as tmp_fasta:
#        tmp_fasta.write(">"+name+"\n")
#        tmp_fasta.write(str(prot_cbs) + "\n")
#        tmp_fasta.write(">BY\n")
#        tmp_fasta.write(str(prot_by))
#    out_file_prot= "tmp3/" + name + ".out.prot.fasta" 
#    muscle_cline = ClustalOmegaCommandline(infile=in_file_prot, outfile=out_file_prot,force=True,outfmt="clu")
#    subprocess.check_call(str(muscle_cline),shell=True)
#    proteins = AlignIO.read(out_file_prot,'clustal',alphabet=IUPAC.protein)
#    prot_start1 = str(proteins[0].seq)
#    prot_start2 = str(proteins[1].seq)
#    dna_start1 = str(seq_cbs.seq)
#    dna_start2 = str(seq_by.seq)
#    #print(len(dna_start1))
#    #print(len(dna_start2))
#    # Identify the indels....
#    prot_start1_idx = [ i for i, x in enumerate(prot_start1) if x != "-" ]
#    prot_start2_idx = [ i for i, x in enumerate(prot_start2) if x != "-" ] 
#    prot_start = [ i for i, (x, y) in enumerate(zip(prot_start1,prot_start2)) if x != "-" and y !="-" ] 
#    prot_start1 = "".join([ prot_start1[i] for i in prot_start])
#    prot_start2 = "".join([ prot_start2[i] for i in prot_start])
#    dna_start1 = "".join([dna_start1[(x*3):(x*3 + 3)] for x in prot_start2_idx]) 
#    dna_start2 = "".join([dna_start2[(x*3):(x*3 + 3)] for x in prot_start1_idx])
#    seq_cbs.seq = Seq.Seq(dna_start1, alphabet=IUPAC.IUPACUnambiguousDNA()) 
#    seq_by.seq = Seq.Seq(dna_start2, alphabet=IUPAC.IUPACUnambiguousDNA()) 
#    if dna_start1.endswith("TAG") or dna_start1.endswith("TAA") or dna_start1.endswith("TGA"):
#        dna_start1 = dna_start1[:(len(dna_start1)-3)]
#    if dna_start2.endswith("TAG") or dna_start2.endswith("TAA") or dna_start2.endswith("TGA"):
#    #    print("HER")
#        dna_start2 = dna_start2[:(len(dna_start2)-3)]
#    in_file = "tmp3/"+ name + ".in.fasta"
#    with open(in_file, "w") as tmp_fasta:
#        tmp_fasta.write(">"+name+"\n")
#        tmp_fasta.write(str(seq_cbs.seq) + "\n")
#        tmp_fasta.write(">BY\n")
#        tmp_fasta.write(str(seq_by.seq))
#    out_file = "tmp3/" +name +".out.fasta" 
#    muscle_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file,force=True,outfmt="clu")
#    subprocess.check_call(str(muscle_cline),shell=True)
#    dna = AlignIO.read(out_file,'clustal',alphabet=IUPAC.IUPACUnambiguousDNA())
#    if(bootstrap):
#        max_dn_ds = 3.23
#        dn_array = [] 
#        #prot_start1 = str(proteins[0].seq)
#        #prot_start2 = str(proteins[1].seq)
#        #dna_start1 = str(dna[0].seq)
#        #dna_start2 = str(dna[1].seq)
#        # Identify the indels....
#        #prot_start = [ i for i, (x, y) in enumerate(zip(prot_start1,prot_start2)) if x != "-" and y !="-" ] 
#        #prot_start1 = "".join([ prot_start1[i] for i in prot_start])
#        #prot_start2 = "".join([ prot_start2[i] for i in prot_start])
#        prot_length = len(prot_start1) 
#        #print(len(prot_start1))
#        #if ("-" in prot_start1):
#        #    print(prot_start.index("-"))
#        #dna_start1 =  "".join([ dna_start1[(i*3):((i*3)+3)] for i in prot_start])
#        #print(dna_start1)
#        #dna_start2 =  "".join([ dna_start2[(i*3):((i*3)+3)] for i in prot_start])
#        #print(len(dna_start2))
#        #prot_start = [  for x, y in zip(prot_start1,prot_start2) if x != "-" and y !="-" ] 
#        for i in range(200):
#            idxs = np.random.choice(prot_length, prot_length)
#            #proteins[0].seq[idxs]
#            #print(idxs)
#            prot1 = "".join([prot_start1[i] for i in idxs])
#            prot2 = "".join([prot_start2[i] for i in idxs])
#            proteins[0].seq = Seq.Seq(prot1, alphabet=IUPAC.IUPACProtein()) 
#            proteins[1].seq = Seq.Seq(prot2, alphabet=IUPAC.IUPACProtein()) 
#            # bootstrap DNA.
#            seq_by.seq = Seq.Seq( "".join([ dna_start2[(i*3):((i*3)+3)] for i in idxs ] ),alphabet=IUPAC.IUPACUnambiguousDNA()) 
#            seq_cbs.seq = Seq.Seq( "".join([ dna_start1[(i*3):((i*3)+3)] for i in idxs ]), alphabet=IUPAC.IUPACUnambiguousDNA()) 
#            #print(proteins)
#            codon_aln = codonalign.build(proteins,[seq_cbs, seq_by], alphabet=codonalign.default_codon_alphabet)
#            #dna_alignment = AlignIO.read(ouput_file,"clustal")
#            #calculator = DistanceCalculator('blosum62')
#            #dm = calculator.get_distance(dna_alignment)
#            dn_ds= (codonalign.codonseq.cal_dn_ds(codon_aln[0],codon_aln[1]))
#            if(dn_ds[1] == -1):
#                dn_array.append((max_dn_ds))
#            else:
#                dn_array.append((dn_ds[1]))
#        print(np.var(dn_array))
#        print(np.median(dn_array))
#        return(dn_array)
#    else:
#        proteins[0].seq = Seq.Seq(prot_start1, alphabet=IUPAC.IUPACProtein()) 
#        proteins[1].seq = Seq.Seq(prot_start2, alphabet=IUPAC.IUPACProtein()) 
#        ouput_file = out_file 
#        print(ouput_file)
#        codon_aln = codonalign.build(proteins,[seq_cbs, seq_by], alphabet=codonalign.default_codon_alphabet)
#        dna_alignment = AlignIO.read(ouput_file,"clustal")
#        calculator = DistanceCalculator('identity')
#        dm = calculator.get_distance(dna_alignment)
#        dn_ds= (codonalign.codonseq.cal_dn_ds(codon_aln[0],codon_aln[1]))
#        return([name, dn_ds[0], dn_ds[1], dm[0,1]])

def get_dn_ds_from_seq(gene,bootstrap,seq,method="NG86",genes_in_fasta=False,all_genes=False, skip_alignment=False):
    #gal_genes_cbs = pyfasta.Fasta(fasta_input)
    #if(all_genes):
    #print(gene)
    #prin#t(seq)
    #    genes = gal_genes_cbs.keys()
    dn_ds = {}
    name = seq[0]
    seq= seq[1]
    if genes_in_fasta:
        tmp_seq = seq
       # print(tmp_seq[0])
      #  print(tmp_seq[1])
        seq = str(tmp_seq[0].seq)
        #print(seq)
     #   print("HERE")
        by_sequence = str(tmp_seq[1].seq)
        #print(by_sequence)
        seq = seq.replace("-","")
        by_sequence=  by_sequence.replace("-","")
        seq = seq.replace("N","")
        by_sequence=  by_sequence.replace("N","")
    else:
        s288c_protein_hash = _get_protein_hash()
        by_sequence = s288c_protein_hash[gene]
    #for gene in genes:
    #name = gene
    seq =str(seq)
    #print(seq)
    #print(seq)
    seq_list = (extend_ambiguous_dna(seq))
    seq = random.choice(seq_list)
    seq_cbs = SeqRecord(Seq.Seq(seq, alphabet=IUPAC.IUPACUnambiguousDNA()), id=name)
    seq_by = SeqRecord(Seq.Seq(by_sequence,alphabet=IUPAC.IUPACUnambiguousDNA()),id='BY')
    prot_cbs = seq_cbs.seq.translate(to_stop=True)
    prot_by = seq_by.seq.translate(to_stop=True)
#    print(prot_cbs)
#    print(prot_by)
    in_file_prot = "tmp3/"+ name + ".in.prot.fasta"
    with open(in_file_prot, "w") as tmp_fasta:
        tmp_fasta.write(">"+name+"\n")
        tmp_fasta.write(str(prot_cbs) + "\n")
        tmp_fasta.write(">BY\n")
        tmp_fasta.write(str(prot_by))
    out_file_prot= "tmp3/" + name + ".out.prot.fasta" 
    #print(out_file_prot)
    muscle_cline = ClustalOmegaCommandline(infile=in_file_prot, outfile=out_file_prot,force=True,outfmt="clu")
    subprocess.check_call(str(muscle_cline),shell=True)
    proteins = AlignIO.read(out_file_prot,'clustal',alphabet=IUPAC.protein)
    prot_start1 = str(proteins[0].seq)
    prot_start2 = str(proteins[1].seq)
    dna_start1 = str(seq_cbs.seq)
    dna_start2 = str(seq_by.seq)
    #print(len(dna_start1))
    #print(len(dna_start2))
    #print(prot_start1)
    #print(prot_start2)
    # Identify the indels....
    #print(prot_start1)
    prot_start1_idx = [ i for i, x in enumerate(prot_start1) if x != "-" ]
    prot_start1_idx = ([ i for i, x in enumerate(prot_start1) if x == "-" ])
    prot_start2_idx = ([ i for i, x in enumerate(prot_start2) if x == "-" ])
    #prot_start2_idx = [ i for i, x in enumerate(prot_start2) if x != "-" ] 
    prot_start = [ i for i, (x, y) in enumerate(zip(prot_start1,prot_start2)) if x != "-" and y !="-" ] 
    prot_start = [ i for i, (x, y) in enumerate(zip(prot_start1,prot_start2)) if x != "-" and y !="-" ] 
    print(len(prot_start))
    prot_start1 = "".join([ prot_start1[i] for i in prot_start])
    #print(prot_start1)
    #print(prot_start2)
    #print(prot_start2[281:295])
    #print(dna_start2[(280*3):(281*3 +6)])
    prot_start2 = "".join([ prot_start2[i] for i in prot_start])
    #print(prot_start2)
    ## Calculate offset....
    #print(prot_start1_idx)
    #print(prot_start2_idx)
    if len(prot_start2_idx) != 0 and len(prot_start1_idx) !=0:
        count_less_than = []
        for pos in prot_start1_idx: 
            aa = sum(i < pos for i in prot_start2_idx)
            count_less_than.append( pos - aa )
        prot_start1_idx = count_less_than
        count_less_than = []
        for pos in prot_start2_idx: 
            aa = sum(i < pos for i in prot_start1_idx)
            count_less_than.append( pos - aa )
     #   print(count_less_than)
        prot_start2_idx = count_less_than

    dna_start1 = "".join([dna_start1[(x*3):(x*3 + 3)] for x in range(len(proteins[0].seq)) if x not in prot_start2_idx]) 
    dna_start2 = "".join([dna_start2[(x*3):(x*3 + 3)] for x in range(len(proteins[1].seq)) if x not in prot_start1_idx])
#    print(len(dna_start2))
    if dna_start1.endswith("TAG") or dna_start1.endswith("TAA") or dna_start1.endswith("TGA"):
        dna_start1 = dna_start1[:(len(dna_start1)-3)]
    if dna_start2.endswith("TAG") or dna_start2.endswith("TAA") or dna_start2.endswith("TGA"):
      #  print("HER")
        dna_start2 = dna_start2[:(len(dna_start2)-3)]
    seq_cbs.seq = (Seq.Seq(dna_start1, alphabet=IUPAC.IUPACUnambiguousDNA())) 
    seq_by.seq   = (Seq.Seq(dna_start2, alphabet=IUPAC.IUPACUnambiguousDNA())) 
    in_file = "tmp3/"+ name + ".in.fasta"
    with open(in_file, "w") as tmp_fasta:
        tmp_fasta.write(">"+name+"\n")
        tmp_fasta.write(str(seq_cbs.seq) + "\n")
        tmp_fasta.write(">BY\n")
        tmp_fasta.write(str(seq_by.seq))
    out_file = "tmp3/" +name +".out.fasta" 
    muscle_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file,force=True,outfmt="clu")
    #print(muscle_cline)
    out = subprocess.Popen(str(muscle_cline),shell=True)
    out.communicate()
    print(muscle_cline)
    #sys.exit(1)
    dna = AlignIO.read(out_file,'clustal',alphabet=IUPAC.IUPACUnambiguousDNA())
    if(bootstrap):
        max_dn_ds = 3.23
        dn_array = [] 
        #prot_start1 = str(proteins[0].seq)
        #prot_start2 = str(proteins[1].seq)
        #dna_start1 = str(dna[0].seq)
        #dna_start2 = str(dna[1].seq)
        # Identify the indels....
        #prot_start = [ i for i, (x, y) in enumerate(zip(prot_start1,prot_start2)) if x != "-" and y !="-" ] 
        #prot_start1 = "".join([ prot_start1[i] for i in prot_start])
        #prot_start2 = "".join([ prot_start2[i] for i in prot_start])
        prot_length = len(prot_start1) 
        #print(len(prot_start1))
        #if ("-" in prot_start1):
        #    print(prot_start.index("-"))
        #dna_start1 =  "".join([ dna_start1[(i*3):((i*3)+3)] for i in prot_start])
        #print(dna_start1)
        #dna_start2 =  "".join([ dna_start2[(i*3):((i*3)+3)] for i in prot_start])
        #print(len(dna_start2))
        #prot_start = [  for x, y in zip(prot_start1,prot_start2) if x != "-" and y !="-" ] 
        for i in range(200):
            idxs = np.random.choice(prot_length, prot_length)
            #proteins[0].seq[idxs]
            #print(idxs)
            prot1 = "".join([prot_start1[i] for i in idxs])
            prot2 = "".join([prot_start2[i] for i in idxs])
            proteins[0].seq = Seq.Seq(prot1, alphabet=IUPAC.IUPACProtein()) 
            proteins[1].seq = Seq.Seq(prot2, alphabet=IUPAC.IUPACProtein()) 
            # bootstrap DNA.
            seq_by.seq = Seq.Seq( "".join([ dna_start2[(i*3):((i*3)+3)] for i in idxs ] ),alphabet=IUPAC.IUPACUnambiguousDNA()) 
            seq_cbs.seq = Seq.Seq( "".join([ dna_start1[(i*3):((i*3)+3)] for i in idxs ]), alphabet=IUPAC.IUPACUnambiguousDNA()) 
            #print(proteins)
            codon_aln = codonalign.build(proteins,[seq_cbs, seq_by], alphabet=codonalign.default_codon_alphabet)
            #dna_alignment = AlignIO.read(ouput_file,"clustal")
            #calculator = DistanceCalculator('blosum62')
            #dm = calculator.get_distance(dna_alignment)
            dn_ds= (codonalign.codonseq.cal_dn_ds(codon_aln[0],codon_aln[1]))
            if(dn_ds[1] == -1):
                dn_array.append((max_dn_ds))
            else:
                dn_array.append((dn_ds[1]))
        #print(np.var(dn_array))
        #print(np.median(dn_array))
        return(dn_array)
    else:
        proteins[0].seq = Seq.Seq(prot_start1, alphabet=IUPAC.IUPACProtein()) 
        proteins[1].seq = Seq.Seq(prot_start2, alphabet=IUPAC.IUPACProtein()) 
        #print(proteins[1].seq)
        ouput_file = out_file 
        #print("BLAH")
        #print(seq_by.seq)
        #print(seq_by.seq.translate())
        #print([(x,y) for x,y in zip(proteins[1].seq,seq_by.seq.translate()) if x!=y])
        #print([(x,y) for x,y in zip(proteins[0].seq,seq_cbs.seq.translate()) if x!=y])
        codon_aln = codonalign.build(proteins,[seq_cbs, seq_by], alphabet=codonalign.default_codon_alphabet)
        dna_alignment = AlignIO.read(ouput_file,"clustal")
        calculator = DistanceCalculator('blosum62')
        #print(dna_alignment)
        dm = calculator.get_distance(dna_alignment)
        #print(dm)
        #print(codon_aln[0])
        #print(codon_aln[1])
        import pickle
        with open("pickle.dump", "wb") as out_pick:
            pickle.dump(codon_aln,out_pick)
        #print(codonalign.codonseq.cal_dn_ds(codon_aln[0],codon_aln[1]))
        dn_ds= (codonalign.codonseq.cal_dn_ds(codon_aln[0],codon_aln[1]))
        dm = (np.array(dm)[0,1])
        #print(dm)
        #print(name)
        #print(dn_ds[0])
        #print(dn_ds[1])
        #print([name, dn_ds[0], dn_ds[1], dm])
        return([name, dn_ds[0], dn_ds[1], dm])

def get_dn_ds(fasta_input, genes,s288c_protein_hash,method="NG86",all_genes=False):
    gal_genes_cbs = pyfasta.Fasta(fasta_input)
    if(all_genes):
        genes = gal_genes_cbs.keys()
    dn_ds = {}

    for gene in genes:
        seq=  str(gal_genes_cbs[gene])
        name = gene
        by_sequence = s288c_protein_hash[gene]
        seq_cbs = SeqRecord(Seq(seq, alphabet=IUPAC.IUPACUnambiguousDNA()),id=' CBS')
        seq_by = SeqRecord(Seq(by_sequence,alphabet=IUPAC.IUPACUnambiguousDNA()),id='BY')
        prot_cbs = seq_cbs.seq.translate(to_stop=True)
        prot_by = seq_by.seq.translate(to_stop=True)
        with open("tmp.in.fasta", "w") as tmp_fasta:
            tmp_fasta.write(">CBS\n")
            tmp_fasta.write(str(seq_cbs.seq) + "\n")
            tmp_fasta.write(">BY\n")
            tmp_fasta.write(str(seq_by.seq))
        in_file = "tmp.in.fasta"
        out_file = "tmp.out.dna.fasta"
        muscle_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file,force=True,outfmt="clu")
        subprocess.check_call(str(muscle_cline),shell=True)
        dna = AlignIO.read("tmp.out.dna.fasta",'clustal',alphabet=IUPAC.IUPACUnambiguousDNA())
        print(dna)
        print(seq_by)
        if(len(prot_cbs) < len(prot_by)*0.9 or len(prot_cbs) > len(prot_by)*1.1):
            print("Skipped gene " + name + " its too short")
            continue
        with open("tmp.fasta", "w") as tmp_fasta:
            tmp_fasta.write(">CBS\n")
            tmp_fasta.write(str(prot_cbs) + "\n")
            tmp_fasta.write(">BY\n")
            tmp_fasta.write(str(prot_by))
        in_file = "tmp.fasta"
        out_file = "tmp.out.fasta"
        muscle_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file,force=True,outfmt="clu")
        subprocess.check_call(str(muscle_cline),shell=True)
        proteins = AlignIO.read("tmp.out.fasta",'clustal',alphabet=IUPAC.protein)
        codon_aln = codonalign.build(proteins,[seq_cbs, seq_by], alphabet=codonalign.default_codon_alphabet)
        ouput_file = "tmp.out.dna.fasta"
        dna_alignment = AlignIO.read(ouput_file,"clustal")
        calculator = DistanceCalculator('blosum62')
        dm = calculator.get_distance(dna_alignment)
        dn_ds[name] = (codonalign.codonseq.cal_dn_ds(codon_aln[0],codon_aln[1]),dm[0,1])
        print(dn_ds)
    return(dn_ds)
                

def get_dn_ds_windowed_alignments(fasta_input_list,window=100,step=50,output_folder="outputs/windowed"):
    step = step * 3 
    window = window * 3
    try:
        os.makedirs(output_folder)
    except:
        pass
    for fasta in fasta_input_list:
        print(fasta)
        gene_fasta = pyfasta.Fasta(fasta)
        in_file = fasta
        out_file = os.path.join(output_folder, os.path.basename(in_file))
        muscle_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file,force=True,outfmt="fasta")
        subprocess.check_call(str(muscle_cline),shell=True)
        dna = AlignIO.read(out_file,'fasta',alphabet=IUPAC.IUPACUnambiguousDNA())
        dna_list = {}
        prot_list = {}
        for dn in dna:
            # Identify codons and remove ones having Rs and Ys.
            codons = str(dn.seq).replace("-","")
            idxs = [] 
            for i in range(0,len(codons),3):
                codon_tmp = codons[i:(i+3)]
                if "R" in codon_tmp or "Y" in codon_tmp:
                    idxs.extend(range(i,i+3))
            codons = list(codons)
            codons = [x for i,x in enumerate(codons) if i not in idxs]
            codons = "".join(codons)
            dna_list[dn.id] = SeqRecord(Seq.Seq(codons,alphabet=IUPAC.IUPACUnambiguousDNA()), id=dn.id)
            prot_list[dn.id] = SeqRecord(dna_list[dn.id].seq.translate(to_stop=True), id=dn.id)
        prot_out = out_file + ".aa" 
        with open(out_file + ".aa" , "w") as tmp_fasta:
            for key, seq in prot_list.items():
                tmp_fasta.write(">" + key + "\n" + str(seq.seq) + "\n")
        i=0
        out_file=out_file +".aligned.fasta"
        muscle_cline = ClustalOmegaCommandline(infile=prot_out, outfile=out_file,force=True,outfmt="clu")
        subprocess.check_call(str(muscle_cline),shell=True)
        proteins = AlignIO.read(out_file,'clustal',alphabet=IUPAC.protein)
        #out_file = "tmp.out.fasta"
        i = i + step
        mapping = { dn.id:dn.id for dn in dna}
        dna_list = [dn for key, dn in dna_list.items()]          
        codon_align = codonalign.build(proteins, dna_list, corr_dict=mapping, alphabet=codonalign.default_codon_alphabet)
        dn_matrix, ds_matrix = codon_align.get_dn_ds_matrix()
        print(dn_matrix)
        print(ds_matrix)
        ouput_file = os.path.join(output_folder, os.path.basename(in_file))
        dna_alignment = AlignIO.read(ouput_file,"fasta")
        calculator = DistanceCalculator('blosum62')
        dm = calculator.get_distance(dna_alignment)
       # with open("output.phyllip","w") as out_w:
        #    ds_matrix.format_phylip(out_w)
        print("PIECEMEAL")
        return(ds_matrix) 
        while i < (len(dna[0].seq) - window):
            dna_list = {}
            prot_list = {}
            out_file= os.path.basename(fasta) 
            print(i)
            for dn in dna:
                dna_list[dn.id] = SeqRecord(Seq(str(dn.seq)[i:(i+window)].replace("-",""), alphabet=IUPAC.IUPACUnambiguousDNA()), id=dn.id)
                prot_list[dn.id] = SeqRecord(dna_list[dn.id].seq.translate(to_stop=True), id=dn.id)
            prot_out = output_folder + "/" + out_file + "_" + str(i) + ".aa"
            with open(output_folder + "/" + out_file + "_" + str(i) , "w") as tmp_fasta:
                for key, seq in dna_list.items():
                    tmp_fasta.write(">" + key + "\n" + str(seq.seq) + "\n")

            with open(prot_out , "w") as tmp_fasta:
                for key, seq in prot_list.items():
                    tmp_fasta.write(">" + key + "\n" + str(seq.seq) + "\n")
            #in_file = "tmp.fasta"
            try:
                in_file = output_folder + "/" + out_file + "_" + str(i) 
                ouput_file = os.path.join(output_folder, os.path.basename(in_file) + ".dna.fasta")
                muscle_cline = ClustalOmegaCommandline(infile=in_file, outfile=ouput_file,force=True,outfmt="fasta")
                subprocess.check_call(str(muscle_cline),shell=True)
                out_file="tmp.out.fasta"
                muscle_cline = ClustalOmegaCommandline(infile=prot_out, outfile=out_file,force=True,outfmt="clu")
                subprocess.check_call(str(muscle_cline),shell=True)
            except:
                #TODO fix
                break 
            i = i + step
            proteins = AlignIO.read("tmp.out.fasta",'clustal',alphabet=IUPAC.protein)
            #out_file = "tmp.out.fasta"
            mapping = { dn.id:dn.id for dn in dna}
            dna_list = [dn for key, dn in dna_list.items()]           
            codon_align = codonalign.build(proteins, dna_list, corr_dict=mapping, alphabet=codonalign.default_codon_alphabet)
            dn_matrix, ds_matrix = codon_align.get_dn_ds_matrix()
            print(dn_matrix)
            print(ds_matrix)
            dna_alignment = AlignIO.read(ouput_file,"fasta")
            calculator = DistanceCalculator('blosum62')
            dm = calculator.get_distance(dna_alignment)
            # Make sure the last step includes the final parts of the alignment. We have a 3' flanking reverse gene here.


def _get_hash_from_fasta(fasta):
    s288c_protein_hash = {}
    with open(fasta) as fasta_in:
        tmp_line = fasta_in.readline()
        i = 1
        while(tmp_line):
            name = tmp_line.split(">")[1].split(" ")[0]
            seq = ""
            tmp_line=fasta_in.readline().strip()
            while(">" not in tmp_line and tmp_line):
                seq = seq + tmp_line
                tmp_line = fasta_in.readline().strip()
            s288c_protein_hash[name] = seq
            i = i  + 1
    return(s288c_protein_hash)


def _get_protein_hash():
    s288c_protein_hash = _get_hash_from_fasta("/u/home/s/smilefre/project-kruglyak/ref/orf_coding_all.fasta") 
    return(s288c_protein_hash)

def _get_bam_hash():
    s288c_protein_hash = {}
    with open("//media/theboocock/data/Dropbox/PHDTHESIS/projects/gal_final_github/popgen_phylo_notebooks/data/fastas/bam_genes.fasta") as fasta_in:
        tmp_line = fasta_in.readline()
        i = 1
        while(tmp_line):
            name = tmp_line.strip().split(">")[1]
            seq = ""
            tmp_line=fasta_in.readline().strip()
            while(">" not in tmp_line and tmp_line):
                seq = seq + tmp_line
                tmp_line = fasta_in.readline().strip()
            s288c_protein_hash[name] = seq
    return(s288c_protein_hash)


def _get_cbs_hash():
    s288c_protein_hash = {}
    with open("//media/theboocock/data/Dropbox/PHDTHESIS/projects/gal_final_github/popgen_phylo_notebooks/data/fastas/cbs.fasta") as fasta_in:
        tmp_line = fasta_in.readline()
        i = 1
        while(tmp_line):
            name = tmp_line.strip().split(">")[1]
            seq = ""
            tmp_line=fasta_in.readline().strip()
            while(">" not in tmp_line and tmp_line):
                seq = seq + tmp_line
                tmp_line = fasta_in.readline().strip()
            s288c_protein_hash[name] = seq
    return(s288c_protein_hash)

import pyfasta

def windowed_divergence(alignment, do_window=False, window=100, fasta_reference="/media/theboocock/data//PHDTHESIS/Single-CellRNASEQ/eqtls/ref_data/yeast/reference_genome/sacCer3.fasta",sample="bam"):
   
    OUT_BLAST="reference_blast"
    blast_cmd="makeblastdb -in {0} -dbtype 'nucl' -out reference_blast".format(fasta_reference)
    #subprocess.check_call(blast_cmd, shell=True)
    fasta_align = pyfasta.Fasta(alignment)
    blast_search = """blastn -outfmt "6 qseqid sseqid pident qstart qend qlen length mismatch gapope evalue bitscore sstart send sstrand" -db reference_blast -query {0} -num_alignments 15 -word_size 10""".format(alignment)
    blast_output = subprocess.check_output(blast_search, shell=True)
    query_length = len(fasta_align[sample])
    end = False
    start = False
    start_ref = None
    end_ref =None
    for line in blast_output.splitlines():
        line_s = line.decode().split()
        percent = float(line_s[2])
        query_end = int(line_s[4])
        query_start = int(line_s[3])
        subject_end = int(line_s[11])
        subject_start = int(line_s[10])
        ref = line_s[1]
        if query_end >= 0.9 * query_length: 
            if subject_end <= subject_start:
                end = subject_start
            else:
                end = subject_end
            end_ref = ref
        if query_start <= 0.05  * query_length: 
            if subject_start >= subject_end:
                start = subject_end
            else:
                start = subject_start
            start_ref = ref
        if start != False and end != False: 
            break
    fasta_reference = pyfasta.Fasta(fasta_reference)
    if start > end: 
        tmp = end 
        end = start
        start = tmp
    if start_ref == end_ref:
        ref_seq_region = fasta_reference[end_ref][start:end]
    with open("tmp.fasta.out","w") as in_f:
        in_f.write(">REF\n")
        in_f.write(ref_seq_region + "\n")
        with open(alignment) as in_align:
            for line in in_align:
                in_f.write(line) 
    muscle_cline = MuscleCommandline("muscle", out="align.tmp.fasta", input="tmp.fasta.out")
    #subprocess.check_call(str(muscle_cline),shell=True) 
    alignments  = AlignIO.read("align.tmp.fasta",'fasta',alphabet=IUPAC.IUPACUnambiguousDNA())
    calculator = DistanceCalculator('blosum62')
    dm = calculator.get_distance(alignments)
    for i in range(int(len(alignments[0])/100)):
        new_alignments =alignments[:,(i*100):((i*100)+100)]
        dm = calculator.get_distance(new_alignments)



def get_dn_ds_from_alignment(alignment,do_window=False, window=100, step=10,method="NG86", these_samples=None, gene_name="YBR018C",cbs_reference=False, distance_only=False,hoffman=False):
    if not cbs_reference:
        s288c_protein_hash = _get_protein_hash()
    else:
        s288c_protein_hash = _get_cbs_hash()
    gal7_sliding_window = pyfasta.Fasta(alignment)
    dn_ds = {}
    step = step *3
    window = window * 3
    try:
        reference_gene = s288c_protein_hash[gene_name] 
    except:
        return None
    try:
        os.mkdir("tmp3/")
    except:
        pass
    import tempfile
    temp = tempfile.NamedTemporaryFile(dir="tmp3/",delete=False)

    #from_input= str(item)
    
    with open(temp.name,"w") as tmp_fasta:
        tmp_fasta.write(">REF\n")
        tmp_fasta.write(reference_gene +"\n")
        for key, item in gal7_sliding_window.items():
            if key in these_samples:
                tmp_fasta.write(">" + key + "\n")
                item = str(item).upper()
                seq_list = (extend_ambiguous_dna((item)))
                seq = random.choice(seq_list)
                tmp_fasta.write(seq + "\n")
    import uuid
    output_filename= os.path.join("tmp3/",str(uuid.uuid4()))
    muscle_cline = ClustalOmegaCommandline(infile=temp.name, outfile=output_filename, force=True, outfmt="clu")
    if hoffman:
        muscle_cline = "/u/project/kruglyak/smilefre/anaconda3/bin/"+ str(muscle_cline)
    subprocess.check_call(str(muscle_cline),shell=True) 
    alignments  = AlignIO.read(output_filename,'clustal',alphabet=IUPAC.IUPACUnambiguousDNA())
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignments)
    #proteins =
    #sys.exit(1)
    os.remove(temp.name)
    os.remove(output_filename)
    if distance_only:
        return(dm[0,1])
    dn_ds_list = {}
    for alignment in alignments:
        if alignment.id == "REF":
            ref_alignment = list(alignment.seq)
            seq_by = alignment 
            ref_translate = SeqRecord(Seq.Seq(str(alignment.seq).replace("-",""), alphabet=IUPAC.IUPACUnambiguousDNA()), id="REF")
            ref_protein = ref_translate.seq.translate(to_stop=True)
            ref_protein = SeqRecord(ref_protein, id="REF")
        else:
            #("".join(ref_alignment))
            #tmp_alignment = "".join([char if char != "-" else ref_alignment[i] for i, char in enumerate(list(alignment.seq))])
            tmp_alignment = str(alignment.seq) 
            #print(tmp_alignment)
            is_full_protein = (len(tmp_alignment)%3)
            if is_full_protein != 0:
                return(None)
            dn_ds_list[alignment.id] = []
            dna_cheese_sample = SeqRecord(Seq.Seq(tmp_alignment, alphabet=IUPAC.IUPACUnambiguousDNA()), id=alignment.id)
            if not str(dna_cheese_sample.seq).startswith("ATG"):
                return(None)
            seq_cbs = dna_cheese_sample
            i=0
            j=0
            test = Seq.Seq(str(dna_cheese_sample.seq).replace("-",""))
            if(len(str(test)) % 3 != 0):
                return None
            prot_cheese_sample = SeqRecord(test.translate(to_stop=True),id=alignment.id)
            #print(test)
            input_filename= os.path.join("tmp3/",str(uuid.uuid4()))
            output_filename= os.path.join("tmp3/",str(uuid.uuid4()))
            with open(input_filename, "w") as tmp_fasta:
                tmp_fasta.write(">CBS\n")
                tmp_fasta.write(str(prot_cheese_sample.seq) + "\n")
                tmp_fasta.write(">BY\n")
                tmp_fasta.write(str(ref_protein.seq))
            in_file = input_filename 
            out_file = output_filename 
            muscle_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file,force=True,outfmt="clu")
            if hoffman:
                muscle_cline = "/u/project/kruglyak/smilefre/anaconda3/bin/"+ str(muscle_cline)
            subprocess.check_call(str(muscle_cline),shell=True)
            proteins = AlignIO.read(out_file,'clustal',alphabet=IUPAC.protein)
            os.remove(input_filename)
            os.remove(output_filename)
            prot_start1 = str(proteins[0].seq)
            prot_start2 = str(proteins[1].seq)
            dna_start1 = str(test)
            dna_start2 = str(ref_translate.seq)
            #print(len(dna_start1))
            #print(len(dna_start2))
            #print(prot_start1)
            #print(prot_start2)
            # Identify the indels....
            #print(prot_start1)
            prot_start1_idx = [ i for i, x in enumerate(prot_start1) if x != "-" ]
            prot_start1_idx = ([ i for i, x in enumerate(prot_start1) if x == "-" ])
            prot_start2_idx = ([ i for i, x in enumerate(prot_start2) if x == "-" ])
            #prot_start2_idx = [ i for i, x in enumerate(prot_start2) if x != "-" ] 
            prot_start = [ i for i, (x, y) in enumerate(zip(prot_start1,prot_start2)) if x != "-" and y !="-" ] 
            prot_start1 = "".join([ prot_start1[i] for i in prot_start])
            #print(prot_start1)
            #print(prot_start2)
            #print(prot_start2[281:295])
            #print(dna_start2[(280*3):(281*3 +6)])
            prot_start2 = "".join([ prot_start2[i] for i in prot_start])
            #print(prot_start2)
            ## Calculate offset....
            #print(prot_start1_idx)
            #print(prot_start2_idx)
            if len(prot_start2_idx) != 0 and len(prot_start1_idx) !=0:
                count_less_than = []
                for pos in prot_start1_idx: 
                    aa = sum(i < pos for i in prot_start2_idx)
                    count_less_than.append( pos - aa )
                prot_start1_idx = count_less_than
                count_less_than = []
                for pos in prot_start2_idx: 
                    aa = sum(i < pos for i in prot_start1_idx)
                    count_less_than.append( pos - aa )
             #   print(count_less_than)
                prot_start2_idx = count_less_than

            dna_start1 = "".join([dna_start1[(x*3):(x*3 + 3)] for x in range(len(proteins[0].seq)) if x not in prot_start2_idx]) 
            dna_start2 = "".join([dna_start2[(x*3):(x*3 + 3)] for x in range(len(proteins[1].seq)) if x not in prot_start1_idx])
        #    print(len(dna_start2))
            if dna_start1.endswith("TAG") or dna_start1.endswith("TAA") or dna_start1.endswith("TGA"):
                dna_start1 = dna_start1[:(len(dna_start1)-3)]
            if dna_start2.endswith("TAG") or dna_start2.endswith("TAA") or dna_start2.endswith("TGA"):
              #  print("HER")
                dna_start2 = dna_start2[:(len(dna_start2)-3)]
            seq_cbs.seq = (Seq.Seq(dna_start1, alphabet=IUPAC.IUPACUnambiguousDNA())) 
            seq_by.seq   = (Seq.Seq(dna_start2, alphabet=IUPAC.IUPACUnambiguousDNA())) 
            in_file = "tmp3/"+ str(uuid.uuid4())
            with open(in_file, "w") as tmp_fasta:
                tmp_fasta.write(">CBS\n")
                tmp_fasta.write(str(seq_cbs.seq) + "\n")
                tmp_fasta.write(">BY\n")
                tmp_fasta.write(str(seq_by.seq))
            out_file= "tmp3/"+ str(uuid.uuid4())
            #out_file = "tmp3/" +name +".out.fasta" 
            muscle_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file,force=True,outfmt="clu")
            if hoffman:
                muscle_cline = "/u/project/kruglyak/smilefre/anaconda3/bin/"+ str(muscle_cline)
            subprocess.check_call(str(muscle_cline),shell=True)
            dna = AlignIO.read(out_file,'clustal',alphabet=IUPAC.IUPACUnambiguousDNA())
            proteins[0].seq = Seq.Seq(prot_start1, alphabet=IUPAC.IUPACProtein()) 
            proteins[1].seq = Seq.Seq(prot_start2, alphabet=IUPAC.IUPACProtein()) 
            try:
                codon_aln = codonalign.build(proteins,[dna[0], dna[1]], alphabet=codonalign.default_codon_alphabet)
            except:
                os.remove(out_file)
                os.remove(in_file)
                return(None)  
            try:
                dn_ds_all =  (codonalign.codonseq.cal_dn_ds(codon_aln[0],codon_aln[1],method=method))
            except:
                return(None)
            os.remove(in_file)
            os.remove(out_file)
            if do_window:
                in_var  = True 
                while i < max(len(ref_protein) * 3 -window, len(prot_cheese_sample)* 3 - window):
                    tmp_seq_cbs = str(dna_cheese_sample.seq)[i:(i + window)].replace("-","")
                    tmp_seq_by = str(ref_translate.seq)[i:(i + window)].replace("-","")
                    seq_cbs = SeqRecord(Seq.Seq(tmp_seq_cbs, alphabet=IUPAC.IUPACUnambiguousDNA()),id=' CBS')
                    seq_by = SeqRecord(Seq.Seq(tmp_seq_by,alphabet=IUPAC.IUPACUnambiguousDNA()),id='BY')
                    if in_var:
                        #print(seq_by.seq)
                        #print(seq_cbs.seq)
                        in_var  = False 
                    prot_cbs = seq_cbs.seq.translate(to_stop=True)
                    prot_by = seq_by.seq.translate(to_stop=True)
                    if(len(prot_cbs) <= 100 or len(prot_by) <= 100):
                        break
                    #print(tmp_seq_cbs)
                    j = j + window * 3
                    i = i + step
                    #print(len(prot_cbs)/window)
                    #print(len(prot_by)/window)
                    input_filename= os.path.join("tmp3/",str(uuid.uuid4()))
                    output_filename= os.path.join("tmp3/",str(uuid.uuid4()))
                    with open(input_filename, "w") as tmp_fasta:
                        tmp_fasta.write(">CBS\n")
                        tmp_fasta.write(str(prot_cbs) + "\n")
                        tmp_fasta.write(">BY\n")
                        tmp_fasta.write(str(prot_by))
                    in_file = input_filename 
                    out_file = output_filename 
                    muscle_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file,force=True,outfmt="clu")
                    if hoffman:
                        muscle_cline = "/u/project/kruglyak/smilefre/anaconda3/bin/"+ str(muscle_cline)
                    subprocess.check_call(str(muscle_cline),shell=True)
                    proteins = AlignIO.read(out_file,'clustal',alphabet=IUPAC.protein)
                    try:
                        codon_aln = codonalign.build(proteins,[seq_cbs, seq_by], alphabet=codonalign.default_codon_alphabet)
                    except:
                        os.remove(out_file)
                        os.remove(in_file)
                        dn_ds_list[alignment.id].append(("NA","NA"))
                        continue
                    try:
                        ds =  (codonalign.codonseq.cal_dn_ds(codon_aln[0],codon_aln[1],method=method))
                    except:
                        dn_ds_list[alignment.id].append(("NA","NA"))
                        continue
                    if ds[1] == -1 or ds[1] > 3.23:
                        ds= (ds[0], 3.23)
                    dn_ds_list[alignment.id].append(ds)
                    os.remove(out_file)
                    os.remove(in_file)
                if i == 0:
                    tmp_seq_cbs = str(dna_cheese_sample.seq)[i:(i + window)].replace("-","")
                    tmp_seq_by = str(ref_translate.seq)[i:(i + window)].replace("-","")
                    seq_cbs = SeqRecord(Seq.Seq(tmp_seq_cbs, alphabet=IUPAC.IUPACUnambiguousDNA()),id=' CBS')
                    seq_by = SeqRecord(Seq.Seq(tmp_seq_by,alphabet=IUPAC.IUPACUnambiguousDNA()),id='BY')
                    if in_var:
                    #    print(seq_by.seq)
                    #    print(seq_cbs.seq)
                        in_var  = False 
                    prot_cbs = seq_cbs.seq.translate(to_stop=True)
                    prot_by = seq_by.seq.translate(to_stop=True)
                    if(len(prot_cbs) == 0 or len(prot_by) == 0):
                        break
                    #print(tmp_seq_cbs)
                    j = j + window * 3
                    i = i + step
                    input_filename= os.path.join("tmp3/",str(uuid.uuid4()))
                    output_filename= os.path.join("tmp3/",str(uuid.uuid4()))
                    #print(len(prot_cbs)/window)
                    #print(len(prot_by)/window)
                    with open(input_filename, "w") as tmp_fasta:
                        tmp_fasta.write(">CBS\n")
                        tmp_fasta.write(str(prot_cbs) + "\n")
                        tmp_fasta.write(">BY\n")
                        tmp_fasta.write(str(prot_by))
                    in_file = input_filename 
                    out_file = output_filename 
                    muscle_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file,force=True,outfmt="clu")
                    if hoffman:
                        muscle_cline = "/u/project/kruglyak/smilefre/anaconda3/bin/"+ str(muscle_cline)
                    subprocess.check_call(str(muscle_cline),shell=True)

                    proteins = AlignIO.read(out_file,'clustal',alphabet=IUPAC.protein)
                    try:
                        codon_aln = codonalign.build(proteins,[seq_cbs, seq_by], alphabet=codonalign.default_codon_alphabet)
                    except ValueError:
                        os.remove(out_file)
                        os.remove(in_file)
                        return(None)    
                    ds =  (codonalign.codonseq.cal_dn_ds(codon_aln[0],codon_aln[1],method=method))
                    if ds[1] == -1 or ds[1] > 3.23:
                        ds= (ds[0], 3.23)
                    dn_ds_list[alignment.id].append(ds)
                    os.remove(out_file)
                    os.remove(in_file)

                if gene_name.endswith("C"):
                    dn_ds_list[alignment.id]= dn_ds_list[alignment.id][::-1]

            else:
                align = MultipleSeqAlignment([prot_cheese_sample, ref_protein])
                #print(ref_protein)

                #print(prot_cheese_sample)
                codon_aln = codonalign.build(align ,[dna_cheese_sample , ref_translate], alphabet=codonalign.default_codon_alphabet)
                dn_ds_list[alignment.id]=  (codonalign.codonseq.cal_dn_ds(codon_aln[0],codon_aln[1],method=method))
                #print(dn_ds_list)
    return(dn_ds_all, dn_ds_list)

def get_dn_ds_window(fasta_input, genes,s288c_protein_hash, window=200,step=50,method="NG86"):
    gal_genes_cbs = pyfasta.Fasta(fasta_input)
    dn_ds = {}
    step = step * 3 
    window = window * 3
    for gene in genes:
        seq=  str(gal_genes_cbs[gene])
        name = gene
        by_sequence = s288c_protein_hash[gene]
        seq_cbs = SeqRecord(Seq(seq, alphabet=IUPAC.IUPACUnambiguousDNA()),id=' CBS')
        seq_by = SeqRecord(Seq(by_sequence,alphabet=IUPAC.IUPACUnambiguousDNA()),id='BY')
        prot_cbs_all= seq_cbs.seq.translate(to_stop=True)
        prot_by_all= seq_by.seq.translate(to_stop=True)
        i=0
        j=0
        if(len(prot_cbs_all) < len(prot_by_all)*0.9 or len(prot_cbs_all) > len(prot_by_all)*1.1):
            print("Skipped gene " + name + " its too short")
            continue
        dn_ds_list = []
        with open("tmp.in.fasta", "w") as tmp_fasta:
            tmp_fasta.write(">CBS\n")
            tmp_fasta.write(str(seq_cbs.seq) + "\n")
            tmp_fasta.write(">BY\n")
            tmp_fasta.write(str(seq_by.seq))
        in_file = "tmp.in.fasta"
        out_file = "tmp.out.dna.fasta"
        muscle_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file,force=True,outfmt="clu")
        subprocess.check_call(str(muscle_cline),shell=True)
        dna = AlignIO.read("tmp.out.dna.fasta",'clustal',alphabet=IUPAC.IUPACUnambiguousDNA())
        while i < max(len(prot_cbs_all) * 3 -window,len(prot_by_all) * 3 - window):
            tmp_seq_cbs = str(dna[0].seq)[i:(i + window)].replace("-","")
            tmp_seq_by = str(dna[1].seq)[i:(i + window)].replace("-","")
            seq_cbs = SeqRecord(Seq(tmp_seq_cbs, alphabet=IUPAC.IUPACUnambiguousDNA()),id=' CBS')
            #print(tmp_seq_by)
            seq_by = SeqRecord(Seq(tmp_seq_by,alphabet=IUPAC.IUPACUnambiguousDNA()),id='BY')
            prot_cbs = seq_cbs.seq.translate(to_stop=True)
            #print(prot_cbs)
            prot_by = seq_by.seq.translate(to_stop=True)
            if(len(prot_cbs) == 0 or len(prot_by) == 0):
                break
            #print(tmp_seq_cbs)
            j = j + 300
            i = i + step
            #print(len(prot_cbs)/window)
            #print(len(prot_by)/window)
            with open("tmp.fasta", "w") as tmp_fasta:
                tmp_fasta.write(">CBS\n")
                tmp_fasta.write(str(prot_cbs) + "\n")
                tmp_fasta.write(">BY\n")
                tmp_fasta.write(str(prot_by))
            in_file = "tmp.fasta"
            out_file = "tmp.out.fasta"
            #if create_tree:
            #    return None
            muscle_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file,force=True,outfmt="clu")
            subprocess.check_call(str(muscle_cline),shell=True)
            proteins = AlignIO.read("tmp.out.fasta",'clustal',alphabet=IUPAC.protein)
            codon_aln = codonalign.build(proteins,[seq_cbs, seq_by], alphabet=codonalign.default_codon_alphabet)
            ds =  (codonalign.codonseq.cal_dn_ds(codon_aln[0],codon_aln[1],method=method))
            dn_ds_list.append(codonalign.codonseq.cal_dn_ds(codon_aln[0],codon_aln[1]))
        if name.endswith("C"):
            dn_ds_list = dn_ds_list[::-1]
        dn_ds[name] = dn_ds_list
    return(dn_ds)


def align_pgm1_promoter_region():
    in_file = "data/pgm1_region.txt"
    out_file = "tmp.out.dna.fasta"
    muscle_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file,force=True,outfmt="clu")
    subprocess.check_call(str(muscle_cline),shell=True)





