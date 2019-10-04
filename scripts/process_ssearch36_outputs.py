import os
import argparse
import pandas
import glob 
import pyfasta
def load_ssearch36_inputs(ssearch_in):
    """
        @date 03 Oct 2019
        @author James Boocock
    """
    print(ssearch_in) 
    pandas_in = pandas.read_csv(ssearch_in,sep="\t", header=None)
    (pandas_in.columns) = ["gene","strain","percentage_match","3","4","5","ref_start","ref_end","start","end","E","11"]
    return(os.path.basename(ssearch_in), pandas_in)


def process_ssearch36_df(name, ssearch_df, fasta_inputs, out_dir):
    """
        Process ssearch36 files second time mapping. 
    """
    fasta_in = [x for x in fasta_inputs if name.split(".fasta.remap.ss")[0] == os.path.basename(x).split(".fasta")[0]][0]
    fasta_file = pyfasta.Fasta(fasta_in)
    print(name) 
    try:
        os.mkdir(out_dir)
    except:
        pass
    with open(os.path.join(out_dir, os.path.basename(fasta_in)),"w") as out_f:
        for gene in (ssearch_df["gene"].unique()):
            ssearch_tmp =  ssearch_df[ssearch_df["gene"] == gene]
            gene_match = (ssearch_tmp.sort_values(by="11",ascending=False).head(n=1)["strain"])
            if (any(gene_match.isin([gene]))):
                out_f.write(">" + gene +"\n")
                out_f.write(str(fasta_file[gene]) + "\n")
            


def main():
    parser = argparse.ArgumentParser(description="process ssearch36 outputs into a pandas dataframe")
    parser.add_argument("-s","--ssearch36_inputs", dest="ssearch_in", help="ssearch36_folder") 
    parser.add_argument("-f","--input-folder", dest="in_folder", help="input folder")
    parser.add_argument("-o","--output-folder", dest="out_folder", help="output folder")
    args = parser.parse_args()
    #ssearch_inputs = glob.glob(os.path.join(args.ssearch_in_folder,"*ss"))
    fasta_inputs = glob.glob(os.path.join(args.in_folder, "*.fasta"))
    (key ,value ) = load_ssearch36_inputs(args.ssearch_in)
    process_ssearch36_df(key, value, fasta_inputs,args.out_folder)

if __name__=="__main__":
    main()
