"""
This script is for splitting gene count files that contain columns for multiple
samples into many files with just one sample per file. The original intent for
this script is for preprocessing GeneLab gene count data for use in the R
package DESeq2, but the function could be useful for other purposes. 
"""

import argparse
import pandas as pd

def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", type=str, help="Enter the name or directory of the gene count file you want to split. File must be in a specific format; each column must have a header, the header for the genes must be 'Geneid' (case sensitive), the other headers should represent each sample")
    parser.add_argument("output", type=str, help="Enter the name or directory of the folder you want the data to be output into")
    return parser.parse_args()

def process_file(fname, output):
    data = pd.read_csv(f'{fname}')
    for i in data:
        if i == 'Geneid':
            continue
        data[i] = pd.to_numeric(data[i], errors='raise').astype(int)
        feat_count = data[['Geneid', i]]
        feat_count.to_csv(f'{output}/{i}.csv', index=False)
        
if __name__ == '__main__':
    args = arg_parser()
    process_file(args.fname, args.output)