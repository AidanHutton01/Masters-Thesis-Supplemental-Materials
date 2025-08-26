import argparse
import pandas as pd
import glob


def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, help="The path of the folder containing the files of differential expression data you want to filter")
    return parser.parse_args()


def extract_ids(path):
    dirs = glob.glob(f'{path}/*')
    cols = ['GeneID', 'Base mean', 'log2(FC)', 'StdErr', 'Wald-Stats', 'P-value', 'P-adj']
    geneids = pd.Series()
    for i in dirs:
        if '.ipynb_checkpoints' in i:
            continue
        data = pd.read_csv(f'{i}', delimiter='\t', names=cols)
        truth = data['P-adj'] < 0.05
        filtered_data = data[truth]
        filtered_data.to_csv(f'{i}_filtered.csv', index=False)
        geneids = pd.concat([geneids, filtered_data['GeneID']])
    geneids.to_csv(f'{path}/GeneIDs.csv', index=False, header=False)
    

if __name__ == '__main__':
    args = arg_parser()
    extract_ids(args.path)